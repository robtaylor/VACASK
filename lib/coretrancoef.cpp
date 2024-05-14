#include "coretrancoef.h"
#include "common.h"

namespace NAMESPACE {

IntegratorCoeffs::IntegratorCoeffs(Method method, Int order, Int historyOffset) 
    : method_(method), order_(order), historyOffset_(historyOffset) {
    setOrder(order);
    b1_ = 0.0;
    err_ = 0.0;

    // // Test solver
    // matrix = {1, 2, 3, 2, 3, -5, -6, -8, 1};
    // rhs = {-7, 9, 22};
    // solve(3);
    // for(auto it : rhs) {
    //     std::cout << it*25 << " ";
    // }
    // std::cout << "\n";
    // // Result should be: -393 217 72
}

bool IntegratorCoeffs::setMethod(Method method, Int order, double xmu) {
    method_ = method;
    order_ = order;
    xmu_ = xmu;
    size();
    return xmu>=0 && xmu<=0.5;
}

bool IntegratorCoeffs::setOrder(Int order) {
    order_ = order;
    size();
    return true;
}

bool IntegratorCoeffs::size() {
    // Compute sizes of arrays
    switch (method_) {
    case Method::AdamsMoulton:
        implicit_ = true;
        numX_ = 1;
        numXdot_ = order_ - 1;
        // numX_ a coeffs + 1 (b_{-1}) + numXdot_ b coeffs
        // Don't need equation for a_0, because a_0 = 1.0
        n_ = numXdot_ + 1; 
        break;
    case Method::BDF:
        implicit_ = true;
        numX_ = order_;
        numXdot_ = 0;
        // numX_ a coeffs + 1 (b_{-1})
        n_ = numX_ + 1;
        break;
    case Method::AdamsBashforth:
        implicit_ = false;
        numX_ = 1;
        numXdot_ = order_;
        // numXdot_ b coeffs
        // Don't need equation for a_0, because a_0 = 1.0
        n_ = numXdot_;
        break;
    case Method::PolynomialExtrapolation:
        implicit_ = false;
        numX_ = order_ + 1;
        numXdot_ = 0;
        // numX_ a coeffs
        n_ = numX_;
        break;
    }
    return true;
}

// bool IntegratorCoeffs::compute(std::vector<double>& pastSteps, double newStep, Status& s) {
bool IntegratorCoeffs::compute(CircularBuffer<double>& pastSteps, double newStep, Status& s) {
    pastSteps_ = &pastSteps;
    // Compute sizes of arrays
    switch (method_) {
    case Method::AdamsMoulton:
        implicit_ = true;
        numX_ = 1;
        numXdot_ = order_ - 1;
        // numX_ a coeffs + 1 (b_{-1}) + numXdot_ b coeffs
        // Don't need equation for a_0, because a_0 = 1.0
        n_ = numXdot_ + 1; 
        break;
    case Method::BDF:
        implicit_ = true;
        numX_ = order_;
        numXdot_ = 0;
        // numX_ a coeffs + 1 (b_{-1})
        n_ = numX_ + 1;
        break;
    case Method::AdamsBashforth:
        implicit_ = false;
        numX_ = 1;
        numXdot_ = order_;
        // numXdot_ b coeffs
        // Don't need equation for a_0, because a_0 = 1.0
        n_ = numXdot_;
        break;
    case Method::PolynomialExtrapolation:
        implicit_ = false;
        numX_ = order_ + 1;
        numXdot_ = 0;
        // numX_ a coeffs
        n_ = numX_;
        break;
    }

    // Unknowns order: a_0, ..., a_{numx-1}
    //                 b_{-1}
    //                 b_0, ..., b_{numxdot-1}
    Int matSize = n_*n_;
    matrix.resize(matSize);
    rhs.resize(n_);

    // Prepare space for coeffs
    a_.resize(numX_);
    b_.resize(numXdot_);
    
    // Coeffs are 0.0 by default
    a_.assign(numX_, 0.0);
    b1_ = 0.0;
    b_.assign(numXdot_, 0.0);

    // All matrix and RHS entries will be set, so there is no need to set them to 0

    /*
    // Check if we have enough past steps
    if (pastSteps.size()<numX_-1 || pastSteps.size()<numXdot_-1) {
        s.set(Status::Internal, "Timestep history is too short.");
        return false;
    }
    
    // Compute past timepoints, index 0 is timepoint 0.0 (last computed solution)
    normalizedTimePoint.resize(pastSteps.size()+1);
    normalizedTimePoint[0] = 0.0;
    for(Int i=0;i<pastSteps.size();i++) {
        normalizedTimePoint[i+1] = normalizedTimePoint[i]-pastSteps[i];
    }

    // Normalize by newStep
    for(Int i=1; i<normalizedTimePoint.size(); i++) {
        normalizedTimePoint[i] /= newStep;
    }
    */
    if (pastSteps.valueCount()+1<numX_ || pastSteps.valueCount()+1<numXdot_) {
        s.set(Status::Internal, "Timestep history is too short.");
        return false;
    }
    
    // Compute past timepoints, index 0 is timepoint 0.0 (last computed solution)
    normalizedTimePoint.resize(pastSteps.valueCount()+1);
    normalizedTimePoint[0] = 0.0;
    for(Int i=0;i<pastSteps.valueCount();i++) {
        normalizedTimePoint[i+1] = normalizedTimePoint[i]-pastSteps.at(i);
    }

    // Normalize by newStep
    for(Int i=1; i<normalizedTimePoint.size(); i++) {
        normalizedTimePoint[i] /= newStep;
    }

    // System of equations
    // 0:
    //   sum_{i=0}^{numX-1} a_i = 1
    // j = 1..n:
    //   sum_{i=0}^{numX-1} (t_{k-i}/h_k)^j a_i + sum_{i=-1}^{numXdot-1} j (t_{k-i}/h_k)^(j-1) b_i = 1
    // after setting t_k = 0
    //   sum_{i=1}^{numX-1} (t_{k-i}/h_k)^j a_i + sum_{i=-1}^{numXdot-1} j (t_{k-i}/h_k)^(j-1) b_i = 1
    // make sure (t_{k}/h_k)^(j-1) for k=0, j=1 is 1
    switch (method_) {
    case Method::AdamsMoulton:
        // Handle single step methods without equations (constant coeffs)
        // First equation is always a_0 = 1.0, so we don't need it
        a_[0] = 1.0;
        switch (order_) {
        case 1:
            // Backward Euler
            b1_ = 1.0;
            break;
        case 2:
            // Trapezoidal
            //
            {
                b_[0] = xmu_; // 0.5;
                b1_ = 1.0-xmu_; // 0.5;

                // xmu=0.5 - pure trapezoidal
                //   b_[0] = 0.5
                //   b1_   = 0.5
                //
                // xmu=0 - pure Euler
                //   b_[0] = 0
                //   b1_   = 1
            }
            break;
        default:
            // Multistep Adams-Moulton methods
            // Fill matrix and RHS
            // t_{k}=0.0, a_0 term vanishes
            // j = 1..order
            //   sum_{i=-1}^{numXdot-1} j (t_{k-i}/h_k)^(j-1) b_i = 1
            // Unknowns: b_{-1}, b_0, b_1, ...
            Int base = 0;
            for(Int j=1; j<=order_; j++, base+=n_) {
                // Manually add b_{-1}
                matrix[base+0] = j;
                // b coeffs
                for(Int i=0; i<numXdot_; i++) {
                    // Treat b_0 differently for j=1
                    if (i==0 && j==1) {
                        matrix[base+1+i] = j;
                    } else {
                        matrix[base+1+i] = j*pow(normalizedTimePoint[i], j-1);
                    }
                }
                rhs[j-1] = 1.0;
            }
            
            // Dump system
            // std::cout << "A: ";
            // for(auto it : matrix) {
            //     std::cout << it << " ";
            // }
            // std::cout << "\nb: ";
            // for(auto it : rhs) {
            //     std::cout << it << " ";
            // }
            // std::cout << "\n";
            
            // Solve 
            solve(n_);

            // Unpack
            b1_ = rhs[0];
            for(Int i=1; i<=numXdot_; i++) {
                b_[i-1] = rhs[i];
            }
            break;
        }
        break;

    case Method::BDF:
        // Handle single step methods without equations (constant coeffs)
        switch (order_) {
        case 1:
            // Backward Euler
            a_[0] = 1.0;
            b1_ = 1.0;
            break;
        default: 
            // Multistep BDF methods
            // Fill matrix and RHS
            // 0:
            //   sum_{i=0}^{numX-1} a_i = 1
            // j = 1..order
            //   sum_{i=1}^{numX-1} (t_{k-i}/h_k)^j a_i + j (t_{k+1}/h_k)^(j-1) b_{-1} = 1
            // Unknowns: a_0, a_1, ..., a_{numX-1}, b_{-1}
            // First equation
            Int base = 0;
            for(Int i=0; i<numX_; i++) {
                matrix[base+i] = 1.0;
            }
            matrix[base+numX_+1] = 0.0;
            rhs[0] = 1.0;
            // Remaining order_ equations
            base += n_;
            for(Int j=1; j<=order_; j++, base+=n_) {
                // a coeffs
                for(Int i=0; i<numX_; i++) {
                    matrix[base+i] = pow(normalizedTimePoint[i], j);
                }
                // b_{-1} coeff
                matrix[base+numX_] = j;
                // RHS
                rhs[j] = 1.0;
            }

            // Solve
            solve(n_);

            // Unpack
            for(Int i=0; i<numX_; i++) {
                a_[i] = rhs[i]; 
            }
            b1_ = rhs[numX_];
            break;
        }
        break;

    case Method::AdamsBashforth:
        // Handle single step methods without equations (constant coeffs)
        // First equation is always a_0 = 1.0, so we don't need it
        a_[0] = 1.0;
        switch (order_) {
        case 1:
            // Forward Euler
            b_[0] = 1.0;
            break;
        default:
            // Multistep Adams-Bashforth methods
            // Fill matrix and RHS
            // t_{k}=0.0, a_0 term vanishes
            // j = 1..order:
            //   sum_{i=0}^{numXdot-1} j (t_{k-i}/h_k)^(j-1) b_i = 1
            // Unknowns: b_0, b_1, ...
            Int base = 0;
            for(Int j=1; j<=order_; j++, base+=n_) {
                // b coeffs
                for(Int i=0; i<numXdot_; i++) {
                    // Treat b_0 differently for j=1
                    if (i==0 && j==1) {
                        matrix[base+i] = j;
                    } else {
                        matrix[base+i] = j*pow(normalizedTimePoint[i], j-1);
                    }
                }
                // RHS
                rhs[j-1] = 1.0;
            }

            // Solve
            solve(n_);

            // Unpack
            for(Int i=0; i<numXdot_; i++) {
                b_[i] = rhs[i];
            }
            break;
        }
        break;
    
    case Method::PolynomialExtrapolation:
        // Multistep Adams-Bashforth methods
        // Fill matrix and RHS
        // 0:
        //   sum_{i=0}^{numX-1} a_i = 1
        // j = 1..order:
        //   sum_{i=1}^{numX-1} (t_{k-i}/h_k)^j a_i = 1
        // Unknowns: a_0, a_1, ...
        // First equation
        Int base = 0;
        for(Int i=0; i<numX_; i++) {
            matrix[base+i] = 1.0;
        }
        matrix[base+numX_+1] = 0.0;
        rhs[0] = 1.0;
        // Remaining order_ equations
        base += n_;
        for(Int j=1; j<=order_; j++, base+=n_) {
            // a coeffs
            for(Int i=0; i<numX_; i++) {
                matrix[base+i] = pow(normalizedTimePoint[i], j);
            }
            // RHS
            rhs[j] = 1.0;
        }

        // Solve
        solve(n_);

        // Unpack
        for(Int i=0; i<numX_; i++) {
            a_[i] = rhs[i];
        }
        break;
    }
    
    // Compute error coefficient without the j! in denominator
    // For a method of order j-1 LTE is determined by coefficient 
    // 
    //        1         numX-1      t_{k-i} - t_k        numXdot-1  t_{k-i} - t_k
    // C_j = --- ( -1 +  sum   a_i (-------------)^j + j   sum     (-------------)^(j-1)
    //       j!          i=0             h_k               i=-1           h_k
    //
    // LTE at t_{k+1} = C_j x^(j)(t_k) h_k^j
    // 
    // x^(j)(t) is the j-th derivative of x(t) wrt. t
    err_ = -1.0;
    for(Int i=0; i<numX_; i++) {
        err_ += a_[i]*std::pow(normalizedTimePoint[i], order_+1);
    }
    for(Int i=0; i<numXdot_; i++) {
        err_ += (order_+1)*b_[i]*std::pow(normalizedTimePoint[i], order_);
    }
    err_ += (order_+1)*b1_;

    return true;
}

bool IntegratorCoeffs::solve(Int n) {
    // Gaussian elimination with partial pivoting
    Int diagRow = 0;
    double* mxptr = matrix.data();
    double* rhsPtr = rhs.data();
    double* diagRowPtr = mxptr;
    for(Int i=0; i<n; i++, diagRowPtr+=n) {
        // Look for pivot
        double* pivRowPtr = diagRowPtr;
        double pivot = std::abs(pivRowPtr[i]);
        Int pivRowNdx = i; // Pivot row index
        double* atRowPtr = pivRowPtr+n;
        for(Int j=i+1; j<n; j++, atRowPtr+=n) {
            if (std::abs(atRowPtr[i])>pivot) {
                pivot = std::abs(atRowPtr[i]);
                // pivRow = atRow;
                pivRowPtr = atRowPtr;
                pivRowNdx = j;
            }
        }
        // Do we need to swap
        if (diagRowPtr!=pivRowPtr) {
            // Swap RHS
            double tmp;
            tmp = rhsPtr[i];
            rhsPtr[i] = rhsPtr[pivRowNdx];
            rhsPtr[pivRowNdx] = tmp;
            // Swap rows
            for(Int j=i; j<n; j++) {
                tmp = diagRowPtr[j];
                diagRowPtr[j] = pivRowPtr[j];
                pivRowPtr[j] = tmp;
            }
        }
        // Eliminate
        double* subDiagRowPtr = diagRowPtr + n;
        for(Int j=i+1; j<n; j++, subDiagRowPtr+=n) {
            double factor = subDiagRowPtr[i]/diagRowPtr[i];
            // Below diagonal we get 0.0, no need to compute it
            for(Int k=i+1; k<n; k++) {
                subDiagRowPtr[k] -= factor*diagRowPtr[k];
            }
            // Subtract RHS
            rhsPtr[j] -= factor*rhsPtr[i];
            // Store factor (LU decomposition) - no need to do this
            // subDiagRowPtr[i] = factor;
            // For actual LU decomposition we also need to store row permutations
        }
    }
    // Back substitution, start with last row
    double* rowPtr = diagRowPtr - n;
    for(Int i=n-1; i>=0; i--, rowPtr -= n) {
        for(Int j=i+1; j<n; j++) {
            rhsPtr[i] -= rhsPtr[j]*rowPtr[j];
        }
        // Divide
        rhsPtr[i] /= rowPtr[i];
    }
    return true;
}

bool IntegratorCoeffs::scaleDifferentiator(double hk, Status& s) {
    if (hk==0.0) {
        s.set(Status::DivZero, "Division by zero due to zero timestep.");
        return false;
    }
    hk_ = hk;
    aScaled_ = a_;
    bScaled_ = b_;
    
    // Scale only implicit algorithm coeffs
    leading_ = 1 / (hk*b1_);
    for(auto& aIt : aScaled_) {
        aIt *= leading_;
    }
    for(auto& bIt : bScaled_) {
        bIt /= b1_;
    }

    return true;
}

bool IntegratorCoeffs::scalePredictor(double hk, Status& s) {
    if (implicit_) {
        s.set(Status::BadArguments, "Predictor cannot be implicit.");
        return false;
    }
    aScaled_ = a_;
    bScaled_ = b_;
    
    // Scale only implicit algorithm coeffs
    for(auto& bIt : bScaled_) {
        bIt *= hk;
    }

    return true;
}
 
void IntegratorCoeffs::dump(std::ostream& os, bool scaled) {
    switch (method_) {
        case Method::AdamsMoulton:
            os << "AM ";
            break;
        case Method::AdamsBashforth:
            os << "AB ";
            break;
        case Method::BDF:
            os << "BDF ";
            break;
        case Method::PolynomialExtrapolation:
            os << "Poly ";
            break;
    }
    for(Int i=0; i<numX_; i++) {
        if (scaled) {
            os << "a" << i << "/(hk b_{-1})=" << aScaled_[i] << " ";
        } else {
            os << "a" << i << "=" << a_[i] << " ";
        }
    }
    if (implicit_) {
        if (scaled) {
            os << "1/(hk b_{-1})=" << leading_ << " ";
        } else {
            os << "b_{-1}=" << b1_ << " ";
        }
    }
    for(Int i=0; i<numXdot_; i++) {
        if (scaled) {
            os << "b" << i << "/b_{-1}=" << bScaled_[i] << " ";
        } else {
            os << "b" << i << "=" << b_[i] << " ";
        }
    }
}

}
