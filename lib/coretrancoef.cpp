#include "coretrancoef.h"
#include "common.h"

namespace NAMESPACE {

IntegratorCoeffs::IntegratorCoeffs(Method method, Int order) 
    : method_(method), order_(order) {
    setOrder(order);
    b1_ = 0.0;
    err_ = 0.0;

    /*
    // Test solver
    matrix = DenseMatrix<double>({1, 2, 3, 2, 3, -5, -6, -8, 1}, 3, 3);
    rhs = {-7, 9, 22};
    solve(3);
    for(auto it : rhs) {
        std::cout << it*25 << " ";
    }
    std::cout << "\n";
    // Result should be: -393 217 72
    */ 
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

bool IntegratorCoeffs::compute(CircularBuffer<double>& pastSteps, double newStep) {
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
    matrix.resize(n_, n_);
    rhs.resize(n_);

    // Prepare space for coeffs
    a_.resize(numX_);
    b_.resize(numXdot_);
    
    // Coeffs are 0.0 by default
    a_.assign(numX_, 0.0);
    b1_ = 0.0;
    b_.assign(numXdot_, 0.0);

    // All matrix and RHS entries will be set, so there is no need to set them to 0
    DBGCHECK(pastSteps.valueCount()+1<numX_ || pastSteps.valueCount()+1<numXdot_, "Timestep history is too short.");
    
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
            for(Int j=1; j<=order_; j++) {
                auto row = matrix.row(j-1);
                // Manually add b_{-1}
                row[0] = j;
                // b coeffs
                for(Int i=0; i<numXdot_; i++) {
                    // Treat b_0 differently for j=1
                    if (i==0 && j==1) {
                        row[1+i] = j;
                    } else {
                        row[1+i] = j*pow(normalizedTimePoint[i], j-1);
                    }
                }
                rhs[j-1] = 1.0;
            }
            
            // Solve 
            if (!solve(n_)) {
                return false;
            }

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
            auto row0 = matrix.row(0);
            for(Int i=0; i<numX_; i++) {
                row0[i] = 1.0;
            }
            row0[numX_] = 0.0;
            rhs[0] = 1.0;
            // Remaining order_ equations
            for(Int j=1; j<=order_; j++) {
                auto row = matrix.row(j);
                // a coeffs
                row[0] = 0;
                for(Int i=1; i<numX_; i++) {
                    row[i] = pow(normalizedTimePoint[i], j);
                }
                // b_{-1} coeff
                row[numX_] = j;
                // RHS
                rhs[j] = 1.0;
            }

            // Solve
            if (!solve(n_)) {
                return false;
            }
            
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
            for(Int j=1; j<=order_; j++) {
                auto row = matrix.row(j-1);
                // b coeffs
                for(Int i=0; i<numXdot_; i++) {
                    // Treat b_0 differently for j=1
                    if (i==0 && j==1) {
                        row[i] = j;
                    } else {
                        row[i] = j*pow(normalizedTimePoint[i], j-1);
                    }
                }
                // RHS
                rhs[j-1] = 1.0;
            }

            // Solve
            if (!solve(n_)) {
                return false;
            }
            
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
        auto row0 = matrix.row(0);
        for(Int i=0; i<numX_; i++) {
            row0[i] = 1.0;
        }
        rhs[0] = 1.0;
        // Remaining order_ equations
        for(Int j=1; j<=order_; j++) {
            auto row = matrix.row(j);
            // a coeffs
            for(Int i=0; i<numX_; i++) {
                row[i] = pow(normalizedTimePoint[i], j);
            }
            // RHS
            rhs[j] = 1.0;
        }

        // Solve
        if (!solve(n_)) {
            return false;
        }

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
    // C_j = --- ( -1 +  sum   a_i (-------------)^j + j   sum     (-------------)^(j-1) )
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
    // auto dmv = DenseMatrixView(matrix.data(), n, n, n, 1);
    auto vv = VectorView(rhs.data(), n, 1);
    return matrix.destructiveSolve(vv);
}

bool IntegratorCoeffs::scaleDifferentiator(double hk) {
    if (hk==0.0) {
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

bool IntegratorCoeffs::scalePredictor(double hk) {
    if (implicit_) {
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

bool IntegratorCoeffs::test() {
    // Object
    IntegratorCoeffs ic;

    // Test outcome
    bool ok = true;

    // Expected values are for order=3
    int order = 3;

    // Past steps
    CircularBuffer<double> pastSteps(order);
    for(int i=0; i<order; i++) {
        pastSteps.add(0.1);
    }
    
    // AM3, uniform step
    ic.setMethod(Method::AdamsMoulton, order);
    if (ic.compute(pastSteps, 0.1)) {
        ic.dump(std::cout);
        auto errExpect = (1.0/24)*ffactorial(order+1);
        std::cout << " C=" << ic.err_ << "\n";
        std::cout << "Expected: " << 1.0 << " " << (5.0/12) << " " << (8.0/12) << " " << (-1.0/12) 
                << " " << errExpect << "\n";
        std::cout << "\n";
        if (std::abs(ic.err_-errExpect)>1e-12) {
            ok = false;
        }
    } else {
        std::cout << "AM failed\n";
        ok = false;
    }

    // BDF3, uniform step
    ic.setMethod(Method::BDF, order);
    if (ic.compute(pastSteps, 0.1)) {
        ic.dump(std::cout);
        auto errExpect = (3.0/22)*ffactorial(order+1);
        std::cout << " C=" << ic.err_ << "\n";
        std::cout << "Expected: " << " " << (18.0/11) << " " << (-9.0/11) << " " << (2.0/11) << " " << (6.0/11) 
                << " " << errExpect << "\n";
        std::cout << "\n";
        if (std::abs(ic.err_-errExpect)>1e-12) {
            ok = false;
        }
    } else {
        std::cout << "BDF failed\n";
        ok = false;
    }
    
    // AB3, uniform step
    ic.setMethod(Method::AdamsBashforth, order);
    if (ic.compute(pastSteps, 0.1)) {
        ic.dump(std::cout);
        auto errExpect = (-3.0/8)*ffactorial(order+1);
        std::cout << " C=" << ic.err_ << "\n";
        std::cout << "Expected: " << " " << (1) << " " << (23.0/12) << " " << (-16.0/12) << " " << (5.0/12) 
                << " " << errExpect << "\n";
        std::cout << "\n";
        if (std::abs(ic.err_-errExpect)>1e-12) {
            ok = false;
        }
    } else {
        std::cout << "AB failed\n";
        ok = false;
    }

    // Polynomial extrapolation
    ic.setMethod(Method::PolynomialExtrapolation, order);
    if (ic.compute(pastSteps, 0.1)) {
        ic.dump(std::cout);
        auto errExpect = (-1.0)*ffactorial(order+1);
        std::cout << " C=" << ic.err_ << "\n";
        std::cout << "Expected: " << " " << (4) << " " << (-6) << " " << (4) << " " << (-1) 
                << " " << errExpect << "\n";
        std::cout << "\n";
        if (std::abs(ic.err_-errExpect)>1e-12) {
            ok = false;
        }
    } else {
        std::cout << "Polynomial extrapolation failed\n";
        ok = false;
    }
    std::cout << "Integrator coeffs test " << (ok ? "OK" : "FAILED") << "\n";
    return ok;
}

}
