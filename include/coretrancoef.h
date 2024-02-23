#ifndef __CORETRANCOEF_DEFINED
#define __CORETRANCOEF_DEFINED

#include "value.h"
#include "ansupport.h"
#include "common.h"
#include <algorithm>


namespace NAMESPACE {

class IntegratorCoeffs {
public:
    enum class Method { AdamsMoulton, BDF, AdamsBashforth, PolynomialExtrapolation };

    IntegratorCoeffs(Method method=Method::AdamsMoulton, Int order=1, Int historyOffset=0); 

    // Method
    Method method() const { return method_; };
    
    // Order
    Int order() const { return order_; };

    // xmu for trapezoidal algorithm
    double xmu() const { return xmu_; }; 

    // Change history offset
    void setHistoryOffset(Int offs) { historyOffset_ = offs; };
    
    // Change method
    bool setMethod(Method method, Int order, double xmu=0.5);

    // Change order
    bool setOrder(Int order);

    // Set xmu for trapezoidal algorithm
    bool setXmu(double xmu=0.5) { xmu_ = xmu; return xmu>=0 && xmu<=0.5; };

    // Compute coefficients
    bool compute(CircularBuffer<double>& pastSteps, double newStep, Status& s=Status::ignore);
    
    // Coefficients for past values
    const std::vector<double>& a() const { return a_; };

    // Coefficients for past derivatives
    const std::vector<double>& b() const { return b_; };

    // Coefficient for new derivative
    double b1() const { return b1_; };

    // Compute scaled coefficients
    bool scaleDifferentiator(double hk, Status& s=Status::ignore);
    bool scalePredictor(double hk, Status& s=Status::ignore);

    // Scaled coefficients
    const std::vector<double> aScaled() const { return aScaled_; };
    const std::vector<double> bScaled() const { return bScaled_; };
    double leadingCoeff() const { return leading_; }; 

    // Differentiate state at tk+hk based on 
    // - value history including future value and
    // - derivative history
    // To be used with implicit integration algorithms
    double differentiate(VectorRepository<double>& states, double futureValue, StateIndex state);

    // Predict value based on value history 
    // To be used with explicit algorithms (predictors)
    // Assumes prediction is a zero vector
    void predict(VectorRepository<double>& history, Vector<double>& prediction) {
        if (numXdot_>0) {
            // Does not handle predictors that use past derivative values
            throw std::length_error("Predictors using past derivatives are not supported.");
        }
        auto n = history.length();
        for(decltype(n) i=0; i<n; i++) {
            for(Int j=0; j<a_.size(); j++) {
                prediction[i] += aScaled_[j] * history.at(historyOffset_+j)[i];
            }
        }
    };

    // Minimal number of past points needed by the predictor
    size_t minimalHistoryForPredictor() {
        return std::max(numX_, numXdot_);
    };

    // Error coefficient multiplied by (order+1)!
    // Note that all return values remain reasonable, 
    // except for polynomial extrapolation where they grow with (order+1)!
    double errorCoeff() const { return err_; };

    // Required history length
    Int pastStatesNeeded() const { return std::max(numX_, numXdot_); }

    void dump(std::ostream& os, bool scaled=false);
    
private:
    bool size();
    bool solve(Int n);

    Method method_;
    Int order_;
    double xmu_;
    Int numX_;
    Int numXdot_;
    bool implicit_;
    bool multistep_;
    Int historyOffset_;

    // Number of equations, matrix (ordered by rows), and RHS
    Int n_;
    std::vector<double> matrix; // row1, row2, ...
    std::vector<double> rhs;
    
    // New timepoint: 
    //   t_{k+1} = h_k 
    // Old timepoints:
    //   0: t_k = 0.0
    //   1: t_{k-i} = h_{k-1} + h_{k-2} + ... + h_{k-i}
    //   ...

    // Old timepoints, normalized by h_k
    //   t_k/h_k, t_{k-1}/h_k, t_{k-2}/h_k, ...
    std::vector<double> normalizedTimePoint; 

    // Computed coefficients
    std::vector<double> a_; // Coeffs for x(t_{k-i}), i>=0
    std::vector<double> b_; // Coeffs for xdot(t_{k-i}), i>=0
    double b1_;             // Coeff for xdot(t_{k+1})
    double err_;            // Error coefficient multiplied by (order+1)!
    std::vector<double> aScaled_; // a_ / (h_k b_{-1})
    std::vector<double> bScaled_; // b_ / b_{-1}
    double hk_;
    double leading_;
    CircularBuffer<double>* pastSteps_;
};

}

#endif
