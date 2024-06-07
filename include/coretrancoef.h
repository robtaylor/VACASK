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

    IntegratorCoeffs(Method method=Method::AdamsMoulton, Int order=1); 

    // Method
    Method method() const { return method_; };
    
    // Order
    Int order() const { return order_; };

    // xmu for trapezoidal algorithm
    double xmu() const { return xmu_; }; 

    // Change method
    bool setMethod(Method method, Int order, double xmu=0.5);

    // Change order
    bool setOrder(Int order);

    // Set xmu for trapezoidal algorithm
    bool setXmu(double xmu=0.5) { xmu_ = xmu; return xmu>=0 && xmu<=0.5; };

    // Compute coefficients
    bool compute(CircularBuffer<double>& pastSteps, double newStep);
    
    // Coefficients for past values
    const std::vector<double>& a() const { return a_; };

    // Coefficients for past derivatives
    const std::vector<double>& b() const { return b_; };

    // Coefficient for new derivative
    double b1() const { return b1_; };

    // Compute scaled coefficients
    bool scaleDifferentiator(double hk);
    bool scalePredictor(double hk);

    // Scaled coefficients
    const std::vector<double> aScaled() const { return aScaled_; };
    const std::vector<double> bScaled() const { return bScaled_; };
    double leadingCoeff() const { return leading_; }; 

    // Minimal number of past points needed by the predictor
    size_t minimalPredictorHistory() {
        return numX_;
    };

    // Minimal number of past points needed by the differentiator
    size_t minimalDifferentiatorHistory() {
        return std::max(numX_, numXdot_);
    };

    // Prepare fast array pointers on which predict() will operate
    bool preparePredictorHistory(VectorRepository<double>& repo, Int historyOffset=1) {
        DBGCHECK(numXdot_>0, "Predictors using past derivatives are not supported.");
        predictorHistory.clear();
        auto n = minimalPredictorHistory();
        for(decltype(n) i=0; i<n; i++) {
            predictorHistory.push_back(repo.data(historyOffset+i));
        }
        return true;
    };

    // Prepare fast array pointers on which differentiate() will operate
    // Also prepare fast pointer to future states
    bool prepareDifferentiatorHistory(VectorRepository<double>& repo, Int historyOffset=1) {
        DBGCHECK(!implicit_, "Explicit algorithms cannot be used for computing future derivative.");
        differentiatorHistory.clear();
        auto n = minimalDifferentiatorHistory();
        for(decltype(n) i=0; i<n; i++) {
            differentiatorHistory.push_back(repo.data(historyOffset+i));
        }
        return true;
    };

    // Differentiate state at tk+hk based on 
    // - future value 
    // - value history (state) and
    // - derivative history (state+1)
    // To be used with implicit integration algorithms
    double differentiate(double futureValue, StateIndex state) {
        // Contribution of future value
        double deriv = leading_ * futureValue;
        // Contribution of past values
        for(Int i=0; i<aScaled_.size(); i++) {
            deriv -= aScaled_[i] * differentiatorHistory[i][state];
        }
        // Contribution of past derivatives
        for(Int i=0; i<bScaled_.size(); i++) {
            deriv -= bScaled_[i] * differentiatorHistory[i][state+1];
        }
        return deriv;
    };

    // Predict value based on value history 
    // To be used with explicit algorithms (predictors)
    // No need to zero prediction before this function is called
    void predict(Vector<double>& prediction) {
        auto n = prediction.size();
        for(decltype(n) i=0; i<n; i++) {
            double pred = 0;
            for(Int j=0; j<a_.size(); j++) {
                pred += aScaled_[j] * predictorHistory[j][i];
            }
            prediction[i] = pred;
        }
    };

    // Error coefficient multiplied by (order+1)!
    // Note that all return values remain reasonable, 
    // except for polynomial extrapolation where they grow with (order+1)!
    double errorCoeff() const { return err_; };

    // Required history length
    Int pastStatesNeeded() const { return std::max(numX_, numXdot_); }

    void dump(std::ostream& os, bool scaled=false);

    static double ffactorial(int n) {
        static std::vector<double> cache = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320 };
        if (n>=cache.size()) {
            for(int i=cache.size(); i<=n; i++) {
                cache.push_back(cache[i-1]*i);
            }
        }
        return cache[n];
    };
    
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

    // History - fast pointers
    std::vector<double*> predictorHistory;
    std::vector<double*> differentiatorHistory;
};

}

#endif
