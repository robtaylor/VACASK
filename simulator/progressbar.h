#ifndef __PROGRESSBAR_DEFINED
#define __PROGRESSBAR_DEFINED

#include <string>
#include <iomanip>
#include "common.h"

namespace NAMESPACE { 

class ProgressBar {
public:
    ProgressBar(int width=20, int prec=0, double norm=100) : width(width), prec(prec), norm(norm) {};

    std::string render(double val) {
        auto rel = val/norm;
        if (rel<0) {
            rel = 0;
        }
        int nchar = std::round(rel*(width-2));
        
        std::stringstream ss;
        ss << std::fixed << std::setprecision(prec);
        ss << rel*100 << "%";

        auto num = ss.str();
        int nnum = num.size();

        int leftSpace = (width - 2 - nnum)/2;
        int rightBegin = leftSpace + nnum;
        int rightSpace = width - 2 - leftSpace - nnum;

        int nleft = nchar;
        if (nleft>leftSpace-1) {
            nleft = leftSpace-1;
        }
        int nright = nchar - leftSpace - 1 - nnum;
        if (nright<0) {
            nright = 0;
        }

        return (
            "|" +
            std::string(nleft, '=') + 
            (leftSpace-nleft>0 ? std::string(leftSpace-1-nleft, ' ') : "") + 
            " " + num + " " + 
            std::string(nright, '=') + 
            (rightSpace-nright>0 ? std::string(rightSpace-1-nright, ' ') : "") 
            + "|"
        );
    }

private:
    int width;
    int prec;
    double norm;
};


}

#endif
