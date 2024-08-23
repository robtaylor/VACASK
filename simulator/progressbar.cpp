#include "progressbar.h"
#include "platform.h"
#include "common.h"


namespace NAMESPACE {

AnalysisProgress::AnalysisProgress(int indent, std::ostream& os, double dt) 
    : indent(indent), os(os), ProgressReporter(dt), lastLen(0) {
    if (Platform::isTty(os)) {
        columns = Platform::ttyColumns(os);
    } else {
        columns = 0;
    }
    barSize = (columns-1)/4 - indent;
    enabled_ = barSize>=10;
}

std::string AnalysisProgress::renderBar(double norm, double val) {
    auto padding = 1;
    auto precision = 0;

    auto padStr = std::string(padding, ' ');
    
    auto rel = val/norm; 
    if (rel<0) {
        rel = 0;
    }
    
    // Format number
    std::stringstream ss;
    ss << std::fixed << std::setprecision(precision) << rel*100 << "%";
    auto num = ss.str();
    int nnum = num.size();

    // Compute components
    int active = barSize - 2;
    int leftSpace = (active - nnum - padding*2)/2;
    int rightBegin = leftSpace + nnum + padding*2;
    int rightSpace = active - rightBegin;

    // At least 5 chars + decimals for number
    // At least 2 visible bar characters
    if (active<2+padding+5+precision || leftSpace<1 || rightSpace<1) {
        // Too small
        return "";
    }

    // Left and right bar
    int nbar = std::round(rel*active);
    int nleft = nbar;
    if (nleft>leftSpace) {
        nleft = leftSpace;
    }
    int nright = nbar - leftSpace - padding*2 - nnum;
    if (nright<0) {
        nright = 0;
    }

    return (
        "|" +
        std::string(nleft, '=') + 
        (leftSpace-nleft>0 ? std::string(leftSpace-nleft, ' ') : "") + 
        padStr + num + padStr + 
        std::string(nright, '=') + 
        (rightSpace-nright>0 ? std::string(rightSpace-nright, ' ') : "") 
        + "|"
    );
}

bool AnalysisProgress::reporter() {
    std::string txt;

    if (!enabled_) {
        return enabled_;
    }

    txt += (indent>0 ? std::string(indent, ' ') : "");
    
    auto barTxt = renderBar(extent_, pos_);
    if (barTxt.size()<=0 || barTxt.size()>columns-1) {
        enabled_ = false;
        return enabled_;
    }
    txt += barTxt;
    
    std::stringstream ss;
    ss << " " << std::fixed << std::setprecision(6) << time() << "s";
    auto timeTxt = ss.str();
    if (txt.size()+timeTxt.size()<columns-1) {
        txt += timeTxt;
        
        if (value_.has_value()) {
            ss.str("");
            ss << " @ ";
            switch (valueFormat_) {
                case ValueFormat::Fixed:
                    ss << std::fixed;
                    break;
                case ValueFormat::Scientific:
                    ss << std::scientific;
                    break;
                case ValueFormat::Default:
                    ss << std::defaultfloat;
                    break;
            }
            ss << std::setprecision(precision_>=0 ? precision_ : 0);
            ss << prefix_ << value_.value() << postfix_;
            auto posStr = ss.str();
            if (txt.size()+posStr.size()<columns+1) {
                txt += posStr;
            }
        }
    }
    
    auto len = txt.size();
    // If new string is shorter, overwrite what is left of the old string
    if (len<lastLen) {
        txt += std::string(lastLen-len, ' ');
    }
    if (enabled_) {
        os << '\r' << txt << std::flush;
    }
    lastStr = std::move(txt);
    lastLen = len;
    return enabled_;
}

}
