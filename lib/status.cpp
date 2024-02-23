#include <fstream>
#include <sstream>
#include "status.h"
#include "filesystem.h"
#include "common.h"


namespace NAMESPACE {

Status Status::ignore = Status(true);

void Status::set(Code c, const std::string& msg) { 
    if (!ignoreFlag) {
        if (message_.size()>0) {
            message_ += "\n";
        }
        message_ += msg; 
        if (c!=OK) {
            error = true;
        }
    }
}

void Status::set(const Status& other) {
    if (!ignoreFlag) {
        // Anything to append
        if (other.message_.size()>0) {
            // Need newline (already have some message)
            if (message_.size()>0) {
                message_ += "\n";
            }
            message_ += other.message_; 
        }
        if (other.error) {
            error = true;
        }
    }
}

void Status::extend(const std::string& msg) { 
    set(OK, msg);
} 


void Status::extend(const Loc& point) {
    if (point) {
        extend(point.toString());
    }
}

void Status::prefix(const std::string& msg) {
    if (!ignoreFlag) {
        if (msg.size()>0) {
            message_ = msg + "\n" + message_;
        } else {
            message_ = msg;
        }
    }
}

const std::string& Status::message() const {
    return message_;
}

void Status::clear() {
    if (!ignoreFlag) {
        message_ = "";
        error = false;
    }
}

Status& Status::operator=(Status&& other) {
    if (!ignoreFlag) {
        message_ = std::move(other.message_);
        error = other.error;
    } else {
        auto dummy1 = std::move(other.message_);
    }
    return *this;
};


}
