#include "outrawfile.h"
#include <ctime>
#include <iomanip>
#include <filesystem>
#include "common.h"


namespace NAMESPACE {

// std::endl is slow because it flushes the stream, use \n

OutputRawfile::OutputRawfile(const std::string& baseName, Output::DescriptorList& descriptors, Output::SourcesList& sources, Flags f) 
    : Output(baseName, descriptors, sources), FlagBase(f), fileName(baseName+".raw") {
    // Delete old file before opening
    if (std::filesystem::exists(fileName)) {
        std::filesystem::remove(fileName);
    }
}

OutputRawfile::~OutputRawfile() {
    outStream.close();
}

bool OutputRawfile::prologue(Status& s) {
    outStream.open(fileName, std::ios::binary | std::ios::out | std::ios::trunc);
    outStream << "Title: " << title_ << "\n";
    if (outStream.fail()) {
        s.set(Status::CreationFailed, "Failed to write file '"+fileName+"'.");
        return false;
    }
    
    auto t = std::time(nullptr);
    outStream << "Date: " << std::asctime(std::localtime(&t)); // asctime() adds a newline

    outStream << "Plotname: " << plotname_ << "\n";
    outStream << "Flags: " << (checkFlags(Flags::Complex) ? "complex" : "real" ) 
              << (checkFlags(Flags::Padded) ? "" : " unpadded") << "\n";
    outStream << "No. Variables: " << std::to_string(srcs.size()) << "\n";
    
    outStream << "No. Points: "; 
    pointCountPos = outStream.tellp();
    outStream << "                \n"; // 16 spaces for point count

    // Add "Dimensions: n,n,n,...\n" for plots with multidimensional vectors

    outStream << "Variables:\n";
    for(size_t i=0; i<descrs.size(); i++) {
        outStream << "\t" << std::to_string(i) << "\t" << std::string(descrs[i].name) << "\t" << "notype"; 
        // TODO: add " dims=n,n,n,...\n" for rawfiles with vectors of different length
        outStream << "\n"; 
    }

    if (checkFlags(Flags::Binary)) {
        outStream << "Binary:\n";
    } else {
        outStream << "Values:\n";
        outStream << std::scientific << std::setprecision(15);
    }

    count=0;

    return true;
}

bool OutputRawfile::addPoint(Status& s) {
    if (checkFlags(Flags::Binary)) {
        if (checkFlags(Flags::Complex)) {
            for(size_t i=0; i<descrs.size(); i++) {
                Complex c = srcs[i].getC();
                outStream.write(reinterpret_cast<char*>(&c), sizeof(Complex));
            }
        } else {
            for(size_t i=0; i<descrs.size(); i++) {
                double r = srcs[i].getR();
                outStream.write(reinterpret_cast<char*>(&r), sizeof(double));
            }
        }
    } else {
        if (checkFlags(Flags::Complex)) {
            for(size_t i=0; i<descrs.size(); i++) {
                if (i==0) {
                    outStream << " " << count;
                }
                Complex c = srcs[i].getC();
                outStream << "\t" << c.real() << "," << c.imag() << "\n";
            }
        } else {
            for(size_t i=0; i<descrs.size(); i++) {
                if (i==0) {
                    outStream << " " << count;
                }
                double r = srcs[i].getR();
                outStream << "\t" << r << "\n";
            }
        }
        outStream << "\n";
    }
    count++;

    return true;
}

bool OutputRawfile::epilogue(Status& s) {
    outStream.seekp(pointCountPos);
    outStream << count;
    outStream.seekp(outStream.end);
    outStream.close(); 

    return true;
}

bool OutputRawfile::remove(Status& s) {
    if (std::filesystem::exists(fileName)) {
        std::filesystem::remove(fileName);
    }

    return true;
}

}
