#ifndef __DFLLOCATION_DEFINED
#define __DFLLOCATION_DEFINED

#include <iostream>
#include <string>
#include <vector>
#include <variant>
#include <sstream>
#include <type_traits>
#include "filestack.h"
#include "sourceloc.h"
#include "common.h"


namespace NAMESPACE {

// Position and Location, Bison style
// We just extend it with FileStack
// Note that locations are stored in parser tables with a different type (Loc)
// which takes less storage. 
// We translate Bison Location to Loc as tables are generated. 

// A point in a source file.
class Position {
public:
    // Type for file name.
    typedef const std::string FilenameType;
    // Type for line and column numbers.
    typedef std::make_signed<size_t>::type CounterType;

    explicit Position(FilenameType* f=nullptr, CounterType l=1, CounterType c=1, CounterType offset=0)
        : filename (f), line (l), column (c), offset(offset), fileStack(nullptr), fileStackPosition(0) {};

    // Set position to given file, line, column
    void initialize(FilenameType* f=nullptr, CounterType l=1, CounterType c=1, CounterType o=0) {
        filename = f;
        line = l;
        column = c;
        offset = o;
        fileStack = nullptr;
        fileStackPosition = 0;
    };

    // Change filename
    void setFilename(FilenameType* f) { filename = f; };

    // Change filestack information
    void setFileStack(FileStack& fstack, FileStackFileIndex fstackPos) { fileStack = &fstack; fileStackPosition = fstackPos; }

    // Advance by count lines, reset columns to 1
    void lines(CounterType count=1) {
        // If count>0 increase line by count, lower bound for result is 1
        if (count) {
            column = 1;
            line = add_(line, count, 1);
        }
    }

    // Advance column by count, lower bound for result is 1
    void columns (CounterType count = 1) {
        column = add_ (column, count, 1);
        offset += count;
    }

    // Convert to simulator location
    Loc loc() const {
        return Loc(fileStack->id(), fileStackPosition, line, offset);
    };

    // File name to which this position refers.
    FilenameType* filename;
    // Current line number.
    CounterType line;
    // Current column number.
    CounterType column;
    // Offset within file
    CounterType offset;
    // File stack
    FileStack *fileStack;
    // File stack position of file 
    FileStackFileIndex fileStackPosition;

private:
    /// Compute max(min, lhs+rhs).
    static CounterType add_(CounterType lhs, CounterType rhs, CounterType min) {
        return lhs + rhs < min ? min : lhs + rhs;
    }
};

/// Add width columns, in place.
inline Position& operator+=(Position& res, Position::CounterType width) {
    res.columns(width);
    return res;
}

/// Add width columns.
inline Position operator+(Position res, Position::CounterType width) {
    return res += width;
}

/// Subtract width columns, in place.
inline Position& operator-=(Position& res, Position::CounterType width) {
    return res += -width;
}

/// Subtract width columns.
inline Position operator-(Position res, Position::CounterType width) {
    return res -= width;
}

// Output positionb to a stream
template <typename YYChar> std::basic_ostream<YYChar>& operator<<(
    std::basic_ostream<YYChar>& ostr, const Position& pos
) {
    if (pos.filename)
        ostr << *pos.filename << ':';
    return ostr << pos.line << '.' << pos.column;
}

// Range in the input defined by two Position objects
class Location {
public:
    // Type for file name.
    typedef Position::FilenameType FilenameType;
    // Type for line and column numbers.
    typedef Position::CounterType CounterType;

    // Construct a location from two Positions
    Location(const Position& b, const Position& e) : begin(b), end(e) {};

    // Construct a 0-width location at p.
    explicit Location(const Position& p=Position()) : begin(p), end(p) {};

    // Construct a 0-width location in f at line l, column c
    explicit Location (FilenameType* f, CounterType l=1, CounterType c=1)
        : begin(f, l, c), end(f, l, c) {};

    // Initialization.
    void initialize (FilenameType* f=nullptr, CounterType l=1, CounterType c=1) {
        begin.initialize(f, l, c);
        end = begin;
    }

    // Reset initial location to final location.
    void step () { begin = end; };

    // Extend the current location to the count next columns.
    void columns(CounterType count=1) { end += count; };  

    // Extend the current location to the COUNT next lines.
    void lines(CounterType count=1) { end.lines (count); }

    // Convert to simulator location
    Loc loc() const {
        return Loc(begin.fileStack->id(), begin.fileStackPosition, begin.line, begin.offset);
    };
    
    // Beginning of the located region.
    Position begin;
    // End of the located region.
    Position end;
};

// Join two locations, in place.
inline Location& operator+=(Location& res, const Location& end) {
    res.end = end.end;
    return res;
};

// Join two locations.
inline Location operator+(Location res, const Location& end) { return res += end; };

// Add width columns to the end position, in place.
inline Location& operator+=(Location& res, Location::CounterType width) { res.columns(width); return res; };

// Add width columns to the end position.
inline Location operator+(Location res, Location::CounterType width) { return res += width; };

// Subtract width columns from the end position, in place.
inline Location& operator-=(Location& res, Location::CounterType width) { return res += -width; };

// Subtract width columns from the end position.
inline Location operator-(Location res, Location::CounterType width) { return res -= width; };

// Output location to a stream
template <typename YYChar> std::basic_ostream<YYChar>&
    operator<<(std::basic_ostream<YYChar>& ostr, const Location& loc) {
    Location::CounterType end_col = loc.end.column>0 ? loc.end.column-1 : 0;
    
    // Output beginning location
    ostr << loc.begin;
    
    if (
        loc.end.filename && 
        (!loc.begin.filename || *loc.begin.filename != *loc.end.filename)
    ) {
        // End filename is given and 
        // is different from begin filename or there is no begin filename
        // Print range with filename, line, and column
        ostr << '-' << loc.end.filename << ':' << loc.end.line << '.' << end_col;
    } else if (loc.begin.line < loc.end.line) {
        // Line of the beginning is smaller than line of the end
        // Print range with line and column
        ostr << '-' << loc.end.line << '.' << end_col;
    } else if (loc.begin.column < end_col) {
        // Print range with column only
        ostr << '-' << end_col;
    }
    return ostr;
}

}

#endif
