#ifndef __COMMON_DEFINED
#define __COMMON_DEFINED

#define NAMESPACE sim

#include <cstdint>
#include <cstdlib>
#include <complex>
#include <type_traits>


namespace NAMESPACE {

// Identifier number type
typedef uint32_t IdentifierIndex;

// String pool block size and block closing free space
const size_t stringPoolBlockSize = 8192; 
const size_t stringPoolGrowthFactor = 2; 
const size_t stringPoolRetries = 4;

// Source line and column number type
typedef uint32_t SourceLineNumber;
typedef uint32_t SourceColumnNumber;

// File stack identifier type
typedef uint32_t FileStackIndex;

// File stack file id type, negative values are reserved for RPN strings
typedef int32_t FileStackFileIndex;


// Location identifier type
typedef uint32_t LocationIndex;


// Value class type used for integers (signed)
typedef int32_t IntegerValue;


// These refer to a single instance
// Index of a parameter
typedef uint32_t ParameterIndex;

// Index of a terminal
typedef uint32_t TerminalIndex;

// Index of a state within the states subvector belonging to a particular instance
typedef uint32_t LocalStateIndex;


// Index of a state within the states vector (holding states of all instances)
typedef size_t StateIndex;

// Index of history entry
typedef uint32_t HistoryDepthIndex;

// Positive (past) and negative (future) history index change
typedef std::make_signed<HistoryDepthIndex>::type HistoryDepthIndexDelta;


// RPN function arity
typedef uint32_t RpnArity;

// RPN jump offset
typedef int32_t RpnJumpOffset;


// Reference count indes
typedef uint32_t RefCountIndex;


// Change these to uint64_t/int64_t to use 64-bit indexing in KLU
// Equation index
typedef uint32_t EquationIndex;

// Unknown and node index
typedef uint32_t UnknownIndex;

// Node index (unknown index)
typedef UnknownIndex NodeIndex;

// Matrix entry index - should be a signed 32- or 64-bit number for KLU
typedef int32_t MatrixEntryIndex;

// Complex numbers used in numerical solvers
typedef std::complex<double> Complex;


// Time relative tolerance
const double timeRelativeTolerance = std::numeric_limits<double>::epsilon()*8;

}

#endif
