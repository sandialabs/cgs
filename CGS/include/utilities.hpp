//*****************************************************************//
//    Classifier-Guided Sampling (CGS) 1.0                         //
//    Copyright 2015 Sandia Corporation                            //
//    This software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level CGS directory     //
//*****************************************************************//

// FILE CONTENTS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
///  \file
///  \brief Contains the declarations of general utility namespace functions. 
///         Definitions of template functions are at the end of this file.
////////////////////////////////////////////////////////////////////////////////


// INCLUDE GUARD
// #############################################################################
#ifndef CGS_INCLUDE_UTILITIES_HPP
#define CGS_INCLUDE_UTILITIES_HPP


// INCLUDES
// #############################################################################
#include <fstream>
#include <vector>
#include <xhash>


// #############################################################################
// BEGIN NAMESPACE
// #############################################################################
namespace utilities {

    // ALGORITHMS NAMESPACE
    // =====================================
    namespace algorithms {

        std::vector< std::vector<int>  >
        generatePermutationsWithRepetition(const std::vector<int>& domainSizes);

        double
        factorial(size_t n);

        double
        divideFactorials(size_t numArg, size_t denomArg);

        double
        logFactorial(size_t n);

        unsigned int
        generateRandomSeed();

    }   // namespace algorithms


    // HASH NAMESPACE
    // =====================================
    namespace hash {
        
        template <typename T>
        size_t 
        hash_combine(const size_t& seed, const T& toHash);

        template <typename FWIT>
        size_t 
        hash_range(FWIT start, const FWIT& end);

    }   // namespace hash


    // READWRITE NAMESPACE
    // =================================
    namespace readwrite {

        template <class T>
        std::ifstream&
        readCsvRow(std::ifstream& ifs, std::vector<T>& row);

        template <class T>
        void
        writeVecToCsv(std::ofstream& ofs, const std::vector<T>& vecToWrite);

    }   // namespace readwrite


// TEMPLATE DEFINITIONS
// #############################################################################
////////////////////////////////////////////////////////////////////////////////
/// \brief  Combines two hash values to create a new hash value.
///
/// \param  seed    Beginning hash value.
/// \param  toHash  Variable to combine with existing hash value.
///
/// \return     The new hash value.
////////////////////////////////////////////////////////////////////////////////
template <typename T>
size_t 
hash::hash_combine(const size_t& seed, const T& toHash)
{
    return seed ^ (
        stdext::hash_value(toHash) + 0x9e3779b9 + (seed<<6) + (seed>>2)
        );

}   // hash_combine


////////////////////////////////////////////////////////////////////////////////
/// \brief  Calculates the combined hash value of the elements in an iterator
///         range.
///
/// \param  start   Tterator to start of the range to combine.
/// \param  end     Iterator to the end of the range to combine.
///
/// \return     The combined hash value.
////////////////////////////////////////////////////////////////////////////////
template <typename FWIT>
size_t 
hash::hash_range(FWIT start, const FWIT& end)
{
    size_t ret = 0;
    for(; start!=end; ++start) hash_combine(ret, *start);
    return ret;

}   // hash_range


////////////////////////////////////////////////////////////////////////////////
/// \brief  Reads a row of elements from a .csv file into a std::vector. Returns
///         a reference to the ifstream, which will be in an error state if the
///         end of the file was reached or it encountered something that could
///         not be read into the vector (e.g. inconsistent types).
///
/// \param  ifs     The ifstream that reads from the file.
/// \param  row     A vector into which the row of the file will be read.
///
/// \return     Reference to the input parameter ifstream.
////////////////////////////////////////////////////////////////////////////////
template <class T>
std::ifstream&
readwrite::readCsvRow(std::ifstream& ifs, std::vector<T>& row)
{
    // check if ifs is in an error state
    if (ifs) {
        // clear contents of row
        row.clear();
        T x;
        // attempt to read from the file until the end of the line is reached
        while (ifs >> x) {
            row.push_back(x);

            // ignore the commas
            while (char c = ifs.peek() == ',') {
                ifs.ignore(1, ',');
            }
            
            // if end of line reached, break
            if (char c = ifs.peek() == '\n') {
                break;
            }
        }
    }
    // ifs will be in an error state if eof was reached or it encountered
    // something that could not be read into the vector

    return ifs;

}   // readCsvRow


////////////////////////////////////////////////////////////////////////////////
/// \brief  Writes each element of a vector to an ofstream object, separated
///         by commas so that it can be saved in .csv format. A new row is 
///         created in the ofstream at the end of the vector.
///
/// \param  ofs     The output file stream that the vector is being written to.
/// \param  v       The vector to write to the output file.
////////////////////////////////////////////////////////////////////////////////
template <class T>
void
readwrite::writeVecToCsv(std::ofstream& ofs, const std::vector<T>& v)
{
    
    // write each element of v to the ofstream separated by commas
    for (auto i = v.begin(); i != v.end(); ++i) {

        ofs << *i << ',';
    }
    
    // return to new row
    ofs << '\n';

}   // writeToCsv


// #############################################################################
// END NAMESPACE
// #############################################################################
}   // namespace utilities


#endif  // CGS_INCLUDE_UTILITIES_HPP