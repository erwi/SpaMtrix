#ifndef WRITER_H
#define WRITER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// SpaMtrix INCLUDES
#include <spamtrix_setup.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_ircmatrix.hpp>

namespace SpaMtrix{
class Writer{
    std::fstream file;
    public:
        Writer();
        ~Writer();
        
        bool writeCSV(const std::string &filename,
                      const Vector &data,
                      idx numC=0);

        static bool writeMatrixMarket(const std::string &filename,
                                      const IRCMatrix &A);

    };
} // end namespace SpaMtrix
#endif
