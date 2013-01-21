#ifndef WRITER_H
#define WRITER_H
// SpaMtrix INCLUDES
#include <spamtrix_setup.hpp>
namespace SpaMtrix{
// forward declarations
class Vector;
class IRCMatrix;

class Writer{
    Writer() {}
public:
    ~Writer();
    static bool writeCSV(const std::string &filename,
                  const Vector &data,
                  idx numC=0);
    static bool writeMatrixMarket(const std::string &filename,
                                  const IRCMatrix &A);
};// end class Writer
} // end namespace SpaMtrix
#endif
