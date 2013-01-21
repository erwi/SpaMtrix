#ifndef READER_H
#define READER_H
#include <fstream>
namespace SpaMtrix{
// Forward declarations
class IRCMatrix;

class Reader{
/*! This is a class containing functionality for reading matrices and/or
* vectors from files located on disk.
*/
    static const char* FILE_OPEN_ERROR_STRING;
    static const char* FILE_FORMAT_ERROR_STRING;
    std::fstream file;
public:
    Reader(){};
    ~Reader(){};
    static SpaMtrix::IRCMatrix readMatrixMarket(const std::string &filename);
};// end class Reader
}// end namespace SpaMtrix

#endif