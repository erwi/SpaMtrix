#ifndef READER_H
#define READER_H
#include <fstream>
namespace SpaMtrix{
class IRCMatrix;

/**
 * @brief File input helpers for sparse matrices.
 */
class Reader{
    static const char* FILE_OPEN_ERROR_STRING;
    static const char* FILE_FORMAT_ERROR_STRING;
    std::fstream file;
public:
    /** @brief Construct a reader helper. */
    Reader(){};
    /** @brief Destroy the reader helper. */
    virtual ~Reader(){};
    /**
     * @brief Read a Matrix Market file into an IRC matrix.
     *
     * @param filename Input file path.
     * @return Parsed sparse matrix.
     */
    static SpaMtrix::IRCMatrix readMatrixMarket(const std::string &filename);
};// end class Reader
}// end namespace SpaMtrix

#endif
