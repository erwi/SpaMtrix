#ifndef WRITER_H
#define WRITER_H
#include <spamtrix_setup.hpp>
namespace SpaMtrix{
class Vector;
class IRCMatrix;

/**
 * @brief File output helpers for vectors and sparse matrices.
 */
class Writer{
    Writer() {}
public:
    /** @brief Destroy the writer helper. */
    virtual ~Writer();
    /**
     * @brief Write a vector to CSV.
     *
     * @param filename Output file path.
     * @param data Vector to write.
     * @param numC Number of columns per row in the CSV output.
     * @return True on success.
     */
    static bool writeCSV(const std::string &filename,
                  const Vector &data,
                  idx numC=0);
    /**
     * @brief Write a sparse matrix in Matrix Market format.
     *
     * @param filename Output file path.
     * @param A Matrix to write.
     * @return True on success.
     */
    static bool writeMatrixMarket(const std::string &filename,
                                  const IRCMatrix &A);
};// end class Writer
} // end namespace SpaMtrix
#endif
