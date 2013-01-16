#ifndef READER_H
#define READER_H
#include <fstream>

class FlexiMatrix;

namespace SpaMtrix{
    class Reader{
        /*! This is a class containing functionality for reading matrices and/or
         * vectors from files located on disk.
         */
        std::fstream file;
    public:
        Reader(){};
        ~Reader(){};
        bool readMatrixMarket(const std::string &filename,      
                              FlexiMatrix &mat);
    };
}// end namespace SpaMtrix

#endif