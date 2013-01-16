#include <reader.h>
#include <iostream>
namespace SpaMtrix{
bool Reader::readMatrixMarket(const std::string &filename,
                              FlexiMatrix &matrix){
    /*!Reads sparse matrix from a mattrix market formated (text) file saving it
     * in 'matrix'
     */
    
    std::ifstream file;
    file.open(filename);
    
    if ( !file.is_open() ) 
        return false;
    std::string line;
    while (file.good() ){
        std::getline( file, line );
        std::cout << line << std::endl;
    }
    
    file.close();
    return true;
}
}