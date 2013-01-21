#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <spamtrix_writer.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_ircmatrix.hpp>
namespace SpaMtrix
{

Writer::~Writer(){
}

bool Writer::writeCSV(const std::string &filename,
		      const Vector &data,
		      idx numC){
/*!
*  Writes the contents of a vector into a comma separated 
*  text file. The line-lengths can be limited by supplying 
*  a non-zero numC value. This can be used e.g. to write 
*  data from a negular 2D FD grid to a file that can be read 
*  and plotted using MATLAB/Octave/spreadsheet
*/
    int numD = data.getLength();
    if (numC==0){
        numC = numD;
    }
    // ATTEMPT TO OPEN FILE FOR WRITIGN. RETURN FALSE IF FAILS  
    std::fstream file;
    file.open(filename.c_str() , std::fstream::out);
    if (!file.is_open()){
        return false;
    }
    // WRITE ALL DATA  
    idx cc = 0;
    for (int i = 0 ; i < numD ; i++ ){
        file << data[i]; 
        cc++;
        if (cc < numC){ // IF NOT END OF ROW ...
            file << ",";
        }
        else{ // ROW END REACHED - NEW LINE
            file << "\n";
            cc = 0;
        }
    }// end for i
    file.close();	
    return true;
} // end writeCSV


bool Writer::writeMatrixMarket(const std::string &filename,
                               const IRCMatrix &A){
/*!
* Writes the sparse matrix A to a file in matrix market format.
* NOTE: Matrix Market format indexing starts from 1, i.e.
* A[1,1] is the first row, column of matrix A, not A[0,0].
*/
    std::fstream file;
    // OPEN A FILE FOR WRITING. RETURN FALSE IF FAILS
    file.open(filename.c_str(), std::fstream::out);
    if (!file.is_open()){
        return false;
    }
    // WRITE HEADER
    file << "% MatrixMarket matrix coordinate real general" << std::endl;
    file << "% NOTE: Matrix market indexing is 1-based!"<< std::endl;
    // WRITE MATRIX SIZE DESCRIPTOR
    file << A.getNumRows() << " " << A.getNumCols() << " " << A.getnnz() << std::endl;
    // WRITE MATRIX DATA
    // TODO: NEED TO IMPLEMENT ITERATORS OVER NON-ZEROS TO MAKE THIS FASTER!!
    for (idx r = 0 ; r < A.getNumRows() ; r++){
        for (idx c = 0 ; c < A.getNumCols() ; c++){
            real val;
            if ( A.isNonZero(r,c, val) ){
                file << r+1 <<"\t" << c+1 << "\t" << val << std::endl;
            }
        }
    }
    return true;
}// end bool writeMatrixMarket
} // end namespace SpaMtrix
