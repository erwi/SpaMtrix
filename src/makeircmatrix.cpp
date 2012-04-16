#include <ircmatrix.h>
#include <algorithm>
#include <iostream>
IRCMatrix& makeIRCMatrix(const idx* t,
                         const idx npt,
                         const idx nt)
{
/*! 
    Makes an IRCMatrix based on the connectivity defined in the array t
*/

// DETEMINE REQUIRED MATRIX SIZE. NUMBER OF ROWS/COLUMNS IS 
// LARGEST VALUE FOUND IN t
idx numCols = *std::max_element(t , t+(npt * nt ) ) + 1;
std::cout<< "numCols " << numCols << std::endl;

// 
IRCMatrix* ircm = new IRCMatrix();
return *ircm;

}