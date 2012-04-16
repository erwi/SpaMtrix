#ifndef IRCMATRIX
#define IRCMATRIX
typedef unsigned int idx;
typedef double real;
struct ColVal{
	idx col;	// COLUMN INDEX
	real val;	// VALUE
};



// Interleaved Row Compressed Matrix
class IRCMatrix
{
	idx* rows;			// ROW COUNTER
	ColVal* cvPairs;	// COLUMN-VALUE PAIRS
	idx nnz;			// NUMBER OF NON-ZEROS
public:
	
	IRCMatrix():rows(0), cvPairs(0), nnz(0) {}
	idx getnnz() {return nnz;}

};

// MAKES SIMPLE IRCMatrix OUT OF AN ARRAY OF CONNECTED INDEXES
IRCMatrix& makeIRCMatrix(const idx* t, const idx npt, const idx nt);

#endif
