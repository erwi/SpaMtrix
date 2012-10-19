#ifndef WRITER_H
#define WRITER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// SpaMtrix INCLUDES
#include <setup.h>
#include <vector.h>
#include <ircmatrix.h>
class Writer
{
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

#endif
