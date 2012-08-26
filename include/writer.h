#ifndef WRITER_H
#define WRITER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// SpaMtrix INCLUDES
#include <setup.h>
#include <vector.h>

class Writer
{
	std::fstream file;
	public:
		Writer();
		~Writer();
		
	bool writeCSV(const std::string &filename,
		      const Vector &data,
		      idx numC=0);

};

#endif
