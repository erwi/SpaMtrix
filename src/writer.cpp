#include <writer.h>

Writer::Writer():
file()
{
}

Writer::~Writer()
{
	if (file.is_open() )
		file.close();

}

bool Writer::writeCSV(const std::string &filename,
		      const Vector &data,
		      idx numC)
{
  /*!
   *  Writes the contents of a vector into a comma separated 
   *  text file. The line-lengths can be limited by supplying 
   *  a non-zero numC value. This can be used e.g. to write 
   *  data from a negular 2D FD grid to a file that can be read 
   *  and plotted using MATLAB/Octave/spreadsheet
   */
    int numD = data.getLength();
    if (numC==0)
      numC = numD;
    
    // ATTEMPT TO OPEN FILE FOR WRITIGN. RETURN FALSE IF FAILS  
    this->file.open(filename.c_str() , std::fstream::out);
      if (!file.is_open())
      	return false;
      
    // WRITE ALL DATA  
      int cc = 0;
      for (int i = 0 ; i < numD ; i++ )
      {
	file << data[i]; 
	cc++;
	if (cc < numC) // IF NOT END OF ROW ...
	  file << ",";
	else // ROW END REACHED - NEW LINE
	{
	  file << "\n";
	  cc = 0;
	}
      }
      
      file.close();	
      return true;
}

