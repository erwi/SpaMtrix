  /*!
    A siple c++11 cross-platform stopwatch template class that mesures 
	time by counting ticks of a specified duration. 
	
	The duration can be one of:
	
		std::chrono::nanoseconds
		std::chrono::microseconds
		std::chrono::milliseconds
		std::chrono::seconds
		std::chrono::minutes
		std::chrono::hours
  */
#include <chrono>

template <class TimeUnit>
class TickCounter
{
    // TIMEPOINTS FOR START AND STOP OF MEASUREMENT
    std::chrono::high_resolution_clock::time_point startTime;
    std::chrono::high_resolution_clock::time_point stopTime; 
    bool isRunning;
  public:
      
    // CONSTRUCTOR
    TickCounter(): 
    isRunning(false)
    {
      this->reset();
    }
    
    void reset() 
    {
	/*!
	* Resets TickCounter by setting start time and stop time to equal.
	*/
		startTime = std::chrono::high_resolution_clock::now();
		stopTime = startTime;
		
		// IF CURRENTLY RUNNING
		/*
		if (isRunning)
		{
			this->stop(); 
			startTime = std::chrono::high_resolution_clock::now();
			stopTime = startTime;
			this->start();
		}	
		else
		{
			// SET START AND STOP TO EQUAL
			startTime = std::chrono::high_resolution_clock::now(); 
			stopTime = std::chrono::high_resolution_clock::now();
		}
		*/
    }
    
    void start() 
    {
	/*!
	* Starts the TickCounter. 
	* If the TickCounter was already running, previous counts are cleared
	*/
		// RESET IF ALREADY RUNNING
		if (isRunning)
		{
			this->stop();
			this->reset();
		}
	
		isRunning = true;
		startTime = std::chrono::high_resolution_clock::now();
    }
    void stop()  
    {
	/*!
	* Stops the TickCounter, if it was running.
	* If the TickCounter was not running, nothing is done
	*/
	
		if (isRunning)
		{
			stopTime = std::chrono::high_resolution_clock::now();
			isRunning = false;
		}
    }

	
	// TODO: USE idx INSTEAD OF size_t
    size_t getElapsed()
    {
	/*!
	* Returns elapsed time as number of ticks of TimeUnits
	* Elapsed time can mean 2 different things:
	* 
	* 1. 	If TickCounter is running, elapsed time is measured between starTime
	* 	and now.
	* 2. 	If TickCounter is not running (i.e. it is stopped), elapsed time
	* 	is measured between startTime and stopTime.
	*/
      
		// DIFFERENCE BETWEEN THE TWO "TIMEPOINTS" MUST BE CAST TO
		// A "DURATION" USING THE std::chrono::duration_cast 
		TimeUnit duration;
	
		// CASE 1 - TickCounter IS RUNNING, USE CURRENT TIME
		if (isRunning)
		{
			duration = std::chrono::duration_cast <TimeUnit> 
			( std::chrono::high_resolution_clock::now() - startTime );
		}
		// CASE 2 - NOT RUNNING, USE stopTime
		else
		{
			duration = std::chrono::duration_cast<TimeUnit> 
			(stopTime - startTime);
		}
		// CAST TO TICK COUNT
		return static_cast <size_t> ( duration.count() );
    }


  };
    
  