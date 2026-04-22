// Copyright (c) 2026 Eero
// SPDX-License-Identifier: MIT
/**
 * @brief Cross-platform stopwatch based on a chosen chrono time unit.
 *
 * @tparam TimeUnit Desired duration unit, such as std::chrono::milliseconds.
 */
#include <chrono>

template <class TimeUnit>
class TickCounter
{
    std::chrono::high_resolution_clock::time_point startTime;
    std::chrono::high_resolution_clock::time_point stopTime; 
    bool isRunning;
  public:
    /** @brief Construct and reset the timer. */
    TickCounter() : isRunning(false)
    {
      this->reset();
    }
    
    /** @brief Reset the timer to the current instant. */
    void reset() 
    {
		startTime = std::chrono::high_resolution_clock::now();
		stopTime = startTime;
    }
    
    /** @brief Start measuring elapsed time. */
    void start() 
    {
		if (isRunning)
		{
			this->stop();
			this->reset();
		}
	
		isRunning = true;
		startTime = std::chrono::high_resolution_clock::now();
    }
    /** @brief Stop measuring elapsed time. */
    void stop()  
    {
		if (isRunning)
		{
			stopTime = std::chrono::high_resolution_clock::now();
			isRunning = false;
		}
    }

	
    /**
     * @brief Return the elapsed time in ticks.
     *
     * If the timer is running, the interval runs from start to now.
     * Otherwise, the interval runs from start to stop.
     *
     * @return Elapsed ticks as a size_t.
     */
    size_t getElapsed()
    {
		TimeUnit duration;
	
		if (isRunning)
		{
			duration = std::chrono::duration_cast <TimeUnit> 
			( std::chrono::high_resolution_clock::now() - startTime );
		}
		else
		{
			duration = std::chrono::duration_cast<TimeUnit> 
			(stopTime - startTime);
		}
		return static_cast <size_t> ( duration.count() );
    }


  };
    
  