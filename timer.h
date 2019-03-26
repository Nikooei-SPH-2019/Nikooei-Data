//****************************************************************************************
//
// SPH Code
//
// created	85/12 By R. Ramazani
// Last change		85/12 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************
// Timer file
//	Calculates runnung time of the program and whole CPU time
//****************************************************************************************

#ifndef _TIMER_H_
#define _TIMER_H_

#include "general.h"
#pragma warning(disable : 4996)

	void PrintTotalTime(std::time_t ts, std::time_t tf);

	class Timer
	{
	private:
		clock_t cks;
		clock_t cke;
		clock_t ckp;
	public:
		void Start();
		void Stop();
		void Pause();
		void Resume();

		double CPUTime();
	};

#endif
