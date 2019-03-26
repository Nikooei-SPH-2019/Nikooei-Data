//****************************************************************************************
//
// SPH Code
//
// created	85/12 By R. Ramazani
// changed		85/12 By Fatehi
//  Last change	97/11 By Nikooei
//****************************************************************************************
// Timer file
//	Calculates runnung time of the program and whole CPU time
//****************************************************************************************

#include "timer.h"

void Timer::Start()
{
	cks=clock();
}

void Timer::Stop()
{
	cke=clock();
}

void Timer::Pause()
{
	ckp=clock();
}

void Timer::Resume()
{
	cks=clock()-(ckp-cks);
}

double Timer::CPUTime()
{
	return double(cke-cks)/CLOCKS_PER_SEC;
}

void PrintTotalTime(std::time_t ts, std::time_t tf)
{
	long int h, m, s;
	//char timebuf[26];
	//ctime_s(timebuf, 26, &ts);
	char* timebuf;			//use in gcc
	timebuf = ctime(&ts);		//use in gcc
	
	std::cout << std::endl;
	std::cout << "-------------------------------------------------------------------------------" << std::endl;
	std::cout << "Computation started at: " << timebuf;
	//ctime_s(timebuf, 26, &tf);
	timebuf = ctime(&tf);		//use in gcc
	std::cout << "and finished at: " << timebuf;
	s=(long int)(tf-ts);
	h=int(s/3600);
	s-=3600*h;
	m=int(s/60);
	s-=60*m;
	std::cout << "Total time of computation is: ";
	std::cout.width(4);
	std::cout << h << " hours, ";
	std::cout.width(2);
	std::cout << m << " minutes and ";
	std::cout.width(2);
	std::cout << s << " seconds ";
	std::cout << std::endl;
}
