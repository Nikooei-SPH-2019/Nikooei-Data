//****************************************************************************************
//
// SPH Code
//
// created	85/12 By R. Ramazani
// Last change	86/01 By Fatehi
// Last change	87/11 By Fatehi
// Last change	88/01 By Fatehi
// Last change	88/04 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************

#ifndef _GENERAL_H_
#define _GENERAL_H_

#define D2
//#define TWOPHASE_PERIODIC

#ifdef D2	
#define DIM 2
#else 
#define DIM 3
#endif


#define WithoutCorr



#ifdef Corr	
#define Corr
#else 
#define WithoutCorr
#endif
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
 

const double pi = 3.1415926535897932384626433832795028841971693993751;
const double e = 2.7182818284590452353602874713526624977572470937000;

typedef std::vector<int> intvec;
typedef std::vector<double> doublevec;

template<class T>
inline T max(const T &a, const T &b)
{
	return b > a ? b : a;
}

template<class T>
inline T min(const T &a, const T &b)
{
	return b < a ? b : a;
}


template<class T>
inline T mean(const T &a, const T &b)
{
	return 0.5*(a + b);
}




template<class T>
inline void Swap(T &a, T &b)
{
	T temp=a;
	a=b;
	b=temp;
}

inline void FatalError(const std::string err)
{
	std::cerr << "The following error has been occured: " << std::endl;
	std::cerr << err << std::endl;
	std::cerr << "The program is going to exit." << std::endl;
	//_getch();
	exit(1);
}

 class TRunError : public std::runtime_error {
 public:
   TRunError(const std::string& s)
     : std::runtime_error(s)
     { }
 };
 inline double ConvertToDouble(const std::string& s)
 {
   std::istringstream i(s);
   double x;
   if (!(i >> x))
     throw TRunError("convertToDouble(\"" + s + "\")");
   return x;
 } 


#endif
