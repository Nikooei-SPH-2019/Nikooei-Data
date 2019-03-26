//****************************************************************************************
//
// SPH Code
//
// created	85/12 By R. Ramazani
// Last change		85/12 By Fatehi
// Last change	97/11 By Nikooei
//****************************************************************************************

#ifndef _KERNEL_H_
#define _KERNEL_H_

#include "general.h"

double Kernel(double d, double h, short KERNELTYPE);
double DKernel(double d, double h, short KERNELTYPE);

double D2Kernel(double d, double h, short KERNELTYPE);
#endif
