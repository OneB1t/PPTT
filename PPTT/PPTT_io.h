/**********************************************************/
/*					    PPTT 0.1						  */
/*				Created by Jakub Mesicek				  */
/*						11/2015							  */
/**********************************************************/
/*														  */
/*	PPTT_io defines input/output functions				  */
/**********************************************************/

#ifndef PPTT_IO_H
#define PPTT_IO_H

#include <iostream>
#include "PPTT_core.h"
#include <fstream>

using namespace std;

void Introduction();
char ChooseSteadyOrTime();
char OpenCLOrCPU();
char OpenCLPlatform();
char OpenCLDevice();
long HowManyPhotons();

void WriteTemperature(Medium * m, Heat * h);

#endif