/**********************************************************/
/*    PPTT 0.1											  */
/*    Created by Jakub Mesicek							  */
/*      10/2015											  */
/**********************************************************/
/*														  */
/* PPTT_main is a main source file calling the whole      */
/* program												  */
/**********************************************************/

#include <iostream>
#include <ctime>
#include <thread>
#include "PPTT_core.h"
#include "PPTT_io.h"

using namespace std;

int e;
long numPhotons = 1000;
const int numThreads = 2;
std::thread myThreads[numThreads];

int main() {

    clock_t start, end;

	Introduction();
	int Time_Selection = ChooseSteadyOrTime();          // 1 for steady state, 2 for time resolved
	numPhotons = HowManyPhotons();

    start = clock();
    Medium * m = new Medium;
    Heat * h = new Heat;

    //  inserting brain layers
    m->CreateCube(0, 0, 0, 10, 10, 10, 0.07, 37.4, 0.977, 1.37);		// adipose tissue @ 700 nm
    //m->CreateBall(1,1,1,1,0.02,9.0,0.89,1.37);
  
    m->CreateCube(1, 1, 0, 8, 8, 8, 0.15, 1.67, 0.7, 1.37);		// AuNR in intralipid
	m->CreateCube(2, 2, 2, 6, 6, 6, 0.045, 29.5, 0.96, 1.37);	// breast carcinoma @ 700nm
    //m->CreateCube(6, 6, 6, 2, 2, 2, 0.08, 40.9, 0.84, 1.37);		// white-matter

    h->AddThermalCoef(m, 2, 3.800, 0.001000, 0.000500, 0);                     // spinus
    h->AddThermalCoef(m, 1, 1.590, 0.001520, 0.000650, 0);                       // bone
    h->AddThermalCoef(m, 3, 3.680, 0.001030, 0.000565, 0);                     // grey-matter 
//	h->AddThermalCoef(m, 4, 3.600, 0.001030, 0.000505);                     // white-matter
    //m->PrintMediumProperties();

    Source * s1 = new Source;
    
    //     for(long i = 0; i < numPhotons; i++)
      //       RunPhotonNew(m, s);

	if (Time_Selection == 2)
	{
		for (long i = 0; i < numThreads; i++)
		{
			myThreads[i] = thread(CreateNewThread_time, m, s1, (long)floor(numPhotons / numThreads));
		}

		for (long i = 0; i < numThreads; i++)
		{
			myThreads[i].join();
		}
		m->RescaleEnergy_Time(numPhotons, time_step);
		WriteAbsorbedEnergyToFile_Time(m);
	}
	else
	{
		for (long i = 0; i < numThreads; i++)
		{
			myThreads[i] = thread(CreateNewThread_steady, m, s1, (long)floor(numPhotons / numThreads));
		}

		for (long i = 0; i < numThreads; i++)
		{
			myThreads[i].join();
		}
		m->RescaleEnergy(numPhotons);
		WriteAbsorbedEnergyToFile(m);

		h->PennesEquation(m,35);
		WriteTemperature(m, h);
	}
	//	Run second pulse
	Source * s2 = new Source;
	/*Prepare_SecondPulse(m, s2, 0.1);
	for (long i = 0; i < numThreads; i++)
	{
		myThreads[i] = thread(CreateNewThread_secondPulse, m, s2, (long)floor(numPhotons / numThreads));
	}

	for (long i = 0; i < numThreads; i++)
	{
		myThreads[i].join();
	}*/

	
	//m->RecordFluence();

	
	// WritePhotonFluenceToFile(m);

	/*m->RescaleEnergy_Time_secondPulse(numPhotons, time_step);
	WriteAbsorbedEnergyToFile_Time_secondPulse(m);*/

    delete m;
	delete s1; delete s2;
    delete h;

    end = clock();

    cout << "Simulation duration was " << (float)(end - start) / CLOCKS_PER_SEC << " seconds." << endl;
    system("pause");
    exit(0);
}
