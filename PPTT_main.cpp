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

using namespace std;

int e;
const long numPhotons = 100000;
const int numThreads = 2;
std::thread myThreads[numThreads];

int main() {

    clock_t start, end;

    start = clock();
    Medium * m = new Medium;
    Heat * h = new Heat;

    //  inserting brain layers
    m->CreateCube(0, 0, 0, 10, 10, 10, 0.019, 7.8, 0.89, 1.37);		// scalp and skull
    //m->CreateBall(1,1,1,1,0.02,9.0,0.89,1.37);
    m->CreateCube(2, 2, 2, 6, 6, 6, 0.004, 0.009, 0.89, 1.37);	// cerebro-spinal fluid
    m->CreateCube(3, 3, 3, 4, 4, 4, 0.02, 9.0, 0.89, 1.37);		// gray-matter
    //m->CreateCube(6, 6, 6, 2, 2, 2, 0.08, 40.9, 0.84, 1.37);		// white-matter

    h->AddThermalCoef(m, 2, 3.800, 0.001000, 0.000500);                     // spinus
    h->AddThermalCoef(m, 1, 1.590, 0.001520, 0.000650);                       // bone
    h->AddThermalCoef(m, 3, 3.680, 0.001030, 0.000565);                     // grey-matter 
//	h->AddThermalCoef(m, 4, 3.600, 0.001030, 0.000505);                     // white-matter
    //m->PrintMediumProperties();

    Source * s = new Source;
    s->Collimated_gaussian_beam(5.0, 5.0, 0.0, 0.5, 0.0, 0.0, 1.0);
    //     for(long i = 0; i < numPhotons; i++)
      //       RunPhotonNew(m, s);

    for(long i = 0; i < numThreads; i++)
    {
        myThreads[i] = thread(CreateNewThread, m, s, (long)floor(numPhotons / numThreads));
    }

    for(long i = 0; i < numThreads; i++)
    {
        myThreads[i].join();
    }

    m->RescaleEnergy_Time(numPhotons, time_step);
    //m->RecordFluence();

    //WriteAbsorbedEnergyToFile_Time(m);
   // WritePhotonFluenceToFile(m);

    delete m;
    delete s;
    delete h;

    end = clock();

    cout << "Simulation duration was " << (float)(end - start) / CLOCKS_PER_SEC << " seconds." << endl;
    system("pause");
    exit(0);

}
