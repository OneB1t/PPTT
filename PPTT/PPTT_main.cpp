/**********************************************************/
/*    PPTT 0.1											  */
/*    Created by Jakub Mesicek							  */
/*      10/2015											  */
/**********************************************************/
/*														  */
/* PPTT_main is a main source file calling the whole      */
/* program												  */
/**********************************************************/
#define VIENNACL_WITH_OPENCL
#include <iostream>
#include <ctime>
#include <thread>
#include "PPTT_core.h"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"


using namespace std;
using namespace viennacl;

int e;
const long numPhotons = 800000;
const int numThreads = 8;
thread myThreads[numThreads];

int main(int argc, char *argv[]) {

    clock_t start, end;

    start = clock();
    Medium * m = new Medium;
    Heat * h = new Heat;

    thread tList[numThreads];

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
    s->Collimated_gaussian_beam(5.0, 5.0, 0.0, 0.5, 0.0, 0.0, 1.0); // this causing crash with big number of photons if used for each of them so moved back to main


    // fully working viennacl
    viennacl::scalar<float> vcl_s1;
    viennacl::scalar<float> vcl_s2 = 1.0;
    viennacl::scalar<float> vcl_s3 = 1.0;

    vcl_s1 = vcl_s2 + vcl_s3;
    std::cout << "GPU scalar vcl_s2: " << vcl_s1 << std::endl;



    for(long i = 0; i < numThreads; i++)
    {
        tList[i] = thread(CreateNewThread, m, s, (long)floor(numPhotons / numThreads));
    }

    for(long i = 0; i < numThreads; i++)
    {
        tList[i].join();
    }

    m->RescaleEnergy_Time(numPhotons, time_step);
    //m->RecordFluence();

    WriteAbsorbedEnergyToFile_Time(m);
    // WritePhotonFluenceToFile(m);

    delete m;
    delete s;
    delete h;

    end = clock();

    cout << "Simulation duration was " << (float)(end - start) / CLOCKS_PER_SEC << " seconds." << endl;
    system("pause");
    exit(0);

}
