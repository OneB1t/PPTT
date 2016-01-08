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
#include "CL\CL.h"
#include <GL\glut.h>
#include "PPTT_io.h"
#include "PPTT_core.h"
#include "GLView.h"
#include "PPTT_OpenCLinit.h"

using namespace std;
int usePlatform = 1; // 1 - openCL 2 - CPU
int openCLPlatform;
int openCLDevice;
long numPhotons = 0;
int timeSelection = 0;
bool debugMode = 1;
clock_t startTime, endTime, simulationStart, simulationEnd;

int main(int argc, char *argv[]) {
    Introduction();
    SelectMode();
    startTime = clock();
    Medium * m = new Medium;
    Heat * h = new Heat;
    Source * s = new Source;

    //  inserting brain layers
    m->CreateCube(0, 0, 0, 15, 15, 15, 0.07, 37.4, 0.977, 1.37);		// adipose tissue @ 700 nm
    m->CreateBall(1, 1, 1, 0.5, 0.02, 9.0, 0.89, 1.37);

    //m->CreateCube(1, 1, 0, 8, 8, 8, 0.15, 1.67, 0.7, 1.37);		// AuNR in intralipid
    //m->CreateCube(1, 2.5, 1, 2.8, 2.8, 2.8, 0.045, 89.5, 0.96, 1.37);	// breast carcinoma @ 700nm
    //m->CreateCube(1, 1, 1, 1.2, 1.2, 1.2, 0.045, 89.5, 0.96, 1.37);	// breast carcinoma @ 700nm
    //m->CreateCube(6, 6, 6, 2, 2, 2, 0.08, 40.9, 0.84, 1.37);		// white-matter

    h->AddThermalCoef(m, 2, 3.800, 0.001000, 0.000500, 0.001);                     // spinus
    h->AddThermalCoef(m, 1, 1.590, 0.001520, 0.000650, 0.001);                       // bone
    h->AddThermalCoef(m, 3, 3.680, 0.001030, 0.000565, 0.001);                     // grey-matter 
	h->AddThermalCoef(m, 4, 3.600, 0.001030, 0.000505, 0.001);                     // white-matter

    s->CollimatedGaussianBeam(0.0, 2.0, 2.5, 5.0, 0.0, 1.0, 0.0); // this causing crash with big number of photons if used for each of them so moved back to main


    switch(usePlatform)
    {
        case 1: // OpenCL
        {
            
            OpenCL *cl = new OpenCL(m,h,s,debugMode,openCLDevice,openCLPlatform,numPhotons,timeSelection);
            cl->DetectOpenCLDevices();
            // compute Photons travel
            cl->InitPhotonCompute();
            cl->CopyIntoOpenCLStructuresPhoton();
            cl->CopyAndExecuteKernelParametersPhoton();
            cl->CopyResultsPhoton();

            // compute Heat Transfer
            cl->InitHeatCompute();
            cl->CopyIntoOpenCLStructuresHeat();       
            cl->CopyAndExecuteKernelParametersHeat(10);
            cl->CopyResultsHeat();

            // cleanup
            cl->ReleaseOpenCLStructures();
        }
        break;
        case 2: // CPU
        {
            switch(timeSelection)
            {
                case 1:
                CreateNewThread_steady(m, s, numPhotons);
                break;
                case 2:
                CreateNewThread_time(m, s, numPhotons);
                break;
                default:
                CreateNewThread_time(m, s, numPhotons);
                break;
            }
            break;
            m->RescaleEnergy(numPhotons);
            m->RescaleEnergy_Time(numPhotons, timeStep);
            h->PennesEquation(m, 36);
        }
    }
    endTime = clock();
    cout << "Simulation duration was " << (float)(endTime - startTime) / CLOCKS_PER_SEC << " seconds." << endl;
    //m->RecordFluence();

    //WriteAbsorbedEnergyToFile_Time(m);
    // WritePhotonFluenceToFile(m);

    GLView * view = new GLView();
    view->SaveMedium(m, h);
    view->Init(argc, argv); //Initialize rendering
    view->Run();
    delete s;
    delete h;

    system("pause");
    exit(0);

}
void SelectMode()
{
    if(!debugMode)
    {
        timeSelection = ChooseSteadyOrTime();          // 1 for steady state, 2 for time resolved
        numPhotons = HowManyPhotons();
        usePlatform = OpenCLOrCPU();
    }
    else
    {
        timeSelection = 1;          // 1 for steady state, 2 for time resolved
        numPhotons = 10000000;
        usePlatform = 1;  // 1 - openCL 2 - CPU
        openCLPlatform = 0;
        openCLDevice = 0;
    }
}