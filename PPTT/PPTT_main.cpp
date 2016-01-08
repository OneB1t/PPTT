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
int viewID = 0;
clock_t startTime, endTime, simulationStart, simulationEnd;
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
        numPhotons = 5000;
        usePlatform = 1;  // 1 - openCL 2 - CPU
        openCLPlatform = 0;
        openCLDevice = 0;
        viewID = 6;
    }
}

int main(int argc, char *argv[]) {
    Introduction();
    SelectMode();
    startTime = clock();
    Medium * m = new Medium;
    Heat * h = new Heat;
    Source * s = new Source;
    s->Collimated_launch(4, 4, 0.1, 0, 0, 1); // this causing crash with big number of photons if used for each of them so moved back to main
    CreateEnviroment(m, h);

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
                m->RescaleEnergy(numPhotons);
                h->PennesEquation(m, 36);
                break;
                case 2:
                CreateNewThread_time(m, s, numPhotons);
                m->RescaleEnergy_Time(numPhotons, timeStep);
                break;
            }
            break;
            
            
            
        }
    }
    endTime = clock();
    cout << "Simulation duration was " << (float)(endTime - startTime) / CLOCKS_PER_SEC << " seconds." << endl;
    //m->RecordFluence();

    //WriteAbsorbedEnergyToFile_Time(m);
    // WritePhotonFluenceToFile(m);

    GLView * view = new GLView();
    view->SaveMedium(m, h,viewID);
    view->Init(argc, argv); //Initialize rendering
    view->Run();
    delete s;
    delete h;

    system("pause");
    exit(0);

}
