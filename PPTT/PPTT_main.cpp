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
int openCLPlatform = 0;
int openCLDevice = 0;
long numPhotons = 0;
int timeSelection = 0;
bool debugMode = 1;
int viewID = 0;
bool saveResultsToFile = false;
float power = 0; // laser power
clock_t startTime, endTime, simulationStart, simulationEnd;
void SelectMode()
{
    if(!debugMode)
    {
        timeSelection = ChooseSteadyOrTime();
        numPhotons = HowManyPhotons();
        usePlatform = OpenCLOrCPU();
        saveResultsToFile = SaveFile();
        power = LaserPower();
    }
    else
    {
        timeSelection = STEADY_STATE;
        numPhotons = 1000000;
        power = 10000; //  now this is bulgarian constant representing power of laser
        usePlatform = OPENCL_MODE;  // 1 - openCL 2 - CPU
        openCLPlatform = 0; // preselect platform 
        openCLDevice = 0; // preselect device
    }
    switch(timeSelection)
    {
        case 1: // steady state
        viewID = 6;
        break;
        case 2: // time resolved
        viewID = 1;
        break;
    }
}

int main(int argc, char *argv[]) {
    Introduction();
    SelectMode();
    startTime = clock();
    Medium * m = new Medium;
    Heat * h = new Heat;
    Source * s = new Source;
    s->Collimated_launch(4, 4, 1.0, 0, 0, 1);
    CreateEnviroment(m, h);

    switch(usePlatform)
    {
        case OPENCL_MODE: // OpenCL
        {

            OpenCL *cl = new OpenCL(m, h, s, debugMode, openCLDevice, openCLPlatform, numPhotons, timeSelection);
            cl->DetectOpenCLDevices();
            // compute Photons travel
            cl->InitPhotonCompute();
            cl->CopyIntoOpenCLStructuresPhoton();
            cl->CopyAndExecuteKernelParametersPhoton();
            cl->CopyResultsPhoton();

            // compute Heat Transfer
            cl->InitHeatCompute();
            cl->CopyIntoOpenCLStructuresHeat(power);
            cl->CopyAndExecuteKernelParametersHeat(1);
            cl->CopyResultsHeat();

            // cleanup
            cl->ReleaseOpenCLStructures();
        }
        break;
        case C_MODE: // CPU
        {
            switch(timeSelection)
            {
                case STEADY_STATE:
                simulationStart = clock();
                CreateNewThread_steady(m, s, numPhotons);
                std::cout << "Photon Simulation duration was " << (float)(clock() - simulationStart) / CLOCKS_PER_SEC << " seconds." << endl;
                m->RescaleEnergy(numPhotons);
                simulationStart = clock();
                h->PennesEquation(m, 36);
                std::cout << "Heat Simulation duration was " << (float)(clock() - simulationStart) / CLOCKS_PER_SEC << " seconds." << endl;
                break;
                case TIME_RESOLVED:
                simulationStart = clock();
                CreateNewThread_time(m, s, numPhotons);
                std::cout << "Photon Simulation duration was " << (float)(clock() - simulationStart) / CLOCKS_PER_SEC << " seconds." << endl;
                m->RescaleEnergy_Time(numPhotons, timeStep);
                break;
            }
            break;



        }
    }
    endTime = clock();
    cout << "Simulation duration was " << (float)(endTime - startTime) / CLOCKS_PER_SEC << " seconds." << endl;

    if(saveResultsToFile)
    {
        m->RecordFluence();
        WriteAbsorbedEnergyToFile_Time(m);
        WritePhotonFluenceToFile(m);
    }

    GLView * view = new GLView();
    view->SaveMedium(m, h, s, viewID, timeSelection);
    view->Init(argc, argv); //Initialize rendering
    view->Run();
    delete s;
    delete h;
    delete m;

    system("pause");
    exit(0);

}
