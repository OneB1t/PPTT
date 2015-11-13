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
#include "simpleCL.h"



using namespace std;

int e;
const long numPhotons = 800;
const int numThreads = 1;
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


                                                                    /**
                                                                    * Initialize OpenCL vectors:
                                                                    **/

    cl_int error;
    cl_platform_id platform;
    cl_device_id device;
    cl_uint platforms, devices;
    char build_c[4096];
    size_t srcsize, worksize;

    /* Fetch the Platforms, we only want one. */
    error = clGetPlatformIDs(1, &platform, &platforms);
    if(error != CL_SUCCESS) {
        printf("\n Error number %d", error);
    }
    /* Fetch the Devices for this platform */
    error = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 1, &device, &devices);
    if(error != CL_SUCCESS) {
        printf("\n Error number %d", error);
    }
    /* Create a memory context for the device we want to use  */
    cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform,0 };
    /* Note that nVidia's OpenCL requires the platform property */
    cl_context context = clCreateContext(properties, 1, &device, NULL, NULL, &error);
    if(error != CL_SUCCESS) {
        printf("\n Error number %d", error);
    }
    /* Create a command queue to communicate with the device */
    cl_command_queue cq = clCreateCommandQueue(context, device, 0, &error);
    if(error != CL_SUCCESS) {
        printf("\n Error number %d", error);
    }


    FILE *fp;
    char fileName[] = "c:/PPTT/PPTT/OpenCL/photoncompute.cl";
    char *source_str;
    size_t source_size;

    /* Load the source code containing the kernel*/
    fp = fopen(fileName, "r");
    if(!fp) {
        fprintf(stderr, "Failed to load kernel.\n");
    }
    
    source_str = (char*)malloc(MAX_SOURCE_SIZE);
    srcsize = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
    fclose(fp);

    const char *srcptr[] = { source_str };
    /* Submit the source code of the kernel to OpenCL, and create a program object with it */
    cl_program prog = clCreateProgramWithSource(context,1, srcptr, &srcsize, &error);
    if(error != CL_SUCCESS) {
        printf("\n Error number %d", error);
    }

    /* Compile the kernel code (after this we could extract the compiled version) */
    error = clBuildProgram(prog, 0, NULL, "", NULL, NULL);
    if(error != CL_SUCCESS) {
        printf("Error on buildProgram ");
        printf("\n Error number %d", error);
        fprintf(stdout, "\nRequestingInfo\n");
        clGetProgramBuildInfo(prog, 0, CL_PROGRAM_BUILD_LOG, 4096, build_c, NULL);
        printf("Build Log for %s_program:\n%s\n", "example", build_c);
    }
    else
    {
        printf("compile sucessfull");
    }
    

    cl_kernel computePhoton = clCreateKernel(prog, "computePhoton", &error);
    if(error != CL_SUCCESS) {
        printf("\n Error number %d", error);
    }

    cl_int status = 0;
    my_struct* ms = new my_struct[5];

    cl_mem mem = clCreateBuffer(context, 0, sizeof(my_struct) * 5, NULL, &status);
    clEnqueueWriteBuffer(cq, mem, CL_TRUE, 0, sizeof(mem) * 5, &ms, 0, NULL, NULL);

    status = clSetKernelArg(computePhoton, 0, sizeof(mem), &mem);

    size_t global[] = { 5 };
    status = clEnqueueNDRangeKernel(cq, computePhoton, 1, NULL, global, NULL, 0, NULL, NULL);

    status = clEnqueueReadBuffer(cq, mem, CL_TRUE, 0, sizeof(my_struct) * 5, ms, 0, NULL, NULL);

    for(int i = 0; i < 5; i++)
        cout << (ms + i)->a << " " << (ms + i)->b << " " << (ms + i)->c << endl;

    cout << ms[0].energy[0][0][0] << endl;
    /* Tell the Device, through the command queue, to execute que Kernel */
    if(error != CL_SUCCESS) {
        printf("\n Error number %d", error);
    }
    /* Await completion of all the above */
    error = clFinish(cq);
    if(error != CL_SUCCESS) {
        printf("\n Error number %d", error);
    }
    /* Finally, output the result */
    system("pause");
    return 0;


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
