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
#include "GLView.h"
#include "CL\CL.h"
#include <GL\glut.h>
#include "clRNG/clRNG.h"
#include "clRNG/mrg31k3p.h"



using namespace std;

int e;
const long numBatches = 1;
const long numPhotons = 32000 * numBatches;
const int numThreads = 1;
thread myThreads[numThreads];


int main(int argc, char *argv[]) {

    clock_t start, end,batch_start,batch_end;

    start = clock();
    Medium * m = new Medium();
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
    cl_event event;
    cl_int status = 0;
    cl_uint platforms, devices;
    clrngMrg31k3pStream *streams = 0;
    size_t streamBufferSize = 0;
    char *clrng_root;
    char include_str[1024];
    char build_c[8192];
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

    /* Make sure CLRNG_ROOT is specified to get library path */
    clrng_root;
    if(clrng_root == NULL) printf("\nSpecify environment variable CLRNG_ROOT as described\n");
    strcpy(include_str, "-I ");
    strcat(include_str, clrng_root);
    strcat(include_str, "/include");

    FILE *fp;
    char fileName[] = "photoncompute.cl";
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
        clGetProgramBuildInfo(prog, 0, CL_PROGRAM_BUILD_LOG, 8192, build_c, NULL);
        printf("Build Log for %s_program:\n%s\n", "example", build_c);
    }
    else
    {
        printf("compile sucessfull \n");
    }
    

    cl_kernel computePhoton = clCreateKernel(prog, "computePhoton", &error);
    if(error != CL_SUCCESS) {
        printf("\n Error number %d", error);
    }

    // this is usable for more than 1 medium and source
    m_str* ms = new m_str[1];
    s_str * ss = new s_str[1];

    // copy all into openCL structures

    ss[0].release_time = s->release_time;
    ss[0].x = s->x;
    ss[0].y = s->y;
    ss[0].z = s->z;
    ss[0].ux = s->ux;
    ss[0].uy = s->uy;
    ss[0].uz = s->uz;

   for(int temp = 0; temp < voxels_x; temp++)			// matrix with photon fluence inicialization
        for(int temp2 = 0; temp2 < voxels_y; temp2++)
            for(int temp3 = 0; temp3 < voxels_z; temp3++)
            {
                ms[0].energy[temp][temp2][temp3] = m->energy[temp][temp2][temp3];
                ms[0].structure[temp][temp2][temp3] = m->structure[temp][temp2][temp3];
                ms[0].fluence[temp][temp2][temp3] = m->fluence[temp][temp2][temp3];
            }


    for(int temp = 0; temp < max_regions; temp++)
    {
        ms[0].ua[temp] = m->ua[temp];
        ms[0].us[temp] = m->us[temp];
        ms[0].inv_albedo[temp] = m->inv_albedo[temp];
        ms[0].g[temp] = m->g[temp];
        ms[0].n[temp] = m->n[temp];
        ms[0].k[temp] = m->k[temp];
        ms[0].rho[temp] = m->rho[temp];
        ms[0].c_h[temp] = m->c_h[temp];
        cout << ms[0].n[temp] << endl;
    }

    // this should go into another structure
    ms[0].time_end = time_end;
    ms[0].time_start = time_start;
    ms[0].pulseDuration = pulseDuration;
    ms[0].time_step = time_step;

    size_t globalWorkItems[] = { numPhotons / numBatches };  // basically number of photons per batch

    // copy memory buffers to GPU
    cl_mem mediumMemoryBlock = clCreateBuffer(context, 0, sizeof(ms[0]), NULL, &status);
    clEnqueueWriteBuffer(cq, mediumMemoryBlock, CL_TRUE, 0, sizeof(ms[0]), &ms[0], 0, NULL, NULL);
    status = clSetKernelArg(computePhoton, 0, sizeof(ms), &mediumMemoryBlock);

    cl_mem structureMemoryBlock = clCreateBuffer(context, 0, sizeof(ss[0]), NULL, &status);
    clEnqueueWriteBuffer(cq, structureMemoryBlock, CL_TRUE,0 , sizeof(ss[0]), &ss[0], 0, NULL, NULL);
    status = clSetKernelArg(computePhoton, 1, sizeof(ss), &structureMemoryBlock);

    /* Create streams for clRNG*/
    streams = clrngMrg31k3pCreateStreams(NULL, numPhotons / numBatches, &streamBufferSize, (clrngStatus *)&status);
    cl_mem randomNumber = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, streamBufferSize, streams, &status);
    status = clSetKernelArg(computePhoton, 2, sizeof(randomNumber), &randomNumber);


    for(int i = 0; i < numBatches;i++)
    {
        batch_start = clock();
        status = clEnqueueNDRangeKernel(cq, computePhoton, 1, NULL, globalWorkItems, NULL, 0, NULL, &event);
        cout << "run number:" << i << " ";
        clWaitForEvents(1, &event);
        batch_end = clock();
        cout << "Batch Simulation duration was " << (float)(batch_end - batch_start) / CLOCKS_PER_SEC << " seconds." << endl;
    }

    status = clEnqueueReadBuffer(cq, mediumMemoryBlock, CL_TRUE, 0, sizeof(ms[0]), ms,0, NULL, NULL);

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
    end = clock();

    cout << "Simulation duration was " << (float)(end - start) / CLOCKS_PER_SEC << " seconds." << endl;

    int num_time_steps = (int)ceil((time_end - time_start) / time_step);

    for(int temp = 0; temp < voxels_x; temp++)
        for(int temp2 = 0; temp2 < voxels_y; temp2++)
            for(int temp3 = 0; temp3 < voxels_z; temp3++)
                for(int temp4 = 0; temp4 < num_time_steps; temp4++)
                    m->energy_t[temp][temp2][temp3][temp4] = ms[0].energy_t[temp][temp2][temp3][temp4];

    m->RescaleEnergy_Time(numPhotons, time_step);
    //m->RecordFluence();

    //WriteAbsorbedEnergyToFile_Time(m);
    // WritePhotonFluenceToFile(m);

    GLView * view = new GLView();
    view->savemedium(m);
    view->init(argc, argv); //Initialize rendering
    view->run();
    delete s;
    delete h;

    system("pause");
    exit(0);

}