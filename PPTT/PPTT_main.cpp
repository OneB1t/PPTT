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
#include "PPTT_io.h"

using namespace std;
const int numThreads = 8;
const int usePlatform = 0; // 0 - openCL 1 - CPU c++
long numPhotons = 0;
thread myThreads[numThreads];

int main(int argc, char *argv[]) {
    Introduction();
    clock_t start, end, batch_start, batch_end;

    int Time_Selection = ChooseSteadyOrTime();          // 1 for steady state, 2 for time resolved
    numPhotons = HowManyPhotons();

    start = clock();
    Medium * m = new Medium;
    Heat * h = new Heat;

    //  inserting brain layers
    m->CreateCube(0, 0, 0, 15, 15, 15, 0.07, 37.4, 0.977, 1.37);		// adipose tissue @ 700 nm
    //m->CreateBall(1,1,1,1,0.02,9.0,0.89,1.37);

    m->CreateCube(1, 1, 0, 8, 8, 8, 0.15, 1.67, 0.7, 1.37);		// AuNR in intralipid
    m->CreateCube(2, 2, 2, 6, 6, 6, 0.045, 29.5, 0.96, 1.37);	// breast carcinoma @ 700nm
    //m->CreateCube(6, 6, 6, 2, 2, 2, 0.08, 40.9, 0.84, 1.37);		// white-matter

    h->AddThermalCoef(m, 2, 3.800, 0.001000, 0.000500, 0);                     // spinus
    h->AddThermalCoef(m, 1, 1.590, 0.001520, 0.000650, 0);                       // bone
    h->AddThermalCoef(m, 3, 3.680, 0.001030, 0.000565, 0);                     // grey-matter 
//	h->AddThermalCoef(m, 4, 3.600, 0.001030, 0.000505);                     // white-matter
    //m->PrintMediumProperties();

    Source * s = new Source;
    s->CollimatedGaussianBeam(5.0, 2.0, 2.5, 2.0, 0.0, 1.0, 0.0); // this causing crash with big number of photons if used for each of them so moved back to main
    cl_int error;
    cl_platform_id platform;
    cl_device_id device;
    cl_event event;
    cl_int status = 0;
    cl_uint platforms, devices;
    size_t streamBufferSize = 0;
    char build_c[8192];
    size_t srcsize, worksize;

    switch(usePlatform)
    {
        case 0: // OpenCL
        {
            clErrorCheck(clGetPlatformIDs(1, &platform, &platforms));
            clErrorCheck(clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, &devices));
            cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform,0 };
            cl_context context = clCreateContext(properties, 1, &device, NULL, NULL, &error);
            clErrorCheck(error);
            cl_command_queue cq = clCreateCommandQueueWithProperties(context, device, 0, &error);
            clErrorCheck(error);

            FILE *fp;
            char fileName[] = "photoncompute.cl";

            /* Load the source code containing the kernel*/
            fp = fopen(fileName, "r");
            if(!fp) {
                fprintf(stderr, "Failed to load kernel.\n");
            }
            char * source_str = (char*)malloc(MAX_SOURCE_SIZE);
            srcsize = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
            fclose(fp);

            const char *srcptr[] = { source_str, };
            /* Submit the source code of the kernel to OpenCL, and create a program object with it */
            cl_program prog = clCreateProgramWithSource(context, 1, srcptr, &srcsize, &error);
            clErrorCheck(error);

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
                printf("OpenCL kernel compile sucessfull. \n");
            }
            cl_kernel computePhoton = clCreateKernel(prog, "computePhoton", &error);
            clErrorCheck(error);

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
                ms[0].w_g[temp] = m->w_g[temp];
            }

            // this should go into another structure
            ms[0].time_end = time_end;
            ms[0].time_start = time_start;
            ms[0].pulseDuration = pulseDuration;
            ms[0].time_step = time_step;
            ms[0].finished = 0;
            size_t globalWorkItems[] = { numPhotons };  // basically number of photons per batch
            size_t localWorkItems[] = { 1 };  // basically number of photons per batch

            // copy memory buffers to GPU
            cl_mem mediumMemoryBlock = clCreateBuffer(context, NULL, sizeof(ms[0]), NULL, &status);
            clEnqueueWriteBuffer(cq, mediumMemoryBlock, CL_TRUE, 0, sizeof(ms[0]), ms, 0, NULL, NULL);
            status = clSetKernelArg(computePhoton, 0, sizeof(ms), &mediumMemoryBlock);

            cl_mem structureMemoryBlock = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(ss[0]), NULL, &status);
            clEnqueueWriteBuffer(cq, structureMemoryBlock, CL_TRUE, 0, sizeof(ss[0]), &ss[0], 0, NULL, NULL);
            status = clSetKernelArg(computePhoton, 1, sizeof(ss), &structureMemoryBlock);

            int random = rand();
            cl_mem randomValue = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
            status = clSetKernelArg(computePhoton, 2, sizeof(int), &randomValue);

            batch_start = clock();
            status = clEnqueueNDRangeKernel(cq, computePhoton, 1, NULL, globalWorkItems, localWorkItems, 0, NULL, NULL);
            random = rand();
            clEnqueueWriteBuffer(cq, randomValue, CL_TRUE, 0, sizeof(int), &random, 0, NULL, NULL);
            batch_end = clock();
            cout << "Batch Simulation duration was " << (float)(batch_end - batch_start) / CLOCKS_PER_SEC << " seconds." << endl;

            status = clEnqueueReadBuffer(cq, mediumMemoryBlock, CL_TRUE, 0, sizeof(ms[0]), ms, 0, NULL, NULL);
            clErrorCheck(error);
            error = clFinish(cq);
            clErrorCheck(error);

            end = clock();
            cout << "Number of unfinished photons:" << ms[0].finished << endl;
            cout << "Simulation duration was " << (float)(end - start) / CLOCKS_PER_SEC << " seconds." << endl;

            int num_time_steps = (int)ceil((time_end - time_start) / time_step);
            for(int temp = 0; temp < voxels_x; temp++)
                for(int temp2 = 0; temp2 < voxels_y; temp2++)
                    for(int temp3 = 0; temp3 < voxels_z; temp3++)
                        for(int temp4 = 0; temp4 < num_time_steps; temp4++)
                            m->energy_t[temp][temp2][temp3][temp4] = ms[0].energy_t[temp][temp2][temp3][temp4];

            clReleaseMemObject(mediumMemoryBlock);
            clReleaseMemObject(structureMemoryBlock);
            if(cq != 0)
                clReleaseCommandQueue(cq);
            if(computePhoton != 0)
                clReleaseKernel(computePhoton);
            if(prog != 0)
                clReleaseProgram(prog);
            if(context != 0)
                clReleaseContext(context);
        }
        break;
        case 1: // CPU
        {
            if(Time_Selection == 2)
            {
                for(long i = 0; i < numThreads; i++)
                    myThreads[i] = thread(CreateNewThread_time, m, s, (long)floor(numPhotons / numThreads));
                for(long i = 0; i < numThreads; i++)
                    myThreads[i].join();
            }
            else
            {
                for(long i = 0; i < numThreads; i++)
                    myThreads[i] = thread(CreateNewThread_time, m, s, (long)floor(numPhotons / numThreads));
                for(long i = 0; i < numThreads; i++)
                    myThreads[i].join();
            }
        }
        break;
    }
    m->RescaleEnergy_Time(numPhotons, time_step);
    h->PennesEquation(m, 35);
    //m->RecordFluence();

    //WriteAbsorbedEnergyToFile_Time(m);
    // WritePhotonFluenceToFile(m);

    GLView * view = new GLView();
    view->savemedium(m, h);
    view->init(argc, argv); //Initialize rendering
    view->run();
    delete s;
    delete h;

    system("pause");
    exit(0);

}

void clErrorCheck(cl_int error)
{
    if(error != CL_SUCCESS) {
        printf("\n Error number %d", error);
    }
}