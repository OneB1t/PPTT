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
#include <vector>
#include <iomanip> 
#include "PPTT_core.h"
#include "GLView.h"
#include "CL\CL.h"
#include <GL\glut.h>
#include "PPTT_io.h"

using namespace std;
const int numThreads = 8;
int usePlatform = 2; // 1 - openCL 2 - CPU c++
int openCLPlatform = 0;
int openCLDevice = 0;
long numPhotons = 0;
thread myThreads[numThreads];

int main(int argc, char *argv[]) {
    Introduction();
    clock_t start, end, batch_start, batch_end;

    int Time_Selection = ChooseSteadyOrTime();          // 1 for steady state, 2 for time resolved
    numPhotons = HowManyPhotons();
    usePlatform = OpenCLOrCPU();

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
    cl_uint devices;
    size_t streamBufferSize = 0;
    char build_c[8192];
    size_t srcsize, worksize;

    char name[256];
    char version[256];
    cl_uint num_platforms;

    switch(usePlatform)
    {
        case 1: // OpenCL
        {
            cout << "Detecting OpenCL platforms" << endl;
            char name[256];
            char version[256];
            cl_uint num_platforms;
            clErrorCheck(clGetPlatformIDs(0, NULL, &num_platforms));
            vector<cl_platform_id> platforms(num_platforms);
            clErrorCheck(clGetPlatformIDs(num_platforms, platforms.data(), NULL));
            vector<cl_uint> num_platform_devices(num_platforms);
            cl_uint num_devices = 0;
            for(cl_uint i = 0; i < num_platforms; ++i)
            {
                const auto platform = platforms[i];
                cout << "ID   Name" << endl;
                clErrorCheck(clGetPlatformInfo(platform, CL_PLATFORM_NAME, sizeof(name), name, NULL));
                clErrorCheck(clGetPlatformInfo(platform, CL_PLATFORM_VERSION, sizeof(version), version, NULL));
                cout << i << ' ' << name << ", " << version << endl;
                clErrorCheck(clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 0, NULL, &num_platform_devices[i]));
                num_devices += num_platform_devices[i];
            }
            platform = platforms[OpenCLPlatform()];
            cout << "Detecting OpenCL devices" << endl;
            if(!num_devices)
            {
                cerr << "No OpenCL devices detected" << endl;
                return 2;
            }
            vector<cl_device_id> devices(num_devices);
            vector<bool> cl12(num_devices);
            vector<cl_bool> host_unified_memory(num_devices);
            cout << "ID                          Name  CL CU GMEM(MB) LMEM(KB) CMEM(KB) LMEMTYPE ECC" << endl;
            const char* local_mem_types[] = { "NONE", "LOCAL", "GLOBAL" };
            for(cl_uint i = 0, dev = 0; i < num_platforms; ++i)
            {
                const auto platform = platforms[i];
                const auto npd = num_platform_devices[i];
                clErrorCheck(clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, npd, &devices[dev], NULL));
                for(int d = 0; d < npd; ++d, ++dev)
                {
                    const auto device = devices[dev];
                    cl_uint max_compute_units;
                    cl_ulong global_mem_size;
                    cl_ulong local_mem_size;
                    cl_ulong max_constant_buffer_size;
                    cl_bool error_correction_support;
                    cl_device_local_mem_type local_mem_type;
                    clErrorCheck(clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(name), name, NULL));
                    clErrorCheck(clGetDeviceInfo(device, CL_DEVICE_OPENCL_C_VERSION, sizeof(version), version, NULL));
                    clErrorCheck(clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(max_compute_units), &max_compute_units, NULL));
                    clErrorCheck(clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(global_mem_size), &global_mem_size, NULL));
                    clErrorCheck(clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(local_mem_size), &local_mem_size, NULL));
                    clErrorCheck(clGetDeviceInfo(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(max_constant_buffer_size), &max_constant_buffer_size, NULL));
                    clErrorCheck(clGetDeviceInfo(device, CL_DEVICE_ERROR_CORRECTION_SUPPORT, sizeof(error_correction_support), &error_correction_support, NULL));
                    clErrorCheck(clGetDeviceInfo(device, CL_DEVICE_HOST_UNIFIED_MEMORY, sizeof(host_unified_memory[dev]), &host_unified_memory[dev], NULL));
                    clErrorCheck(clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_TYPE, sizeof(local_mem_type), &local_mem_type, NULL));
                    cl12[dev] = version[9] > '1' || version[11] >= '2';
                    const string name_str(name);
                    const size_t name_len = name_str.size();
                    cout << dev << setw(31) << (name_len <= 30 ? name_str : name_str.substr(name_str.find(' ', name_len - 31) + 1)) << ' ' << version[9] << '.' << version[11] << setw(3) << max_compute_units << setw(9) << global_mem_size / 1048576 << setw(9) << local_mem_size / 1024 << setw(9) << max_constant_buffer_size / 1024 << setw(9) << local_mem_types[local_mem_type] << setw(4) << error_correction_support << endl;
                }
            }
            cl_device_id device = devices[OpenCLDevice()];

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
            status = clEnqueueNDRangeKernel(cq, computePhoton, 1, NULL, globalWorkItems, NULL, 0, NULL, NULL);
            clEnqueueWriteBuffer(cq, randomValue, CL_TRUE, 0, sizeof(int), &random, 0, NULL, NULL);
            batch_end = clock();
            cout << "Batch Simulation duration was " << (float)(batch_end - batch_start) / CLOCKS_PER_SEC << " seconds." << endl;

            status = clEnqueueReadBuffer(cq, mediumMemoryBlock, CL_TRUE, 0, sizeof(ms[0]), ms, 0, NULL, NULL);
            clErrorCheck(error);
            error = clFinish(cq);
            clErrorCheck(error);

            int num_time_steps = (int)ceil((time_end - time_start) / time_step);
            for(int temp = 0; temp < voxels_x; temp++)
                for(int temp2 = 0; temp2 < voxels_y; temp2++)
                    for(int temp3 = 0; temp3 < voxels_z; temp3++)
                        for(int temp4 = 0; temp4 < num_time_steps; temp4++)
                            m->energy_t[temp][temp2][temp3][temp4] = ms[0].energy_t[temp][temp2][temp3][temp4];

            for(int side = 0; side < 2; side++)
            {
                for(int temp1 = 0; temp1 < voxels_x; temp1++)
                {
                    for(int temp2 = 0; temp2 < voxels_y; temp2++)// this will stop working if medium is not square
                    {
                        m->surrounding_x[temp1][temp2][side] = ms[0].surrounding_x[temp1][temp2][side];
                        m->surrounding_y[temp1][temp2][side] = ms[0].surrounding_y[temp1][temp2][side];
                        m->surrounding_z[temp1][temp2][side] = ms[0].surrounding_z[temp1][temp2][side];
                    }
                }
            }


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
        case 2: // CPU
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
    end = clock();
    cout << "Simulation duration was " << (float)(end - start) / CLOCKS_PER_SEC << " seconds." << endl;
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