/**********************************************************/
/*    PPTT 0.1											  */
/*    Created by Jakub Mesicek							  */
/*      10/2015											  */
/**********************************************************/
/*														  */
/* PPTT_main is a main source file calling the whole      */
/* program												  */
/**********************************************************/
#include <iomanip> 
#include <iostream>
#include <vector>
#include <string>
#include "CL\CL.h"
#include "PPTT_OpenCLinit.h"
#include "PPTT_io.h"

using namespace std;

OpenCL::OpenCL(Medium * m, Heat * h, Source * s, bool debugMode, int openCLDevice, int openCLPlatform, long numPhotons, int timeSelection)
{
    this->m = m;
    this->h = h;
    this->s = s;
    this->debugMode = debugMode;
    this->openCLDevice = openCLDevice;
    this->openCLPlatform = openCLPlatform;
    this->numPhotons = numPhotons;
    this->timeSelection = timeSelection;
}


void OpenCL::DetectOpenCLDevices()
{
    std::cout << "Detecting OpenCL platforms" << endl;
    char name[256];
    char version[256];
    cl_uint num_platforms;
    ClErrorCheck(clGetPlatformIDs(0, NULL, &num_platforms));
    vector<cl_platform_id> platforms(num_platforms);
    ClErrorCheck(clGetPlatformIDs(num_platforms, platforms.data(), NULL));
    vector<cl_uint> num_platform_devices(num_platforms);
    cl_uint num_devices = 0;
    for(cl_uint i = 0; i < num_platforms; ++i)
    {
        const auto platform = platforms[i];
        std::cout << "ID   Name" << endl;
        ClErrorCheck(clGetPlatformInfo(platform, CL_PLATFORM_NAME, sizeof(name), name, NULL));
        ClErrorCheck(clGetPlatformInfo(platform, CL_PLATFORM_VERSION, sizeof(version), version, NULL));
        std::cout << i << ' ' << name << ", " << version << endl;
        ClErrorCheck(clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 0, NULL, &num_platform_devices[i]));
        num_devices += num_platform_devices[i];
    }
    if(!debugMode)
    {
        platform = platforms[OpenCLPlatform()];
    }

    else
    {
        platform = platforms[openCLPlatform];
    }
    std::cout << "Detecting OpenCL devices" << endl;
    if(!num_devices)
    {
        cerr << "No OpenCL devices detected" << endl;
    }
    vector<cl_device_id> devices(num_devices);
    vector<bool> cl12(num_devices);
    vector<cl_bool> host_unified_memory(num_devices);
    std::cout << "ID                          Name  CL CU GMEM(MB) LMEM(KB) CMEM(KB) LMEMTYPE ECC" << endl;
    const char* local_mem_types[] = { "NONE", "LOCAL", "GLOBAL" };
    for(cl_uint i = 0, dev = 0; i < num_platforms; ++i)
    {
        const auto platform = platforms[i];
        const auto npd = num_platform_devices[i];
        ClErrorCheck(clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, npd, &devices[dev], NULL));
        for(cl_uint d = 0; d < npd; ++d, ++dev)
        {
            const auto device = devices[dev];
            cl_uint max_compute_units;
            cl_ulong global_mem_size;
            cl_ulong local_mem_size;
            cl_ulong max_constant_buffer_size;
            cl_bool error_correction_support;
            cl_device_local_mem_type local_mem_type;
            ClErrorCheck(clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(name), name, NULL));
            ClErrorCheck(clGetDeviceInfo(device, CL_DEVICE_OPENCL_C_VERSION, sizeof(version), version, NULL));
            ClErrorCheck(clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(max_compute_units), &max_compute_units, NULL));
            ClErrorCheck(clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(global_mem_size), &global_mem_size, NULL));
            ClErrorCheck(clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(local_mem_size), &local_mem_size, NULL));
            ClErrorCheck(clGetDeviceInfo(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(max_constant_buffer_size), &max_constant_buffer_size, NULL));
            ClErrorCheck(clGetDeviceInfo(device, CL_DEVICE_ERROR_CORRECTION_SUPPORT, sizeof(error_correction_support), &error_correction_support, NULL));
            ClErrorCheck(clGetDeviceInfo(device, CL_DEVICE_HOST_UNIFIED_MEMORY, sizeof(host_unified_memory[dev]), &host_unified_memory[dev], NULL));
            ClErrorCheck(clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_TYPE, sizeof(local_mem_type), &local_mem_type, NULL));
            cl12[dev] = version[9] > '1' || version[11] >= '2';
            const string name_str(name);
            const size_t name_len = name_str.size();
            cout << dev << setw(31) << (name_len <= 30 ? name_str : name_str.substr(name_str.find(' ', name_len - 31) + 1)) << ' ' << version[9] << '.' << version[11] << setw(3) << max_compute_units << setw(9) << global_mem_size / 1048576 << setw(9) << local_mem_size / 1024 << setw(9) << max_constant_buffer_size / 1024 << setw(9) << local_mem_types[local_mem_type] << setw(4) << error_correction_support << endl;
        }
    }
    if(!debugMode)
    {
        device = devices[OpenCLDevice()];
    }

    else
    {
        device = devices[openCLDevice];
    }
}
void OpenCL::InitPhotonCompute()
{
    cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform,0 };
    context = clCreateContext(properties, 1, &device, NULL, NULL, &error);
    ClErrorCheck(error);
    cq = clCreateCommandQueue(context, device, 0, &error);
    ClErrorCheck(error);

    char fileName[] = "photoncompute.cl";
    FILE *fp = fopen(fileName, "r");
    if(!fp) {
        fprintf(stderr, "Failed to load kernel.\n");
    }
    char * sourceString = (char*)malloc(MAX_SOURCE_SIZE);
    srcsize = fread(sourceString, 1, MAX_SOURCE_SIZE, fp);
    fclose(fp);

    const char *srcptr[] = { sourceString, };
    /* Submit the source code of the kernel to OpenCL, and create a program object with it */
    prog = clCreateProgramWithSource(context, 1, srcptr, &srcsize, &error);
    ClErrorCheck(error);

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
        printf("OpenCL code for Photon propagation compile sucessfull. \n");
    }
    kernel = clCreateKernel(prog, "computePhoton", &error);
    ClErrorCheck(error);
}

void OpenCL::InitHeatCompute()
{
    char fileNameH[] = "heatcompute.cl";
    FILE *fph = fopen(fileNameH, "r");
    if(!fph) {
        fprintf(stderr, "Failed to load kernel.\n");
    }
    char * sourceStringH = (char*)malloc(MAX_SOURCE_SIZE);
    srcsize = fread(sourceStringH, 1, MAX_SOURCE_SIZE, fph);
    fclose(fph);

    const char *srcptrh[] = { sourceStringH };
    /* Submit the source code of the kernel to OpenCL, and create a program object with it */
    prog = clCreateProgramWithSource(context, 1, srcptrh, &srcsize, &error);
    ClErrorCheck(error);

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
        printf("OpenCL code for Heat transfer compile sucessfull. \n");
    }
    kernel = clCreateKernel(prog, "PennesEquation", &error);
    ClErrorCheck(error);
}

void OpenCL::CopyIntoOpenCLStructuresPhoton()
{
    // copy all into openCL structures
    ss[0].releaseTime = s->releaseTime;
    ss[0].simulationType = timeSelection;
    ss[0].x = s->x;
    ss[0].y = s->y;
    ss[0].z = s->z;
    ss[0].ux = s->ux;
    ss[0].uy = s->uy;
    ss[0].uz = s->uz;
    ss[0].freq = s->freq;
    ss[0].power = s->power;
    ss[0].energy = s->energy;

    for(int temp = 0; temp < voxelsX; temp++)			// matrix with fluence inicialization
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
            {
                ms[0].energy[temp][temp2][temp3] = m->energy[temp][temp2][temp3];
                ms[0].structure[temp][temp2][temp3] = m->structure[temp][temp2][temp3];
                ms[0].fluence[temp][temp2][temp3] = m->fluence[temp][temp2][temp3];
                mh[0].structure[temp][temp2][temp3] = ms[0].structure[temp][temp2][temp3];
            }

    for(int temp = 0; temp < maxRegions; temp++)
    {
        ms[0].ua[temp] = m->ua[temp];
        ms[0].us[temp] = m->us[temp];
        ms[0].inv_albedo[temp] = m->invAlbedo[temp];
        ms[0].g[temp] = m->g[temp];
        ms[0].n[temp] = m->n[temp];
        ms[0].k[temp] = m->k[temp];
        ms[0].rho[temp] = m->rho[temp];
        ms[0].c_h[temp] = m->c_h[temp];
        ms[0].w_g[temp] = m->w_g[temp];
    }

    // this should go into another structure
    ms[0].time_end = timeEnd;
    ms[0].time_start = timeStart;
    ms[0].pulseDuration = pulseDuration;
    ms[0].time_step = timeStep;
    ms[0].finished = 0;


}

void OpenCL::CopyIntoOpenCLStructuresHeat(float power)
{

    for(int temp = 0; temp < voxelsX; temp++)			// matrix with fluence inicialization
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
            {
                mh[0].structure[temp][temp2][temp3] = m->structure[temp][temp2][temp3];
                mh[0].energy[temp][temp2][temp3] = m->energy[temp][temp2][temp3];
                mh[0].temperature[temp][temp2][temp3] = 36.5f;
                mh[0].temperature_help[temp][temp2][temp3] = 36.5f;
            }

    for(int temp = 0; temp < maxRegions; temp++)
    {
        mh[0].k[temp] = m->k[temp];
        mh[0].w_g[temp] = m->w_g[temp];
    }

    mh[0].arterial_temperature = 36.5f;
    mh[0].timeSelection = timeSelection;
    mh[0].power = power;


}

void OpenCL::CopyResultsPhoton()
{
    // energy and energy_t back from openCL
    int num_time_steps = (int)ceil((timeEnd - timeStart) / timeStep);
    for(int temp = 0; temp < voxelsX; temp++)
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
            {
                m->energy[temp][temp2][temp3] = ms[0].energy[temp][temp2][temp3];
                for(int temp4 = 0; temp4 < num_time_steps; temp4++)
                {
                    m->energy_t[temp][temp2][temp3][temp4] = ms[0].energy_t[temp][temp2][temp3][temp4];
                }

            }
    switch(timeSelection)
    {
        case STEADY_STATE:
        m->RescaleEnergy(numPhotons);
        break;
        case TIME_RESOLVED:
        m->RescaleEnergy_Time(numPhotons, timeStep);
        break;
    }



    // results for surrounding matrix
    for(int side = 0; side < 2; side++)
    {
        for(int temp1 = 0; temp1 < voxelsY; temp1++)
        {
            for(int temp2 = 0; temp2 < voxelsZ; temp2++)// this will stop working if medium is not square
            {
                m->surroundingX[temp1][temp2][side] = ms[0].surrounding_x[temp1][temp2][side];

            }
        }
        for(int temp1 = 0; temp1 < voxelsX; temp1++)
        {
            for(int temp2 = 0; temp2 < voxelsZ; temp2++)// this will stop working if medium is not square
            {
                m->surroundingY[temp1][temp2][side] = ms[0].surrounding_y[temp1][temp2][side];

            }
        }
        for(int temp1 = 0; temp1 < voxelsY; temp1++)
        {
            for(int temp2 = 0; temp2 < voxelsX; temp2++)// this will stop working if medium is not square
            {
                m->surroundingZ[temp1][temp2][side] = ms[0].surrounding_z[temp1][temp2][side];

            }
        }
    }
}

void OpenCL::CopyAndExecuteKernelParametersPhoton()
{
    clock_t simulationStart = clock();
    size_t globalWorkItems[] = { (size_t)numPhotons };

    // copy memory buffers to GPU
    cl_mem mediumMemoryBlock = clCreateBuffer(context, NULL, sizeof(ms[0]), NULL, &status);
    clEnqueueWriteBuffer(cq, mediumMemoryBlock, CL_TRUE, 0, sizeof(ms[0]), ms, 0, NULL, NULL);
    status = clSetKernelArg(kernel, 0, sizeof(ms), &mediumMemoryBlock);

    cl_mem structureMemoryBlock = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(ss[0]), NULL, &status);
    clEnqueueWriteBuffer(cq, structureMemoryBlock, CL_TRUE, 0, sizeof(ss[0]), &ss[0], 0, NULL, NULL);
    status = clSetKernelArg(kernel, 1, sizeof(ss), &structureMemoryBlock);

    int random = rand();
    cl_mem randomValue = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
    clEnqueueWriteBuffer(cq, randomValue, CL_TRUE, 0, sizeof(int), &random, 0, NULL, NULL);
    status = clSetKernelArg(kernel, 2, sizeof(int), &randomValue);

    status = clEnqueueNDRangeKernel(cq, kernel, 1, NULL, globalWorkItems, NULL, 0, NULL, NULL);
    status = clEnqueueReadBuffer(cq, mediumMemoryBlock, CL_TRUE, 0, sizeof(ms[0]), ms, 0, NULL, NULL);

    ClErrorCheck(error);
    error = clFinish(cq);
    ClErrorCheck(error);
    std::cout << "Photon Simulation duration was " << (float)(clock() - simulationStart) / CLOCKS_PER_SEC << " seconds." << endl;

}

void OpenCL::CopyAndExecuteKernelParametersHeat(int iteration)
{

    for(int temp = 0; temp < voxelsX; temp++)
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
            {
                mh[0].temperature[temp][temp2][temp3] = 36.0;
            }


    // copy memory buffers to GPU
    size_t globalWorkItemsHeat[] = { (size_t)voxelsX - 2, (size_t)voxelsY - 2,(size_t)voxelsZ - 2 };
    cl_mem heatMemoryBlock = clCreateBuffer(context, NULL, sizeof(mh[0]), NULL, &status);
    clEnqueueWriteBuffer(cq, heatMemoryBlock, CL_TRUE, 0, sizeof(mh[0]), mh, 0, NULL, NULL);
    status = clSetKernelArg(kernel, 0, sizeof(mh), &heatMemoryBlock);
    clock_t simulationStart = clock();

    for(int i = 0; i < iteration; i++)
    {
        status = clEnqueueNDRangeKernel(cq, kernel, 3, NULL, globalWorkItemsHeat, NULL, 0, NULL, NULL);
        status = clEnqueueReadBuffer(cq, heatMemoryBlock, CL_TRUE, 0, sizeof(mh[0]), mh, 0, NULL, NULL);
        for(int temp = 0; temp < voxelsX; temp++)
            for(int temp2 = 0; temp2 < voxelsY; temp2++)
                for(int temp3 = 0; temp3 < voxelsZ; temp3++)
                {
                    mh[0].temperature[temp][temp2][temp3] = mh[0].temperature_help[temp][temp2][temp3];
                }
        clEnqueueWriteBuffer(cq, heatMemoryBlock, CL_TRUE, 0, sizeof(mh[0]), mh, 0, NULL, NULL);
    }
    status = clEnqueueReadBuffer(cq, heatMemoryBlock, CL_TRUE, 0, sizeof(mh[0]), mh, 0, NULL, NULL);
    std::cout << "Heat Simulation duration was " << (float)(clock() - simulationStart) / CLOCKS_PER_SEC << " seconds." << endl;
    ClErrorCheck(error);
    error = clFinish(cq);
    ClErrorCheck(error);
}

void OpenCL::CopyResultsHeat()
{
    // results for heat transfer
    for(int temp = 0; temp < voxelsX; temp++)
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
            {
                h->temperature[temp][temp2][temp3] = mh[0].temperature[temp][temp2][temp3];
            }
}

void OpenCL::ReleaseOpenCLStructures()
{
    if(cq != 0)
        clReleaseCommandQueue(cq);
    if(kernel != 0)
        clReleaseKernel(kernel);
    if(prog != 0)
        clReleaseProgram(prog);
    if(context != 0)
        clReleaseContext(context);
}


void OpenCL::ClErrorCheck(cl_int error)
{
    if(error != CL_SUCCESS) {
        printf("\n Error number %d", error);
    }
}