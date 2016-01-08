#pragma once

#include "PPTT_core.h"

class OpenCL
{
public:
    OpenCL(Medium *m,Heat *h,Source *s,bool debugMode,int openCLDevice, int openCLPlatform, long numPhotons,int timeSelection);
    ~OpenCL();

    Medium *m;
    Heat *h;
    Source *s;
    bool debugMode;
    int openCLDevice;
    int openCLPlatform;
    long numPhotons;
    int timeSelection;

    cl_int error = 0;
    cl_platform_id platform = 0;
    cl_device_id device = 0;
    cl_event event = 0;
    cl_int status = 0;
    cl_uint devices = 0;
    char build_c[8192];
    size_t srcsize, streamBufferSize = 0;
    cl_context context;
    cl_command_queue cq;
    cl_program prog;
    cl_kernel kernel;
    m_str* ms = new m_str[1];
    s_str * ss = new s_str[1];
    m_str_heat* mh = new m_str_heat[1];

    void InitPhotonCompute();
    void InitHeatCompute();
    void DetectOpenCLDevices();
    void ClErrorCheck(cl_int error);
    void CopyIntoOpenCLStructuresPhoton();
    void CopyIntoOpenCLStructuresHeat();
    void CopyResultsPhoton();
    void CopyAndExecuteKernelParametersPhoton();
    void CopyAndExecuteKernelParametersHeat(int iteration);
    void ReleaseOpenCLStructures();
    void CopyResultsHeat();
};
