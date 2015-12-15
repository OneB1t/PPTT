/**********************************************************/
/*					    PPTT 0.1						  */
/*				Created by Jakub Mesicek				  */
/*						10/2015							  */
/**********************************************************/
/*														  */
/*	PPTT_core defines basic objects such as Photon etc.   */
/**********************************************************/

/* TODO list: -withdraw Region class, replace with matrix in Space region */

#ifndef PPTT_CORE_H
#define	PPTT_CORE_H

#include <cmath>
#include "CL\CL.h"
#define PHOTON_DEATH	0.0001
#define MAX_SOURCE_SIZE (0x100000)

const float PI = 3.14159265358979323846f;
const float light_speed = 299.792458f;  // mm per ns
const int voxels_x = 150;
const int voxels_y = 150;
const int voxels_z = 150;
const int max_regions = 16;
const int units = 10;               // voxels per mm


const float time_start = 0;
const float time_step = 0.00125f;	// in ns
const float time_end = 0.12f;
const float pulseDuration = 0.0025f;
const int timeSegments = 96;

typedef struct tag_my_struct {
    float time_start;
    float time_step;
    float time_end;
    float pulseDuration;
    int finished;
    int num_time_steps;
    // medium struct
    int structure[voxels_x][voxels_y][voxels_z];	//	matrix with id of every media, air = 0
    float energy[voxels_x][voxels_y][voxels_z];		//	matrix with absorbed energy
    float fluence[voxels_x][voxels_z][voxels_z];            //      matrix with photon fluence
    float ua[max_regions];							//	absorption coef by id
    float us[max_regions];							//  scattering coef by id
    float inv_albedo[max_regions];					//  1 / (ua + us)
    float g[max_regions];							//  anisotropy parameter
    float n[max_regions];                                               //	refractive index
    float k[max_regions];                                           // heat conduction coeficient
    float rho[max_regions];                                         // tissue density
    float c_h[max_regions];                                         // specific heat of tissue
    float w_g[max_regions];
    float energy_t[voxels_x][voxels_y][voxels_z][timeSegments];
    float surrounding_x[voxels_x][voxels_y][2];
    float surrounding_y[voxels_y][voxels_z][2];
    float surrounding_z[voxels_x][voxels_z][2];
}m_str;

typedef struct source_struct
{
    float x, y, z;	// entering position
    float ux, uy, uz; // direction
    float release_time;
}s_str;

class Source;
class Photon;
class Medium;

class Source
{
public:
    Source();
    ~Source();

    enum beam_type {
        collimated_launch, isotropic_point_source, collimated_gaussian_beam, focused_gaussian_beam, circular_flat_field
    }; // beam type from Optical-Thermal Response of Laser-irradiated tissue
    float x, y, z;	// entering position
    float ux, uy, uz; // direction
    float release_time;
    float freq;			// laser repetition rate
    float power, energy;	// average power and energy in a pulse

    void Collimated_launch(float x, float y, float z, float ux, float uy, float uz);
    void Isotropic_point_source(float x, float y, float z);
    void CollimatedGaussianBeam(float x, float y, float z, float radius, float ux, float uy, float uz); // radius 1/e
    void Focused_gaussian_beam(float x, float y, float z, float radius, float focus_z, float ux, float uy, float uz);
    void Circular_flat_beam(float x, float y, float z, float radius, float ux, float uy, float uz);
    void TimeProfile_infiniteSharp();
    void TimeProfile_flat(float pulse_duration); // in ns
    void TimeProfile_gaussian(float pulse_duration);
    void TimeProfile_sech(float pulse_duration);
    void Set_RepetitionRate(float repetition_rate);
};

class Photon
{
public:
    Photon();			// photon inicialization & termination
    ~Photon();

    float x, y, z;					// photon position
    int round_x, round_y, round_z, prev_round_x, prev_round_y, prev_round_z;	//  rounded position for faster operation with matrix indices 
    float ux, uy, uz;				// photon direction
    float w;						// photon weight
    int regId, lastRegId;						// regionId position of the photon
    float step, remStep, stepToNextVoxel;			// photon generated step and remaing step
    float cosTheta, phi;
    float time_of_flight;				// in nanoseconds
    int timeId;

    void GetSourceParameters(Source * s);
    float GenStep(float invAlbedo);		// generate step size
    void UpdatePos(Medium * m);			// update position and direction
    float FindEdgeDistance();
    void Move(Medium * m);              // update position to another voxel
    float GetTOF(Medium * m, float step_size);			// calculate time of flight in medium
    int GetTimeId();					// returns Id of time gate
    int CheckTOF(float end);			// check if TOF doesn't exceed time_end
    void GetVoxId();					// get voxel ID where the photon is located
    void PosAndDir();					// print current position and direction 
    void RoundPosition();				// round position and store to int variables
    int CheckBoundaries(Medium * m);				// check whether the photon is in medium
    void GenDir(float g);				// generate cos(theta) a phi according to HG phase function
    void UpdateDir(Medium * m);			// update direction
    int CheckRefIndexMismatch(Medium * m);		// returns 0 if refractive index has changed
    float GetReflectionCoef(Medium * m);	// calculate reflection coefficient on the interface
    void Reflect(Medium * m);
    void Transmis(Medium * m);
};

class Medium
{
public:
    Medium();
    ~Medium();

    int ***structure;	//	matrix with id of every media, air = 0
    int number_of_regions;							//  number of defined regions
    float ***energy;								//	matrix with absorbed energy
    float ***fluence;            //      matrix with photon fluence
    float ***surrounding_x;
    float ***surrounding_y;
    float ***surrounding_z;
    float ua[max_regions];							//	absorption coef by id
    float us[max_regions];							//  scattering coef by id
    float inv_albedo[max_regions];					//  1 / (ua + us)
    float g[max_regions];							//  anisotropy parameter
    float n[max_regions];                                               //	refractive index
    float k[max_regions];                                           // heat conduction coeficient
    float rho[max_regions];                                         // tissue density
    float c_h[max_regions];                                         // specific heat of tissue
    float w_g[max_regions];							// blood perfusivity

    /////////////////////////////////////////////////////////
    //   Timing variables
    int num_time_steps;
    float time;
    float ****energy_t;		//  score energy in first run
    float ****energy_next;	//	score energy in second run
    /////////////////////////////////////////////////////////

    void PrintMediumProperties();					//  print properties of all regions within medium
    int RetRegId(Photon * p);						//  returns id stored in structure matrix
    void AbsorbEnergy(Photon * p);					//  absorb energy of the photon in the voxel
    void AbsorbEnergyBeer(Photon * p);			//  absorb energy exponentially
    void AbsorbEnergyBeer_Time(Photon * p);		//  absorb energy with time resolution
    void AbsorbEnergyBeer_Time_secondPulse(Photon * p);		//  absorb energy for second pulse
    void RescaleEnergy(long num_photons);                                  //  rescale absorbed energy 
    void RescaleEnergy_Time(long num_photons, float time_min_step);                                  //  rescale absorbed energy 
    void RescaleEnergy_Time_secondPulse(long num_photons, float time_min_step);                      //  rescale absorbed energy for second pulse 
    void RecordFluence();                                  //  record photon fluence in a voxel
    void CreateCube(int start_x, int start_y, int start_z, int dim_x, int dim_y, int dim_z, float ua, float us, float g, float n);	// create cube with specfied properties in mm	
    void CreateBall(int center_x, int center_y, int center_z, int radius, float ua, float us, float g, float n);
    void CreateLine(float start_x, float start_y, float start_z, float dir_x, float dir_y, float dir_z, float length, float set_ua, float set_us, float set_g, float set_n);
};

class Heat          // applied heat transfer equation
{
public:
    Heat();
    ~Heat();

    float time_start, time_step, time_end;
    float ***temperature;
    void AddThermalCoef(Medium * m, int mediumId, float specific_heat, float density, float conduction, float blood_perfusivity);

    /* Numerical methods for solving heat transfer */
    void LaplaceOperator();
    void PoissonEquation(Medium * m);

    void PennesEquation(Medium * m, float arterial_temperature);
    float AproximateBloodPerfusivity(float omega0, float omega1, float omega2, float omega3, float temperature);  // in 
};

// run one photon cycle
void RunPhoton_steady(Medium * m, Source * s);
void RunPhoton_time(Medium * m, Source * s);

// print absorbed energy
void PrintAbsEnergy(Medium * m);

// write results in file
void WriteAbsorbedEnergyToFile(Medium * m);
void WritePhotonFluenceToFile(Medium * m);
void WriteAbsorbedEnergyToFile_Time(Medium * m);
void WriteAbsorbedEnergyToFile_Time_secondPulse(Medium * m);
void CreateNewThread_time(Medium * m, Source * s, long numPhotons);
void CreateNewThread_steady(Medium * m, Source * s, long numPhotons);
void clErrorCheck(cl_int error);

// second pulse
void Prepare_SecondPulse(Medium * m, Source * s, float delay);
void CreateNewThread_secondPulse(Medium * m, Source * s, long numPhotons);
void RunPhotonNew_secondPulse(Medium * m, Source * s);

float RandomNumber();

float RandomNumber(); // from 0 to 1
inline int int_floor(float x);
#endif	/* PPTT_CORE_H */

