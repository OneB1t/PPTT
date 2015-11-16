#define PHOTON_DEATH	0.0001
const float PI = 3.14159265358979323846;
const float light_speed = 299.792458;  // mm per ns



typedef struct medium_struct{
    float time_start;
    float time_step;
    float time_end;
    float pulseDuration;

  // medium struct
  int structure[100][100][100];	//	matrix with id of every media, air = 0
  float energy[100][100][100];		//	matrix with absorbed energy
  float fluence[100][100][100];            //      matrix with photon fluence
    float ua[16];							//	absorption coef by id
    float us[16];							//  scattering coef by id
    float inv_albedo[16];					//  1 / (ua + us)
    float g[16];							//  anisotropy parameter
    float n[16];                                               //	refractive index
    float k[16];                                           // heat conduction coeficient
    float rho[16];                                         // tissue density
    float c_h[16];                                         // specific heat of tissue
	//   Timing variables
	int num_time_steps;
	float time;

}m_str;

typedef struct photon_struct
{
	float x, y, z;					// photon position
	int round_x, round_y, round_z;	//  rounded position for faster operation with matrix indices 
    float ux, uy, uz;				// photon direction
	float w;						// photon weight
	int regId, lastRegId;						// regionId position of the photon
	float step, remStep, stepToNextVoxel;			// photon generated step and remaing step
	float cosTheta, phi;
	float time_of_flight;				// in nanoseconds
	int timeId;
}p_str;

typedef struct source_struct
{
	float x, y, z;	// entering position
	float ux, uy, uz; // direction
	float release_time;
}s_str;


float timeProfile_flat(float pulse_duration)
{
	return pulse_duration;  // hacked this need random
}

__kernel void computePhoton(__global m_str *myStruct)
{
    local p_str photon; 
    local s_str source; 

    // create new photon
    photon.x = 10;
    photon.y = 10;
    photon.z = 0;
    photon.ux = 0;
    photon.uy = 0;
    photon.uz = 1;
    photon.w = 1;
    photon.round_x = floor(photon.x);
    photon.round_y = floor(photon.y);
    photon.round_z = floor(photon.z);
    photon.time_of_flight = 0;
    photon.timeId = 0;

    // fill photon variables for photon
    source.release_time = timeProfile_flat(myStruct[0].pulseDuration);
    photon.x = source.x;
    photon.y = source.y;
    photon.z = source.z;    
    photon.ux = source.ux;
    photon.uy = source.uy;
    photon.uz = source.uz;

	photon.time_of_flight = source.release_time;
    photon.regId = myStruct[0].structure[(int)floor(photon.x)][(int)floor(photon.y)][(int)floor(photon.z)];
    photon.lastRegId = photon.regId;
    photon.remStep = -0.01; // this is unimplemented

    int stop = 0;
    while(photon.w > PHOTON_DEATH) 
    {
        // this is hack to test photon movement
        photon.x = photon.x + photon.ux * photon.step;
        photon.y = photon.y + photon.uy * photon.step;
        photon.z = photon.z + photon.uz * photon.step;

        if((photon.z > 0 && photon.z < 100) || (photon.y > 0 && photon.y < 100) || (photon.x > 0 && photon.x < 100))
        {
        }
        else 
        {
            break;
        }
        
        if(photon.time_of_flight < myStruct[0].time_end)
        {
            
        }
        else
        {
            break;
        }

        photon.timeId = floor(photon.time_of_flight / myStruct[0].time_step);
        photon.lastRegId = photon.regId;
        photon.regId = myStruct[0].structure[(int)floor(photon.x)][(int)floor(photon.y)][(int)floor(photon.z)];

	    float temp = photon.w * (1 - (myStruct[0].ua[photon.regId] * photon.step) + (myStruct[0].ua[photon.regId] * myStruct[0].ua[photon.regId] * photon.step * photon.step / 2)); // Taylor expansion series of Lambert-Beer law
	    myStruct[0].energy[photon.round_x][photon.round_y][photon.round_z] += (photon.w - temp); // this is not finished
	    photon.w = temp;
        stop++;
        if(stop > 100)
            break;
        
    }
}

