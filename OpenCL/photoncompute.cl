typedef struct medium_struct{
  int a;
  int b;
  int c;
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

__kernel void computePhoton(__global m_str *myStruct)
{
    local p_str photon; 
    photon.x = 9;
    int gid = get_global_id(0);
    myStruct[0].a = photon.x;
}