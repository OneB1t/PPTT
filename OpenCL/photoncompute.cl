#define PHOTON_DEATH	0.0001
static const float PI     =         3.14159265358979323846;
static const float light_speed =    299.792458;  // mm per ns
static const int units = 10;               // voxels per mm
static const int voxels_x = 100;
static const int voxels_y = 100;
static const int voxels_z = 100;
static const int max_regions =  16;


typedef struct medium_struct{
    float time_start;
    float time_step;
    float time_end;
    float pulseDuration;

  // medium struct
  int structure[voxels_x][voxels_y][voxels_z];	//	matrix with id of every media, air = 0
  float energy[voxels_x][voxels_y][voxels_z];		//	matrix with absorbed energy
  float fluence[voxels_x][voxels_y][voxels_z];            //      matrix with photon fluence
    float ua[max_regions];							//	absorption coef by id
    float us[max_regions];							//  scattering coef by id
    float inv_albedo[max_regions];					//  1 / (ua + us)
    float g[max_regions];							//  anisotropy parameter
    float n[max_regions];                                               //	refractive index
    float k[max_regions];                                           // heat conduction coeficient
    float rho[max_regions];                                         // tissue density
    float c_h[max_regions];                                         // specific heat of tissue
	//   Timing variables
	int num_time_steps;
	float time;
    float energy_t[voxels_x][voxels_y][voxels_z][6];

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

float FindEdgeDistance(p_str photon)
{
    float temp_x;
    float temp_y;
    float temp_z;

    //  Calculate distances
    if(photon.ux > 0.0)
        temp_x = ((photon.round_x + 1.0f - photon.x) / photon.ux);
    else
    {
        if(photon.ux < 0.0)
        {
            if(photon.round_x != photon.x) temp_x = ((photon.round_x - photon.x) / photon.ux);
            else temp_x = (-1 / photon.ux);
        }
        else temp_x = 0.0;
    }
    if(photon.uy > 0.0)
        temp_y = ((photon.round_y + 1.0f - photon.y) / photon.uy);
    else
    {
        if(photon.uy < 0.0)
        {
            if(photon.round_y != photon.y) temp_y = ((photon.round_y - photon.y) / photon.uy);
            else temp_y = (-1 / photon.uy);
        }
        else temp_y = 0.0;
    }
    if(photon.uz > 0.0)
        temp_z = ((photon.round_z + 1.0f - photon.z) / photon.uz);
    else
    {
        if(photon.uz < 0.0)
        {
            if(photon.round_z != photon.z) temp_z = ((photon.round_z - photon.z) / photon.uz);
            else temp_z = (-1 / photon.uz);
        }
        else temp_z = 0.0;
    }

    float temp_first = 0, temp_second = 0, temp_third = 0;

    if(temp_x == 0 && temp_y == 0 && temp_z == 0)
        return 1.0;
    else
    {
        if(temp_x < temp_y)
        {
            if(temp_z < temp_x)
            {
                temp_first = temp_z;
                temp_second = temp_x;
                temp_third = temp_y;
            }
            else
            {
                temp_first = temp_x;
                if(temp_z < temp_y)
                {
                    temp_second = temp_z;
                    temp_third = temp_y;
                }
                else
                {
                    temp_second = temp_y;
                    temp_third = temp_z;
                }
            }
        }
        else
        {
            if(temp_z < temp_y)
            {
                temp_first = temp_z;
                temp_second = temp_y;
                temp_third = temp_x;
            }
            else
            {
                temp_first = temp_y;
                if(temp_z < temp_x)
                {
                    temp_second = temp_z;
                    temp_third = temp_x;
                }
                else
                {
                    temp_second = temp_x;
                    temp_third = temp_z;
                }
            }
        }
    }
    if(temp_first != 0)
        return temp_first;
    else
    {
        if(temp_second != 0)
            return temp_second;
        else
            return temp_third;
    }
}
float RandomNumber()
{
    return 0.5f;
}

void GenDir(float g,p_str photon)
{
    if(g != 0)
        photon.cosTheta = 1.0 / (2.0 * g) * (1 + g*g - pow((1 - g*g) / (1 - g + 2 * g * RandomNumber()), 2));
    else
        photon.cosTheta = 2 * RandomNumber() - 1;
    photon.phi = 2 * PI * RandomNumber();
}

float GetTOF(__global m_str *myStruct, p_str photon, float step_size)
{
    return step_size * myStruct[0].n[photon.regId] / (light_speed * units);
}

float GenStep(float invAlbedo)
{
    float temp = -log(RandomNumber()) * invAlbedo; //  random number hack
    return temp;
}

void UpdateDir(__global m_str *myStruct, p_str photon)
{
	GenDir(myStruct[0].g[photon.regId],photon);
	float temp_ux, temp_uy, temp_uz;
	float sinTheta = sin(acos(photon.cosTheta)); //cout << "sinTheta " << sinTheta << endl;
	float temp_sqrt = 1.0 / sqrt(1 - photon.uz*photon.uz); //cout << "temp sqrt " << temp_sqrt << endl;
	if (fabs(photon.uz) < 0.99999)
	{
		temp_ux = sinTheta * (photon.ux * photon.uz * cos(photon.phi) - photon.uy * sin(photon.phi)) * temp_sqrt + photon.ux * photon.cosTheta;
		temp_uy = sinTheta * (photon.uy * photon.uz * cos(photon.phi) + photon.ux * sin(photon.phi)) * temp_sqrt + photon.uy * photon.cosTheta;
		temp_uz = -sinTheta * cos(photon.phi) / temp_sqrt + photon.uz * photon.cosTheta;
		photon.ux = temp_ux;
		photon.uy = temp_uy;
		photon.uz = temp_uz;
	}
	else
	{
		photon.ux = sinTheta * cos(photon.phi);
		photon.uy = sinTheta * sin(photon.phi);
		photon.uz = photon.uz * photon.cosTheta / fabs(photon.uz);
	}
}

void RoundPosition(p_str photon)
{
    photon.round_x = floor(photon.x);
    photon.round_y = floor(photon.y);
    photon.round_z = floor(photon.z);
}

void Move(__global m_str *myStruct, p_str photon)
{
    float temp_step = FindEdgeDistance(photon);
    photon.remStep = photon.remStep * (myStruct[0].us[photon.lastRegId] / myStruct[0].us[photon.regId]);

    if(temp_step > photon.remStep)
    {
        photon.step = photon.remStep;
 
        photon.x = photon.x + photon.ux * photon.step;
        photon.y = photon.y + photon.uy * photon.step;
        photon.z = photon.z + photon.uz * photon.step;
        photon.time_of_flight += GetTOF(myStruct, photon, photon.step);
        photon.remStep = GenStep(myStruct[0].inv_albedo[photon.regId]);
        UpdateDir(myStruct,photon);
        RoundPosition(photon);
    }
    else
    {
        photon.step = temp_step;
        //cout << "Making " << step << " step" << endl;
        photon.x = photon.x + photon.ux * photon.step;
        photon.y = photon.y + photon.uy * photon.step;
        photon.z = photon.z + photon.uz * photon.step;
        photon.time_of_flight += GetTOF(myStruct, photon, photon.step);
        photon.remStep -= temp_step;
        if(photon.remStep > 0)
            RoundPosition(photon);
        else
        {
            RoundPosition(photon);
            photon.remStep = GenStep(myStruct[0].inv_albedo[photon.regId]);
            UpdateDir(myStruct,photon);
        }
    }
    // cout << "Moving to position x: " << x << " y: " << y << " z: " << z << endl;
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

    while(photon.w > PHOTON_DEATH) 
    {

        Move(myStruct ,photon);

        if(photon.z < 0 && photon.z > voxels_z)
            break;
        if(photon.y < 0 && photon.y > voxels_y)
            break;
        if(photon.x < 0 && photon.x > voxels_x)
            break;        


        if(photon.time_of_flight >= myStruct[0].time_end)
        {
            break;
        }

        photon.timeId = floor(photon.time_of_flight / myStruct[0].time_step);
        photon.lastRegId = photon.regId;
        photon.regId = myStruct[0].structure[(int)floor(photon.x)][(int)floor(photon.y)][(int)floor(photon.z)];

	    float temp = photon.w * (1 - (myStruct[0].ua[photon.regId] * photon.step) + (myStruct[0].ua[photon.regId] * myStruct[0].ua[photon.regId] * photon.step * photon.step / 2)); // Taylor expansion series of Lambert-Beer law
	    myStruct[0].energy_t[photon.round_x][photon.round_y][photon.round_z][photon.timeId] += (photon.w - temp); // this is not finished
	    photon.w = temp;
        
    }
}

