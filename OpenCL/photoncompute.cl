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

float RandomNumber()
{
    return 0.8f;
}

float timeProfile_flat(float pulse_duration)
{
	return pulse_duration * RandomNumber();  // hacked this need random
}

float FindEdgeDistance(p_str *photon)
{
    float temp_x;
    float temp_y;
    float temp_z;

    //  Calculate distances
    if((*photon).ux > 0.0)
        temp_x = (((*photon).round_x + 1.0f - (*photon).x) / (*photon).ux);
    else
    {
        if((*photon).ux < 0.0)
        {
            if((*photon).round_x != (*photon).x) temp_x = (((*photon).round_x - (*photon).x) / (*photon).ux);
            else temp_x = (-1 / (*photon).ux);
        }
        else temp_x = 0.0;
    }
    if((*photon).uy > 0.0)
        temp_y = (((*photon).round_y + 1.0f - (*photon).y) / (*photon).uy);
    else
    {
        if((*photon).uy < 0.0)
        {
            if((*photon).round_y != (*photon).y) temp_y = (((*photon).round_y - (*photon).y) / (*photon).uy);
            else temp_y = (-1 / (*photon).uy);
        }
        else temp_y = 0.0;
    }
    if((*photon).uz > 0.0)
        temp_z = (((*photon).round_z + 1.0f - (*photon).z) / (*photon).uz);
    else
    {
        if((*photon).uz < 0.0)
        {
            if((*photon).round_z != (*photon).z) temp_z = (((*photon).round_z - (*photon).z) / (*photon).uz);
            else temp_z = (-1 / (*photon).uz);
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

void GenDir(float g,p_str *photon)
{
    if(g != 0)
        (*photon).cosTheta = 1.0 / (2.0 * g) * (1 + g*g - pow((1 - g*g) / (1 - g + 2 * g * RandomNumber()), 2));
    else
        (*photon).cosTheta = 2 * RandomNumber() - 1;
    (*photon).phi = 2 * PI * RandomNumber();
}

float GetTOF(__global m_str *m_str,p_str *photon, float step_size)
{
    return step_size * m_str[0].n[(*photon).regId] / (2997.92458);
}

float GenStep(float invAlbedo)
{
    float temp = -log(RandomNumber()) * invAlbedo;
    return temp;
}

void UpdateDir(__global m_str *m_str,p_str *photon)
{
	GenDir(m_str[0].g[(*photon).regId],photon);
	float temp_ux, temp_uy, temp_uz;
	float sinTheta = sin(acos((*photon).cosTheta)); //cout << "sinTheta " << sinTheta << endl;

	if (fabs((*photon).uz) < 0.99999)
	{
		float temp_sqrt = 1.0 / sqrt(1 - (*photon).uz*(*photon).uz);
		(*photon).ux = sinTheta * ((*photon).ux * (*photon).uz * cos((*photon).phi) - (*photon).uy * sin((*photon).phi)) * temp_sqrt + (*photon).ux * (*photon).cosTheta;
		(*photon).uy = sinTheta * ((*photon).uy * (*photon).uz * cos((*photon).phi) + (*photon).ux * sin((*photon).phi)) * temp_sqrt + (*photon).uy * (*photon).cosTheta;
		(*photon).uz = -sinTheta * cos((*photon).phi) / temp_sqrt + (*photon).uz * (*photon).cosTheta;
	}
	else
	{
		(*photon).ux = sinTheta * cos((*photon).phi);
		(*photon).uy = sinTheta * sin((*photon).phi);
		(*photon).uz = (*photon).uz * (*photon).cosTheta / fabs((*photon).uz);
	}
}

void RoundPosition(p_str *photon)
{
    (*photon).round_x = floor((*photon).x);
    (*photon).round_y = floor((*photon).y);
    (*photon).round_z = floor((*photon).z);
}

void Move(__global m_str *m_str,p_str *photon)
{
    float temp_step = FindEdgeDistance(photon);
    (*photon).remStep = (*photon).remStep * (m_str[0].us[(*photon).lastRegId] / m_str[0].us[(*photon).regId]);

    if(temp_step > (*photon).remStep)
    {
        (*photon).step = (*photon).remStep;
 
        (*photon).x = (*photon).x + (*photon).ux * (*photon).step;
        (*photon).y = (*photon).y + (*photon).uy * (*photon).step;
        (*photon).z = (*photon).z + (*photon).uz * (*photon).step;
        (*photon).time_of_flight += GetTOF(m_str, photon, (*photon).step);
        (*photon).remStep = GenStep(m_str[0].inv_albedo[(*photon).regId]);
        UpdateDir(m_str,photon);
        RoundPosition(photon);
    }
    else
    {
        (*photon).step = temp_step;
        //cout << "Making " << step << " step" << endl;
        (*photon).x = (*photon).x + (*photon).ux * (*photon).step;
        (*photon).y = (*photon).y + (*photon).uy * (*photon).step;
        (*photon).z = (*photon).z + (*photon).uz * (*photon).step;
        (*photon).time_of_flight += GetTOF(m_str, photon, (*photon).step);
        (*photon).remStep -= temp_step;
        if((*photon).remStep > 0)
            RoundPosition(photon);
        else
        {
            RoundPosition(photon);
            (*photon).remStep = GenStep(m_str[0].inv_albedo[(*photon).regId]);
            UpdateDir(m_str,photon);
        }
    }
}


__kernel void computePhoton(__global m_str *m_str,__global s_str *source)
{
    p_str photon; 

		// create new photon
    photon.w = 1;
    photon.timeId = 0;

    // fill photon variables for photon
    //source[0].release_time = timeProfile_flat(m_str[0].pulseDuration); send random from host
    photon.x = source[0].x;
    photon.y = source[0].y;
    photon.z = source[0].z;    
    photon.ux = source[0].ux;
    photon.uy = source[0].uy;
    photon.uz = source[0].uz;
    photon.round_x = floor(photon.x);
    photon.round_y = floor(photon.y);
    photon.round_z = floor(photon.z);

		photon.time_of_flight = source[0].release_time;
    photon.regId = m_str[0].structure[(int)floor(photon.x)][(int)floor(photon.y)][(int)floor(photon.z)];
    photon.lastRegId = photon.regId;
    photon.remStep = GenStep(m_str[0].inv_albedo[photon.regId]);

    while(photon.w > PHOTON_DEATH) 
    {

        Move(m_str , &photon);

        if(photon.z < 0.0 && photon.z > voxels_z)
            break;
        if(photon.y < 0.0 && photon.y > voxels_y)
            break;
        if(photon.x < 0.0 && photon.x > voxels_x)
            break;        


        if(photon.time_of_flight >= m_str[0].time_end)
        {
            break;
        }

        photon.timeId = floor(photon.time_of_flight / m_str[0].time_step);
        photon.lastRegId = photon.regId;
        photon.regId = m_str[0].structure[(int)floor(photon.x)][(int)floor(photon.y)][(int)floor(photon.z)];

	    float temp = photon.w * (1 - (m_str[0].ua[photon.regId] * photon.step) + (m_str[0].ua[photon.regId] * m_str[0].ua[photon.regId] * photon.step * photon.step / 2)); // Taylor expansion series of Lambert-Beer law
	    m_str[0].energy_t[photon.round_x][photon.round_y][photon.round_z][photon.timeId] += (photon.w - temp);
	    photon.w = temp;
        
    }
}

