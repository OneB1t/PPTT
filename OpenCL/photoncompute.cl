#define PHOTON_DEATH	0.0001

static const float PI     =         3.14159265358979323846;
static const float light_speed =    299.792458;  // mm per ns
static const int units = 10;               // voxels per mm
static const int voxels_x = 150;
static const int voxels_y = 150;
static const int voxels_z = 150;
static const int timeSegment = 96;
static const int max_regions =  16;

typedef struct { uint x; uint c; } mwc64x_state_t;


typedef struct medium_struct{
    float time_start;
    float time_step;
    float time_end;
    float pulseDuration;
    int finished;
	int num_time_steps;
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
    float w_g[max_regions];
    float energy_t[voxels_x][voxels_y][voxels_z][timeSegment];


}m_str;

typedef struct photon_struct
{
    float4 position;  // this also implements photon weight
    int3 roundposition;
    int3 prevroundposition;
    float3 vector;
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

ulong AddMod64(ulong a, ulong b, ulong M)
{
	ulong v = a + b;
	if (v >= M || v < a) v -= M;
	return v;
}
ulong MulMod64(ulong a, ulong b, ulong M)
{
	ulong r = 0;
	while (a)
	{
		if (a & 1) r = AddMod64(r, b, M);
		b = AddMod64(b, b, M);
		a = a >> 1;
	}
	return r;
}
ulong PowMod64(ulong a, ulong e, ulong M)
{
	ulong sqr = a, acc = 1;
	while (e)
	{
		if (e & 1) acc = MulMod64(acc, sqr, M);
		sqr = MulMod64(sqr, sqr, M);
		e = e >> 1;
	}
	return acc;
}

enum { A = 4294883355U };
enum { M = 18446383549859758079UL };
enum { B = 4077358422479273989UL };

void skip(mwc64x_state_t *s, ulong d)
{
	ulong m = PowMod64(A, d, M);
	ulong x = MulMod64(s->x * (ulong)A + s->c, m, M);
	s->x = x / A;
	s->c = x % A;
}

void seed(mwc64x_state_t *s, ulong baseOffset, ulong perStreamOffset)
{
	ulong d = baseOffset + get_global_id(0) * perStreamOffset;
	ulong m = PowMod64(A, d, M);
	ulong x = MulMod64(B, m, M);
	s->x = x / A;
	s->c = x % A;
}

uint next(mwc64x_state_t *s)
{
	uint X = s->x;
	uint C = s->c;
	uint r = X ^ C;
	uint Xn = A * X + C;
	uint carry = Xn < C;
	uint Cn = mad_hi(A, X, carry);
	s->x = Xn;
	s->c = Cn;
	return r;
}


float RandomNumber(mwc64x_state_t *rng)
{
    return ((float)next(rng) / 0xffffffff);
}

float timeProfile_flat(float pulse_duration,mwc64x_state_t *rng)
{
	return pulse_duration * RandomNumber(rng);  // hacked this need random
}

float FindEdgeDistance(p_str *photon,mwc64x_state_t *rng)
{
		float3 temp;
	
    //  Calculate distances
    if((*photon).vector.x > 0.0)
        temp.x = (((*photon).roundposition.x + 1.0f - (*photon).position.x) / (*photon).vector.x);
    else
    {
        if((*photon).vector.x < 0.0)
        {
            if((*photon).roundposition.x != (*photon).position.x) 
							temp.x = (((*photon).roundposition.x - (*photon).position.x) / (*photon).vector.x);
            else 
							temp.x = (-1 / (*photon).vector.x);
        }
        else temp.x = 0.0;
    }
    if((*photon).vector.y > 0.0)
        temp.y = (((*photon).roundposition.y + 1.0f - (*photon).position.y) / (*photon).vector.y);
    else
    {
        if((*photon).vector.y < 0.0)
        {
            if((*photon).roundposition.y != (*photon).position.y) 
							temp.y = (((*photon).roundposition.y - (*photon).position.y) / (*photon).vector.y);
            else 
							temp.y = (-1 / (*photon).vector.y);
        }
        else temp.y = 0.0;
    }
    if((*photon).vector.z > 0.0)
        temp.z = (((*photon).roundposition.z + 1.0f - (*photon).position.z) / (*photon).vector.z);
    else
    {
        if((*photon).vector.z < 0.0)
        {
            if((*photon).roundposition.z != (*photon).position.z) 
							temp.z = (((*photon).roundposition.z - (*photon).position.z) / (*photon).vector.z);
            else 
							temp.z = (-1 / (*photon).vector.z);
        }
        else temp.z = 0.0;
    }

		float3 temp2 = (0,0,0);
		
        if(temp.x == 0 && temp.y == 0 && temp.z == 0)
        return 1.0;
    else
    {
        if(temp.x < temp.y)
        {
            if(temp.z < temp.x)
            {
                temp2.x = temp.z;
                temp2.y = temp.x;
                temp2.z = temp.y;
            }
            else
            {
                temp2.x = temp.x;
                if(temp.z < temp.y)
                {
                    temp2.y = temp.z;
                    temp2.z = temp.y;
                }
                else
                {
                    temp2.y = temp.y;
                    temp2.z = temp.z;
                }
            }
        }
        else
        {
            if(temp.z < temp.y)
            {
                temp2.x = temp.z;
                temp2.y = temp.y;
                temp2.z = temp.x;
            }
            else
            {
                temp2.x = temp.y;
                if(temp.z < temp.x)
                {
                    temp2.y = temp.z;
                    temp2.z = temp.x;
                }
                else
                {
                    temp2.y = temp.x;
                    temp2.z = temp.z;
                }
            }
        }
    }
    if(temp2.x != 0)
        return temp2.x;
    else
    {
        if(temp2.y != 0)
            return temp2.y;
        else
            return temp2.z;
    }
}

void GenDir(float g,p_str *photon,mwc64x_state_t *rng)
{
    if(g != 0)
        (*photon).cosTheta = 1.0 / (2.0 * g) * (1 + g*g - pow((1 - g*g) / (1 - g + 2 * g * RandomNumber(rng)), 2));
    else
        (*photon).cosTheta = 2 * RandomNumber(rng) - 1;
    (*photon).phi = 2 * PI * RandomNumber(rng);
}

float GetTOF(__global m_str *m_str,p_str *photon, float step_size)
{
    return step_size * m_str[0].n[(*photon).regId] / (2997.92458);
}

float GenStep(float invAlbedo, mwc64x_state_t *rng)
{
    float temp = -log(RandomNumber(rng)) * invAlbedo;
    return temp;
}

float GetReflectionCoef(__global m_str *m_str,p_str *photon)
{
	return ((m_str[0].n[(*photon).lastRegId] - m_str[0].n[(*photon).regId])*(m_str[0].n[(*photon).lastRegId] - m_str[0].n[(*photon).regId])) / ((m_str[0].n[(*photon).lastRegId] + m_str[0].n[(*photon).regId])*(m_str[0].n[(*photon).lastRegId] + m_str[0].n[(*photon).regId]));
}

void Reflect(__global m_str *m_str,p_str *photon)
{
	if ((*photon).prevroundposition.z != (*photon).roundposition.z)
		(*photon).vector.z = -(*photon).vector.z;
  if ((*photon).prevroundposition.y != (*photon).roundposition.y)
		(*photon).vector.y = -(*photon).vector.y;
	if ((*photon).prevroundposition.x != (*photon).roundposition.x)
		(*photon).vector.x = -(*photon).vector.x;
}

void Transmis(__global m_str *m_str,p_str *photon)
{
  if ((*photon).prevroundposition.z != (*photon).roundposition.z)
	{
		float alpha_i = acos(fabs((*photon).vector.z));
		float alpha_t = asin(sin(alpha_i) * m_str[0].n[(*photon).lastRegId] / m_str[0].n[(*photon).regId]);
		(*photon).vector.x *= m_str[0].n[(*photon).lastRegId] / m_str[0].n[(*photon).regId];
		(*photon).vector.y *= m_str[0].n[(*photon).lastRegId] / m_str[0].n[(*photon).regId];
		(*photon).vector.z *= sin(alpha_t) / fabs((*photon).vector.z);
	}
	if ((*photon).prevroundposition.y != (*photon).roundposition.y)
	{
		float alpha_i = acos(fabs((*photon).vector.y));
		float alpha_t = asin(sin(alpha_i) * m_str[0].n[(*photon).lastRegId] / m_str[0].n[(*photon).regId]);
		(*photon).vector.x *= m_str[0].n[(*photon).lastRegId] / m_str[0].n[(*photon).regId];
		(*photon).vector.z *= m_str[0].n[(*photon).lastRegId] / m_str[0].n[(*photon).regId];
		(*photon).vector.y *= sin(alpha_t) / fabs((*photon).vector.y);
	}
	if ((*photon).prevroundposition.x != (*photon).roundposition.x)
	{
		
		float alpha_i = acos(fabs((*photon).vector.x));
		float alpha_t = asin(sin(alpha_i) * m_str[0].n[(*photon).lastRegId] / m_str[0].n[(*photon).regId]);
		(*photon).vector.y *= m_str[0].n[(*photon).lastRegId] / m_str[0].n[(*photon).regId];
		(*photon).vector.z *= m_str[0].n[(*photon).lastRegId] / m_str[0].n[(*photon).regId];
		(*photon).vector.x *= sin(alpha_t) / fabs((*photon).vector.x);
	}
}

int CheckRefIndexMismatch(__global m_str *m_str,p_str *photon)
{
	if (m_str[0].n[(*photon).regId] - m_str[0].n[(*photon).lastRegId])
		return 0;
	else return 1;
}

void UpdateDir(__global m_str *m_str,p_str *photon,mwc64x_state_t *rng)
{
	GenDir(m_str[0].g[(*photon).regId],photon,rng);
	float sinTheta = sin(acos((*photon).cosTheta));
	float sinPhotonPhi = sin((*photon).phi);
	float cosPhotonPhi = cos((*photon).phi);
	if (fabs((*photon).vector.z) < 0.99999)
	{
		float temp = 1.0 / sqrt(1 - (*photon).vector.z*(*photon).vector.z);
		(*photon).vector.x = sinTheta * ((*photon).vector.x * (*photon).vector.z * cosPhotonPhi - (*photon).vector.y * sinPhotonPhi) * temp + (*photon).vector.x * (*photon).cosTheta;
		(*photon).vector.y = sinTheta * ((*photon).vector.y * (*photon).vector.z * cosPhotonPhi + (*photon).vector.x * sinPhotonPhi) * temp + (*photon).vector.y * (*photon).cosTheta;
		(*photon).vector.z = -sinTheta * cosPhotonPhi / temp + (*photon).vector.z * (*photon).cosTheta;
	}
	else
	{
		(*photon).vector.x = sinTheta * cosPhotonPhi;
		(*photon).vector.y = sinTheta * sinPhotonPhi;
		(*photon).vector.z = (*photon).vector.z * (*photon).cosTheta / fabs((*photon).vector.z);
	}
}

void RoundPosition(p_str *photon)
{
    (*photon).prevroundposition = (*photon).roundposition;
    (*photon).roundposition.x = (int)floor((*photon).position.x);
    (*photon).roundposition.y = (int)floor((*photon).position.y);
    (*photon).roundposition.z = (int)floor((*photon).position.z);
}

void Move(__global m_str *m_str,p_str *photon,mwc64x_state_t *rng)
{
    float temp_step = FindEdgeDistance(photon,rng);
    (*photon).remStep = (*photon).remStep * (m_str[0].us[(*photon).lastRegId] / m_str[0].us[(*photon).regId]);

    if(temp_step > (*photon).remStep)
    {
			  (*photon).step = (*photon).remStep;
			  (*photon).position.xyz = (*photon).position.xyz + (*photon).vector.xyz * (*photon).step;

        (*photon).time_of_flight += GetTOF(m_str, photon, (*photon).step);
        (*photon).remStep = GenStep(m_str[0].inv_albedo[(*photon).regId],rng);
        UpdateDir(m_str,photon,rng);
        RoundPosition(photon);
    }
    else
    {
        (*photon).step = temp_step;
			  (*photon).position.xyz = (*photon).position.xyz + (*photon).vector.xyz * (*photon).step;
        (*photon).time_of_flight += GetTOF(m_str, photon, (*photon).step);
        (*photon).remStep -= temp_step;
        if((*photon).remStep > 0)
            RoundPosition(photon);
        else
        {
            RoundPosition(photon);
            (*photon).remStep = GenStep(m_str[0].inv_albedo[(*photon).regId],rng);
            UpdateDir(m_str,photon,rng);
        }
    }
}


__kernel void computePhoton(__global m_str *m_str,__global s_str *source,int random)
{
    // random number generator
    mwc64x_state_t rng;
	seed(&rng, 0, random);
	// create new photon
    p_str photon;     
    photon.timeId = 0;

    photon.position.x = source[0].x;
    photon.position.y = source[0].y;
    photon.position.z = source[0].z;
    photon.position.w = 1;
    photon.vector.x = source[0].ux;
    photon.vector.y = source[0].uy;
    photon.vector.z = source[0].uz;
    photon.roundposition.x = floor(photon.position.x);
	photon.roundposition.y = floor(photon.position.y);
	photon.roundposition.z = floor(photon.position.z);

	photon.time_of_flight = source[0].release_time;
    photon.regId = m_str[0].structure[photon.roundposition.x][photon.roundposition.y][photon.roundposition.z];
    photon.lastRegId = photon.regId;
    photon.remStep = GenStep(m_str[0].inv_albedo[photon.regId],&rng);
    //int i= 0;
    while(photon.position.w > PHOTON_DEATH) 
    {
        if (CheckRefIndexMismatch(m_str,&photon))
        {
            Move(m_str , &photon, &rng);
        }
        else
        {
            if (RandomNumber(&rng) < GetReflectionCoef(m_str,&photon))
            {
                Reflect(m_str,&photon);
            }
            else
            {
                Transmis(m_str,&photon);
            }
            Move(m_str , &photon, &rng);
        }
        

        if(photon.position.x < 1.0 || photon.position.x > voxels_x - 2)
        {
            m_str[0].energy_t[photon.roundposition.x][photon.roundposition.y][photon.roundposition.z][photon.timeId] += photon.position.w;
            break;        
        }
        if(photon.position.y < 1.0 || photon.position.y > voxels_y - 2)
        {
            m_str[0].energy_t[photon.roundposition.x][photon.roundposition.y][photon.roundposition.z][photon.timeId] += photon.position.w;
            break;        
        }
        if(photon.position.z < 1.0 || photon.position.z > voxels_z - 2)
        {
            m_str[0].energy_t[photon.roundposition.x][photon.roundposition.y][photon.roundposition.z][photon.timeId] += photon.position.w;
            break;        
        }
        if(photon.time_of_flight >= m_str[0].time_end)
            break;


        photon.timeId = floor(photon.time_of_flight / m_str[0].time_step);
        photon.lastRegId = photon.regId;
        photon.regId = m_str[0].structure[photon.roundposition.x][photon.roundposition.y][photon.roundposition.z];

	    float temp = photon.position.w * (1 - (m_str[0].ua[photon.regId] * photon.step) + (m_str[0].ua[photon.regId] * m_str[0].ua[photon.regId] * photon.step * photon.step / 2)); // Taylor expansion series of Lambert-Beer law
	    m_str[0].energy_t[photon.roundposition.x][photon.roundposition.y][photon.roundposition.z][photon.timeId] += (photon.position.w - temp);
	    photon.position.w = temp;
        //i++;
        // failsafe
        //if(i > 10000)
        //{
        //    m_str[0].energy_t[photon.roundposition.x][photon.roundposition.y][photon.roundposition.z][photon.timeId] += photon.position.w;
        //    m_str[0].finished = m_str[0].finished + 1; // count number of time ended photons
        //   break;
        //}
        
    }

}

