/**********************************************************/
/*        PPTT 0.1					  */
/*	Created by Jakub Mesicek			  */
/*      10/2015						  */
/**********************************************************/
/*							  */
/*  PPTT_core defines basic objects such as Photon etc.   */
/**********************************************************/

#include "PPTT_core.h"
#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <sstream>

using namespace std;

/////////////////////////////
// Source class
/////////////////////////////

Source::Source()
{
    x = y = z = 0;
    ux = uy = uz = 0;
	release_time = 0.0;
	freq = 1.0;
}

Source::~Source()
{
}

void Source::Collimated_launch(float set_x, float set_y, float set_z, float set_ux, float set_uy, float set_uz) // parameters in mm
{
    float temp = sqrt(set_ux*set_ux + set_uy*set_uy + set_uz*set_uz);
    x = set_x * units;
    y = set_y * units;
    z = set_z * units;
    ux = set_ux / temp;
    uy = set_uy / temp;
    uz = set_uz / temp;    
}

void Source::Isotropic_point_source(float set_x, float set_y, float set_z)  // parameters in mm
{
    x = set_x * units;
    y = set_y * units;
    z = set_z * units;
    
    float cosTheta = 2 * ((float)rand() / RAND_MAX) - 1;
    float sinTheta = sqrt(1 - cosTheta*cosTheta);
    
    float phi = 2 * PI * ((float)rand() / RAND_MAX);
    float sinPhi, cosPhi = cos(phi);
    
    if(phi < PI)
        sinPhi = sqrt(1 - cosPhi * cosPhi);
    else 
        sinPhi = -sqrt(1 - cosPhi * cosPhi);
    
    ux = sinTheta * sinPhi;
    uy = sinTheta * cosPhi;
    uz = cosTheta;
}

void Source::Collimated_gaussian_beam(float set_x, float set_y, float set_z, float radius, float set_ux, float set_uy, float set_uz)
{
    set_x = set_x * units;
    set_y = set_y * units;
    set_z = set_z * units;    
    radius = radius * units;
    
    float temp = sqrt(set_ux*set_ux + set_uy*set_uy + set_uz*set_uz);
    ux = set_ux / temp;
    uy = set_uy / temp;
    uz = set_uz / temp;
    
    float phi = 2*PI*((float)rand() / RAND_MAX);
    float r = radius * sqrt(-log((float)rand() / RAND_MAX));
    
    x = set_x + r * cos(phi);
    y = set_y + r * sin(phi);
    z = set_z;
}

void Source::Circular_flat_beam(float set_x, float set_y, float set_z, float radius, float set_ux, float set_uy, float set_uz)
{
    set_x *= units;
    set_y *= units;
    set_z *= units;
    radius *= units;
    
    float phi = 2 * PI * ((float)rand() / RAND_MAX);
    float r = radius * sqrt((float)rand() / RAND_MAX);
    
    x = set_x + r * cos(phi);
    y = set_y + r * sin(phi);
    z = set_z;
}

void Source::TimeProfile_infiniteSharp()
{
	release_time = 0;
}

void Source::TimeProfile_flat(float pulse_duration)
{
	release_time = ((float)rand() / RAND_MAX) * pulse_duration;
}

void Source::TimeProfile_gaussian(float pulse_duration)
{
	release_time = pulse_duration * sqrt(-log(GenerateRandomNumber())) + (3 * pulse_duration);
}

void Source::TimeProfile_sech(float pulse_duration)
{
	release_time = 1 / cosh(GenerateRandomNumber());
}

void Source::Set_RepetitionRate(float repetition_rate)
{
	freq = repetition_rate;
}
/////////////////////////////
// Photon class
/////////////////////////////

Photon::Photon() 
{
	x = y = 1;
	z = 0;
	ux = uy = 0;
	uz = 1;
	w = 1;
    round_x = (int)floor(x);
    round_y = (int)floor(y);
    round_z = (int)floor(z);
	prev_round_x = round_x;
	prev_round_y = round_y;
	prev_round_z = round_z;
	time_of_flight = 0;
	timeId = 0;
}

Photon::~Photon()
{
}

void Photon::GetSourceParameters(Source * s)
{
    x = s->x;
    y = s->y;
    z = s->z;
    
	round_x = (int)floor(x);
	round_y = (int)floor(y);
	round_z = (int)floor(z);
	
    ux = s->ux;
    uy = s->uy;
    uz = s->uz;

	time_of_flight = s->release_time;
}

float Photon::GenStep(float invAlbedo)
{
	float temp = -log((float)rand() / RAND_MAX) * invAlbedo;
	return temp;
}

void Photon::UpdatePos(Medium * m)
{
	float temp = GenStep(m->inv_albedo[regId]);
	
	x = x + ux * temp;
	y = y + uy * temp;
	z = z + uz * temp;
}

float Photon::FindEdgeDistance()
{
    float temp_x, temp_y, temp_z;
    //  Calculate distances
    if (ux > 0.0)
            temp_x = ((round_x + 1 - x) / ux);
    else 
    {
        if (ux < 0.0)
        {
            if(round_x != x) temp_x = ((round_x - x) / ux);
            else temp_x = (-1 / ux);
        }
        else temp_x = 0.0;
    }
    if (uy > 0.0)
            temp_y = ((round_y + 1 - y) / uy);
    else 
    {
        if (uy < 0.0)
        {
            if(round_y != y) temp_y = ((round_y - y) / uy);
            else temp_y = (-1 / uy);
        }
        else temp_y = 0.0;
    }
    if (uz > 0.0)
            temp_z = ((round_z + 1 - z) / uz);
    else 
    {
        if (uz < 0.0)
        {
            if(round_z != z) temp_z = ((round_z - z)/ uz);
            else temp_z = (-1 / uz);
        }
        else temp_z = 0.0;
    }
       
    float temp_first = 0, temp_second = 0, temp_third = 0;
    
    if (temp_x == 0 && temp_y == 0 && temp_z == 0)
        return 1.0;
    else
    {
        if (temp_x < temp_y)
        {
            if (temp_z < temp_x)
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
            if (temp_z < temp_y)
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

void Photon::Move(Medium * m)
{
    float temp_step = FindEdgeDistance();
    remStep *= (m->us[lastRegId] / m->us[regId]);
     
    if (temp_step > remStep)
    {
        step = remStep;
        x = x + ux * step;
		y = y + uy * step;
		z = z + uz * step;
		time_of_flight += GetTOF(m, step);
        remStep = GenStep(m->inv_albedo[regId]);
        UpdateDir(m);
        RoundPosition();
    }
    else
    {
        step = temp_step;
        x = x + ux * step;
		y = y + uy * step;
		z = z + uz * step;
		time_of_flight += GetTOF(m, step);
        remStep -= temp_step;
        if(remStep > 0)
            RoundPosition();
        else
        {
            RoundPosition();
            remStep = GenStep(m->inv_albedo[regId]);
            UpdateDir(m);
        }
    }
}

float Photon::GetTOF(Medium * m, float step_size)
{
	return step_size * m->n[regId] / (light_speed * units);
}

int Photon::GetTimeId()
{
	return (int)floor(time_of_flight / time_step);
}

int Photon::CheckTOF(float end)
{
	if (time_of_flight < end) return 0;
	else return 1;
}
void Photon::PosAndDir() 
{
	cout << "Current position is (" << x << " " << y << " "<< z << ")" << endl;	
	cout << "Propagation direction is (" << ux << " " << uy << " " << uz << ")" << endl;
}

void Photon::RoundPosition()
{
	prev_round_x = round_x;
	prev_round_y = round_y;
	prev_round_z = round_z;
	
	round_x = (int)floor(x);
	round_y = (int)floor(y);
	round_z = (int)floor(z);
}

int Photon::CheckBoundaries(Medium * m)
{
	if (round_z < 0)
	{
		if (round_x < 0)
			round_x = 0;
		if (round_y < 0)
			round_y = 0;

		if (round_x >= voxels_x)
			round_x = voxels_x - 1;
		if (round_y <= voxels_y)
			round_y = voxels_y - 1;

		m->surrounding_z[round_x][round_y][0] += w;
		return 1;
	}
	if (round_z >= voxels_z)
	{
		if (round_x < 0)
			round_x = 0;
		if (round_y < 0)
			round_y = 0;

		if (round_x >= voxels_x)
			round_x = voxels_x - 1;
		if (round_y <= voxels_y)
			round_y = voxels_y - 1;

		m->surrounding_z[round_x][round_y][1] += w;
		return 1;
	}
	if (round_y < 0)
	{
		if (round_x < 0)
			round_x = 0;
		if (round_z < 0)
			round_z = 0;

		if (round_x >= voxels_x)
			round_x = voxels_x - 1;
		if (round_z <= voxels_z)
			round_z = voxels_z - 1;

		m->surrounding_y[round_x][round_z][0] += w;
		return 1;
	}
	if (round_y >= voxels_y)
	{
		if (round_x < 0)
			round_x = 0;
		if (round_z < 0)
			round_z = 0;

		if (round_x >= voxels_x)
			round_x = voxels_x - 1;
		if (round_z <= voxels_z)
			round_z = voxels_z - 1;

		m->surrounding_y[round_x][round_z][1] += w;
		return 1;
	}
	if (round_x < 0)
	{
		if (round_y < 0)
			round_y = 0;
		if (round_z < 0)
			round_z = 0;

		if (round_y >= voxels_y)
			round_y = voxels_y - 1;
		if (round_z <= voxels_z)
			round_z = voxels_z - 1;

		m->surrounding_x[round_y][round_z][0] += w;
		return 1;
	}
	if (round_x >= voxels_x)
	{
		if (round_y < 0)
			round_y = 0;
		if (round_z < 0)
			round_z = 0;

		if (round_y >= voxels_y)
			round_y = voxels_y - 1;
		if (round_z <= voxels_z)
			round_z = voxels_z - 1;

		m->surrounding_x[round_y][round_z][1] += w;
		return 1;
	}
	return 0;
}
void Photon::GenDir(float g) // generate cos(theta) and phi based on anisotropy parameter g
{
	if (g != 0)
		cosTheta = 1.0 / (2.0 * g) * (1 + g*g - pow((1 - g*g) / (1 - g + 2 * g * ((float)rand() / RAND_MAX)), 2));
	else
		cosTheta = 2 * ((float)rand() / RAND_MAX) - 1;
	phi = 2 * PI * ((float)rand() / RAND_MAX);
}

void Photon::UpdateDir(Medium * m)
{
	GenDir(m->g[regId]);
	float temp_ux, temp_uy, temp_uz;
	float sinTheta = sin(acos(cosTheta)); 
	float temp_sqrt = 1.0 / sqrt(1 - uz*uz); 
	if (abs(uz) < 0.999)
	{
		temp_ux = sinTheta * (ux * uz * cos(phi) - uy * sin(phi)) * temp_sqrt + ux * cosTheta;
		temp_uy = sinTheta * (uy * uz * cos(phi) + ux * sin(phi)) * temp_sqrt + uy * cosTheta;
		temp_uz = -sinTheta * cos(phi) / temp_sqrt + uz * cosTheta;
		ux = temp_ux;
		uy = temp_uy;
		uz = temp_uz;
	}
	else
	{
		ux = sinTheta * cos(phi);
		uy = sinTheta * sin(phi);
		uz = uz * cosTheta / abs(uz);
	}
}

int Photon::CheckRefIndexMismatch(Medium * m)
{
	if (m->n[regId] - m->n[lastRegId])
		return 0;
	else return 1;
}

float Photon::GetReflectionCoef(Medium * m)
{
	return ((m->n[lastRegId] - m->n[regId])*(m->n[lastRegId] - m->n[regId])) / ((m->n[lastRegId] + m->n[regId])*(m->n[lastRegId] + m->n[regId]));
}

void Photon::Reflect(Medium * m)
{
	if (prev_round_z != round_z)
		uz = -uz;
	if (prev_round_y != round_y)
		uy = -uy;
	if (prev_round_x != round_x)
		ux = -ux;
}

void Photon::Transmis(Medium * m)
{
	if (prev_round_z != round_z)
	{
		float alpha_i = acos(abs(uz));
		float alpha_t = asin(sin(alpha_i) * m->n[lastRegId] / m->n[regId]);
		ux *= m->n[lastRegId] / m->n[regId];
		uy *= m->n[lastRegId] / m->n[regId];
		uz *= sin(alpha_t) / abs(uz);
	}
	if (prev_round_y != round_y)
	{
		float alpha_i = acos(abs(uy));
		float alpha_t = asin(sin(alpha_i) * m->n[lastRegId] / m->n[regId]);
		ux *= m->n[lastRegId] / m->n[regId];
		uz *= m->n[lastRegId] / m->n[regId];
		uy *= sin(alpha_t) / abs(uy);
	}
	if (prev_round_x != round_x)
	{
		float alpha_i = acos(abs(ux));
		float alpha_t = asin(sin(alpha_i) * m->n[lastRegId] / m->n[regId]);
		uz *= m->n[lastRegId] / m->n[regId];
		uy *= m->n[lastRegId] / m->n[regId];
		ux *= sin(alpha_t) / abs(ux);
	}
}

//////////////////////////////////////////////////
//			Medium class
//////////////////////////////////////////////////
Medium::Medium()
{
	number_of_regions = 1;								// number of regions inicialization
	
	structure = new int**[voxels_x];						// dynamic allocation 
	for (int temp = 0; temp < voxels_x; temp++)
		structure[temp] = new int*[voxels_y];
	for (int temp1 = 0; temp1 < voxels_x; temp1++)
		for (int temp2 = 0; temp2 < voxels_x; temp2++)
			structure[temp1][temp2] = new int[voxels_z];
	
	energy = new float**[voxels_x];						// dynamic allocation 
	for (int temp = 0; temp < voxels_x; temp++)
		energy[temp] = new float*[voxels_y];
	for (int temp1 = 0; temp1 < voxels_x; temp1++)
		for (int temp2 = 0; temp2 < voxels_x; temp2++)
			energy[temp1][temp2] = new float[voxels_z];

	fluence = new float**[voxels_x];						// dynamic allocation 
	for (int temp = 0; temp < voxels_x; temp++)
		fluence[temp] = new float*[voxels_y];
	for (int temp1 = 0; temp1 < voxels_x; temp1++)
		for (int temp2 = 0; temp2 < voxels_x; temp2++)
			fluence[temp1][temp2] = new float[voxels_z];

	surrounding_x = new float**[voxels_y];					// allocation of surrounding matrix 
	for (int temp = 0; temp < voxels_y; temp++)				// this one for front and rear x plane
		surrounding_x[temp] = new float*[voxels_z];
	for (int temp1 = 0; temp1 < voxels_y; temp1++)
		for (int temp2 = 0; temp2 < voxels_z; temp2++)
			surrounding_x[temp1][temp2] = new float[2];
	for (int temp = 0; temp < voxels_y; temp++)			
		for (int temp2 = 0; temp2 < voxels_z; temp2++)
			for (int temp3 = 0; temp3 < 2; temp3++)
				surrounding_x[temp][temp2][temp3] = 0;

	surrounding_y = new float**[voxels_x];					// allocation of surrounding matrix 
	for (int temp = 0; temp < voxels_x; temp++)				// this one for front and rear y plane
		surrounding_y[temp] = new float*[voxels_z];
	for (int temp1 = 0; temp1 < voxels_x; temp1++)
		for (int temp2 = 0; temp2 < voxels_z; temp2++)
			surrounding_y[temp1][temp2] = new float[2];
	for (int temp = 0; temp < voxels_x; temp++)
		for (int temp2 = 0; temp2 < voxels_z; temp2++)
			for (int temp3 = 0; temp3 < 2; temp3++)
				surrounding_y[temp][temp2][temp3] = 0;

	surrounding_z = new float**[voxels_x];					// allocation of surrounding matrix 
	for (int temp = 0; temp < voxels_x; temp++)				// this one for front and rear z plane
		surrounding_z[temp] = new float*[voxels_y];
	for (int temp1 = 0; temp1 < voxels_x; temp1++)
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			surrounding_z[temp1][temp2] = new float[2];
	for (int temp = 0; temp < voxels_x; temp++)
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < 2; temp3++)
				surrounding_z[temp][temp2][temp3] = 0;

	for (int temp = 0; temp < voxels_x; temp++)			// structure Ids inicialization
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
				structure[temp][temp2][temp3] = 0;
	for (int temp = 0; temp < voxels_x; temp++)			// matrix with absorbed energy inicialization
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
				energy[temp][temp2][temp3] = 0;
    for (int temp = 0; temp < voxels_x; temp++)			// matrix with photon fluence inicialization
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
				fluence[temp][temp2][temp3] = 0;
	for (int temp = 0; temp < max_regions; temp++)
		ua[temp] = 0;
	for (int temp = 0; temp < max_regions; temp++)
		us[temp] = 0;
	for (int temp = 0; temp < max_regions; temp++)
		inv_albedo[temp] = 0;
	for (int temp = 0; temp < max_regions; temp++)
		g[temp] = 1;
	for (int temp = 0; temp < max_regions; temp++)
		n[temp] = 1;
    for (int temp = 0; temp < max_regions; temp++)
		k[temp] = 0;
    for (int temp = 0; temp < max_regions; temp++)
		rho[temp] = 0;
    for (int temp = 0; temp < max_regions; temp++)
		c_h[temp] = 0;
	for (int temp = 0; temp < max_regions; temp++)
		w_g[temp] = 0;

	// Preparing array for time-resolved simulations
	num_time_steps = (int)ceil((time_end - time_start) / time_step);

	energy_t = new float***[voxels_x];						// dynamic allocation 
	for (int temp = 0; temp < voxels_x; temp++)
		energy_t[temp] = new float**[voxels_y];
	for (int temp1 = 0; temp1 < voxels_x; temp1++)
		for (int temp2 = 0; temp2 < voxels_x; temp2++)
			energy_t[temp1][temp2] = new float*[voxels_z];
	for (int temp = 0; temp < voxels_x; temp++)
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
				energy_t[temp][temp2][temp3] = new float[num_time_steps];
    
	for (int temp = 0; temp < voxels_x; temp++)
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
                            for (int temp4 = 0; temp4 < num_time_steps; temp4++)
                                energy_t[temp][temp2][temp3][temp4] = 0;
	
	energy_next = new float***[voxels_x];						// dynamic allocation 
	for (int temp = 0; temp < voxels_x; temp++)
		energy_next[temp] = new float**[voxels_y];
	for (int temp1 = 0; temp1 < voxels_x; temp1++)
		for (int temp2 = 0; temp2 < voxels_x; temp2++)
			energy_next[temp1][temp2] = new float*[voxels_z];
	for (int temp = 0; temp < voxels_x; temp++)
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
				energy_next[temp][temp2][temp3] = new float[num_time_steps];
	
	for (int temp = 0; temp < voxels_x; temp++)
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
				for (int temp4 = 0; temp4 < num_time_steps; temp4++)
					energy_next[temp][temp2][temp3][temp4] = 0;
}

Medium::~Medium()
{
}

void Medium::PrintMediumProperties()
{
	cout << "Structure: " << endl;
	for (int temp = 0; temp < voxels_x; temp++)
	{
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
		{
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
				cout << "  " << structure[temp][temp2][temp3];
			cout << endl;
		}
		cout << endl;
	}
	cout << endl;
	cout << "Absorption coefficients: " << endl;
	for (int temp = 0; temp < max_regions; temp++)
		cout << "  " << ua[temp];
	cout << endl;
	cout << "Scattering coefficients: " << endl;
	for (int temp = 0; temp < max_regions; temp++)
		cout << "  " << us[temp];
	cout << endl;	
	cout << "Inverse Albedo values: " << endl;
	for (int temp = 0; temp < max_regions; temp++)
		cout << "  " << inv_albedo[temp];
	cout << endl;
	cout << "Anisotropy parameter: " << endl;
	for (int temp = 0; temp < max_regions; temp++)
		cout << "  " << g[temp];
	cout << endl;
	cout << "Refractive indices: " << endl;
	for (int temp = 0; temp < max_regions; temp++)
		cout << "  " << n[temp];
	cout << endl;
}

int Medium::RetRegId(Photon * p)
{
	return structure[p->round_x][p->round_y][p->round_z];
}

void Medium::CreateCube(int start_x, int start_y, int start_z, int dim_x, int dim_y, int dim_z, float set_ua, float set_us, float set_g, float set_n)
{
    start_x *= units;
    start_y *= units;
    start_z *= units;
    dim_x *= units;
    dim_y *= units;
    dim_z *= units;
    set_ua /= (float)units;
    set_us /= (float)units;
    
	for (int temp1 = 0; temp1 < dim_x; temp1++)
		for (int temp2 = 0; temp2 < dim_y; temp2++)
			for (int temp3 = 0; temp3 < dim_z; temp3++)
				structure[start_x + temp1][start_y + temp2][start_z + temp3] = number_of_regions;
	ua[number_of_regions] = set_ua;
	us[number_of_regions] = set_us;
	inv_albedo[number_of_regions] = 1 / (set_ua + set_us);
	g[number_of_regions] = set_g;
	n[number_of_regions] = set_n;
	number_of_regions++;
}

void Medium::CreateBall(int center_x, int center_y, int center_z, int radius, float set_ua, float set_us, float set_g, float set_n)
{
    center_x *= units;
    center_y *= units;
    center_z *= units;
    radius *= units;
    set_ua /= (float)units;
    set_us /= (float)units;
    
    for(int temp1 = -radius; temp1 < radius; temp1++)
        for(int temp2 = -radius; temp2 < radius; temp2++)
            for(int temp3 = -radius; temp3 < radius; temp3++)
                if((temp1*temp1 + temp2*temp2 + temp3*temp3) < radius*radius)
                    structure[center_x + temp1][center_y + temp2][center_z + temp3] = number_of_regions;
    ua[number_of_regions] = set_ua;
	us[number_of_regions] = set_us;
	inv_albedo[number_of_regions] = 1 / (set_ua + set_us);
	g[number_of_regions] = set_g;
	n[number_of_regions] = set_n;
	number_of_regions++;        
}

void Medium::CreateLine(float start_x, float start_y, float start_z, float dir_x, float dir_y, float dir_z, float length, float set_ua, float set_us, float set_g, float set_n)
{
	float temp_norm = sqrt(dir_x*dir_x + dir_y*dir_y + dir_z*dir_z);
	start_x *= units;
	start_y *= units;
	start_z *= units;
	length *= units;
	dir_x /= temp_norm;
	dir_y /= temp_norm;
	dir_z /= temp_norm;
	set_ua /= (float)units;
	set_us /= (float)units;
	
	for (int temp = 0; temp < (int)round(length); temp++)
		structure[(int)round(start_x + dir_x*temp)][(int)round(start_y + dir_y*temp)][(int)round(start_z + dir_z*temp)] = number_of_regions;

	ua[number_of_regions] = set_ua;
	us[number_of_regions] = set_us;
	inv_albedo[number_of_regions] = 1 / (set_ua + set_us);
	g[number_of_regions] = set_g;
	n[number_of_regions] = set_n;
	number_of_regions++;
}

void Medium::AbsorbEnergy(Photon * p)
{
    float temp = p->w * (ua[p->regId] * inv_albedo[p->regId]);
    energy[p->prev_round_x][p->prev_round_y][p->prev_round_z] += temp;
    p->w -= temp;
}

void Medium::AbsorbEnergyBeer(Photon * p)
{
	float temp = p->w * exp(-ua[p->regId] * p->step);							// Lambert-Beer law
    energy[p->prev_round_x][p->prev_round_y][p->prev_round_z] += (p->w - temp);
    p->w = temp;
}

void Medium::AbsorbEnergyBeer_Time(Photon * p)
{
	float temp = p->w * exp(-ua[p->regId] * p->step);							// Lambert-Beer law
	energy_t[p->prev_round_x][p->prev_round_y][p->prev_round_z][p->timeId] += (p->w - temp);
	p->w = temp;
}

void Medium::AbsorbEnergyBeer_Time_secondPulse(Photon * p)
{
	float temp = p->w * exp(-ua[p->regId] * p->step);							// Lambert-Beer law
	energy_next[p->prev_round_x][p->prev_round_y][p->prev_round_z][p->timeId] += (p->w - temp);
	p->w = temp;
}
void Medium::RescaleEnergy(long num_photons)
{
    for (int temp = 0; temp < voxels_x; temp++)			// structure Ids inicialization
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
                            energy[temp][temp2][temp3] /= (num_photons / pow(units,3));
}

void Medium::RescaleEnergy_Time(long num_photons, float time_min_step)
{
	for (int temp = 0; temp < voxels_x; temp++)			// structure Ids inicialization
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
				for (int temp4 = 0; temp4 < num_time_steps; temp4++)
				energy_t[temp][temp2][temp3][temp4] /= (time_min_step * num_photons / pow(units, 3));
}

void Medium::RescaleEnergy_Time_secondPulse(long num_photons, float time_min_step)
{
	for (int temp = 0; temp < voxels_x; temp++)			// structure Ids inicialization
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
				for (int temp4 = 0; temp4 < num_time_steps; temp4++)
					energy_next[temp][temp2][temp3][temp4] /= (time_min_step * num_photons / pow(units, 3));
}

void Medium::RecordFluence()
{
    for (int temp = 0; temp < voxels_x; temp++)			// structure Ids inicialization
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
				fluence[temp][temp2][temp3] = energy[temp][temp2][temp3] / ua[structure[temp][temp2][temp3]];
}

//////////////////////////////////////////////////////
//          Heat class
//////////////////////////////////////////////////////

Heat::Heat()
{
	temperature = new float**[voxels_x];						// dynamic allocation 
	for (int temp = 0; temp < voxels_x; temp++)
		temperature[temp] = new float*[voxels_y];
	for (int temp1 = 0; temp1 < voxels_x; temp1++)
		for (int temp2 = 0; temp2 < voxels_x; temp2++)
			temperature[temp1][temp2] = new float[voxels_z];
	
	for (int temp = 0; temp < voxels_x; temp++)					// set human body temperature
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
					temperature[temp][temp2][temp3]= 36.5;		
}

Heat::~Heat()
{
}

void Heat::AddThermalCoef(Medium * m, int mediumId, float specific_heat, float density, float conduction, float blood_perfusivity) // in g, mm, K
{
    m->c_h[mediumId] = specific_heat;
    m->rho[mediumId] = density;
    m->k[mediumId] = conduction;
	m->w_g[mediumId] = blood_perfusivity;
}

void Heat::LaplaceOperator()
{
	float ***temperature_help = new float**[voxels_x];						// dynamic allocation 
	for (int temp = 0; temp < voxels_x; temp++)
		temperature_help[temp] = new float*[voxels_y];
	for (int temp1 = 0; temp1 < voxels_x; temp1++)
		for (int temp2 = 0; temp2 < voxels_x; temp2++)
			temperature_help[temp1][temp2] = new float[voxels_z];

	for (int temp = 0; temp < voxels_x; temp++)					// set human body temperature
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
				temperature_help[temp][temp2][temp3] = 36.5;

	for (int temp = 1; temp < voxels_x - 1; temp++)					// calculate operator
		for (int temp2 = 1; temp2 < voxels_y - 1; temp2++)
			for (int temp3 = 1; temp3 < voxels_z - 1; temp3++)
				temperature_help[temp][temp2][temp3] = (temperature[temp + 1][temp2][temp3] + temperature[temp - 1][temp2][temp3] + temperature[temp][temp2 + 1][temp3] + temperature[temp][temp2 - 1][temp3] + temperature[temp][temp2][temp3 + 1] + temperature[temp][temp2][temp3 - 1]) / 6.0;

	for (int temp = 1; temp < voxels_x - 1; temp++)					// give back into main matrix
		for (int temp2 = 1; temp2 < voxels_y - 1; temp2++)
			for (int temp3 = 1; temp3 < voxels_z - 1; temp3++)
				temperature[temp][temp2][temp3] = temperature_help[temp][temp2][temp3];
}

void Heat::PoissonEquation(Medium * m)
{
	float ***temperature_help = new float**[voxels_x];						// dynamic allocation 
	for (int temp = 0; temp < voxels_x; temp++)
		temperature_help[temp] = new float*[voxels_y];
	for (int temp1 = 0; temp1 < voxels_x; temp1++)
		for (int temp2 = 0; temp2 < voxels_x; temp2++)
			temperature_help[temp1][temp2] = new float[voxels_z];

	for (int temp = 0; temp < voxels_x; temp++)					// set human body temperature
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
				temperature_help[temp][temp2][temp3] = 36.5;

	for (int temp = 1; temp < voxels_x - 1; temp++)					// calculate operator
		for (int temp2 = 1; temp2 < voxels_y - 1; temp2++)
			for (int temp3 = 1; temp3 < voxels_z - 1; temp3++)
			{
				int tempId = m->structure[temp][temp2][temp3];
				temperature_help[temp][temp2][temp3] = (temperature[temp + 1][temp2][temp3] + temperature[temp - 1][temp2][temp3] + temperature[temp][temp2 + 1][temp3] + temperature[temp][temp2 - 1][temp3] + temperature[temp][temp2][temp3 + 1] + temperature[temp][temp2][temp3 - 1] + (m->energy[temp][temp2][temp3] / m->k[tempId])) / 6.0;
			}

	for (int temp = 1; temp < voxels_x - 1; temp++)					// give back into main matrix
		for (int temp2 = 1; temp2 < voxels_y - 1; temp2++)
			for (int temp3 = 1; temp3 < voxels_z - 1; temp3++)
				temperature[temp][temp2][temp3] = temperature_help[temp][temp2][temp3];

}

void Heat::PennesEquation(Medium * m, float arterial_temperature)
{
	float ***temperature_help = new float**[voxels_x];						// dynamic allocation 
	for (int temp = 0; temp < voxels_x; temp++)
		temperature_help[temp] = new float*[voxels_y];
	for (int temp1 = 0; temp1 < voxels_x; temp1++)
		for (int temp2 = 0; temp2 < voxels_x; temp2++)
			temperature_help[temp1][temp2] = new float[voxels_z];

	for (int temp = 0; temp < voxels_x; temp++)					// set human body temperature
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
				temperature_help[temp][temp2][temp3] = 36.5;

	for (int temp = 1; temp < voxels_x - 1; temp++)					// calculate operator
		for (int temp2 = 1; temp2 < voxels_y - 1; temp2++)
			for (int temp3 = 1; temp3 < voxels_z - 1; temp3++)
			{
				int tempId = m->structure[temp][temp2][temp3];
				temperature_help[temp][temp2][temp3] = (temperature[temp + 1][temp2][temp3] + temperature[temp - 1][temp2][temp3] + temperature[temp][temp2 + 1][temp3] + temperature[temp][temp2 - 1][temp3] + temperature[temp][temp2][temp3 + 1] + temperature[temp][temp2][temp3 - 1] + (m->energy[temp][temp2][temp3] / m->k[tempId]) - AproximateBloodPerfusivity(0.255, -0.137, 2.3589, 315, temperature[temp][temp2][temp3] + 273) * m->c_h[tempId] * (temperature[temp][temp2][temp3] - arterial_temperature) / m->k[tempId]) / 6.0;
			}

	for (int temp = 1; temp < voxels_x - 1; temp++)					// give back into main matrix
		for (int temp2 = 1; temp2 < voxels_y - 1; temp2++)
			for (int temp3 = 1; temp3 < voxels_z - 1; temp3++)
				temperature[temp][temp2][temp3] = temperature_help[temp][temp2][temp3];
}

float Heat::AproximateBloodPerfusivity(float omega0, float omega1, float omega2, float omega3, float temperature)
{
	return (omega0 + omega1 * (atan(omega2*(temperature - omega3)))) / 1e6;		// conversion to milimeters and grams
}

//////////////////////////////////////////////////////
//				Nonclass functions
//////////////////////////////////////////////////////
void RunPhoton_steady(Medium * m, Source * s)
{
	Photon * p = new Photon;
	s->Collimated_gaussian_beam(5.0, 5.0, 0.0, 0.5, 0.0, 0.0, 1.0);
	s->TimeProfile_infiniteSharp();
	p->GetSourceParameters(s);
	p->regId = m->RetRegId(p);
	p->lastRegId = p->regId;
	p->remStep = p->GenStep(m->inv_albedo[p->regId]);


	while (p->w > PHOTON_DEATH)
	{
		if (p->CheckRefIndexMismatch(m))
			p->Move(m);
		else
		{
			if (GenerateRandomNumber() < p->GetReflectionCoef(m))
				p->Reflect(m);
			else
				p->Transmis(m);
			p->Move(m);
		}
		if (p->CheckBoundaries(m)) break;
		p->lastRegId = p->regId;
		p->regId = m->RetRegId(p);
		m->AbsorbEnergyBeer(p);
	}

	delete p;
}

void RunPhoton_time(Medium * m, Source * s)
{
    Photon * p = new Photon;
	s->Collimated_gaussian_beam(5.0, 5.0, 0.0, 0.5, 0.0, 0.0, 1.0);
	s->TimeProfile_flat(pulseDuration);
    p->GetSourceParameters(s);
    p->regId = m->RetRegId(p);
    p->lastRegId = p->regId;
    p->remStep = p->GenStep(m->inv_albedo[p->regId]);
        
        
	while(p->w > PHOTON_DEATH) 
	{
		if (p->CheckRefIndexMismatch(m))
			p->Move(m);
		else
		{
			if (GenerateRandomNumber() < p->GetReflectionCoef(m))
				p->Reflect(m);
			else
				p->Transmis(m);
			p->Move(m);
		}
        if (p->CheckBoundaries(m)) break;
		if (p->CheckTOF(time_end)) break;
		p->timeId = p->GetTimeId();
		p->lastRegId = p->regId;
		p->regId = m->RetRegId(p);
		m->AbsorbEnergyBeer_Time(p);		
	}
        
        delete p;
}

void RunPhotonNew_secondPulse(Medium * m, Source * s)
{
	Photon * p = new Photon;
	s->Collimated_gaussian_beam(5.0, 8.0, 0.0, 0.5, 0.0, 0.0, 1.0);
	s->TimeProfile_flat(pulseDuration);
	p->GetSourceParameters(s);
	p->regId = m->RetRegId(p);
	p->lastRegId = p->regId;
	p->remStep = p->GenStep(m->inv_albedo[p->regId]);

	while (p->w > PHOTON_DEATH)
	{
		p->Move(m);
		if (p->CheckBoundaries(m)) break;
		if (p->CheckTOF(time_end)) break;
		p->timeId = p->GetTimeId();
		p->lastRegId = p->regId;
		p->regId = m->RetRegId(p);
		m->AbsorbEnergyBeer_Time_secondPulse(p);
	}

	delete p;
}

void PrintAbsEnergy(Medium * m)
{
	for(int temp1 = 0; temp1 < voxels_x; temp1++)
	{
		for(int temp2 = 0; temp2 < voxels_y; temp2++)
		{
			cout << endl;
			for(int temp3 = 0; temp3 < voxels_z; temp3++)
				cout << m->energy[temp1][temp2][temp3] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void WriteAbsorbedEnergyToFile(Medium * m)
{
	ofstream myFile;
	myFile.open("Results_absorbedEnergy.txt");
	for(int temp1 = 0; temp1 < voxels_x; temp1++)
	{
		for(int temp2 = 0; temp2 < voxels_y; temp2++)
			{
				for(int temp3 = 0; temp3 < voxels_z; temp3++)
					myFile << m->energy[temp1][temp2][temp3] << " ";
				myFile << endl;
		}
		myFile << endl;
	}
	myFile.close();
}

void WriteAbsorbedEnergyToFile_Time(Medium * m)
{
	for (int temp4 = 0; temp4 < m->num_time_steps; temp4++)
	{
		stringstream ss;
		ss << "ResultsAbsorbedEnergyTime_" << temp4 << ".txt";
		ofstream myFile(ss.str().c_str());
		if(myFile.is_open())
			{
                for (int temp1 = 0; temp1 < voxels_x; temp1++)
				{
					for (int temp2 = 0; temp2 < voxels_y; temp2++)
					{
						for (int temp3 = 0; temp3 < voxels_z; temp3++)
							myFile << m->energy_t[temp1][temp2][temp3][temp4] << " ";
						myFile << endl;
					}
						myFile << endl;
				}
				myFile.close();
			}
            else
				cout << "File could not be opened." << endl;
	}
}

void WriteAbsorbedEnergyToFile_Time_secondPulse(Medium * m)
{
	for (int temp4 = 0; temp4 < m->num_time_steps; temp4++)
	{
		stringstream ss;
		ss << "ResultsAbsorbedEnergyTime_SecondPulse_" << temp4 << ".txt";
		ofstream myFile(ss.str().c_str());
		if (myFile.is_open())
		{
			for (int temp1 = 0; temp1 < voxels_x; temp1++)
			{
				for (int temp2 = 0; temp2 < voxels_y; temp2++)
				{
					for (int temp3 = 0; temp3 < voxels_z; temp3++)
						myFile << m->energy_next[temp1][temp2][temp3][temp4] << " ";
					myFile << endl;
				}
				myFile << endl;
			}
			myFile.close();
		}
		else
			cout << "File could not be opened." << endl;
	}
}

void WritePhotonFluenceToFile(Medium * m)
{
	ofstream myFile;
	myFile.open("Results_photonFluence.txt");
	for(int temp1 = 0; temp1 < voxels_x; temp1++)
	{
		for(int temp2 = 0; temp2 < voxels_y; temp2++)
		{
			for(int temp3 = 0; temp3 < voxels_z; temp3++)
				myFile << m->fluence[temp1][temp2][temp3] << " ";
			myFile << endl;
		}
		myFile << endl;
	}
	myFile.close();
}


void CreateNewThread_time(Medium * m, Source * s, long numPhotons)		// Time-resolved
{
    for(long i = 0; i < numPhotons; i++)
    {
        RunPhoton_time(m,s);
    }
}

void CreateNewThread_steady(Medium * m, Source * s, long numPhotons)
{
	for (long i = 0; i < numPhotons; i++)
	{
		RunPhoton_steady(m, s);
	}
}

void CreateNewThread_secondPulse(Medium * m, Source * s, long numPhotons)
{
	for (long i = 0; i < numPhotons; i++)
	{
		RunPhotonNew_secondPulse(m, s);
	}
}

float GenerateRandomNumber()
{
	return (float)rand() / RAND_MAX;
}

void Prepare_SecondPulse(Medium * m, Source * s, float delay)
{
	float period = 1e9 / s->freq;		//	in ns
	int pulseTimeId = (int)floor(delay / time_step);

	//	copy last matrix
	for (int temp = 0; temp < voxels_x; temp++)
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
				for (int temp4 = 0; temp4 < (int)floor((time_end - delay) / time_step); temp4++)
					m->energy_next[temp][temp2][temp3][temp4] += m->energy_t[temp][temp2][temp3][pulseTimeId + temp4];
}