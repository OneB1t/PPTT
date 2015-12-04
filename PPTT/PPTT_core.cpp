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
	release_time = 0;
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

    float cosTheta = 2 * RandomNumber() - 1;
    float sinTheta = sqrt(1 - cosTheta*cosTheta);

    float phi = 2 * PI * RandomNumber();
    float sinPhi, cosPhi = cos(phi);

    if(phi < PI)
        sinPhi = sqrt(1 - cosPhi * cosPhi);
    else
        sinPhi = -sqrt(1 - cosPhi * cosPhi);

    ux = sinTheta * sinPhi;
    uy = sinTheta * cosPhi;
    uz = cosTheta;
}

void Source::CollimatedGaussianBeam(float set_x, float set_y, float set_z, float radius, float set_ux, float set_uy, float set_uz)
{
    set_x = set_x * units;
    set_y = set_y * units;
    set_z = set_z * units;
    radius = radius * units;

    float temp = sqrt(set_ux*set_ux + set_uy*set_uy + set_uz*set_uz);
    ux = set_ux / temp;
    uy = set_uy / temp;
    uz = set_uz / temp;

    float phi = 2 * PI*(RandomNumber());
    float r = radius * sqrt(-log(RandomNumber()));

    x = set_x + r * cos(phi);
    y = set_y + r * sin(phi);
    z = set_z;

    //cout << "Parameters of the source are: " << endl;
    //cout << "x = " << x << " y = " << y << " z = " << z << endl;
}

void Source::Circular_flat_beam(float set_x, float set_y, float set_z, float radius, float set_ux, float set_uy, float set_uz)
{
    set_x *= units;
    set_y *= units;
    set_z *= units;
    radius *= units;

    float phi = 2 * PI * (RandomNumber());
    float r = radius * sqrt(RandomNumber());

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

/////////////////////////////
// Photon class
/////////////////////////////

Photon::Photon()
{
    // cout << "Photon inicialization..." << endl;
    x = y = 10;
    z = 0;
    ux = uy = 0;
    uz = 1;
    w = 1;
    round_x = floor(x);
    round_y = floor(y);
    round_z = floor(z);
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
    
    ux = s->ux;
    uy = s->uy;
    uz = s->uz;

	time_of_flight = s->release_time;
}

float Photon::GenStep(float invAlbedo)
{
    float temp = -log(RandomNumber()) * invAlbedo;
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
    float temp_x;
    float temp_y;
    float temp_z;

    //  Calculate distances
    if(ux > 0.0)
        temp_x = ((round_x + 1.0f - x) / ux);
    else
    {
        if(ux < 0.0)
        {
            if(round_x != x) temp_x = ((round_x - x) / ux);
            else temp_x = (-1 / ux);
        }
        else temp_x = 0.0;
    }
    if(uy > 0.0)
        temp_y = ((round_y + 1.0f - y) / uy);
    else
    {
        if(uy < 0.0)
        {
            if(round_y != y) temp_y = ((round_y - y) / uy);
            else temp_y = (-1 / uy);
        }
        else temp_y = 0.0;
    }
    if(uz > 0.0)
        temp_z = ((round_z + 1.0f - z) / uz);
    else
    {
        if(uz < 0.0)
        {
            if(round_z != z) temp_z = ((round_z - z) / uz);
            else temp_z = (-1 / uz);
        }
        else temp_z = 0.0;
    }

    /*cout << "Direction vectors and distances to the edges are: " << endl;
    cout << "   ux = " << ux << "   dx = " << temp_x << endl;
    cout << "   uy = " << uy << "   dy = " << temp_y << endl;
    cout << "   uz = " << uz << "   dz = " << temp_z << endl;*/

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
    //cout << "Temp_first " << temp_first << " Temp_second " << temp_second << " Temp _third " << temp_third << endl;
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
    //cout << "Old scattering " << m->us[lastRegId] << " New scattering " << m->us[regId] << endl;
    //cout << "Remaining step size " << remStep << endl;

    if(temp_step > remStep)
    {
        step = remStep;
        // cout << "Making " << step << " step" << endl;
        x = x + ux * step;
        y = y + uy * step;
        z = z + uz * step;
        time_of_flight += GetTOF(m, step);
        //   cout << "Time of flight is " << time_of_flight << endl;
        remStep = GenStep(m->inv_albedo[regId]);
        UpdateDir(m);
        RoundPosition();
    }
    else
    {
        step = temp_step;
        //cout << "Making " << step << " step" << endl;
        x = x + ux * step;
        y = y + uy * step;
        z = z + uz * step;
        time_of_flight += GetTOF(m, step);
        //      cout << "Time of flight is " << time_of_flight << endl;
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
    // cout << "Moving to position x: " << x << " y: " << y << " z: " << z << endl;
}

float Photon::GetTOF(Medium * m, float step_size)
{
    return step_size * m->n[regId] / (light_speed * units);
}

int Photon::GetTimeId()
{
	return int_floor(time_of_flight / time_step);
}

int Photon::CheckTOF(float end)
{
    if(time_of_flight < end) return 0;
    else return 1;
}
void Photon::PosAndDir()
{
    cout << "Current position is (" << x << " " << y << " " << z << ")" << endl;
    cout << "Propagation direction is (" << ux << " " << uy << " " << uz << ")" << endl;
}

void Photon::RoundPosition()
{
    round_x = int_floor(x);
    round_y = int_floor(y);
    round_z = int_floor(z);
}

int Photon::CheckBoundaries()
{
    if(z > 0 && z < voxels_z)
        if(y > 0 && y < voxels_y)
            if(x > 0 && x < voxels_x)
                return 0;
            else return 1;
        else return 1;
    else return 1;
}
void Photon::GenDir(float g) // generate cos(theta) and phi based on anisotropy parameter g
{
    if(g != 0)
    {
        cosTheta = 1.0 / (2.0 * g) * (1 + g*g - pow((1 - g*g) / (1 - g + 2 * g * RandomNumber()), 2));
    }
    else
        cosTheta = 2 * RandomNumber() - 1;
    phi = 2 * PI * RandomNumber();
}

void Photon::UpdateDir(Medium * m)
{
	GenDir(m->g[regId]);
	// cout << "anisotropy parameter " << s->anisotropyParamReg[regId] << endl << "cosTheta " << cosTheta << endl << "phi " << phi << endl;
	float temp_ux, temp_uy, temp_uz;
	float sinTheta = sin(acos(cosTheta)); //cout << "sinTheta " << sinTheta << endl;
	float temp_sqrt = 1.0 / sqrt(1 - uz*uz); //cout << "temp sqrt " << temp_sqrt << endl;
	if (abs(uz) < 0.99999)
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
	//cout << "X direction " << ux << endl;
	//cout << "Y direction " << uy << endl;
	//cout << "Z direction " << uz << endl;
}

//////////////////////////////////////////////////
//			Medium class
//////////////////////////////////////////////////
Medium::Medium()
{
	number_of_regions = 1;								// number of regions inicialization
				
    for (int temp = 0; temp < voxels_x; temp++)			// matrix with photon fluence inicialization
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
            {
                
				fluence[temp][temp2][temp3] = 0;
                structure[temp][temp2][temp3] = 0;
                energy[temp][temp2][temp3] = 0;
            }
	for (int temp = 0; temp < max_regions; temp++)
    {
		ua[temp] = 0;
        us[temp] = 0;
        inv_albedo[temp] = 0;
        g[temp] = 1;
        n[temp] = 1;
        k[temp] = 0;
        rho[temp] = 0;
        c_h[temp] = 0;
    }

	// Preparing array for time-resolved simulations
	num_time_steps = (int)ceil((time_end - time_start) / time_step);

	for (int temp = 0; temp < voxels_x; temp++)
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
				energy_t[temp][temp2][temp3] = new float[num_time_steps];
    for (int temp = 0; temp < voxels_x; temp++)
		for (int temp2 = 0; temp2 < voxels_y; temp2++)
			for (int temp3 = 0; temp3 < voxels_z; temp3++)
                            for (int temp4 = 0; temp4 < num_time_steps; temp4++)
                                energy_t[temp][temp2][temp3][temp4] = 0;
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
	return structure[int_floor(p->x)][int_floor(p->y)][int_floor(p->z)];
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

void Medium::AbsorbEnergy(Photon * p)
{
    float temp = p->w * (ua[p->regId] * inv_albedo[p->regId]);
    energy[p->round_x][p->round_y][p->round_z] += temp;
    p->w -= temp;
}

void Medium::AbsorbEnergyBeer(Photon * p)
{
    //cout << "New weight of the photon is " << p->w << endl;
    //cout << "Absorption coef: " << ua[p->regId] << " Step size: " << p->step << endl;
    float temp = p->w * (1 - (ua[p->regId] * p->step) + (ua[p->regId] * ua[p->regId] * p->step * p->step / 2)); // Taylor expansion series of Lambert-Beer law
    energy[p->round_x][p->round_y][p->round_z] += (p->w - temp);
    p->w = temp;
}

void Medium::AbsorbEnergyBeer_Time(Photon * p)
{
	//cout << "New weight of the photon is " << p->w << endl;
	//cout << "Absorption coef: " << ua[p->regId] << " Step size: " << p->step << endl;
	float temp = p->w * (1 - (ua[p->regId] * p->step) + (ua[p->regId] * ua[p->regId] * p->step * p->step / 2)); // Taylor expansion series of Lambert-Beer law
	energy_t[p->round_x][p->round_y][p->round_z][p->timeId] += (p->w - temp);
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
}

Heat::~Heat()
{
}

void Heat::AddThermalCoef(Medium * m, int mediumId, float specific_heat, float density, float conduction) // in g, mm, K
{
    m->c_h[mediumId] = specific_heat;
    m->rho[mediumId] = density;
    m->k[mediumId] = conduction;
}
//////////////////////////////////////////////////////
//				Nonclass functions
//////////////////////////////////////////////////////
void RunPhoton(Medium * m)
{
	Photon * p = new Photon;
	p->regId = m->RetRegId(p);
        //cout << "RegId of the photon is " << p->regId << endl;
        
	while(p->w > PHOTON_DEATH) 
	{
		p->UpdatePos(m);
		if (p->CheckBoundaries()) break;
		p->RoundPosition();
		p->regId = m->RetRegId(p);
		m->AbsorbEnergy(p);
		p->UpdateDir(m);
	}
        
        delete p;
}

void RunPhotonNew(Medium * m, Source * s)
{
    Photon * p = new Photon;
	//s->Collimated_gaussian_beam(5.0, 5.0, 0.0, 0.5, 0.0, 0.0, 1.0);
	s->TimeProfile_flat(pulseDuration);
    p->GetSourceParameters(s);
        p->regId = m->RetRegId(p);
        p->lastRegId = p->regId;
        p->remStep = p->GenStep(m->inv_albedo[p->regId]);
        //cout << "Photon step " << p->remStep << endl;
        
	while(p->w > PHOTON_DEATH) 
	{
            //cout << "================" << endl;
            //cout << "New step" << endl;
            //cout << "================" << endl;
		p->Move(m);
                //p->PosAndDir();
        if (p->CheckBoundaries()) break;
		if (p->CheckTOF(time_end)) break;
		p->timeId = p->GetTimeId();
                // cout << "Time Id is " << p->timeId << endl;
		p->lastRegId = p->regId;
		p->regId = m->RetRegId(p);
                //p->w /= 2;
		m->AbsorbEnergyBeer_Time(p);		
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
		ss << "ResultsAbsorbedEnergyTime" << temp4 << ".txt";
		ofstream myFile(ss.str().c_str());
		// myFile.open("Results_absorbedEnergy.txt");
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


void CreateNewThread(Medium * m, Source * s, long numPhotons)
{
    for(long i = 0; i < numPhotons; i++)
    {
        RunPhotonNew(m,s);
    }
}

inline int int_floor(float x)
{
    int i = (int)x; /* truncate */
    return i - (i > x); /* convert trunc to floor */
}

float RandomNumber()
{
    union {
        uint32_t d;
        float f;
    } u;
    u.d = (((uint32_t)rand() & 0x7fff) << 8) | 0x3f800000;
    return u.f - 1.0f;
}