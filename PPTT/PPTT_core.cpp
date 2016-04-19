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
    releaseTime = 0.0;
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
    releaseTime = 0;
}

void Source::TimeProfile_flat(float pulse_duration)
{
    releaseTime = ((float)rand() / RAND_MAX) * pulse_duration;
}

void Source::TimeProfile_gaussian(float pulse_duration)
{
    releaseTime = pulse_duration * sqrt(-log(RandomNumber())) + (3 * pulse_duration);
}

void Source::TimeProfile_sech(float pulse_duration)
{
    releaseTime = 1 / cosh(RandomNumber());
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
    roundX = (int)floor(x);
    roundY = (int)floor(y);
    roundZ = (int)floor(z);
    prevRoundX = roundX;
    prevRoundY = roundY;
    prevRoundZ = roundZ;
    timeOfFlight = 0;
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

    roundX = (int)floor(x);
    roundY = (int)floor(y);
    roundZ = (int)floor(z);

    ux = s->ux;
    uy = s->uy;
    uz = s->uz;

    timeOfFlight = s->releaseTime;
}

float Photon::GenStep(float invAlbedo)
{
    float temp = -log(RandomNumber()) * invAlbedo;
    return temp;
}

void Photon::UpdatePos(Medium * m)
{
    float temp = GenStep(m->invAlbedo[regId]);

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
        temp_x = ((roundX + 1.0f - x) / ux);
    else
    {
        if(ux < 0.0)
        {
            if(roundX != x) temp_x = ((roundX - x) / ux);
            else temp_x = (-1 / ux);
        }
        else temp_x = 0.0;
    }
    if(uy > 0.0)
        temp_y = ((roundY + 1.0f - y) / uy);
    else
    {
        if(uy < 0.0)
        {
            if(roundY != y) temp_y = ((roundY - y) / uy);
            else temp_y = (-1 / uy);
        }
        else temp_y = 0.0;
    }
    if(uz > 0.0)
        temp_z = ((roundZ + 1.0f - z) / uz);
    else
    {
        if(uz < 0.0)
        {
            if(roundZ != z) temp_z = ((roundZ - z) / uz);
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
        timeOfFlight += GetTOF(m, step);
        //   cout << "Time of flight is " << time_of_flight << endl;
        remStep = GenStep(m->invAlbedo[regId]);
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
        timeOfFlight += GetTOF(m, step);
        //      cout << "Time of flight is " << time_of_flight << endl;
        remStep -= temp_step;
        if(remStep > 0)
            RoundPosition();
        else
        {
            RoundPosition();
            remStep = GenStep(m->invAlbedo[regId]);
            UpdateDir(m);
        }
    }
    // cout << "Moving to position x: " << x << " y: " << y << " z: " << z << endl;
}

float Photon::GetTOF(Medium * m, float step_size)
{
    return step_size * m->n[regId] / (lightSpeed * units);
}

int Photon::GetTimeId()
{
    return IntFloor(timeOfFlight / timeStep);
}

int Photon::CheckTOF(float end)
{
    if(timeOfFlight < end) return 0;
    else return 1;
}
void Photon::PosAndDir()
{
    cout << "Current position is (" << x << " " << y << " " << z << ")" << endl;
    cout << "Propagation direction is (" << ux << " " << uy << " " << uz << ")" << endl;
}

void Photon::RoundPosition()
{
    prevRoundX = roundX;
    prevRoundY = roundY;
    prevRoundZ = roundZ;

    roundX = (int)floor(x);
    roundY = (int)floor(y);
    roundZ = (int)floor(z);
}

int Photon::CheckBoundaries(Medium * m)
{
    if(roundZ < 0)
    {
        if(roundX < 0)
            roundX = 0;
        if(roundY < 0)
            roundY = 0;

        if(roundX >= voxelsX)
            roundX = voxelsX - 1;
        if(roundY <= voxelsY)
            roundY = voxelsY - 1;

        m->surroundingZ[roundX][roundY][0] += w;
        return 1;
    }
    if(roundZ >= voxelsZ)
    {
        if(roundX < 0)
            roundX = 0;
        if(roundY < 0)
            roundY = 0;

        if(roundX >= voxelsX)
            roundX = voxelsX - 1;
        if(roundY <= voxelsY)
            roundY = voxelsY - 1;

        m->surroundingZ[roundX][roundY][1] += w;
        return 1;
    }
    if(roundY < 0)
    {
        if(roundX < 0)
            roundX = 0;
        if(roundZ < 0)
            roundZ = 0;

        if(roundX >= voxelsX)
            roundX = voxelsX - 1;
        if(roundZ <= voxelsZ)
            roundZ = voxelsZ - 1;

        m->surroundingY[roundX][roundZ][0] += w;
        return 1;
    }
    if(roundY >= voxelsY)
    {
        if(roundX < 0)
            roundX = 0;
        if(roundZ < 0)
            roundZ = 0;

        if(roundX >= voxelsX)
            roundX = voxelsX - 1;
        if(roundZ <= voxelsZ)
            roundZ = voxelsZ - 1;

        m->surroundingY[roundX][roundZ][1] += w;
        return 1;
    }
    if(roundX < 0)
    {
        if(roundY < 0)
            roundY = 0;
        if(roundZ < 0)
            roundZ = 0;

        if(roundY >= voxelsY)
            roundY = voxelsY - 1;
        if(roundZ <= voxelsZ)
            roundZ = voxelsZ - 1;

        m->surroundingX[roundY][roundZ][0] += w;
        return 1;
    }
    if(roundX >= voxelsX)
    {
        if(roundY < 0)
            roundY = 0;
        if(roundZ < 0)
            roundZ = 0;

        if(roundY >= voxelsY)
            roundY = voxelsY - 1;
        if(roundZ <= voxelsZ)
            roundZ = voxelsZ - 1;

        m->surroundingX[roundY][roundZ][1] += w;
        return 1;
    }
    return 0;
}
void Photon::GenDir(float g) // generate cos(theta) and phi based on anisotropy parameter g
{
    if(g != 0)
    {
        cosTheta = 1.0f / (2.0f * g) * (1 + g*g - powf((1 - g*g) / (1 - g + 2 * g * RandomNumber()), 2));
    }
    else
        cosTheta = 2 * RandomNumber() - 1;
    phi = 2 * PI * RandomNumber();
}

void Photon::UpdateDir(Medium * m)
{
    GenDir(m->g[regId]);
    float temp_ux, temp_uy, temp_uz;
    float sinTheta = sin(acos(cosTheta));
    float temp_sqrt = 1.0f / sqrt(1 - uz*uz);
    if(abs(uz) < 0.999)
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
    if(m->n[regId] - m->n[lastRegId])
        return 0;
    else return 1;
}

float Photon::GetReflectionCoef(Medium * m)
{
    return ((m->n[lastRegId] - m->n[regId])*(m->n[lastRegId] - m->n[regId])) / ((m->n[lastRegId] + m->n[regId])*(m->n[lastRegId] + m->n[regId]));
}

void Photon::Reflect(Medium * m)
{
    if(prevRoundZ != roundZ)
        uz = -uz;
    if(prevRoundY != roundY)
        uy = -uy;
    if(prevRoundX != roundX)
        ux = -ux;
}

void Photon::Transmis(Medium * m)
{
    if(prevRoundZ != roundZ)
    {
        float alpha_i = acos(abs(uz));
        float alpha_t = asin(sin(alpha_i) * m->n[lastRegId] / m->n[regId]);
        ux *= m->n[lastRegId] / m->n[regId];
        uy *= m->n[lastRegId] / m->n[regId];
        uz *= sin(alpha_t) / abs(uz);
    }
    if(prevRoundY != roundY)
    {
        float alpha_i = acos(abs(uy));
        float alpha_t = asin(sin(alpha_i) * m->n[lastRegId] / m->n[regId]);
        ux *= m->n[lastRegId] / m->n[regId];
        uz *= m->n[lastRegId] / m->n[regId];
        uy *= sin(alpha_t) / abs(uy);
    }
    if(prevRoundX != roundX)
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
    numberOfRegions = 1;								// number of regions inicialization

    structure = new int**[voxelsX];						// dynamic allocation 
    for(int temp = 0; temp < voxelsX; temp++)
        structure[temp] = new int*[voxelsY];
    for(int temp1 = 0; temp1 < voxelsX; temp1++)
        for(int temp2 = 0; temp2 < voxelsX; temp2++)
            structure[temp1][temp2] = new int[voxelsZ];

    energy = new float**[voxelsX];						// dynamic allocation 
    for(int temp = 0; temp < voxelsX; temp++)
        energy[temp] = new float*[voxelsY];
    for(int temp1 = 0; temp1 < voxelsX; temp1++)
        for(int temp2 = 0; temp2 < voxelsX; temp2++)
            energy[temp1][temp2] = new float[voxelsZ];

    fluence = new float**[voxelsX];						// dynamic allocation 
    for(int temp = 0; temp < voxelsX; temp++)
        fluence[temp] = new float*[voxelsY];
    for(int temp1 = 0; temp1 < voxelsX; temp1++)
        for(int temp2 = 0; temp2 < voxelsX; temp2++)
            fluence[temp1][temp2] = new float[voxelsZ];

    surroundingX = new float**[voxelsY];					// allocation of surrounding matrix 
    for(int temp = 0; temp < voxelsY; temp++)				// this one for front and rear x plane
        surroundingX[temp] = new float*[voxelsZ];
    for(int temp1 = 0; temp1 < voxelsY; temp1++)
        for(int temp2 = 0; temp2 < voxelsZ; temp2++)
            surroundingX[temp1][temp2] = new float[2];
    for(int temp = 0; temp < voxelsY; temp++)
        for(int temp2 = 0; temp2 < voxelsZ; temp2++)
            for(int temp3 = 0; temp3 < 2; temp3++)
                surroundingX[temp][temp2][temp3] = 0;

    surroundingY = new float**[voxelsX];					// allocation of surrounding matrix 
    for(int temp = 0; temp < voxelsX; temp++)				// this one for front and rear y plane
        surroundingY[temp] = new float*[voxelsZ];
    for(int temp1 = 0; temp1 < voxelsX; temp1++)
        for(int temp2 = 0; temp2 < voxelsZ; temp2++)
            surroundingY[temp1][temp2] = new float[2];
    for(int temp = 0; temp < voxelsX; temp++)
        for(int temp2 = 0; temp2 < voxelsZ; temp2++)
            for(int temp3 = 0; temp3 < 2; temp3++)
                surroundingY[temp][temp2][temp3] = 0;

    surroundingZ = new float**[voxelsX];					// allocation of surrounding matrix 
    for(int temp = 0; temp < voxelsX; temp++)				// this one for front and rear z plane
        surroundingZ[temp] = new float*[voxelsY];
    for(int temp1 = 0; temp1 < voxelsX; temp1++)
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            surroundingZ[temp1][temp2] = new float[2];
    for(int temp = 0; temp < voxelsX; temp++)
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            for(int temp3 = 0; temp3 < 2; temp3++)
                surroundingZ[temp][temp2][temp3] = 0;

    for(int temp = 0; temp < voxelsX; temp++)			// structure Ids inicialization
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
                structure[temp][temp2][temp3] = 0;
    for(int temp = 0; temp < voxelsX; temp++)			// matrix with absorbed energy inicialization
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
                energy[temp][temp2][temp3] = 0;
    for(int temp = 0; temp < voxelsX; temp++)			// matrix with photon fluence inicialization
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
                fluence[temp][temp2][temp3] = 0;
    for(int temp = 0; temp < maxRegions; temp++)
        ua[temp] = 0;
    for(int temp = 0; temp < maxRegions; temp++)
        us[temp] = 0;
    for(int temp = 0; temp < maxRegions; temp++)
        invAlbedo[temp] = 0;
    for(int temp = 0; temp < maxRegions; temp++)
        g[temp] = 1;
    for(int temp = 0; temp < maxRegions; temp++)
        n[temp] = 1;
    for(int temp = 0; temp < maxRegions; temp++)
        k[temp] = 0;
    for(int temp = 0; temp < maxRegions; temp++)
        rho[temp] = 0;
    for(int temp = 0; temp < maxRegions; temp++)
        c_h[temp] = 0;
    for(int temp = 0; temp < maxRegions; temp++)
        w_g[temp] = 0;
	blood_density = BLOOD_DENSITY / powf(units, 3);

    // Preparing array for time-resolved simulations
    num_time_steps = (int)ceil((timeEnd - timeStart) / timeStep);

    energy_t = new float***[voxelsX];						// dynamic allocation 
    for(int temp = 0; temp < voxelsX; temp++)
        energy_t[temp] = new float**[voxelsY];
    for(int temp1 = 0; temp1 < voxelsX; temp1++)
        for(int temp2 = 0; temp2 < voxelsX; temp2++)
            energy_t[temp1][temp2] = new float*[voxelsZ];
    for(int temp = 0; temp < voxelsX; temp++)
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
                energy_t[temp][temp2][temp3] = new float[num_time_steps];
    for(int temp = 0; temp < voxelsX; temp++)
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
                for(int temp4 = 0; temp4 < num_time_steps; temp4++)
                    energy_t[temp][temp2][temp3][temp4] = 0;

    energy_next = new float***[voxelsX];						// dynamic allocation 
    for(int temp = 0; temp < voxelsX; temp++)
        energy_next[temp] = new float**[voxelsY];
    for(int temp1 = 0; temp1 < voxelsX; temp1++)
        for(int temp2 = 0; temp2 < voxelsX; temp2++)
            energy_next[temp1][temp2] = new float*[voxelsZ];
    for(int temp = 0; temp < voxelsX; temp++)
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
                energy_next[temp][temp2][temp3] = new float[num_time_steps];

    for(int temp = 0; temp < voxelsX; temp++)
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
                for(int temp4 = 0; temp4 < num_time_steps; temp4++)
                    energy_next[temp][temp2][temp3][temp4] = 0;
}

Medium::~Medium()
{
}

void Medium::PrintMediumProperties()
{
    cout << "Structure: " << endl;
    for(int temp = 0; temp < voxelsX; temp++)
    {
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
        {
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
                cout << "  " << structure[temp][temp2][temp3];
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;
    cout << "Absorption coefficients: " << endl;
    for(int temp = 0; temp < maxRegions; temp++)
        cout << "  " << ua[temp];
    cout << endl;
    cout << "Scattering coefficients: " << endl;
    for(int temp = 0; temp < maxRegions; temp++)
        cout << "  " << us[temp];
    cout << endl;
    cout << "Inverse Albedo values: " << endl;
    for(int temp = 0; temp < maxRegions; temp++)
        cout << "  " << invAlbedo[temp];
    cout << endl;
    cout << "Anisotropy parameter: " << endl;
    for(int temp = 0; temp < maxRegions; temp++)
        cout << "  " << g[temp];
    cout << endl;
    cout << "Refractive indices: " << endl;
    for(int temp = 0; temp < maxRegions; temp++)
        cout << "  " << n[temp];
    cout << endl;
}

int Medium::RetRegId(Photon * p)
{
    return structure[p->roundX][p->roundY][p->roundZ];
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

    for(int temp1 = 0; temp1 < dim_x; temp1++)
        for(int temp2 = 0; temp2 < dim_y; temp2++)
            for(int temp3 = 0; temp3 < dim_z; temp3++)
                structure[start_x + temp1][start_y + temp2][start_z + temp3] = numberOfRegions;
    ua[numberOfRegions] = set_ua;
    us[numberOfRegions] = set_us;
    invAlbedo[numberOfRegions] = 1 / (set_ua + set_us);
    g[numberOfRegions] = set_g;
    n[numberOfRegions] = set_n;
    numberOfRegions++;
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
                    structure[center_x + temp1][center_y + temp2][center_z + temp3] = numberOfRegions;
    ua[numberOfRegions] = set_ua;
    us[numberOfRegions] = set_us;
    invAlbedo[numberOfRegions] = 1 / (set_ua + set_us);
    g[numberOfRegions] = set_g;
    n[numberOfRegions] = set_n;
    numberOfRegions++;
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

    for(int temp = 0; temp < (int)round(length); temp++)
        structure[(int)round(start_x + dir_x*temp)][(int)round(start_y + dir_y*temp)][(int)round(start_z + dir_z*temp)] = numberOfRegions;

    ua[numberOfRegions] = set_ua;
    us[numberOfRegions] = set_us;
    invAlbedo[numberOfRegions] = 1 / (set_ua + set_us);
    g[numberOfRegions] = set_g;
    n[numberOfRegions] = set_n;
    numberOfRegions++;
}

void Medium::AbsorbEnergy(Photon * p)
{
    float temp = p->w * (ua[p->regId] * invAlbedo[p->regId]);
    energy[p->prevRoundX][p->prevRoundY][p->prevRoundZ] += temp;
    p->w -= temp;
}

void Medium::AbsorbEnergyBeer(Photon * p)
{
    float temp = p->w * exp(-ua[p->regId] * p->step);							// Lambert-Beer law
    energy[p->prevRoundX][p->prevRoundY][p->prevRoundZ] += (p->w - temp);
    p->w = temp;
}

void Medium::AbsorbEnergyBeer_Time(Photon * p)
{
    float temp = p->w * exp(-ua[p->regId] * p->step);							// Lambert-Beer law
    energy_t[p->prevRoundX][p->prevRoundY][p->prevRoundZ][p->timeId] += (p->w - temp);
    p->w = temp;
}

void Medium::AbsorbEnergyBeer_Time_secondPulse(Photon * p)
{
    float temp = p->w * exp(-ua[p->regId] * p->step);							// Lambert-Beer law
    energy_next[p->prevRoundX][p->prevRoundY][p->prevRoundZ][p->timeId] += (p->w - temp);
    p->w = temp;
}
void Medium::RescaleEnergy(long num_photons)
{
    for(int temp = 0; temp < voxelsX; temp++)			// structure Ids inicialization
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
                energy[temp][temp2][temp3] /= (num_photons / powf((float)units, 3));
}

void Medium::RescaleEnergy_Time(long num_photons, float time_min_step)
{
    for(int temp = 0; temp < voxelsX; temp++)			// structure Ids inicialization
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
                for(int temp4 = 0; temp4 < num_time_steps; temp4++)
                    energy_t[temp][temp2][temp3][temp4] /= (time_min_step * num_photons / powf((float)units, 3));
}

void Medium::RecordFluence()
{
    for(int temp = 0; temp < voxelsX; temp++)			// structure Ids inicialization
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
                fluence[temp][temp2][temp3] = energy[temp][temp2][temp3] / ua[structure[temp][temp2][temp3]];
}

//////////////////////////////////////////////////////
//          Heat class
//////////////////////////////////////////////////////

Heat::Heat()
{
	h_timeStart = 0;
	h_timeEnd = H_TIME_STEP;
	h_timeStep = H_TIME_END;

	h_num_time_steps = (int)ceil((h_timeEnd - h_timeStart) / h_timeStep);

	temperature = new float**[voxelsX];						// dynamic allocation 
    for(int temp = 0; temp < voxelsX; temp++)
        temperature[temp] = new float*[voxelsY];
    for(int temp1 = 0; temp1 < voxelsX; temp1++)
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            temperature[temp1][temp2] = new float[voxelsZ];

	temperature_time = new float***[voxelsX];
	for (int temp = 0; temp < voxelsX; temp++)
		temperature_time[temp] = new float**[voxelsY];
	for (int temp1 = 0; temp1 < voxelsX; temp1++)
		for (int temp2 = 0; temp2 < voxelsY; temp2++)
			temperature_time[temp1][temp2] = new float*[voxelsZ];
	for (int temp1 = 0; temp1 < voxelsX; temp1++)
		for (int temp2 = 0; temp2 < voxelsY; temp2++)
			for (int temp3 = 0; temp3 < voxelsZ; temp3++)
				temperature_time[temp1][temp2][temp3] = new float[h_num_time_steps];

	for (int temp1 = 0; temp1 < voxelsX; temp1++)
		for (int temp2 = 0; temp2 < voxelsY; temp2++)
			for (int temp3 = 0; temp3 < voxelsZ; temp3++)
				for (int temp4 = 0; temp4 < h_num_time_steps; temp4++)
				temperature_time[temp1][temp2][temp3][temp4] = 36.5f;

    for(int temp = 0; temp < voxelsX; temp++)					// set human body temperature
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
                temperature[temp][temp2][temp3] = 36.5;
}

Heat::~Heat()
{
}

void Heat::AddThermalCoef(Medium * m, int mediumId, float specific_heat, float density, float conduction, float blood_perfusivity) 
{
    m->c_h[mediumId] = specific_heat / 1000;				// input in [J/kg/K]
    m->rho[mediumId] = density / (1000000 * UNITS);					// [kg/m3]
    m->k[mediumId] = conduction / (1000 * UNITS);		// [W/m/K]
    m->w_g[mediumId] = blood_perfusivity;			// [ml/s/ml]
}

void Heat::PennesEquation(Medium * m, float arterial_temperature)
{
	float ***temperature_help = new float**[voxelsX];						// dynamic allocation 
	for (int temp = 0; temp < voxelsX; temp++)
		temperature_help[temp] = new float*[voxelsY];
	for (int temp1 = 0; temp1 < voxelsX; temp1++)
		for (int temp2 = 0; temp2 < voxelsY; temp2++)
			temperature_help[temp1][temp2] = new float[voxelsZ];

	float ***epsilon = new float**[voxelsX];					// score the difference between new and old temperature
	for (int temp = 0; temp < voxelsX; temp++)
		epsilon[temp] = new float*[voxelsY];
	for (int temp1 = 0; temp1 < voxelsX; temp1++)
		for (int temp2 = 0; temp2 < voxelsY; temp2++)
			epsilon[temp1][temp2] = new float[voxelsZ];

	for (int temp = 0; temp < voxelsX; temp++)					// set human body temperature
		for (int temp2 = 0; temp2 < voxelsY; temp2++)
			for (int temp3 = 0; temp3 < voxelsZ; temp3++)
				temperature_help[temp][temp2][temp3] = 36.5f;

	float iterative_condition = JACOBI_ITERATIVE + 0.1f; 
	while (iterative_condition > JACOBI_ITERATIVE)
	{
		float iterative_condition_help = 0;
		for (int temp = 1; temp < voxelsX - 1; temp++)					// calculate operator
			for (int temp2 = 1; temp2 < voxelsY - 1; temp2++)
				for (int temp3 = 1; temp3 < voxelsZ - 1; temp3++)
				{
					int tempId = m->structure[temp][temp2][temp3];
					temperature_help[temp][temp2][temp3] = (temperature[temp + 1][temp2][temp3] + temperature[temp - 1][temp2][temp3] + 
						temperature[temp][temp2 + 1][temp3] + temperature[temp][temp2 - 1][temp3] + temperature[temp][temp2][temp3 + 1] + 
						temperature[temp][temp2][temp3 - 1] + (m->energy[temp][temp2][temp3] / m->k[tempId]) + 
						m->w_g[tempId] * BLOOD_CAPACITY * m->blood_density * (arterial_temperature - temperature[temp][temp2][temp3]) / m->k[tempId]) / 6.0f;
				}

		for (int temp = 1; temp < voxelsX - 1; temp++)					// calculate epsilon and
			for (int temp2 = 1; temp2 < voxelsY - 1; temp2++)			// give back into main matrix
				for (int temp3 = 1; temp3 < voxelsZ - 1; temp3++)
				{
					epsilon[temp][temp2][temp3] = abs((temperature[temp][temp2][temp3] - temperature_help[temp][temp2][temp3]) / temperature_help[temp][temp2][temp3]);
					temperature[temp][temp2][temp3] = temperature_help[temp][temp2][temp3];
					if (iterative_condition_help < epsilon[temp][temp2][temp3])
						iterative_condition_help = epsilon[temp][temp2][temp3];
				}
		iterative_condition = iterative_condition_help;
	}
}

void Heat::PennesEquation_time(Medium * m, float arterial_temperature)
{
	for (int temp4 = 0; temp4 < (h_num_time_steps - 1); temp4++)
		for (int temp = 0; temp < voxelsX; temp++)					// set human body temperature
			for (int temp2 = 0; temp2 < voxelsY; temp2++)
				for (int temp3 = 0; temp3 < voxelsZ; temp3++)
				{
					int tempId = m->structure[temp][temp2][temp3];
					temperature_time[temp][temp2][temp3][temp4 + 1] = temperature_time[temp][temp2][temp3][temp4] +
						(h_timeStep / (m->rho[tempId] * m->c_h[tempId]))*(m->k[tempId] * (temperature_time[temp + 1][temp2][temp3][temp4] + temperature_time[temp - 1][temp2][temp3][temp4] +
							temperature_time[temp][temp2 + 1][temp3][temp4] + temperature_time[temp][temp2 - 1][temp3][temp4] + temperature_time[temp][temp2][temp3 + 1][temp4] +
							temperature_time[temp][temp2][temp3 - 1][temp4] - 6.0f * temperature_time[temp][temp2][temp3][temp4]) +
							m->w_g[tempId] * BLOOD_CAPACITY * m->blood_density * (arterial_temperature - temperature_time[temp][temp2][temp3][temp4]) + m->energy_t[temp][temp2][temp3][temp4]);
				}
}

float Heat::AproximateBloodPerfusivity(float omega0, float omega1, float omega2, float omega3, float temperature)
{
    return (omega0 + omega1 * (atan(omega2*(temperature - omega3)))) / 1e6f;		// conversion to milimeters and grams
}

//////////////////////////////////////////////////////
//				Nonclass functions
//////////////////////////////////////////////////////
void RunPhoton_steady(Medium * m, Source * s)
{
    Photon * p = new Photon;
    //s->CollimatedGaussianBeam(5.0, 5.0, 0.0, 0.5, 0.0, 0.0, 1.0);
    s->TimeProfile_infiniteSharp();
    p->GetSourceParameters(s);
    p->regId = m->RetRegId(p);
    p->lastRegId = p->regId;
    p->remStep = p->GenStep(m->invAlbedo[p->regId]);


    while(p->w > PHOTON_DEATH)
    {
        if(p->CheckRefIndexMismatch(m))
            p->Move(m);
        else
        {
            if(RandomNumber() < p->GetReflectionCoef(m))
                p->Reflect(m);
            else
                p->Transmis(m);
            p->Move(m);
        }
        if(p->CheckBoundaries(m)) break;
        p->lastRegId = p->regId;
        p->regId = m->RetRegId(p);
        m->AbsorbEnergyBeer(p);
    }

    delete p;
}

void RunPhoton_time(Medium * m, Source * s)
{
    Photon * p = new Photon;
    //s->CollimatedGaussianBeam(5.0, 5.0, 0.0, 0.5, 0.0, 0.0, 1.0);
    s->TimeProfile_flat(pulseDuration);
    p->GetSourceParameters(s);
    p->regId = m->RetRegId(p);
    p->lastRegId = p->regId;
    p->remStep = p->GenStep(m->invAlbedo[p->regId]);


    while(p->w > PHOTON_DEATH)
    {
        if(p->CheckRefIndexMismatch(m))
            p->Move(m);
        else
        {
            if(RandomNumber() < p->GetReflectionCoef(m))
                p->Reflect(m);
            else
                p->Transmis(m);
            p->Move(m);
        }
        if(p->CheckBoundaries(m)) break;
        if(p->CheckTOF(timeEnd)) break;
        p->timeId = p->GetTimeId();
        p->lastRegId = p->regId;
        p->regId = m->RetRegId(p);
        m->AbsorbEnergyBeer_Time(p);
    }

    delete p;
}

void PrintAbsEnergy(Medium * m)
{
    for(int temp1 = 0; temp1 < voxelsX; temp1++)
    {
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
        {
            cout << endl;
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
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
    for(int temp1 = 0; temp1 < voxelsX; temp1++)
    {
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
        {
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
                myFile << m->energy[temp1][temp2][temp3] << " ";
            myFile << endl;
        }
        myFile << endl;
    }
    myFile.close();
}

void WriteAbsorbedEnergyToFile_Time(Medium * m)
{
    for(int temp4 = 0; temp4 < m->num_time_steps; temp4++)
    {
        stringstream ss;
        ss << "ResultsAbsorbedEnergyTime_" << temp4 << ".txt";
        ofstream myFile(ss.str().c_str());
        if(myFile.is_open())
        {
            for(int temp1 = 0; temp1 < voxelsX; temp1++)
            {
                for(int temp2 = 0; temp2 < voxelsY; temp2++)
                {
                    for(int temp3 = 0; temp3 < voxelsZ; temp3++)
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
    for(int temp4 = 0; temp4 < m->num_time_steps; temp4++)
    {
        stringstream ss;
        ss << "ResultsAbsorbedEnergyTime_SecondPulse_" << temp4 << ".txt";
        ofstream myFile(ss.str().c_str());
        if(myFile.is_open())
        {
            for(int temp1 = 0; temp1 < voxelsX; temp1++)
            {
                for(int temp2 = 0; temp2 < voxelsY; temp2++)
                {
                    for(int temp3 = 0; temp3 < voxelsZ; temp3++)
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
    for(int temp1 = 0; temp1 < voxelsX; temp1++)
    {
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
        {
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
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
        RunPhoton_time(m, s);
    }
}

inline int IntFloor(float x)
{
    int i = (int)x; /* truncate */
    return i - (i > x); /* convert trunc to floor */
}
void CreateNewThread_steady(Medium * m, Source * s, long numPhotons)
{
    for(long i = 0; i < numPhotons; i++)
    {
        RunPhoton_steady(m, s);
    }
}

float RandomNumber()
{
    // fast way to generate random number 
    union {
        uint32_t d;
        float f;
    } u;
    u.d = (((uint32_t)rand() & 0x7fff) << 8) | 0x3f800000;
    return u.f - 1.0f;
}

void CreateEnviroment(Medium * m, Heat * h)
{
	m->CreateCube(0, 0, 0, voxelsX / units, voxelsY / units, voxelsZ / units, 0.001f, 0.001f, 1.0f, 1.0f);		// background for avoiding errors
	
	// Validation with MCML
	m->CreateCube(0, 0, 0, 10, 10, 1, 1, 100, 0.9f, 1.37f);		// Layer 1
	m->CreateCube(0, 0, 1, 10, 10, 1, 1, 10, 0.0f, 1.37f);		// Layer 2	
	m->CreateCube(0, 0, 2, 10, 10, 2, 5, 10, 0.7f, 1.37f);		// Layer 3
	
	h->AddThermalCoef(m, 2, 2000, 1200, 0.3f, 0.0f);				// Layer 1
	h->AddThermalCoef(m, 3, 3600, 1050, 0.5f, 0.5f);				// Layer 2
	h->AddThermalCoef(m, 4, 2350, 900, 1.2f, 0.012f);				// Layer 3
}