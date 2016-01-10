#define PHOTON_DEATH	0.0001
#define VOXELS_X 100
#define VOXELS_Y 100
#define VOXELS_Z 40
#define MAX_REGIONS 16
#define TIME_SEGMENTS 12
#define UNITS 10
#define PULSE_DURATION 0.0025f
#define PI 3.14159265358979323846	
#define JACOBI_ITERATIVE	0.01
#define H_TIME_STEP			0.1			// in seconds
#define H_TIME_END			1.2	
#define BLOOD_DENSITY		0.000106	// g/mm3	[wiki]
#define BLOOD_CAPACITY		3.617		// J/g?C	[http://www.itis.ethz.ch/virtual-population/tissue-properties/database/heat-capacity/]

typedef struct medium_struct_heat{
    float energy[VOXELS_X][VOXELS_Y][VOXELS_Z];		//	matrix with absorbed energy
    float structure[VOXELS_X][VOXELS_Y][VOXELS_Z];
    float k[MAX_REGIONS];					// heat conduction coeficient
	float w_g[MAX_REGIONS];    	
    float temperature[VOXELS_X][VOXELS_Y][VOXELS_Z];
	float arterial_temperature;
    int timeSelection;
}m_str_h;

__kernel void PennesEquation(__global m_str_h *m_str_h)
{
   
	int temp = get_global_id(0) + 1;
	int temp2 = get_global_id(1) + 1;
	int temp3 = get_global_id(2) + 1;

	int tempId = m_str_h[0].structure[temp][temp2][temp3];

		m_str_h[0].temperature[temp][temp2][temp3] = (m_str_h[0].temperature[temp + 1][temp2][temp3] + m_str_h[0].temperature[temp - 1][temp2][temp3] + 
			m_str_h[0].temperature[temp][temp2 + 1][temp3] + m_str_h[0].temperature[temp][temp2 - 1][temp3] + m_str_h[0].temperature[temp][temp2][temp3 + 1] + 
			m_str_h[0].temperature[temp][temp2][temp3 - 1] + (m_str_h[0].energy[temp][temp2][temp3] / m_str_h[0].k[tempId]) + 
			m_str_h[0].w_g[tempId] * BLOOD_CAPACITY * BLOOD_DENSITY * (m_str_h[0].arterial_temperature - m_str_h[0].temperature[temp][temp2][temp3]) / m_str_h[0].k[tempId]) / 6.0f;
}

__kernel void PennesEquation_Time(__global m_str_h *m_str_h)  // just stub
{
   
	int temp = get_global_id(0) + 1;
	int temp2 = get_global_id(1) + 1;
	int temp3 = get_global_id(2) + 1;

	int tempId = m_str_h[0].structure[temp][temp2][temp3];

		m_str_h[0].temperature[temp][temp2][temp3] = (m_str_h[0].temperature[temp + 1][temp2][temp3] + m_str_h[0].temperature[temp - 1][temp2][temp3] + 
			m_str_h[0].temperature[temp][temp2 + 1][temp3] + m_str_h[0].temperature[temp][temp2 - 1][temp3] + m_str_h[0].temperature[temp][temp2][temp3 + 1] + 
			m_str_h[0].temperature[temp][temp2][temp3 - 1] + (m_str_h[0].energy[temp][temp2][temp3] / m_str_h[0].k[tempId]) + 
			m_str_h[0].w_g[tempId] * BLOOD_CAPACITY * BLOOD_DENSITY * (m_str_h[0].arterial_temperature - m_str_h[0].temperature[temp][temp2][temp3]) / m_str_h[0].k[tempId]) / 6.0f;
}