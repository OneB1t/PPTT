const int voxels_x = 100;
const int voxels_y = 100;
const int voxels_z = 100;
const int max_regions = 16;

typedef struct tag_my_struct{
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
}my_struct;

__kernel void computePhoton(__global my_struct *myStruct)
{
    int gid = get_global_id(0);
    myStruct[gid].a = gid + 100;
    myStruct[gid].b = gid + 2;
    myStruct[gid].c = gid + 200;
    myStruct[gid].energy[0][0][0] = 60;
    myStruct[gid].n[0] = myStruct[gid].n[0] + 1.1;
}
