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

}my_struct;

__kernel void computePhoton(__global my_struct *myStruct)
{
    int gid = get_global_id(0);
    myStruct[gid].a = gid;
    myStruct[gid].b = gid + 100;
    myStruct[gid].c = gid + 2;
    myStruct[gid].energy[0][0][0] = 60;
}
