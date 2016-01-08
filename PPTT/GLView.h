/**********************************************************/
/*					    PPTT 0.1						  */
/*				Created by Jakub Mesicek				  */
/*						10/2015							  */
/**********************************************************/
/*														  */
/*	GLView defines Result viewer                          */
/**********************************************************/


#ifndef GLVIEW_H
#define	GLVIEW_H
#include <cmath>
#include <string>
#include "PPTT_core.h"
#include "Camera.h"
#include <GL\glut.h>

static int stepCounter = 0; // time
static int fullscreen = 0;
static bool showboundary = true;
static float adjustSize = 1;
static int viewSelector = 0; // selected view
static int xOrigin = -1;
static bool viewChanged = true;
static GLuint voxelsList;


// slice view parameters
static int sliceX = 0;
static int sliceY = 0;
static int sliceZ = 0;

static Medium * m_draw;
static Heat * h_draw;
static Camera camera;

class GLView;
class GLView
{
public:
    void Run();
    void Init(int argc, char ** argv);
    void SaveMedium(Medium *m, Heat *h,int viewID);
};


void Draw();
void DrawHelp(std::string s, float x, float y);
void DrawColorScale();
void DrawAxis();
void HandleResize(int w, int h);
void ProcessSpecialKeys(int key, int xx, int yy);
void ProcessNormalKeys(unsigned char key, int x, int y);
void MouseButton(int button, int state, int x, int y);
void MouseMove(int x, int y);
void GetColor(float t);
void CreateDisplayList();
void DrawDisplayList();
#endif	/* GLVIEW_H */