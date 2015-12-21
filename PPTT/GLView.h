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
    void run();
    void init(int argc, char ** argv);
    void savemedium(Medium *m, Heat *h);
};


void draw();
void drawHelp(std::string s, float x, float y);
void drawAxis();
void colorPick(float intensity);
void handleResize(int w, int h);
void processSpecialKeys(int key, int xx, int yy);
void processNormalKeys(unsigned char key, int x, int y);
void mouseButton(int button, int state, int x, int y);
void mouseMove(int x, int y);
void getColor(float t);
void createDisplayList();
void drawDisplayList();
#endif	/* GLVIEW_H */