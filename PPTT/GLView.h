/**********************************************************/
/*					    PPTT 0.1						  */
/*				Created by Jakub Mesicek				  */
/*						10/2015							  */
/**********************************************************/
/*														  */
/*	GLView defines Result viewer                          */
/**********************************************************/

/* TODO list: -withdraw Region class, replace with matrix in Space region */

#ifndef GLVIEW_H
#define	GLVIEW_H
#include <cmath>
#include <string>
#include "PPTT_core.h"

static int counter = 0; // time
static bool showboundary = false;
static float adjustsize = 1;
static int selector = 0; // selected view
static float angle = 0.0;
static float lx = 0.0f, lz = -1.0f, ly = 0.0f;
static float x = 50.0f, z = 150.0f, y = 50.0f;
static float deltaAngle = 0.0f;
static int xOrigin = -1;
static Medium * m_draw;
class GLView;
class GLView
{
public:
    void run();
    void init(int argc, char ** argv);
    void savemedium(Medium *m);
};
void draw();
void drawHelp(std::string s,float x, float y, float z);
void colorPick(float intensity);
void handleResize(int w, int h);
void processSpecialKeys(int key, int xx, int yy);
void processNormalKeys(unsigned char key, int x, int y);
void mouseButton(int button, int state, int x, int y);
void mouseMove(int x, int y);
#endif	/* GLVIEW_H */