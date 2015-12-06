#include "GLview.h"
#include "PPTT_core.h"
#include <GL\glut.h>
#include <sstream>
#include <iostream>
#include <string>

void GLView::run()
{
    glutDisplayFunc(draw);
    glutReshapeFunc(handleResize);
    glutIdleFunc(draw);
    glutKeyboardFunc(processNormalKeys);
    glutSpecialFunc(processSpecialKeys);
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMove);
    glutMainLoop(); //Start the main loop. glutMainLoop doesn't return.
}

void GLView::init(int argc, char** argv)
{
    //Initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glShadeModel(GL_SMOOTH);
    glutInitWindowSize(1200, 800); //Window size
    glutCreateWindow("PPTT - matrix view"); //Create a window
    glEnable(GL_DEPTH_TEST); //Make sure 3D drawing works when one object is in front of another
}

void GLView::savemedium(Medium * m)
{
    m_draw = m;
}

//Called when the window is resized
void handleResize(int w, int h)
{
    //Tell OpenGL how to convert from coordinates to pixel values
    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION); //Switch to setting the camera perspective

                                 //Set the camera perspective
    glLoadIdentity(); //Reset the camera
    gluPerspective(45.0,				  //The camera angle
        (double)w / (double)h, //The width-to-height ratio
        1.0,				   //The near z clipping coordinate
        2000.0);				//The far z clipping coordinate
}

void processSpecialKeys(int key, int xx, int yy) 
{
    float fraction = 10.0f;

    switch(key) {
        case GLUT_KEY_LEFT:
        angle -= 0.05f;
        lx = sin(angle);
        lz = -cos(angle);
        break;
        case GLUT_KEY_RIGHT:
        angle += 0.05f;
        lx = sin(angle);
        lz = -cos(angle);
        break;
        case GLUT_KEY_UP:
        x += lx * fraction;
        z += lz * fraction;
        break;
        case GLUT_KEY_DOWN:
        x -= lx * fraction;
        z -= lz * fraction;
        break;
        case GLUT_KEY_F1:
        stepCounter = stepCounter + 1;
        if(stepCounter > timeSegments)
            stepCounter = 0;
        break;
        case GLUT_KEY_F2:
        stepCounter = stepCounter - 1;
        if(stepCounter < 0)
            stepCounter = timeSegments;
        if(stepCounter > timeSegments)
            stepCounter = 0;
        break;
        case GLUT_KEY_F3:
        adjustSize += 10;
        break;
        case GLUT_KEY_F4:
        {

            if(adjustSize <= 1)
            {
                adjustSize -= 0.1;
            }
            else
            {
                adjustSize -= 0.5;
            }
            if(adjustSize <= 0)
                adjustSize = 0;
            break;
        }
        case GLUT_KEY_F5:
        adjustSize = 1;
        break;
        case GLUT_KEY_F6:
        selector++;
        if(selector > 1)
            selector = 0;
        break;
        case GLUT_KEY_F7:
        if(showboundary)
               showboundary = false;
        else
            showboundary = true;
        break;

    }
}

void mouseButton(int button, int state, int x, int y) 
{
    // only start motion if the left button is pressed
    if(button == GLUT_LEFT_BUTTON) {

        // when the button is released
        if(state == GLUT_UP) {
            angle += deltaAngle;
            xOrigin = -1;
        }
        else {// state = GLUT_DOWN
            xOrigin = x;
        }
    }
}

void mouseMove(int x, int y)
{

    // this will only be true when the left button is down
    if(xOrigin >= 0) {

        // update deltaAngle
        deltaAngle = (x - xOrigin) * 0.001f;

        // update camera's direction
        lx = sin(angle + deltaAngle);
        lz = -cos(angle + deltaAngle);
    }
}

void processNormalKeys(unsigned char key, int x, int y)
{
    float fraction = 10.0f;
    switch(key) {
        case 'a':
        angle -= 0.05f;
        lx = sin(angle);
        lz = -cos(angle);
        break;
        case 'd':
        angle += 0.05f;
        lx = sin(angle);
        lz = -cos(angle);
        break;
        case GLUT_KEY_UP:
        case 'w':
        x += lx * fraction;
        z += lz * fraction;
        break;
        case 's':
        x -= lx * fraction;
        z -= lz * fraction;
        break;
}
}

//Draws the 3D scene
void draw()
{
    //Clear screen
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW); //Switch to the drawing perspective
    glLoadIdentity(); //Reset the drawing perspective
    gluLookAt(x, y, z, x + lx, y + ly, z + lz, 0.0f, 1.0f, 0.0f);
    glPushMatrix();
    glColor3ub(255, 0, 0);
    glTranslatef(0, 0, 0);
    glutSolidCube(1);
    glTranslatef(voxels_x / 2 , voxels_y / 2, voxels_z / 2);
    drawBounds(voxels_x / 2);
    glPopMatrix();

    switch(selector)
    {
        case 0:
        for(int temp1 = 0 + showboundary; temp1 < voxels_x - showboundary; temp1++)
        {
            for(int temp2 = 0 + showboundary; temp2 < voxels_y - showboundary; temp2++)
            {
                for(int temp3 = 0 + showboundary; temp3 < voxels_z - showboundary; temp3++)
                {
                    float size = m_draw->energy_t[temp1][temp2][temp3][stepCounter] * adjustSize;

                    if(size != 0)
                    {
                        if(size > 20)
                            size = 20;
                        glPushMatrix();
                        glColor3ub(128 + size * 50, 128 + size * 50, 128);
                        glTranslatef(temp1, temp2, temp3);
                        glutSolidCube(size / 5);
                        glPopMatrix();
                    }


                }
            }
        }
        break;
        case 1:
        for(int temp1 = 0; temp1 < voxels_x; temp1++)
        {
            for(int temp2 = 0; temp2 < voxels_y; temp2++)
            {
                for(int temp3 = 0; temp3 < voxels_z; temp3++)
                {
                    float color2 = m_draw->structure[temp1][temp2][temp3];
                    if(color2 > 1.0)
                    {
                        glPushMatrix();
                        if(color2 > 2.0)
                            glColor3ub(255, 0,0);
                        else
                        {
                            glColor3ub(0, 255, 0);
                        }
                        glTranslatef(temp1, temp2, temp3);
                        glutSolidCube(color2 / 5);
                        glPopMatrix();
                    }
                }
            }
        }
        break;
    }
    drawHelp("Camera control: WSAD + Mouse", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT)-30);
    drawHelp("Simulation control:", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 50);
    drawHelp("F1 - forward", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 70);
    drawHelp("F2 - back", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 90);
    drawHelp("F3 - increase energy size", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 110);
    drawHelp("F4 - decrease energy size", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 130);
    drawHelp("F5 - reset energy size", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 150);
    drawHelp("F6 - show medium properties", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 170);
    drawHelp("F7 - show boundary energy matrix", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 190);
    drawHelp("Simulation step counter: " + std::to_string(stepCounter) + "/" + std::to_string(timeSegments), glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 250);
    drawHelp("Voxel Energy multiplicator: " + std::to_string(adjustSize), glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 270);
    drawHelp("Time per step: " + std::to_string(timeStep) + "ms", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 290);
    glEnd();
    glutSwapBuffers(); //Send scene to the screen to be shown
}


void drawHelp(std::string s, float x, float y) 
{
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0.0, glutGet(GLUT_WINDOW_WIDTH), 0.0, glutGet(GLUT_WINDOW_HEIGHT));
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glColor3f(1.0f, 0.0f, 0.0f);//needs to be called before RasterPos
    glRasterPos2i(x - 300, y);
    for(int i = 0; i < s.length(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, s[i]);
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glEnable(GL_TEXTURE_2D);

    glutPostRedisplay();
}

void drawBounds(int myVar) 
{

    glutWireCube(myVar);
}