#include "GLview.h"
#include "PPTT_core.h"
#include "CL\CL.h"
#include <GL\glut.h>
#include <sstream>
#include <iostream>
#include <string>

void GLView::Run()
{
    glutDisplayFunc(Draw);
    glutReshapeFunc(HandleResize);
    glutIdleFunc(Draw);
    glutKeyboardFunc(ProcessNormalKeys);
    glutSpecialFunc(ProcessSpecialKeys);
    glutMouseFunc(MouseButton);
    glutMotionFunc(MouseMove);
    glutMainLoop();
}

void GLView::Init(int argc, char** argv)
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
    camera = Camera();
    DrawAxis();
}

void GLView::SaveMedium(Medium * m, Heat * h)
{
    m_draw = m;
    h_draw = h;
}

//Called when the window is resized
void HandleResize(int w, int h)
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

void ProcessSpecialKeys(int key, int xx, int yy)
{
    float fraction = 10.0f;

    switch(key) {
        case GLUT_KEY_F1:
        stepCounter = stepCounter + 1;
        if(stepCounter > timeSegments - 1)
            stepCounter = 0;
        break;
        case GLUT_KEY_F2:
        stepCounter = stepCounter - 1;
        if(stepCounter < 0)
            stepCounter = timeSegments - 1;
        if(stepCounter > timeSegments - 1)
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
        viewSelector++;
        if(viewSelector > 7)
            viewSelector = 0;
        break;
        case GLUT_KEY_F7:
        if(showboundary)
            showboundary = false;
        else
            showboundary = true;
        break;
        case GLUT_KEY_F8:
        fullscreen++;
        if(fullscreen > 1)
            fullscreen = 0;
        if(fullscreen)
            glutFullScreen();
        else
        {
            glutReshapeWindow(1200, 800);        /* Restore us */
            glutPositionWindow(0, 30);
        }
        break;

    }
    viewChanged = true;
}

void MouseButton(int button, int state, int x, int y)
{
    // only start motion if the left button is pressed
    if(button == GLUT_LEFT_BUTTON) {
        // when the button is released
        if(state == GLUT_UP) {
            xOrigin = -1;
        }
        else {// state = GLUT_DOWN
            xOrigin = x;
            camera.ox = x;
            camera.oy = y;
        }
    }
}

void MouseMove(int x, int y)
{

    // this will only be true when the left button is down
    if(xOrigin >= 0) {
        camera.addAzimuth(0.1 * PI * (camera.ox - x) / GLUT_WINDOW_WIDTH);
        camera.addZenith(-0.1 * PI * (y - camera.oy) / GLUT_WINDOW_WIDTH);
        camera.ox = x;
        camera.oy = y;
    }
}

void ProcessNormalKeys(unsigned char key, int x, int y)
{
    float fraction = 10.0f;
    switch(key) {
        case 'a':
        camera.left(5);
        break;
        case 'd':
        camera.right(5);
        break;
        case 'w':
        camera.move(camera.to);
        break;
        case 's':
        camera.moveback(camera.to);
        break;
        case 'q':
        camera.down(5);
        break;
        case 'e':
        camera.up(5);
        break;
        case 'u':
        sliceX++;
        if(sliceX > voxelsX - 1)
            sliceX = 0;
        viewChanged = true;
        break;
        case 'i':
        sliceX--;
        if(sliceX < 0)
            sliceX = voxelsX - 1;
        viewChanged = true;
        break;
        case 'j':
        sliceY++;
        if(sliceY > voxelsY - 1)
            sliceY = 0;
        viewChanged = true;
        break;
        case 'k':
        sliceY--;
        if(sliceY < 0)
            sliceY = voxelsY - 1;
        viewChanged = true;
        break;
        case 'n':
        sliceZ++;
        if(sliceZ > voxelsZ - 1)
            sliceZ = 0;
        viewChanged = true;
        break;
        case 'm':
        sliceZ--;
        if(sliceZ < 0)
            sliceZ = voxelsZ - 1;
        viewChanged = true;
        break;
    }
}

//Draws the 3D scene
void Draw()
{
    //Clear screen
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW); //Switch to the drawing perspective
    glLoadIdentity(); //Reset the drawing perspective   
    gluLookAt(camera.from.x, camera.from.y, camera.from.z, camera.from.x + camera.to.x, camera.from.y + camera.to.y, camera.from.z + camera.to.z, camera.upV.x, camera.upV.y, camera.upV.z);
    DrawAxis();
    if(viewChanged)
    {
        CreateDisplayList();
        viewChanged = false;
    }
    else
    {
        DrawDisplayList();
    }

    DrawHelp("Camera control: WSAD + Mouse", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 30);
    DrawHelp("Simulation control:", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 50);
    DrawHelp("F1 - forward", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 70);
    DrawHelp("F2 - back", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 90);
    DrawHelp("F3 - increase energy size", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 110);
    DrawHelp("F4 - decrease energy size", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 130);
    DrawHelp("F5 - reset energy size", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 150);
    DrawHelp("F6 - show medium properties", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 170);
    DrawHelp("F7 - show boundary energy matrix", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 190);
    DrawHelp("Simulation step counter: " + std::to_string(stepCounter + 1) + "/" + std::to_string(timeSegments), glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 250);
    DrawHelp("Voxel Energy multiplicator: " + std::to_string(adjustSize), glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 270);
    DrawHelp("Time per step: " + std::to_string(timeStep) + "ms", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 290);
    switch(viewSelector)
    {
        case 0:
        DrawHelp("View mode: Basic", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 310);
        break;
        case 1:
        DrawHelp("View mode: Medium properties", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 310);
        break;
        case 2:
        DrawHelp("View mode: X axis slice", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 310);
        DrawHelp("sliceX: " + std::to_string(sliceX), glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 330);
        break;
        case 3:
        DrawHelp("View mode: Y axis slice", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 310);
        DrawHelp("sliceY: " + std::to_string(sliceY), glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 330);
        break;
        case 4:
        DrawHelp("View mode: Z axis slice", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 310);
        DrawHelp("sliceZ: " + std::to_string(sliceZ), glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 330);
        break;
        case 5:
        DrawHelp("View mode: Temperature - Z slice", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 310);
        DrawHelp("sliceX: " + std::to_string(sliceX), glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 330);
        DrawHelp("sliceY: " + std::to_string(sliceY), glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 350);
        DrawHelp("sliceZ: " + std::to_string(sliceZ), glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 370);
        break;
        case 6:
        DrawHelp("View mode: Energy for whole simulation", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 310);
        break;
        case 7:
        DrawHelp("View mode: Escaped energy", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 310);
        break;
    }
    DrawColorScale("", glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) - 500);
    glEnd();
    glutSwapBuffers(); //Send scene to the screen to be shown
}


void DrawHelp(std::string s, float x, float y)
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

void DrawColorScale(std::string s, float x, float y)
{
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0.0, glutGet(GLUT_WINDOW_WIDTH), 0.0, glutGet(GLUT_WINDOW_HEIGHT));
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    int colorScale = glutGet(GLUT_WINDOW_HEIGHT) / 255;
    for(int i = 0; i < glutGet(GLUT_WINDOW_HEIGHT); i++)
    {
        GetColor((i * colorScale * 10) / glutGet(GLUT_WINDOW_HEIGHT));
        glBegin(GL_LINES);
        for(int k = 0; k < colorScale; k++)
        {
            glVertex2i(0, k + i);
            glVertex2i(100, k + i);
        }

        glEnd();
    }

    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glEnable(GL_TEXTURE_2D);
    glutPostRedisplay();
}


void DrawAxis()
{
    /* Create a display list for drawing axes */
    GLuint axes_list = glGenLists(1);
    glNewList(axes_list, GL_COMPILE_AND_EXECUTE);

    glColor4ub(0, 0, 255, 255);
    glBegin(GL_LINE_STRIP);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(1.0f, 0.0f, 0.0f);
    glVertex3f(0.75f, 0.25f, 0.0f);
    glVertex3f(0.75f, -0.25f, 0.0f);
    glVertex3f(1.0f, 0.0f, 0.0f);
    glVertex3f(0.75f, 0.0f, 0.25f);
    glVertex3f(0.75f, 0.0f, -0.25f);
    glVertex3f(1.0f, 0.0f, 0.0f);
    glEnd();
    glBegin(GL_LINE_STRIP);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 1.0f, 0.0f);
    glVertex3f(0.0f, 0.75f, 0.25f);
    glVertex3f(0.0f, 0.75f, -0.25f);
    glVertex3f(0.0f, 1.0f, 0.0f);
    glVertex3f(0.25f, 0.75f, 0.0f);
    glVertex3f(-0.25f, 0.75f, 0.0f);
    glVertex3f(0.0f, 1.0f, 0.0f);
    glEnd();
    glBegin(GL_LINE_STRIP);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.25f, 0.0f, 0.75f);
    glVertex3f(-0.25f, 0.0f, 0.75f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.25f, 0.75f);
    glVertex3f(0.0f, -0.25f, 0.75f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glEnd();

    glColor4ub(255, 255, 0, 255);
    glRasterPos3f(1.1f, 0.0f, 0.0f);

    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 'x');
    glRasterPos3f(0.0f, 1.1f, 0.0f);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 'y');
    glRasterPos3f(0.0f, 0.0f, 1.1f);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 'z');

    glEndList();

}

void GetColor(float t)
{
    float invMaxH = 1.0f / (255);
    float zRel = (t * 10) * invMaxH;
    float cR = 0, cG = 0, cB = 0;
    if(0 <= zRel && zRel < 0.2f)
    {
        cB = 1.0f;
        cG = zRel * 5.0f;
    }
    else if(0.2f <= zRel && zRel < 0.4f)
    {
        cG = 1.0f;
        cB = 1.0f - (zRel - 0.2f) * 5.0f;
    }
    else if(0.4f <= zRel && zRel < 0.6f)
    {
        cG = 1.0f;
        cR = (zRel - 0.4f) * 5.0f;
    }
    else if(0.6f <= zRel && zRel < 0.8f)
    {
        cR = 1.0f;
        cG = 1.0f - (zRel - 0.6f) * 5.0f;
    }
    else
    {
        cR = 1.0f;
        cG = (zRel - 0.8f) * 5.0f;
        cB = cG;
    }
    glColor3f(cR, cG, cB);
}

void CreateDisplayList()
{
    /* Create a display list for drawing axes */
    voxelsList = glGenLists(1);
    glNewList(voxelsList, GL_COMPILE_AND_EXECUTE);
    switch(viewSelector)
    {
        case 0:
        for(int x = 0; x < voxelsX; x++)
        {
            for(int y = 0; y < voxelsY; y++)
            {
                for(int z = 0; z < voxelsZ; z++)
                {
                    float size = m_draw->energy_t[x][y][z][stepCounter] * adjustSize;

                    if(size >= 0.1f)
                    {
                        glBegin(GL_POINTS);
                        if(size > 20)
                            size = 20;
                        GetColor(size);
                        glPointSize(size);
                        glVertex3f(x, y, z);
                        glEnd();
                    }
                }
            }
        }

        break;
        case 1:

        for(int temp1 = 0; temp1 < voxelsX; temp1++)
        {
            for(int temp2 = 0; temp2 < voxelsY; temp2++)
            {
                for(int temp3 = 0; temp3 < voxelsZ; temp3++)
                {
                    float color2 = m_draw->structure[temp1][temp2][temp3];
                    if(color2 > 5.0)
                    {
                        glBegin(GL_POINTS);
                        GetColor(color2);
                        glPointSize(15);
                        glVertex3f(temp1, temp2, temp3);
                        glEnd();
                    }
                }
            }
        }

        break;
        case 2: // sliceX

        for(int temp2 = 0; temp2 < voxelsY; temp2++)
        {
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
            {
                float size = m_draw->energy_t[sliceX][temp2][temp3][stepCounter] * adjustSize;

                if(size != 0)
                {
                    glBegin(GL_POINTS);
                    GetColor(size);
                    glPointSize(size * 5);
                    glVertex3f(sliceX, temp2, temp3);
                    glEnd();
                }


            }
        }
        break;
        case 3: // sliceY

        for(int temp1 = 0; temp1 < voxelsX; temp1++)
        {
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
            {
                float size = m_draw->energy_t[temp1][sliceY][temp3][stepCounter] * adjustSize;

                if(size != 0)
                {
                    glBegin(GL_POINTS);
                    GetColor(size);
                    glPointSize(size * 5);
                    glVertex3f(temp1, sliceY, temp3);
                    glEnd();
                }
            }
        }

        break;
        case 4: // sliceZ

        for(int temp1 = 0; temp1 < voxelsX; temp1++)
        {
            for(int temp2 = 0; temp2 < voxelsY; temp2++)
            {
                float size = m_draw->energy_t[temp1][temp2][sliceZ][stepCounter] * adjustSize;

                if(size != 0)
                {
                    glBegin(GL_POINTS);
                    GetColor(size);
                    glPointSize(size * 5);
                    glVertex3f(temp1, temp2, sliceZ);
                    glEnd();
                }
            }
        }

        break;
        case 5: // temperature
        for(int temp1 = 0; temp1 < voxelsX; temp1++)
        {
            for(int temp2 = 0; temp2 < voxelsY; temp2++)
            {
                float size = h_draw->temperature[temp1][temp2][sliceZ] * adjustSize - 36;
                float size2 = h_draw->temperature[temp1][sliceY][temp2] * adjustSize - 36;
                float size3 = h_draw->temperature[sliceX][temp1][temp2] * adjustSize - 36;

                if(size >= 0.1f)
                {
                    glBegin(GL_POINTS);
                    if(size > 20)
                        size = 20;
                    GetColor(size);
                    glPointSize(15);

                    glVertex3f(temp1, temp2, sliceZ);
                    GetColor(size2);
                    glVertex3f(temp1, sliceY, temp2);
                    GetColor(size3);
                    glVertex3f(sliceX, temp1, temp2);
                    glEnd();
                }


            }
        }
        break;
        case 6: // draw energy for steady simulation

        for(int x = 0; x < voxelsX; x++)
        {
            for(int y = 0; y < voxelsY; y++)
            {
                for(int z = 0; z < voxelsZ; z++)
                {
                    float size = m_draw->energy[x][y][z] * adjustSize;

                    if(size >= 0.1f)
                    {
                        glBegin(GL_POINTS);
                        if(size > 20)
                            size = 20;
                        GetColor(size);
                        glPointSize(size);

                        glVertex3f(x, y, z);
                        glEnd();
                    }


                }
            }
        }
        break;
        case 7: // boundary matrix

        for(int side = 0; side < 2; side++)
        {
            for(int temp1 = 0; temp1 < voxelsX; temp1++)
            {
                for(int temp2 = 0; temp2 < voxelsY; temp2++)
                {
                    float size = m_draw->surroundingX[temp1][temp2][side] * adjustSize;
                    glBegin(GL_POINTS);
                    if(size >= 0.1)
                    {
                        GetColor(size);
                        glPointSize(15);
                        glVertex3f(temp1, temp2, side * voxelsX);
                    }
                    size = m_draw->surroundingY[temp1][temp2][side] * adjustSize;

                    if(size >= 0.1)
                    {
                        GetColor(size);
                        glVertex3f(temp1, side * voxelsY, temp2);

                    }
                    size = m_draw->surroundingZ[temp1][temp2][side] * adjustSize;

                    if(size >= 0.1)
                    {
                        GetColor(size);
                        glVertex3f(side * voxelsZ, temp1, temp2);
                    }
                    glEnd();

                }
            }
        }

        break;
    }
    glEndList();
}

void DrawDisplayList()
{
    glCallList(voxelsList);
}