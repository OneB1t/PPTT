#include "GLview.h"
#include "PPTT_core.h"
#include <GL\glut.h>

void GLView::run()
{
    glutDisplayFunc(draw);
    glutReshapeFunc(handleResize);
    glutIdleFunc(draw);
    glutKeyboardFunc(processNormalKeys);
    glutSpecialFunc(processSpecialKeys);
    glutMainLoop(); //Start the main loop. glutMainLoop doesn't return.
}

void GLView::init(int argc, char** argv)
{
    //Initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(800, 800); //Window size
    glutCreateWindow("Introduction to OpenGL"); //Create a window
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
        200.0);				//The far z clipping coordinate
}

void processSpecialKeys(int key, int xx, int yy) {

    float fraction = 0.1f;

    switch(key) {
        case GLUT_KEY_LEFT:
        angle -= 0.01f;
        lx = sin(angle);
        lz = -cos(angle);
        break;
        case GLUT_KEY_RIGHT:
        angle += 0.01f;
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
    }
}

void processNormalKeys(unsigned char key, int x, int y) 
{
}

//Draws the 3D scene
void draw()
{
        if(counter > 6)
        counter = 0;
    //Clear screen
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW); //Switch to the drawing perspective
    glLoadIdentity(); //Reset the drawing perspective
    gluLookAt(x, 1.0f, z,
        x + lx, 1.0f, z + lz,
        0.0f, 1.0f, 0.0f);


    glBegin(GL_POINTS); //Begin drawing points

    for(int temp1 = 0; temp1 < voxels_x; temp1++)
    {
        for(int temp2 = 0; temp2 < voxels_y; temp2++)
        {
            for(int temp3 = 0; temp3 < voxels_z; temp3++)
            {
                float color = m_draw->energy_t[temp1][temp2][temp3][0] * 10;
                float color2 = m_draw->structure[temp1][temp2][temp3] * 1000;
                glColor3f(color, color, color);
                glVertex3f((-temp1 + 50), (-temp2 + 50), (-temp3 - 10));
                glColor3f(color2, color2, color2);
                glVertex3f((-temp1 + 50), (-temp2 + 50), (-temp3 - 10));

            }
        }
    }

    glEnd();
    glutSwapBuffers(); //Send scene to the screen to be shown
    counter++;
}