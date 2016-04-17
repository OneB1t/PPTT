/**********************************************************/
/*					    PPTT 0.1						  */
/*				Created by Jakub Mesicek				  */
/*						11/2015							  */
/**********************************************************/
/*														  */
/*	PPTT_io defines input/output functions				  */
/**********************************************************/

#include "PPTT_io.h"

void Introduction()
{
    cout << "PPPPPPPPPPPPPPPPP   PPPPPPPPPPPPPPPPP   TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" << endl;
    cout << "P::::::::::::::::P  P::::::::::::::::P  T:::::::::::::::::::::TT:::::::::::::::::::::T" << endl;
    cout << "P::::::PPPPPP:::::P P::::::PPPPPP:::::P T:::::::::::::::::::::TT:::::::::::::::::::::T" << endl;
    cout << "PP:::::P     P:::::PPP:::::P     P:::::PT:::::TT:::::::TT:::::TT:::::TT:::::::TT:::::T" << endl;
    cout << "  P::::P     P:::::P  P::::P     P:::::PTTTTTT  T:::::T  TTTTTTTTTTTT  T:::::T  TTTTTT" << endl;
    cout << "  P::::P     P:::::P  P::::P     P:::::P        T:::::T                T:::::T        " << endl;
    cout << "  P::::PPPPPP:::::P   P::::PPPPPP:::::P         T:::::T                T:::::T        " << endl;
    cout << "  P:::::::::::::PP    P:::::::::::::PP          T:::::T                T:::::T        " << endl;
    cout << "  P::::PPPPPPPPP      P::::PPPPPPPPP            T:::::T                T:::::T        " << endl;
    cout << "  P::::P              P::::P                    T:::::T                T:::::T        " << endl;
    cout << "PP::::::PP          PP::::::PP                TT:::::::TT            TT:::::::TT      " << endl;
    cout << "P::::::::P          P::::::::P                T:::::::::T            T:::::::::T      " << endl;
    cout << "P::::::::P          P::::::::P                T:::::::::T            T:::::::::T      " << endl;
    cout << "PPPPPPPPPP          PPPPPPPPPP                TTTTTTTTTTT            TTTTTTTTTTT       " << endl;
    cout << " " << endl;

    cout << "============================================================" << endl;
    cout << "          Welcome to Planning Photothermal Therapy          " << endl;
    cout << "                     simulation software		             " << endl;
    cout << "                                                            " << endl;
    cout << "                          Brought by	                     " << endl;
    cout << "            Faculty of Informatics and Management           " << endl;
    cout << "                 University of Hradec Kralove               " << endl;
    cout << "                                                            " << endl;
    cout << "                   jakub.mesicek(at)uhk.cz					 " << endl;
    cout << "                     jan(at)zdarsky.org   					 " << endl;
    cout << "============================================================" << endl;
    cout << "                                                            " << endl;
}

char ChooseSteadyOrTime()
{
    cout << "Select type of simulation and press enter: " << endl;
    cout << "	1 - Steady-state simulation" << endl;
    cout << "	2 - Time-resolved simulation" << endl;
    cout << "You have selected ";
    int selection;
    cin >> selection;
    cout << endl;
    return selection;
}

char OpenCLOrCPU()
{
    cout << "Select type of simulation and press enter: " << endl;
    cout << "	1 - OpenCL" << endl;
    cout << "	2 - CPU" << endl;
    cout << "You have selected ";
    int selection;
    cin >> selection;
    cout << endl;
    return selection;
}

char OpenCLPlatform()
{
    cout << "Select OpenCL platform ID: " << endl;
    cout << "You have selected ";
    int selection;
    cin >> selection;
    cout << endl;
    return selection;
}

char OpenCLDevice()
{
    cout << "Select OpenCL device ID: " << endl;
    cout << "You have selected ";
    int selection;
    cin >> selection;
    cout << endl;
    return selection;
}


long HowManyPhotons()
{
    cout << "How many photons should participate in simulation? " << endl;
    long selection;
    cin >> selection;
    return selection;
}

void WriteTemperature(Medium * m, Heat * h)
{
    ofstream myFile;
    myFile.open("Temperature.txt");
    for(int temp1 = 0; temp1 < voxelsX; temp1++)
    {
        for(int temp2 = 0; temp2 < voxelsY; temp2++)
        {
            for(int temp3 = 0; temp3 < voxelsZ; temp3++)
                myFile << h->temperature[temp1][temp2][temp3] << " ";
            myFile << endl;
        }
        myFile << endl;
    }
    myFile.close();
}