#define _MAIN_

#ifdef _MAIN_

#include "main.h"

//#include <headers/include/zCore.h>
//#include <headers/include/zGeometry.h>
//#include <headers/include/zDisplay.h>
//#include <headers/include/zData.h>
//#include <headers/include/zIO.h> 
//
//using namespace zSpace;

////////////////////////////////////////////////////////////////////////// General

// global variables ?
int gridSize = 25;// declare variables globall (outside setup, draw, keyPress etc) .. if you want to change or use it anywhere in the rest fo the code
//int stores numbers without decimals i.e integers
float decimalNUmber = 22.25;// float stores numbers with decimals
double morePrecise_decimalNumber = 22.2535353535353535353;// double stores numbers with decimals
string studentName = "aschawit"; // string stores text
bool compute = true; // bool stores true or false
char c = 'o';

string textToPrint;
////////////////////////////////////////////////////////////////////////// zSpace Objects
// vector field


void setup() // 
{

	// code inside setup runs once, when you press the green run button ;
	////////////////////////////////////////////////////////////////////////// Enable smooth display

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_POINT_SMOOTH);

	////////////////////////////////////////////////////////////////////////// Sliders


}

void update(int value) // code inside here runs 100 times a second
{

}

void draw() // code inside here runs 100 times a second
{
	// put only drawing code / commands here 

	backGround(0.8);
	drawGrid(gridSize);

	Alice::vec constantDir(0.15, 0.15, 0);
	Alice::vec AttractorPt(10, 10, 0); 


	// draw a grid of points 
	for (int i = 0; i < gridSize; i++)
	{
		for (int j = 0; j < gridSize; j++)
		{
			Alice::vec pt(i, j, 0);
			glPointSize(5);
			drawPoint(pt);

			Alice::vec dir = AttractorPt - pt;
			dir.normalise();

			drawLine(pt, pt + dir);
		}
	}

}

void keyPress(unsigned char k, int xm, int ym)
{
	// write code here, if you want it to run when you press a key


}



void mousePress(int b, int s, int x, int y)
{

}

void mouseMotion(int x, int y)
{

}



#endif // _MAIN_
