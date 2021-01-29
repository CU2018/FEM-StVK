#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <string>
#include <fstream> 
#include <iomanip> 
#include <iostream> 

#include "MathHeaders.h"
#include "Mesh.h"
#include "MeshLoader.h"
#include "FEMSimObj.h"


// simulation settings
int maxFrames = 240;  // 240
int maxSubstep = 10;  // 10
float FPS = 24.0f;  // 24.0f
int gdIteration = 10;  // uum of gradient descent iteration

// init FEMSimObj
float timeStep = 1.0f / (FPS*maxSubstep); //1.0/240f;
float dampingCoef = 0.9f;  // damping coefficient 0.9f
float gravityConst = 9.8f;
float restitutionCoef = 1.0;
float frictionCoef = 0.1;
float lsAlpha = 0.03;
float lsBeta = 0.5;
unsigned int iterationNum = 100;

const int meshTotalVertices = 6;  // 2886 for the pig head

TetMesh* gMesh;
FEMSimObj * gFEMSimObj;

// read input file
std::string filePath = "oneTet/oneTetRest.tet";

int main()
{
	// MeshLoader meshLoader(filePath);
	gMesh = new TetMesh();
	gMesh->init(filePath);  // load the mesh info from the file
	gFEMSimObj = new FEMSimObj();
	// init FEMSimObj: pass in necessary variables
	gFEMSimObj->init(gMesh, timeStep, gravityConst, dampingCoef, restitutionCoef, frictionCoef, 
					lsAlpha, lsBeta, iterationNum, maxSubstep);

	for (int frameNum = 1; frameNum <= maxFrames; ++frameNum)
	{
		// update
		gFEMSimObj->update();
		if(gFEMSimObj->saveTetAsHDA(frameNum))
			printf("saving %d sucess!\n", frameNum);
	}
	return 0;
}
