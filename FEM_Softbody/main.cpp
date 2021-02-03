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
int maxFrames = 16;  // 240
int maxSubstep = 2;  // 10
float FPS = 24.0f;  // 24.0f
int gdIteration = 500;  // num of gradient descent iteration

// init FEMSimObj
float timeStep = 1.0f / (FPS*maxSubstep); //1.0/240f;
float dampingCoef = 0.9f;  // damping coefficient 0.9f
float gravityConst = 9.8f;
float restitutionCoef = 1.0f;
float frictionCoef = 0.2f;
float lsAlpha = 0.03f;
float lsBeta = 0.5f;
unsigned int iterationNum = 100;

// material info
MaterialType matType = MATERIAL_TYPE_StVK;
ScalarType matMu = 1e4; // Lame's second param - resistance for shearing
ScalarType matLambda = 4e4; // Lame's first param - resistance for change of volume

TetMesh* gMesh;
FEMSimObj * gFEMSimObj;

// read input files
std::string oneTetFilePath = "oneTet/oneTet.tet";
std::string manyTetsFilePath = "oneTet/manyTets.tet";
std::string pigHeadFilePath = "oneTet/pigHead.tet";

int main()
{	
	// load the mesh info from the file
	gMesh = new TetMesh();
	gMesh->init(oneTetFilePath); 

	// init FEMSimObj: pass in necessary variables
	gFEMSimObj = new FEMSimObj();
	gFEMSimObj->init(gMesh, timeStep, gravityConst, dampingCoef, restitutionCoef, frictionCoef,
		lsAlpha, lsBeta, iterationNum, maxSubstep);
	gFEMSimObj->initMatInfo(matType, matMu, matLambda);

	// simulation loop
	for (int frameNum = 1; frameNum <= maxFrames; ++frameNum)
	{
		// update
		gFEMSimObj->update();
		if(gFEMSimObj->saveTetAsHDA(frameNum))
			printf("Successfully saved frame %d!\n", frameNum);
	}
	return 0;
}
