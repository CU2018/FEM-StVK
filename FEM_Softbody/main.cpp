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
int maxFrames = 100;  // 240
int maxSubstep = 10;  // 10
float FPS = 24.0f;  // 24.0f
bool debug = false;
bool enableStaticAnal = false;

// init FEMSimObj
float deltaTime = 1.0f / (FPS*maxSubstep); //1.0/240f;
float dampingCoef = 0.9f;  // damping coefficient 0.9f
float gravityConst = 9.8f;
float restitutionCoef = 1.0f;
float frictionCoef = 0.2f;
float lsAlpha = 0.03f;
float lsBeta = 0.5f;
unsigned int iterationNum = 200; // num of gradient descent iteration

// material info
MaterialType matType = MATERIAL_TYPE_StVK;
ScalarType matMu = 2e6; // Lame's second param - resistance for shearing
ScalarType matLambda = 8e6; // Lame's first param - resistance for change of volume
ScalarType unitMass = 0.01f;

TetMesh* gMesh;
TetMesh* predMesh;
FEMSimObj * gFEMSimObj;

// read input files
std::string oneTetFilePath = "oneTet/oneTet.tet";
std::string manyTetsFilePath = "oneTet/manyTets.tet";
std::string pigHeadFilePath = "oneTet/pigHead.tet";
std::string staticRestFilePath = "oneTet/staticRest.tet";
std::string staticPredFilePath = "oneTet/staticPred.tet";

int main()
{	
	// load the mesh info from the file
	gMesh = new TetMesh();
	gMesh->init(pigHeadFilePath, unitMass);

	if (enableStaticAnal)
	{
		predMesh = new TetMesh();
		predMesh->init(staticPredFilePath, unitMass);
	}
	printf("---------------------------------------------\n");
	// init FEMSimObj: pass in necessary variables
	gFEMSimObj = new FEMSimObj();
	gFEMSimObj->init(gMesh, deltaTime, gravityConst, dampingCoef, restitutionCoef, frictionCoef,
		lsAlpha, lsBeta, iterationNum, maxSubstep);
	gFEMSimObj->initMatInfo(matType, matMu, matLambda);
	// for debugging
	gFEMSimObj->setDebug(debug);
	gFEMSimObj->staticAnalysisLog = enableStaticAnal;
	if (enableStaticAnal) gFEMSimObj->setPredMeshPos(predMesh);
	
	// simulation loop
	for (int frameNum = 1; frameNum <= maxFrames; ++frameNum)
	{
		if (frameNum == 1) gFEMSimObj->isFirstFrame = true;
		else gFEMSimObj->isFirstFrame = false;
		// update
		gFEMSimObj->update();
		if(gFEMSimObj->saveTetAsHDA(frameNum))
			printf("Successfully saved frame %d!\n", frameNum);
	}
	return 0;
}
