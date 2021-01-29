#include "Mesh.h"

Mesh::~Mesh()
{
	edgesList.clear();
}

bool TetMesh::init(std::string fileName)
{
	loadedMesh = new MeshLoader(fileName);
	if (loadedMesh->info() == false)
	{
		std::cout << "Load mesh error. Using regular Mesh." << std::endl;
		delete loadedMesh;

		return false;
	}
	this->totalMass = 4.0f;
	generateParticleList();
	return true;
}

bool TetMesh::exportToFile(int frameNum)
{
	std::ofstream outfile;
	std::string fileName = "tet." + std::to_string(frameNum) + ".tet";
	outfile.open(fileName, std::ifstream::out);
	if (outfile.is_open())
	{
		outfile << "plist:";
		for (unsigned int i = 0; i < verticesNum; i++)
		{
			// save positions
			outfile << currPos.block_vector(i)[0] << "," << currPos.block_vector(i)[1] << "," << currPos.block_vector(i)[2] << ",";
		}
		outfile << "\nidxlist:";
		for (unsigned int i = 0; i < loadedMesh->tetsList.size(); i++)
		{
			// save positions
			outfile << loadedMesh->tetsList[i].id1 << "," << loadedMesh->tetsList[i].id2 << "," << loadedMesh->tetsList[i].id3 << "," << loadedMesh->tetsList[i].id4 << ",";
		}
		outfile.close();
		return true;
	}
	return false;
}

void TetMesh::generateParticleList()
{
	verticesNum = loadedMesh->verticesList.size();
	systemDimension = verticesNum * 3;
	ScalarType unitMass = totalMass / systemDimension;
	// ScalarType unitMass = 0.05f;

	// Assign initial position, velocity and mass to all the vertices.
	restposePos.resize(systemDimension);
	currPos.resize(systemDimension);
	currVel.resize(systemDimension);
	massMat.resize(systemDimension, systemDimension);
	invMassMat.resize(systemDimension, systemDimension);

	//massMat1d.resize(verticesNum, verticesNum);
	//invMassMat1d.resize(verticesNum, verticesNum);

	// Assign initial position to all the vertices.
	restposePos.setZero();
	unsigned int index;
	for (index = 0; index < verticesNum; ++index)
	{
		restposePos.block_vector(index) = GLM2Eigen(loadedMesh->verticesList[index]);
	}
	assert(restposePos.size() == loadedMesh->verticesList.size() * 3);
	// initialize currVel, currPos, prevPos, prevVel
	currVel.setZero();
	currPos = restposePos;
	prevPos = restposePos;
	prevVel = currVel;

	// initialize mass matrix and an equally sized identity matrix
	std::vector<SparseMatrixTriplet> iTriplets;
	std::vector<SparseMatrixTriplet> mTriplets;
	std::vector<SparseMatrixTriplet> mInvTriplets;
	iTriplets.clear();
	mTriplets.clear();
	ScalarType invUnitMass = 1.0 / unitMass;
	for (index = 0; index < systemDimension; index++)
	{
		iTriplets.push_back(SparseMatrixTriplet(index, index, 1));
		mTriplets.push_back(SparseMatrixTriplet(index, index, unitMass));
		mInvTriplets.push_back(SparseMatrixTriplet(index, index, invUnitMass));
	}
	massMat.setFromTriplets(mTriplets.begin(), mTriplets.end());
	invMassMat.setFromTriplets(mInvTriplets.begin(), mInvTriplets.end());
	// 1d matrices
	//mTriplets.clear();
	//mInvTriplets.clear();
	//iTriplets.clear();
	//for (index = 0; index < verticesNum; index++)
	//{
	//	iTriplets.push_back(SparseMatrixTriplet(index, index, 1));
	//	mTriplets.push_back(SparseMatrixTriplet(index, index, unitMass));
	//	mInvTriplets.push_back(SparseMatrixTriplet(index, index, invUnitMass));
	//}
	//massMat1d.setFromTriplets(mTriplets.begin(), mTriplets.end());
	//invMassMat1d.setFromTriplets(mInvTriplets.begin(), mInvTriplets.end());
	printf("Generate particle list done!\n");
}

void TetMesh::generateEdgeList()
{
	SparseMatrix edgeMatrix(verticesNum, verticesNum);
	edgeMatrix.setZero();

	unsigned int i1, i2;
	for (std::vector<MeshLoader::Tet>::iterator iter = loadedMesh->tetsList.begin(); iter != loadedMesh->tetsList.end(); ++iter)
	{
		// 1-2
		i1 = iter->id1;
		i2 = iter->id2;
		if (edgeMatrix.coeff(i1, i2) < EPSILON)
		{
			edgeMatrix.coeffRef(i1, i2) = 1;
			edgeMatrix.coeffRef(i2, i1) = 1;
			Edge newEdge;
			newEdge.v1 = i1;
			newEdge.v2 = i2;
			edgesList.push_back(newEdge);
		}
		// 1-3
		i1 = iter->id1;
		i2 = iter->id3;
		if (edgeMatrix.coeff(i1, i2) < EPSILON)
		{
			edgeMatrix.coeffRef(i1, i2) = 1;
			edgeMatrix.coeffRef(i2, i1) = 1;
			Edge newEdge;
			newEdge.v1 = i1;
			newEdge.v2 = i2;
			edgesList.push_back(newEdge);
		}
		// 1-4
		i1 = iter->id1;
		i2 = iter->id4;
		if (edgeMatrix.coeff(i1, i2) < EPSILON)
		{
			edgeMatrix.coeffRef(i1, i2) = 1;
			edgeMatrix.coeffRef(i2, i1) = 1;
			Edge newEdge;
			newEdge.v1 = i1;
			newEdge.v2 = i2;
			edgesList.push_back(newEdge);
		}
		// 2-3
		i1 = iter->id2;
		i2 = iter->id3;
		if (edgeMatrix.coeff(i1, i2) < EPSILON)
		{	
			edgeMatrix.coeffRef(i1, i2) = 1;
			edgeMatrix.coeffRef(i2, i1) = 1;
			Edge newEdge;
			newEdge.v1 = i1;
			newEdge.v2 = i2;
			edgesList.push_back(newEdge);
		}
		// 2-4
		i1 = iter->id2;
		i2 = iter->id4;
		if (edgeMatrix.coeff(i1, i2) < EPSILON)
		{	
			edgeMatrix.coeffRef(i1, i2) = 1;
			edgeMatrix.coeffRef(i2, i1) = 1;
			Edge newEdge;
			newEdge.v1 = i1;
			newEdge.v2 = i2;
			edgesList.push_back(newEdge);
		}
		// 3-4
		i1 = iter->id3;
		i2 = iter->id4;
		if (edgeMatrix.coeff(i1, i2) < EPSILON)
		{	
			edgeMatrix.coeffRef(i1, i2) = 1;
			edgeMatrix.coeffRef(i2, i1) = 1;
			Edge newEdge;
			newEdge.v1 = i1;
			newEdge.v2 = i2;
			edgesList.push_back(newEdge);
		}
	}
}