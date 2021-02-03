#ifndef _MESH_H_
#define _MESH_H_

#include "MathHeaders.h"
#include <iostream>
#include "MeshLoader.h"

typedef enum
{
	MESH_TYPE_CLOTH,
	MESH_TYPE_TET,

	MESH_TYPE_TOTAL_NUM
} MeshType;

struct Edge
{
	unsigned int v1, v2; // indices of endpoint vertices
	unsigned int tri1, tri2; // indices of adjacent faces
};

class Mesh
{
	friend class FEMSimObj;

public:
	Mesh() : meshType() {}
	Mesh(MeshType meshType) : meshType(meshType) {}
	virtual ~Mesh() {};
	virtual bool init() { std::cout << "Warning: reach base class virtual init function." << std::endl; return false; }

public:
	MeshType meshType;
	unsigned int verticesNum; // m
	unsigned int systemDimension; // 3m

	// vertices positions/prev positions/ mass
	VectorX restposePos; // 1x3m
	VectorX currPos; // 1x3m
	VectorX currVel; // 1x3m
	VectorX prevPos; // 1x3m
	VectorX prevVel; // 1x3m
	SparseMatrix massMat; // 3mx3m
	SparseMatrix invMassMat; // 3mx3m

	// mass
	ScalarType totalMass;

	// for generating constraints.
	std::vector<Edge> edgesList;

protected:
	// initialize every particle pos / vel / mass / color.
	virtual void generateParticleList() { std::cout << "Warning: reach base class virtual function." << std::endl; }
	// generate edge list from the geometry representation.
	virtual void generateEdgeList() { std::cout << "Warning: reach base class virtual function." << std::endl; }
};


class TetMesh : public Mesh
{
	friend class FEMSimObj;

public:
	TetMesh() : Mesh(MESH_TYPE_TET), loadedMesh(NULL) {}
	virtual ~TetMesh() { if (loadedMesh) { delete loadedMesh; } }

	virtual bool init(std::string fileName);
	bool exportToFile(int frameNum);

protected:
	// tet mesh if loaded from mesh file
	MeshLoader *loadedMesh;

protected:
	// initialize every particle pos / vel / mass / color.
	virtual void generateParticleList();
	// generate edge list from the geometry representation.
	virtual void generateEdgeList();
};
#endif