#ifndef _CONSTRAINT_H_
#define _CONSTRAINT_H_

#include <vector>
#include <iostream>
#include "MathHeaders.h"

typedef enum
{
	MATERIAL_TYPE_StVK,
	MATERIAL_TYPE_TOTAL_NUM
} MaterialType;

typedef enum
{
	CONSTRAINT_TYPE_COLLISION,
	CONSTRAINT_TYPE_TET,
	CONSTRAINT_TYPE_NULL,
	CONSTRAINT_TYPE_TOTAL_NUM
} ConstraintType;

class Constraint
{
public:
	Constraint();
	Constraint(ConstraintType type);
	Constraint(ConstraintType type, ScalarType stiffness);
	Constraint(const Constraint& other);
	virtual ~Constraint();

	virtual void getMaterialProperty(ScalarType& stiffness) { stiffness = constrStiffness; }
	virtual void getMaterialProperty(MaterialType& type, ScalarType& mu, ScalarType& lambda) { std::cout << "Warning: reach <Constraint::getMaterialProperty> base class virtual function." << std::endl; }
	virtual void setMaterialProperty(ScalarType stiffness) { constrStiffness = stiffness; }

	virtual bool vertexIncluded(unsigned int vi) { return false; }

	virtual ScalarType  evalEnergy(const VectorX& pos) { std::cout << "Warning: reach <Constraint::evalPotentialEnergy> base class virtual function." << std::endl; return 0; }
	virtual ScalarType  getEnergy() { std::cout << "Warning: reach <Constraint::GetPotentialEnergy> base class virtual function." << std::endl; return 0; }
	
	virtual void  evalGradient(const VectorX& pos, VectorX& gradient) { std::cout << "Warning: reach <Constraint::EvalGradient> base class virtual function." << std::endl; }
	virtual void  getGradient(VectorX& gradient) { std::cout << "Warning: reach <Constraint::getGradient> base class virtual function." << std::endl; }

	// inline
	const ConstraintType& Type() { return constrType; }

protected:
	ConstraintType constrType;
	ScalarType constrStiffness;  // constraint stiffness

	// saved energy
	ScalarType energy;

	// for visualization and selection
public:
	virtual void writeToFileOBJ(std::ofstream& outfile, int& existing_vertices) { /*do nothing*/ }
	virtual void writeToFileOBJHead(std::ofstream& outfile) { /*do nothing*/ }
	virtual void writeToFileOBJTet(std::ofstream& outfile) { /*do nothing*/ }
};

class CollisionConstraint : public Constraint
{
public:
	CollisionConstraint(unsigned int p0, const EigenVector3& fixedPoint, const EigenVector3& normal);
	CollisionConstraint(ScalarType stiffness, unsigned int p0, const EigenVector3& fixedPoint, const EigenVector3& normal);
	CollisionConstraint(const CollisionConstraint& other);
	virtual ~CollisionConstraint();

	bool isActive(const VectorX& pos);
	virtual ScalarType  evalEnergy(const VectorX& pos);
	virtual ScalarType  getEnergy();
	virtual void  evalGradient(const VectorX& pos, VectorX& gradient);
	virtual void  getGradient(VectorX& gradient);

protected:
	unsigned int p0;  // collided point
	EigenVector3 constrGradient;
	EigenVector3 fixedPoint;  // the position that we want the point p0 to be set
	EigenVector3 normal;
};

class TetConstraint : public Constraint
{
public:
	TetConstraint(unsigned int p1, unsigned int p2, unsigned int p3, unsigned int p4, VectorX& pos);
	TetConstraint(const TetConstraint& other);
	virtual ~TetConstraint();

	virtual bool vertexIncluded(unsigned int vi) { for (unsigned int i = 0; i != 4; i++) { return vi == indicesList[i]; } return false; }

	virtual void getMaterialProperty(MaterialType& type, ScalarType& mu, ScalarType& lambda);
	virtual void setMaterialProperty(MaterialType type, ScalarType mu, ScalarType lambda);

	virtual ScalarType  evalEnergy(const VectorX& pos);
	virtual ScalarType  getEnergy();
	virtual void  evalGradient(const VectorX& pos, VectorX& gradient);
	virtual void  getGradient(VectorX& gradient);

	// set mass matrix
	ScalarType setMassMatrix(std::vector<SparseMatrixTriplet>& massMat);

//public:
//	virtual void WriteToFileOBJTet(std::ofstream& outfile, const VectorX& pos);

private:
	void getMatrixDs(EigenMatrix3& Ds, const VectorX& pos);  // Ds
	void getDeformationGradient(EigenMatrix3& F, const VectorX& pos);  // F
	void getStressTensor(EigenMatrix3& P, const EigenMatrix3& F);  // first Piola-Kirchhoff stress tensor 

protected:
	// material properties
	MaterialType matType;
	ScalarType mu;  // Lame's second parameter - resistance for shearing change
	ScalarType lambda; // Lame's first parameter - resistance for volume change

	unsigned int indicesList[4]; // indices of four vertices
	EigenVector3 gradientsList[4];  // gradient
	EigenMatrix3 Dm; // [x1-x4|x2-x4|x3-x4] reference shape matrix
	EigenMatrix3 DmInv; // inverse of Dm
	ScalarType w; // 1/6 det(Dr); rest-pose volume of the element
};

#endif