#ifndef _FEMSIMOBJ_H_
#define _FEMSIMOBJ_H_

#include "Mesh.h"
#include "Constraint.h"

class Mesh;

class FEMSimObj
{
public:
	FEMSimObj();
	virtual ~FEMSimObj();

	void init(TetMesh* mesh, ScalarType h,
		ScalarType gravityConst, ScalarType dampingCoef, 
		ScalarType restitutionCoef, ScalarType frictionCoef,
		ScalarType lsAlpha,	ScalarType lsBeta,
		unsigned int iterationNum, unsigned int maxSubstep);
	void initMatInfo(MaterialType matType, ScalarType matMu, ScalarType matLambda);
	void update();
	bool saveTetAsHDA(int frameNum);
	bool saveGradientAsOBJ();

protected:
	// simulation constants
	ScalarType h;  // time Step
	ScalarType gravityConst;
	ScalarType dampingCoef;  // damping coefficient
	ScalarType restitutionCoef;  // restituation coefficent
	ScalarType frictionCoef;  // friction coefficient
	// material property
	MaterialType matType;
	ScalarType matMu;
	ScalarType matLambda;

	TetMesh *mesh; // main object

	// constant term in optimization
	VectorX y;

	// external force
	VectorX externalForce;

	// iteration for optimization method (gradient descent)
	unsigned int iterationNum;
	unsigned int maxSubstep;

	// line search (backtracking line search)
	ScalarType lsAlpha;
	ScalarType lsBeta;
	ScalarType lsStepSize;

	// constraint 
	std::vector<Constraint*> constraintsList;

	// collision constraints
	std::vector<CollisionConstraint> collisionConstrList;

	// hard coded collision plane for demo
	bool processCollision;

private:
	// update helper functions
	void clearConstraints(); // cleanup all constraints
	void setupConstraints();   // initialize tet constraints
	void dampVelocity();  // damp velocity at the end of each iteration
	void calculateExternalForce();  // gravity force ONLY
	VectorX collisionDetectionPostProcessing(const VectorX& pos); // detect collision and return a vector of penetration
	void collisionDectection(const VectorX& pos);
	void collisionResolution(const VectorX& penetration, VectorX& pos, VectorX& currVel);
	
	void integrateImplicitProjection();  // implicit projection
	
	bool performGradientDescentOneIter(VectorX& pos); // gradient descent

	// key initializations and constants computations
	void computeConstantVectorY();
	void updatePosAndVel(const VectorX& newPos);

	// evaluate energy
	ScalarType evalEnergy(const VectorX& pos); 
	ScalarType evalEnergyPureConstraint(const VectorX& pos);

	// evaluate gradient	
	void evalGradient(const VectorX& pos, VectorX& gradient);	
	void evalGradientPureConstraint(const VectorX& pos, VectorX& gradient);

	// collisions
	ScalarType evalEnergyCollision(const VectorX& pos);
	void evalGradientCollision(const VectorX& pos, VectorX& gradient);

	// energy conservation
	ScalarType evalPotentialEnergy(const VectorX& pos);
	ScalarType evalKineticEnergy(const VectorX& vel);
	ScalarType evalTotalEnergy(const VectorX& pos, const VectorX& vel);

	// line search
	ScalarType lineSearch(const VectorX& pos, const VectorX& gradient, const VectorX& descentDir);

	// check the correctness of the gradient
	void checkGradient(VectorX& pos, VectorX& gradient); 

	// pre allocate
	VectorX predPos;
};

#endif