#include "FEMSimObj.h"

#define COLLISION_EPSILON 1e-3

// preset plane at the center of the plane
bool planeIntersectionTest(const EigenVector3& p, EigenVector3& normal, ScalarType& dist)
{
	ScalarType height = 0.0f;  // the plane is at the center of the scene with 0.0f height
	dist = p(1) - height - COLLISION_EPSILON;
	normal = EigenVector3(0.0f, 1.0f, 0.0f);  // facing up

	if (dist < 0)
		return true;
	else
		return false;
}

FEMSimObj::FEMSimObj()
{
	processCollision = true;
}

FEMSimObj::~FEMSimObj()
{
	clearConstraints();
}

// TODO: pass in necessary variables and initialize them
void FEMSimObj::init(TetMesh* mesh, ScalarType h,
	ScalarType gravityConst, ScalarType dampingCoef,
	ScalarType restitutionCoef, ScalarType frictionCoef,
	ScalarType lsAlpha, ScalarType lsBeta,
	unsigned int iterationNum, unsigned int maxSubstep)
{
	(this->mesh) = mesh;
	this->h = h;
	this->gravityConst = gravityConst;
	this->dampingCoef = dampingCoef;
	this->restitutionCoef = restitutionCoef;
	this->frictionCoef = frictionCoef;
	this->lsAlpha = lsAlpha;
	this->lsBeta = lsBeta;
	this->iterationNum = iterationNum;
	this->maxSubstep = maxSubstep;
	this->processCollision = true;

	setupConstraints();
}

void FEMSimObj::update()
{
	// compute external force
	calculateExternalForce();

	for (unsigned int currSubstep = 0; currSubstep < maxSubstep; ++currSubstep)
	{
		// compute inertia term
		computeConstantVectorY();
		// perform gradient descent
		integrateImplicitProjection();
		// damping
		dampVelocity();
	}
}


bool FEMSimObj::saveTetAsHDA(int frameNum)
{
	return mesh->exportToFile(frameNum);
}

bool FEMSimObj::saveGradientAsOBJ()
{
	return false;
}

void FEMSimObj::clearConstraints()
{
	for (unsigned int i = 0; i < constraintsList.size(); ++i)
	{
		delete constraintsList[i];
	}
	constraintsList.clear();
}

void FEMSimObj::setupConstraints()
{
	clearConstraints();
	
	stiffnessHigh = 1e5;

	// reset mass matrix for tet simulation:
	ScalarType totalVolume = 0;
	std::vector<SparseMatrixTriplet> massTriplets;
	// std::vector<SparseMatrixTriplet> mass1dTriplets;
	massTriplets.clear();
	//mass1dTriplets.clear();

	VectorX& pos = mesh->currPos;
	TetMesh* tetMesh = dynamic_cast<TetMesh*>(mesh);

	for (unsigned int i = 0; i < tetMesh->loadedMesh->tetsList.size(); ++i)
	{
		MeshLoader::Tet& tet = tetMesh->loadedMesh->tetsList[i];
		TetConstraint *tc;
		tc = new TetConstraint(tet.id1, tet.id2, tet.id3, tet.id4, pos);
		constraintsList.push_back(tc);

		// totalVolume += tc->setMassMatrix(massTriplets, mass1dTriplets);
		totalVolume += tc->setMassMatrix(massTriplets);
		/*mesh->expandedSysDim += 9; 
		mesh->expandedSysDim1d += 3;*/
	}

	mesh->massMat.setFromTriplets(massTriplets.begin(), massTriplets.end());
	//mesh->massMat1d.setFromTriplets(mass1dTriplets.begin(), mass1dTriplets.end());

	mesh->massMat = mesh->massMat * (mesh->totalMass / totalVolume);
	//mesh->massMat1d = mesh->massMat1d * (mesh->totalMass / totalVolume);

	std::vector<SparseMatrixTriplet> massInvTriplets;
	massInvTriplets.clear();
	/*std::vector<SparseMatrixTriplet> massInv1dTriplets;
	massInv1dTriplets.clear();*/
	for (unsigned int i = 0; i != mesh->massMat.rows(); i++)
	{
		ScalarType mi = mesh->massMat.coeff(i, i);
		ScalarType miInv;
		if (std::abs(mi) > 1e-12)
		{
			miInv = 1.0 / mi;
		}
		else
		{
			mesh->massMat.coeffRef(i, i) = 1e-12;
			miInv = 1e12;
		}
		massInvTriplets.push_back(SparseMatrixTriplet(i, i, miInv));
	}
	//for (unsigned int i = 0; i != mesh->massmat1d.rows(); i++)
	//{
	//	scalartype mi = mesh->massmat1d.coeff(i, i);
	//	scalartype miinv;
	//	if (std::abs(mi) > 1e-12)
	//	{
	//		miinv = 1.0 / mi;
	//	}
	//	else
	//	{
	//		mesh->massmat1d.coeffref(i, i) = 1e-12;
	//		miinv = 1e12;
	//	}
	//	massinv1dtriplets.push_back(sparsematrixtriplet(i, i, miinv));
	//}

	mesh->invMassMat.setFromTriplets(massInvTriplets.begin(), massInvTriplets.end());
	// mesh->invMassMat1d.setFromTriplets(massInv1dTriplets.begin(), massInv1dTriplets.end());
}

void FEMSimObj::dampVelocity()
{
	if (std::abs(dampingCoef) < EPSILON)
		return;
	mesh->currVel *= 1 - dampingCoef;
}

void FEMSimObj::calculateExternalForce()
{
	externalForce.resize(mesh->systemDimension);
	externalForce.setZero();
	// gravity
	for (unsigned int i = 0; i < mesh->verticesNum; ++i)
	{
		externalForce[3 * i + 1] += -gravityConst; // y value of a vec3f
	}
	externalForce = mesh->massMat * externalForce;
}

VectorX FEMSimObj::collisionDetectionPostProcessing(const VectorX& pos)
{
	VectorX penetration(mesh->systemDimension);
	penetration.setZero();
	EigenVector3 normal;
	ScalarType dist;

	for (unsigned int i = 0; i != mesh->verticesNum; ++i)
	{
		EigenVector3 onePos = pos.block_vector(i);

		if (planeIntersectionTest(onePos, normal, dist))
		{
			penetration.block_vector(i) += (dist)* normal;
		}
	}

	return penetration;
}

void FEMSimObj::collisionDectection(const VectorX& pos)
{
	collisionConstrList.clear();
	EigenVector3 surfacePoint;
	EigenVector3 normal;
	ScalarType dist;

	for (unsigned int i = 0; i != mesh->verticesNum; ++i)
	{
		EigenVector3 onePos = pos.block_vector(i);

		if (planeIntersectionTest(onePos, normal, dist))
		{
			surfacePoint = onePos - normal * dist; // dist is negative...
			printf("collided with the floor! surface point: %f\n", surfacePoint[1]);
			collisionConstrList.push_back(CollisionSpringConstraint(1e3, i, surfacePoint, normal));
		}
	}

	VectorX penetration;
	if (collisionConstrList.size() != 0)
	{
		penetration = collisionDetectionPostProcessing(pos);
		collisionResolution(penetration, mesh->currPos, mesh->currVel);
	}
}

void FEMSimObj::collisionResolution(const VectorX& penetration, VectorX& pos, VectorX& vel)
{
	EigenVector3 currVPos, currVVel, currVPene, currVNorm;
	EigenVector3 velin, velit;
	for (unsigned int i = 0; i != mesh->verticesNum; ++i)
	{
		currVPos = pos.block_vector(i);
		currVVel = vel.block_vector(i);
		currVPene = penetration.block_vector(i);

		ScalarType dist = currVPene.norm();
		if (dist > EPSILON) // there is collision
		{
			currVNorm = -currVPene / dist; // normalize
			currVPos -= currVPene;
			velin = currVVel.dot(currVNorm)*currVNorm;
			velit = currVVel - velin;
			currVVel = -(restitutionCoef) * velin + (1 - frictionCoef) * velit;
			pos.block_vector(i) = currVPos;
			vel.block_vector(i) = currVVel;
		}
	}
	//for (int i = 0; i < 12; ++i)
	//{
	//	printf(" %f", mesh->currPos[i]);
	//}
}

void FEMSimObj::integrateImplicitProjection()
{
	// take a initial guess: constant y
	VectorX predPos = y; 

	// while loop until converge or exceeds maximum iterations
	bool converge = false;
	for (unsigned int currIter = 0; !converge && currIter < iterationNum; ++currIter)
	{
		if (processCollision)
		{
			// Collision Detection every iteration
			collisionDectection(predPos);
		}

		converge = performGradientDescentOneIter(predPos);
		// if (converge) printf("converge\n");
	}

	// update constants
	updatePosAndVel(predPos);
}

bool FEMSimObj::performGradientDescentOneIter(VectorX& pos)
{
	// calcuate gradient direction
	VectorX gradient;
	evalGradient(pos, gradient);

	if (gradient.norm() < EPSILON)
		return true;

	VectorX descentDir = -gradient;  // descent direction is the opposite of the gradient direction
	
	// line search
	ScalarType stepSize = lineSearch(pos, gradient, descentDir);

	// update pos
	pos = pos + descentDir * stepSize;

	// check convergence
	if (stepSize < EPSILON)
		return true;
	else
		return false;
}

void FEMSimObj::computeConstantVectorY()
{
	// implicit euler
	y = mesh->currPos + mesh->currVel * h + h * h * mesh->invMassMat*externalForce;
}

void FEMSimObj::updatePosAndVel(const VectorX& newPos)
{
	mesh->prevVel = mesh->currVel;
	mesh->prevPos = mesh->currPos;
	mesh->currVel = (newPos - mesh->currPos) / h;
	mesh->currPos = newPos;
	for (int i = 0; i < 12; ++i)
	{
		printf(" %f", mesh->currPos[i]);
	}
}

ScalarType FEMSimObj::evalEnergy(const VectorX& pos)
{
	ScalarType energyPureConstraints, rtnEnergy;
	ScalarType inertia = 0.5 * (pos - y).transpose() * mesh->massMat * (pos - y);
	ScalarType hSquare = h * h;

	energyPureConstraints = evalEnergyPureConstraint(pos);
	rtnEnergy = inertia + hSquare * energyPureConstraints;

	return rtnEnergy;
}

ScalarType FEMSimObj::evalEnergyPureConstraint(const VectorX& pos)
{
	ScalarType rtnEnergy = 0.0;

	for (std::vector<Constraint*>::iterator it = constraintsList.begin(); it != constraintsList.end(); ++it)
	{
		rtnEnergy += (*it)->evalEnergy(pos);
	}

	// hardcoded collision plane
	if (processCollision)
	{
		rtnEnergy += evalEnergyCollision(pos);
	}

	return rtnEnergy;
}

void FEMSimObj::evalGradient(const VectorX& pos, VectorX& gradient)
{
	ScalarType hSquare = h * h;

	evalGradientPureConstraint(pos, gradient);
	gradient = mesh->massMat * (pos - y) + hSquare * gradient;
}

void FEMSimObj::evalGradientPureConstraint(const VectorX& pos, VectorX& gradient)
{
	gradient.resize(mesh->systemDimension);
	gradient.setZero();

	// constraints single thread
	for (std::vector<Constraint*>::iterator iter = constraintsList.begin(); iter != constraintsList.end(); ++iter)
	{
		(*iter)->evalGradient(pos, gradient);
	}

	// hardcoded collision plane
	if (processCollision)
	{
		VectorX gc;

		evalGradientCollision(pos, gc);

		gradient += gc;
	}
}

ScalarType FEMSimObj::evalEnergyCollision(const VectorX& pos)
{
	ScalarType rtnEnergy = 0.0;

	for (std::vector<CollisionSpringConstraint>::iterator iter = collisionConstrList.begin(); iter != collisionConstrList.end(); ++iter)
	{
		rtnEnergy += iter->evalEnergy(pos);
	}

	return rtnEnergy;
}

void FEMSimObj::evalGradientCollision(const VectorX& pos, VectorX& gradient)
{
	gradient.resize(mesh->systemDimension);
	gradient.setZero();

	for (std::vector<CollisionSpringConstraint>::iterator iter = collisionConstrList.begin(); iter != collisionConstrList.end(); ++iter)
	{
		iter->evalGradient(pos, gradient);
	}
}

ScalarType FEMSimObj::lineSearch(const VectorX& pos, const VectorX& gradient, const VectorX& descentDir)
{
	VectorX newPos(mesh->systemDimension);
	ScalarType attemptStepSize = 1.0 / lsBeta;
	ScalarType lhs, rhs;

	ScalarType currObjectiveVal;
	try
	{
		currObjectiveVal = evalEnergy(pos);
	}
	catch (const std::exception& e)
	{
		std::cout << e.what() << std::endl;
	}
	do
	{
		attemptStepSize *= lsBeta;
		newPos = pos + attemptStepSize * descentDir;

		lhs = 1e15;
		rhs = 0;
		try
		{
			lhs = evalEnergy(newPos);
		}
		catch (const std::exception&)
		{
			continue;  // not a good candidate for step size
		}
		rhs = currObjectiveVal + lsAlpha * attemptStepSize * (gradient.transpose() * descentDir)(0);
		if (lhs >= rhs)
			continue;  // not a good candidate for step size
		break; // exit line search looping

	} while (attemptStepSize > 1e-5);

	if (attemptStepSize < 1e-5)
	{
		attemptStepSize = 0.0f;
	}
	lsStepSize = attemptStepSize;
	return lsStepSize;
}

ScalarType FEMSimObj::evalPotentialEnergy(const VectorX& pos)
{
	ScalarType rtnEnergy = evalEnergyPureConstraint(pos);
	rtnEnergy -= externalForce.dot(pos);

	return rtnEnergy;
}

ScalarType FEMSimObj::evalKineticEnergy(const VectorX& vel) 
{
	return (0.5*vel.transpose()*mesh->massMat*vel);
}

ScalarType FEMSimObj::evalTotalEnergy(const VectorX& pos, const VectorX& vel)
{
	return (evalPotentialEnergy(pos) + evalKineticEnergy(vel));
}