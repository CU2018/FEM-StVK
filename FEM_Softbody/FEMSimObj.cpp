#include "FEMSimObj.h"

#define COLLISION_EPSILON 1e-3

// preset plane at the center of the plane
bool planeIntersectionTest(const EigenVector3& p, EigenVector3& normal, ScalarType& dist)
{
	ScalarType height = 0.0f;  // the plane is at the center of the scene with 0.0f height
	dist = p(1) - height - COLLISION_EPSILON;
	normal = EigenVector3(0.0f, 1.0f, 0.0f);  // facing up

	return dist < 0 ? true : false;
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
	y.resize(mesh->systemDimension);
	externalForce.resize(mesh->systemDimension);
	collisionConstrList.clear();
}

void  FEMSimObj::initMatInfo(MaterialType matType, ScalarType matMu, ScalarType matKappa)
{
	this->matType = matType;
	this->matMu = matMu;
	this->matLambda = matLambda;
	setupConstraints();
}

void FEMSimObj::update()
{
	// compute external force
	calculateExternalForce();

	printf("---------------------------------------------\n");
	for (unsigned int currSubstep = 0; currSubstep < maxSubstep; ++currSubstep)
	{
		// compute inertia term
		computeConstantVectorY();
		// perform gradient descent
		integrateImplicitProjection();
		
		// damping
		dampVelocity();
	}
	printf("Frame CurrPos:\n ");
	std::cout << mesh->currPos << std::endl;
	printf("\n");
	// output total energy;
	ScalarType K = evalKineticEnergy(mesh->currVel);
	ScalarType W = evalPotentialEnergy(mesh->currPos);
	ScalarType KPlusW = K + W;
	std::cout << "Kinetic Energy = " << K << std::endl;
	std::cout << "Potential Energy = " << W << std::endl;
	std::cout << "Total Energy = " << KPlusW << std::endl;
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
	if (constraintsList.size() != 0) clearConstraints();

	// reset mass matrix for tet simulation:
	ScalarType totalVolume = 0;
	std::vector<SparseMatrixTriplet> massTriplets;
	massTriplets.clear();

	VectorX& pos = mesh->currPos;
	TetMesh* tetMesh = dynamic_cast<TetMesh*>(mesh);

	for (unsigned int i = 0; i < tetMesh->loadedMesh->tetsList.size(); ++i)
	{
		MeshLoader::Tet& tet = tetMesh->loadedMesh->tetsList[i];
		TetConstraint *tc = new TetConstraint(tet.id1, tet.id2, tet.id3, tet.id4, pos);
		tc->setMaterialProperty(MATERIAL_TYPE_StVK, matMu, matLambda);
		constraintsList.push_back(tc);
		totalVolume += tc->setMassMatrix(massTriplets);
	}

	mesh->massMat.setFromTriplets(massTriplets.begin(), massTriplets.end());
	// mesh->massMat = mesh->massMat * (mesh->totalMass / totalVolume);

	std::vector<SparseMatrixTriplet> massInvTriplets;
	massInvTriplets.clear();
	for (unsigned int i = 0; i != mesh->massMat.rows(); ++i)
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
	mesh->invMassMat.setFromTriplets(massInvTriplets.begin(), massInvTriplets.end());
}

void FEMSimObj::dampVelocity()
{
	if (std::abs(dampingCoef) < EPSILON)
		return;
	// mesh->currVel *= (1 - dampingCoef);

	mesh->currVel *= pow((1 - dampingCoef), h);
}

void FEMSimObj::calculateExternalForce()
{
	externalForce.setZero();
	// gravity
	for (unsigned int i = 0; i < mesh->verticesNum; ++i)
	{
		externalForce[3 * i + 1] += -gravityConst; // y value of a vec3f
	}
	// std::cout << mesh->massMat << std::endl;
	externalForce = mesh->massMat * externalForce;
	// std::cout << externalForce << std::endl;
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
	// std::cout << penetration << std::endl;
	return penetration;
}

void FEMSimObj::collisionDectection(const VectorX& pos)
{
	if(collisionConstrList.size() != 0 ) collisionConstrList.clear();
	EigenVector3 surfacePoint;
	EigenVector3 normal;
	ScalarType dist;

	for (unsigned int i = 0; i != mesh->verticesNum; ++i)
	{
		EigenVector3 onePos = pos.block_vector(i);
		if (planeIntersectionTest(onePos, normal, dist))
		{
			surfacePoint = EigenVector3(onePos(0), 0.0f, onePos(2));
			collisionConstrList.push_back(CollisionConstraint(1e3, i, surfacePoint, normal));
		}
		//if (planeIntersectionTest(onePos, normal, dist))
		//{
		//	surfacePoint = onePos - normal * dist; // dist is negative...
		//	// printf("collided with the floor! surface point: %f\n", surfacePoint[1]);
		//	collisionConstrList.push_back(CollisionConstraint(1e3, i, surfacePoint, normal));
		//}
	}

	VectorX penetration;
	if (collisionConstrList.size() != 0)
	{
		penetration = collisionDetectionPostProcessing(pos);

		// std::cout << penetration << std::endl;
		collisionResolution(penetration, mesh->currPos, mesh->currVel);
	}
}

void FEMSimObj::collisionResolution(const VectorX& penetration, VectorX& pos, VectorX& vel)
{
	//EigenVector3 currVPos, currVVel, currVPene, currVNorm;
	//EigenVector3 velin, velit;
	//for (unsigned int i = 0; i != mesh->verticesNum; ++i)
	//{
	//	currVPos = pos.block_vector(i);
	//	currVVel = vel.block_vector(i);
	//	currVPene = penetration.block_vector(i);

	//	ScalarType dist = currVPene.norm();
	//	if (dist > EPSILON) // there is collision
	//	{
	//		currVNorm = -currVPene / dist; // normalize
	//		currVPos -= currVPene;
	//		velin = currVVel.dot(currVNorm)*currVNorm;
	//		velit = currVVel - velin;
	//		currVVel = -(restitutionCoef) * velin + (1 - frictionCoef) * velit;
	//		pos.block_vector(i) = currVPos;
	//		vel.block_vector(i) = currVVel;
	//	}
	//}
	EigenVector3 xi, vi, pi, ni;
	EigenVector3 vin, vit;
	for (unsigned int i = 0; i != mesh->verticesNum; ++i)
	{
		xi = pos.block_vector(i);
		vi = vel.block_vector(i);
		pi = penetration.block_vector(i);

		ScalarType dist = pi.norm();
		if (dist > EPSILON) // there is collision
		{
			ni = -pi / dist; // normalize
			xi -= pi;
			vin = vi.dot(ni)*ni;
			vit = vi - vin;
			vi = -(restitutionCoef)*vin + (1 - frictionCoef) * vit;
			pos.block_vector(i) = xi;
			vel.block_vector(i) = vi;
		}
	}
	/*printf("Collision resolution:\n");
	std::cout << pos << std::endl;*/
}

void FEMSimObj::integrateImplicitProjection()
{
	// take a initial guess: constant y
	predPos = y; 

	// while loop until converge or exceeds maximum iterations
	bool converge = false;
	unsigned int currIter = 0;
	for (currIter = 0; !converge && currIter < iterationNum; ++currIter)
	{
		if (processCollision)
		{
			// Collision Detection every iteration
			collisionDectection(predPos);
		}

		converge = performGradientDescentOneIter(predPos);
		// if (converge) printf("converge\n");
	}
	printf("converge iter: %d\n", currIter);
	// update constants
	updatePosAndVel(predPos);
}

void FEMSimObj::checkGradient(VectorX& pos, VectorX& gradient)
{
	VectorX smallShiftX(pos.size());
	for (unsigned int i = 0; i < 12; ++i)
	{
		VectorX posCopy = pos;
		ScalarType xAdd	= pos(i) + COLLISION_EPSILON;
		posCopy(i) = xAdd;
		ScalarType energyAdd = evalEnergy(posCopy);

		posCopy = pos;
		ScalarType xSub = pos(i) - COLLISION_EPSILON;
		posCopy(i) = xSub;
		ScalarType energySub = evalEnergy(posCopy);

		smallShiftX(i) = (energyAdd - energySub) / (2 * COLLISION_EPSILON);
	}
	/*printf("Check Gradient:\n");
	std::cout << smallShiftX << std::endl;*/
}

bool FEMSimObj::performGradientDescentOneIter(VectorX& pos)
{
	// calcuate gradient direction
	VectorX gradient;
	// printf("Eval Gradient:\n");
	evalGradient(pos, gradient);
	// std::cout << gradient << std::endl;

	checkGradient(pos ,gradient);
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
	/*for (int i = 0; i < 12; ++i)
	{
		printf(" %f", mesh->currPos[i]);
		printf("\n");
	}*/
}

ScalarType FEMSimObj::evalEnergy(const VectorX& pos)
{
	ScalarType energyPureConstraints, rtnEnergy;
	ScalarType inertia = 0.5 * (pos - y).transpose() * mesh->massMat * (pos - y);
	ScalarType hSquare = h * h;

	energyPureConstraints = evalEnergyPureConstraint(pos);
	rtnEnergy = inertia + hSquare * energyPureConstraints;  // g(x)

	return rtnEnergy;
}

ScalarType FEMSimObj::evalEnergyPureConstraint(const VectorX& pos)
{
	ScalarType rtnEnergy = 0.0;
	// tet constraints
	for (std::vector<Constraint*>::iterator it = constraintsList.begin(); it != constraintsList.end(); ++it)
	{
		rtnEnergy += (*it)->evalEnergy(pos);
	}

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
	gradient = mesh->massMat * (pos - y) + hSquare * gradient; // M*(x-y) + h^2*grad
}

void FEMSimObj::evalGradientPureConstraint(const VectorX& pos, VectorX& gradient)
{
	gradient.resize(mesh->systemDimension);
	gradient.setZero();

	// tet contraints
	for (std::vector<Constraint*>::iterator iter = constraintsList.begin(); iter != constraintsList.end(); ++iter)
	{
		(*iter)->evalGradient(pos, gradient);
	}

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
	// Collision Constraint
	for (std::vector<CollisionConstraint>::iterator iter = collisionConstrList.begin(); iter != collisionConstrList.end(); ++iter)
	{
		rtnEnergy += iter->evalEnergy(pos);
	}

	return rtnEnergy;
}

void FEMSimObj::evalGradientCollision(const VectorX& pos, VectorX& gradient)
{
	gradient.resize(mesh->systemDimension);
	gradient.setZero();

	for (std::vector<CollisionConstraint>::iterator iter = collisionConstrList.begin(); iter != collisionConstrList.end(); ++iter)
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