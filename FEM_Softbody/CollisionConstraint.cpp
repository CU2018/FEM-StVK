#include "Constraint.h"

CollisionConstraint::CollisionConstraint(unsigned int p0, const EigenVector3& fixedPoint, const EigenVector3& normal) :
	Constraint(CONSTRAINT_TYPE_COLLISION),
	p0(p0),
	fixedPoint(fixedPoint),
	normal(normal)
{

}

CollisionConstraint::CollisionConstraint(ScalarType stiffness, unsigned int p0, const EigenVector3& fixedPoint, const EigenVector3& normal) :
	Constraint(CONSTRAINT_TYPE_COLLISION, stiffness),
	p0(p0),
	fixedPoint(fixedPoint),
	normal(normal)
{

}

CollisionConstraint::CollisionConstraint(const CollisionConstraint& other) :
	Constraint(other),
	p0(other.p0),
	fixedPoint(other.fixedPoint),
	normal(other.normal)
{

}

CollisionConstraint::~CollisionConstraint()
{

}

bool CollisionConstraint::isActive(const VectorX& pos)
{
	EigenVector3 pos0 = pos.block_vector(p0);  // position of p0

	if ((pos0 - fixedPoint).dot(normal) > 0)
		return false;
	else
		return true;
}

// 0.5*k*(current_length)^2
ScalarType CollisionConstraint::evalEnergy(const VectorX& pos)
{
	if (isActive(pos))
	{
		ScalarType rtnEnergy = 0.5*(constrStiffness)*(pos.block_vector(p0) - fixedPoint).squaredNorm();
		energy = rtnEnergy;
		return rtnEnergy;
	}
	else
	{
		energy = 0;
		return 0;
	}
}

ScalarType CollisionConstraint::getEnergy()
{
	return energy;
}

void CollisionConstraint::evalGradient(const VectorX& pos, VectorX& gradient)
{
	if (isActive(pos))
	{
		EigenVector3 currVGradient = (constrStiffness)*(pos.block_vector(p0) - fixedPoint);
		gradient.block_vector(p0) += currVGradient;
	}
}

void CollisionConstraint::getGradient(VectorX& gradient)
{
	gradient.block_vector(p0) += constrGradient;
}
