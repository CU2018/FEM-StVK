#include "Constraint.h"

CollisionSpringConstraint::CollisionSpringConstraint(unsigned int p0, const EigenVector3& fixedPoint, const EigenVector3& normal) :
	Constraint(CONSTRAINT_TYPE_COLLISION),
	p0(p0),
	fixedPoint(fixedPoint),
	normal(normal)
{

}

CollisionSpringConstraint::CollisionSpringConstraint(ScalarType stiffness, unsigned int p0, const EigenVector3& fixedPoint, const EigenVector3& normal) :
	Constraint(CONSTRAINT_TYPE_COLLISION, stiffness),
	p0(p0),
	fixedPoint(fixedPoint),
	normal(normal)
{

}

CollisionSpringConstraint::CollisionSpringConstraint(const CollisionSpringConstraint& other) :
	Constraint(other),
	p0(other.p0),
	fixedPoint(other.fixedPoint),
	normal(other.normal)
{

}

CollisionSpringConstraint::~CollisionSpringConstraint()
{

}

bool CollisionSpringConstraint::isActive(const VectorX& pos)
{
	EigenVector3 pos0 = pos.block_vector(p0);  // position of p0

	if ((pos0 - fixedPoint).dot(normal) > 0)
		return false;
	else
		return true;
}

// 0.5*k*(current_length)^2
ScalarType CollisionSpringConstraint::evalEnergy(const VectorX& pos)
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

ScalarType CollisionSpringConstraint::getEnergy()
{
	return energy;
}

void CollisionSpringConstraint::evalGradient(const VectorX& pos, VectorX& gradient)
{
	if (isActive(pos))
	{
		EigenVector3 currVGradient = (constrStiffness)*(pos.block_vector(p0) - fixedPoint);
		gradient.block_vector(p0) += currVGradient;
	}
}

void CollisionSpringConstraint::evalGradient(const VectorX& pos)
{
	if (isActive(pos))
	{
		constrGradient = (constrStiffness)*(pos.block_vector(p0) - fixedPoint);
	}
	else
	{
		constrGradient.setZero();
	}
}

void CollisionSpringConstraint::getGradient(VectorX& gradient)
{
	gradient.block_vector(p0) += constrGradient;
}

ScalarType CollisionSpringConstraint::evalEnergyAndGradient(const VectorX& pos, VectorX& gradient)
{
	evalEnergyAndGradient(pos);
	gradient.block_vector(p0) += constrGradient;

	return energy;
}

ScalarType CollisionSpringConstraint::evalEnergyAndGradient(const VectorX& pos)
{
	if (isActive(pos))
	{
		// energy
		energy = 0.5*(constrStiffness)*(pos.block_vector(p0) - fixedPoint).squaredNorm();
		// gradient
		constrGradient = (constrStiffness)*(pos.block_vector(p0) - fixedPoint);
	}
	else
	{
		energy = 0;
		constrGradient.setZero();
	}

	return energy;
}
ScalarType CollisionSpringConstraint::getEnergyAndGradient(VectorX& gradient)
{
	gradient.block_vector(p0) += constrGradient;
	return energy;
}