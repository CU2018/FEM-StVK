#include "Constraint.h"

//----------Constraint Class----------
Constraint::Constraint()
{
	constrType = CONSTRAINT_TYPE_NULL;
}

Constraint::Constraint(ConstraintType type) :
constrType(type)
{

}

Constraint::Constraint(ConstraintType type, ScalarType stiffness) :
constrType(type),
constrStiffness(stiffness)
{

}

Constraint::Constraint(const Constraint& other) :
constrType(other.constrType),
constrStiffness(other.constrStiffness)
{
}

Constraint::~Constraint()
{
}