#include "Constraint.h"

// combine three column vectors to a matrix
void ThreeVector3ToMatrix3(EigenMatrix3& m, EigenVector3& v1, EigenVector3& v2, EigenVector3& v3)
{
	m.block<3, 1>(0, 0) = v1;
	m.block<3, 1>(0, 1) = v2;
	m.block<3, 1>(0, 2) = v3;
}

//----------TetConstraint Class----------
TetConstraint::TetConstraint(unsigned int p1, unsigned int p2, unsigned int p3, unsigned int p4, VectorX& verticesList) :
	Constraint(CONSTRAINT_TYPE_TET)
{
	indicesList[0] = p1;
	indicesList[1] = p2;
	indicesList[2] = p3;
	indicesList[3] = p4;

	EigenVector3 v1 = verticesList.block_vector(p1) - verticesList.block_vector(p4); // p1 - p4
	EigenVector3 v2 = verticesList.block_vector(p2) - verticesList.block_vector(p4); // p2 - p4
	EigenVector3 v3 = verticesList.block_vector(p3) - verticesList.block_vector(p4); // p3 - p4
	
	ThreeVector3ToMatrix3(Dm, v1, v2, v3);

	w = Dm.determinant();

	w = 1.0 / 6.0 * std::abs(w);  // ???????????

	DmInv = Dm.inverse();

	/*
	      | 1  0  0 -1 |  
		  | 0  1  0 -1 |
		  | 0  0  1 -1 |
 	*/
	/*Eigen::Matrix<ScalarType, 3, 4> IND;
	IND.block<3, 3>(0, 0) = EigenMatrix3::Identity();
	IND.block<3, 1>(0, 3) = EigenVector3(-1, -1, -1);

	G = DmInv.transpose() * IND;*/
}

TetConstraint::TetConstraint(const TetConstraint& other) :
	Constraint(other)
{
	indicesList[0] = other.indicesList[0];
	indicesList[1] = other.indicesList[1];
	indicesList[2] = other.indicesList[2];
	indicesList[3] = other.indicesList[3];

	Dm = other.Dm;
	DmInv = other.DmInv;

	w = other.w;
	// G = other.G;

	matType = other.matType;
	mu = other.mu;
	lambda = other.lambda;
	kappa = other.kappa;
}

TetConstraint::~TetConstraint()
{
}

// get current shape matrix
void TetConstraint::getMatrixDs(EigenMatrix3& Ds, const VectorX& pos)
{
	EigenVector3 v1 = pos.block_vector(indicesList[0]) - pos.block_vector(indicesList[3]);
	EigenVector3 v2 = pos.block_vector(indicesList[1]) - pos.block_vector(indicesList[3]);
	EigenVector3 v3 = pos.block_vector(indicesList[2]) - pos.block_vector(indicesList[3]);

	ThreeVector3ToMatrix3(Ds, v1, v2, v3);
}

void TetConstraint::getMaterialProperty(MaterialType& type, ScalarType& mu, ScalarType& lambda, ScalarType& kappa)
{
	type = matType;
	mu = this->mu;
	lambda = this->lambda;
	kappa = this->kappa;
}

void TetConstraint::setMaterialProperty(MaterialType type, ScalarType mu, ScalarType lambda, ScalarType kappa)
{
	matType = type;
	this->mu = mu;
	this->lambda = lambda;
	this->kappa = kappa;
}

// energy = 0.5*mu*||strainTensor||^2 + 0.5*lambda*(tr(strainTensor))^2
ScalarType TetConstraint::evalEnergy(const VectorX& pos)
{
	EigenMatrix3 F;  // deformation gradient
	getDeformationGradient(F, pos);

	ScalarType rtnEnergy = 0;
	switch (matType)
	{
		case MATERIAL_TYPE_StVK:
		{
			EigenMatrix3 I = EigenMatrix3::Identity();  // identity matrix
			EigenMatrix3 E = 0.5*(F.transpose()*F - I);  // Green strain tensor
			rtnEnergy = mu * E.squaredNorm() + 0.5*lambda*std::pow(E.trace(), 2); // relationship between energy and strain
			ScalarType J = F.determinant();
			if (J < 1)
			{
				rtnEnergy += kappa / 12 * std::pow((1 - J) / 6, 3);
			}
		}
		break;
	}

	rtnEnergy *= w;

	energy = rtnEnergy;

	return rtnEnergy;
}

ScalarType TetConstraint::getEnergy()
{
	return energy;
}

void TetConstraint::evalGradient(const VectorX& pos, VectorX& gradient)
{
	EigenMatrix3 F;
	getDeformationGradient(F, pos);

	EigenMatrix3 P;
	getStressTensor(P, F);

	EigenMatrix3 H = w * P * DmInv.transpose();  // the gradient of elastic potential

	EigenVector3 gradientsList[4];
	gradientsList[0] = H.block<3, 1>(0, 0);
	gradientsList[1] = H.block<3, 1>(0, 1);
	gradientsList[2] = H.block<3, 1>(0, 2);
	gradientsList[3] = -gradientsList[0] - gradientsList[1] - gradientsList[2];

	for (unsigned int i = 0; i < 4; i++)
	{
		gradient.block_vector(indicesList[i]) += gradientsList[i];
	}
}

void TetConstraint::evalGradient(const VectorX& pos)
{
	EigenMatrix3 F;  // deformation gradient
	getDeformationGradient(F, pos);

	EigenMatrix3 P;
	getStressTensor(P, F);

	EigenMatrix3 H = w * P * DmInv.transpose(); // the gradient of elastic potential

	gradientsList[0] = H.block<3, 1>(0, 0);
	gradientsList[1] = H.block<3, 1>(0, 1);
	gradientsList[2] = H.block<3, 1>(0, 2);
	gradientsList[3] = -gradientsList[0] - gradientsList[1] - gradientsList[2];
}

void TetConstraint::getGradient(VectorX& gradient)
{
	for (unsigned int i = 0; i < 4; i++)
	{
		gradient.block_vector(indicesList[i]) += gradientsList[i];
	}
}

ScalarType TetConstraint::evalEnergyAndGradient(const VectorX& pos, VectorX& gradient)
{
	evalEnergyAndGradient(pos);

	for (unsigned int i = 0; i < 4; i++)
	{
		gradient.block_vector(indicesList[i]) += gradientsList[i];
	}
	return energy;
}

ScalarType TetConstraint::evalEnergyAndGradient(const VectorX& pos)
{
	EigenMatrix3 F;
	getDeformationGradient(F, pos);

	EigenMatrix3 P; // first Piola-Kirchhoff stress tensor
	ScalarType rtnEnergy = getStressTensorAndEnergyDensity(P, F);
	energy = rtnEnergy * w;

	EigenMatrix3 H = w * P * DmInv.transpose();  // the gradient of elastic potential

	gradientsList[0] = H.block<3, 1>(0, 0);
	gradientsList[1] = H.block<3, 1>(0, 1);
	gradientsList[2] = H.block<3, 1>(0, 2);
	gradientsList[3] = -gradientsList[0] - gradientsList[1] - gradientsList[2];

	return energy;
}

ScalarType TetConstraint::getEnergyAndGradient(VectorX& gradient)
{
	for (unsigned int i = 0; i < 4; i++)
	{
		gradient.block_vector(indicesList[i]) += gradientsList[i];
	}
	return energy;
}

void TetConstraint::getDeformationGradient(EigenMatrix3& F, const VectorX& pos)
{
	EigenMatrix3 Ds;
	getMatrixDs(Ds, pos);
	F = Ds * DmInv;  // F = Ds * Dm^(-1)
}

void TetConstraint::getStressTensor(EigenMatrix3& P, const EigenMatrix3& F)
{
	switch (matType)
	{
		case MATERIAL_TYPE_StVK:
		{
			EigenMatrix3 I = EigenMatrix3::Identity();
			EigenMatrix3 E = 0.5*(F.transpose()*F - I); // Green strain tensor
			P = F * (2 * mu*E + lambda * E.trace() * I);  // first Piola-Kirchhoff stress tensor
			ScalarType J = F.determinant();
			if (J < 1)
			{
				P += -kappa / 24 * std::pow((1 - J) / 6, 2) * J * F.inverse().transpose();
			}
		}
		break;
		default:
			break;
	}
}

ScalarType TetConstraint::getStressTensorAndEnergyDensity(EigenMatrix3& P, const EigenMatrix3& F)
{
	ScalarType rtnEnergy = 0;
	switch (matType)
	{
		case MATERIAL_TYPE_StVK:
		{
			EigenMatrix3 I = EigenMatrix3::Identity();
			EigenMatrix3 E = 0.5*(F.transpose()*F - I); // Green strain tensor
			P = F * (2 * mu * E + lambda * E.trace() * I); // first Piola-Kirchhoff stress tensor
			rtnEnergy = mu * E.squaredNorm() + 0.5*lambda*std::pow(E.trace(), 2); // relationship between energy and strain
			ScalarType J = F.determinant();
			if (J < 1)
			{
				P += -kappa / 24 * std::pow((1 - J) / 6, 2) * J * F.inverse().transpose();
				rtnEnergy += kappa / 12 * std::pow((1 - J) / 6, 3);
			}
		}
		break;
		default:
			break;
	}

	return rtnEnergy;
}

ScalarType TetConstraint::setMassMatrix(std::vector<SparseMatrixTriplet>& massMat)
{
	ScalarType wInv = 1.0 / w;
	for (unsigned i = 0; i != 4; ++i)
	{
		// mass
		massMat.push_back(SparseMatrixTriplet(3 * indicesList[i], 3 * indicesList[i], 0.25*w));
		massMat.push_back(SparseMatrixTriplet(3 * indicesList[i] + 1, 3 * indicesList[i] + 1, 0.25*w));
		massMat.push_back(SparseMatrixTriplet(3 * indicesList[i] + 2, 3 * indicesList[i] + 2, 0.25*w));

		// mass_1d
		// m_1d.push_back(SparseMatrixTriplet(indicesList[i], indicesList[i], 0.25*w));
	}
	return w;
}

ScalarType TetConstraint::setMassMatrix(std::vector<SparseMatrixTriplet>& m, std::vector<SparseMatrixTriplet>& m_1d)
{
	ScalarType wInv = 1.0 / w;
	for (unsigned i = 0; i != 4; i++)
	{
		// mass
		m.push_back(SparseMatrixTriplet(3 * indicesList[i], 3 * indicesList[i], 0.25*w));
		m.push_back(SparseMatrixTriplet(3 * indicesList[i] + 1, 3 * indicesList[i] + 1, 0.25*w));
		m.push_back(SparseMatrixTriplet(3 * indicesList[i] + 2, 3 * indicesList[i] + 2, 0.25*w));

		// mass_1d
		m_1d.push_back(SparseMatrixTriplet(indicesList[i], indicesList[i], 0.25*w));
	}

	return w;
}