#include "mathHeaders.h"
#include <vector>

glm::vec3 Eigen2GLM(const EigenVector3& eigen_vector)
{
	return glm::vec3(eigen_vector[0], eigen_vector[1], eigen_vector[2]);
}
EigenVector3 GLM2Eigen(const glm::vec3& glm_vector)
{
	return EigenVector3(glm_vector[0], glm_vector[1], glm_vector[2]);
}
void Eigen2GLM(const VectorX& eigen_vector, std::vector<glm::vec3>& glm_vector)
{
	int size = glm_vector.size();
	assert(eigen_vector.size() == glm_vector.size() * 3);
	memcpy(glm_vector.data(), eigen_vector.data(), sizeof(float)*size * 3);
}
void GLM2Eigen(const std::vector<glm::vec3>& glm_vector, VectorX& eigen_vector)
{
	int size = glm_vector.size();
	assert(eigen_vector.size() == glm_vector.size() * 3);
	memcpy(eigen_vector.data(), glm_vector.data(), sizeof(float)*size * 3);
}

void EigenMakeSparseIdentityMatrix(unsigned int rows, unsigned int cols, SparseMatrix& I)
{
	assert(rows == cols);
	std::vector<SparseMatrixTriplet> triplets;
	for (unsigned int i = 0; i != rows; i++)
	{
		triplets.push_back(SparseMatrixTriplet(i, i, 1));
	}
	I.resize(rows, cols);
	I.setFromTriplets(triplets.begin(), triplets.end());
}

void EigenSparseMatrixToTriplets(const SparseMatrix& A, std::vector<SparseMatrixTriplet>& A_triplets)
{
	A_triplets.clear();
	A_triplets.reserve(A.nonZeros());
	for (unsigned int col = 0; col != A.outerSize(); ++col)
	{
		for (SparseMatrix::InnerIterator it(A, col); it; ++it)
		{
			A_triplets.push_back(SparseMatrixTriplet(it.row(), col, it.value()));
		}
	}
}

void EigenExtractDiagonalOffDiagonal(const SparseMatrix& A, VectorX& D, SparseMatrix& OD)
{
	unsigned rows = A.rows();
	unsigned cols = A.cols();
	assert(rows == cols);

	std::vector<SparseMatrixTriplet> A_triplets; A_triplets.clear();
	//std::vector<SparseMatrixTriplet> M_triplets; M_triplets.clear();
	std::vector<SparseMatrixTriplet> N_triplets; N_triplets.clear();

	D.resize(rows);
	OD.resize(rows, cols);

	D.setZero();

	EigenSparseMatrixToTriplets(A, A_triplets);

	ScalarType ind, val;
	for (std::vector<SparseMatrixTriplet>::iterator it = A_triplets.begin(); it != A_triplets.end(); it++)
	{
		if ((ind = it->col()) == it->row()) // digonal
		{
			val = it->value();
			D(ind) += val;
		}
		else
		{
			N_triplets.push_back((*it));
		}
	}

	OD.setFromTriplets(N_triplets.begin(), N_triplets.end());
}
void EigenExtractTriangular(const SparseMatrix& A, SparseMatrix& DL, SparseMatrix& U)
{
	unsigned rows = A.rows();
	unsigned cols = A.cols();
	assert(rows == cols);

	std::vector<SparseMatrixTriplet> A_triplets; A_triplets.clear();
	std::vector<SparseMatrixTriplet> M_triplets; M_triplets.clear();
	std::vector<SparseMatrixTriplet> N_triplets; N_triplets.clear();

	DL.resize(rows, cols);
	U.resize(rows, cols);

	EigenSparseMatrixToTriplets(A, A_triplets);

	ScalarType ind, val;
	for (std::vector<SparseMatrixTriplet>::iterator it = A_triplets.begin(); it != A_triplets.end(); it++)
	{
		if ((ind = it->col()) <= it->row()) // digonal
		{
			M_triplets.push_back((*it));
		}
		else
		{
			N_triplets.push_back((*it));
		}
	}

	DL.setFromTriplets(M_triplets.begin(), M_triplets.end());
	U.setFromTriplets(N_triplets.begin(), N_triplets.end());
}

// for subspace usage
void EigenExtractCols(Matrix& dst, const Matrix& src, const std::vector<unsigned int>& indices)
{
	unsigned int size = indices.size();
	assert(size > 0);

	dst.resize(src.rows(), size);

	for (unsigned int i = 0; i != size; ++i)
	{
		dst.col(i) = src.col(indices[i]);
	}
}

void EigenVectorExtractElements(VectorX& dst, const VectorX& src, const std::vector<unsigned int>& indices)
{
	unsigned int size = indices.size();
	assert(size > 0);

	dst.resize(size);

	for (unsigned int i = 0; i != size; ++i)
	{
		dst[i] = src[indices[i]];
	}
}

void EigenExtractCols(Matrix& dst, const Matrix& src, unsigned int d, bool from_front_end)
{
	std::vector<unsigned int> indices; indices.clear();
	for (unsigned int i = 0; i != d; ++i)
	{
		if (from_front_end)
		{
			indices.push_back(i);
		}
		else
		{
			indices.push_back(src.rows() - i - 1);
		}
	}

	EigenExtractCols(dst, src, indices);
}
void EigenVectorExtractElements(VectorX& dst, const VectorX& src, unsigned int d, bool from_front_end)
{
	std::vector<unsigned int> indices; indices.clear();
	for (unsigned int i = 0; i != d; ++i)
	{
		if (from_front_end)
		{
			indices.push_back(i);
		}
		else
		{
			indices.push_back(src.size() - i - 1);
		}
	}

	EigenVectorExtractElements(dst, src, indices);
}