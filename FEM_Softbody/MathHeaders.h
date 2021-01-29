#ifndef _MATHHEADERS_H_
#define _MATHHEADERS_H_

// eigen
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"

// glm
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>


typedef float ScalarType;
#define TW_TYPE_SCALAR_TYPE TW_TYPE_FLOAT
#define EPSILON 1e-6
#define EPSILON_SQUARE 1e-12
#define LARGER_EPSILON 1e-4

// eigen vectors and matrices
typedef int IndexType;
typedef Eigen::Matrix<ScalarType, 12, 12, 0, 12, 12> EigenMatrix12;
typedef Eigen::Matrix<ScalarType, 12, 1, 0, 12, 1> EigenVector12;
typedef Eigen::Matrix<ScalarType, 4, 4, 0, 4, 4> EigenMatrix4;
typedef Eigen::Matrix<ScalarType, 4, 1, 0, 4, 1> EigenVector4;
typedef Eigen::Matrix<ScalarType, 3, 3, 0, 3, 3> EigenMatrix3;
typedef Eigen::Matrix<ScalarType, 3, 1, 0, 3, 1> EigenVector3;
typedef Eigen::Matrix<ScalarType, 2, 2, 0, 2, 2> EigenMatrix2;
typedef Eigen::Matrix<ScalarType, 2, 1, 0, 2, 1> EigenVector2;
typedef Eigen::Matrix<ScalarType, -1, 3, 0, -1, 3> EigenMatrixx3;
typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> VectorX;
typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::SparseMatrix<ScalarType> SparseMatrix;
typedef Eigen::Triplet<ScalarType, IndexType> SparseMatrixTriplet;

// eigen quaternions
typedef Eigen::AngleAxis<ScalarType> EigenAngleAxis;
typedef Eigen::Quaternion<ScalarType, Eigen::DontAlign> EigenQuaternion;

// eigen vector accessor
#define block_vector(a) block<3,1>(3*(a), 0)

// eigen 2 glm, glm 2 eigen
glm::vec3 Eigen2GLM(const EigenVector3& eigen_vector);
EigenVector3 GLM2Eigen(const glm::vec3& glm_vector);
void Eigen2GLM(const VectorX& eigen_vector, std::vector<glm::vec3>& glm_vector);
void GLM2Eigen(const std::vector<glm::vec3>& glm_vector, VectorX& eigen_vector);

// eigen make sparse identity
void EigenMakeSparseIdentityMatrix(unsigned int rows, unsigned int cols, SparseMatrix& I);

// eigen sparse matrix to triplets
void EigenSparseMatrixToTriplets(const SparseMatrix& A, std::vector<SparseMatrixTriplet>& A_triplets);

void EigenExtractDiagonalOffDiagonal(const SparseMatrix& A, VectorX& D, SparseMatrix& OD);
void EigenExtractTriangular(const SparseMatrix& A, SparseMatrix& DL, SparseMatrix& U);

// for subspace usage
void EigenExtractCols(Matrix& dst, const Matrix& src, const std::vector<unsigned int>& indices);
void EigenVectorExtractElements(VectorX& dst, const VectorX& src, const std::vector<unsigned int>& indices);
void EigenExtractCols(Matrix& dst, const Matrix& src, unsigned int d, bool from_front_end = true);
void EigenVectorExtractElements(VectorX& dst, const VectorX& src, unsigned int d, bool from_front_end = true);


#endif
