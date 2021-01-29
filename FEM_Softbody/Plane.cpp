#include "Plane.h"

#define COLLISION_EPSILON 1e-3

//----------Plane Class----------//
void Plane::init()
{
	positionsList.clear();
	normalsList.clear();
	indicesList.clear();

	glm::vec3 center(0.0, 0.0, 0.0);
	glm::vec3 localX, localZ;
	localX = glm::cross(normal, glm::vec3(0.0, 0.0, 1.0));
	if (glm::length(localX) < 0.00001f)
		localX = glm::cross(normal, glm::vec3(1.0, 0.0, 0.0));
	localX = glm::normalize(localX);
	localZ = glm::normalize(glm::cross(localX, normal));

	glm::vec3 mat_color(1.0);
	unsigned int slice = 24;

	glm::vec3 vertex(center);
	positionsList.push_back(center);
	normalsList.push_back(normal);

	float delta = 360.0 / slice;
	float radius = 100.0;
	glm::vec3 localPos;
	for (float theta = 0.0; theta < 359.99; theta += delta)
	{
		localPos.x = radius * cos(glm::radians(theta));
		localPos.z = radius * sin(glm::radians(theta));

		vertex = localPos.x * localX - localPos.z * localZ + center;

		positionsList.push_back(vertex);
		normalsList.push_back(normal);
	}
	for (unsigned int i = 0; i < slice - 1; ++i)
	{
		indicesList.push_back(0);
		indicesList.push_back(i + 1);
		indicesList.push_back(i + 2);
	}
	indicesList.push_back(0);
	indicesList.push_back(slice);
	indicesList.push_back(1);
}

bool Plane::intersectionTest(const EigenVector3& p, EigenVector3& normal, ScalarType& dist)
{
	ScalarType height = centerPos[1];
	dist = p(1) - height - COLLISION_EPSILON;
	normal = EigenVector3(normal[0], normal[1], normal[2]);

	if (dist < 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}
