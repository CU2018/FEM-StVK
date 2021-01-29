#ifndef _PLANE_H_
#define _PLANE_H_

#include "MathHeaders.h"


class Primitive
{
public:
	Primitive(glm::vec3 pos) : m_pos(pos){};
	Primitive(const Primitive& other) :
		positionsList(other.positionsList), normalsList(other.normalsList), indicesList(other.indicesList)
	{
	}
	virtual ~Primitive()
	{
		positionsList.clear();
		normalsList.clear();
		indicesList.clear();
	}

	virtual void moveTo(const glm::vec3& target) { m_pos = target; }
	virtual bool intersectionTest(const EigenVector3& p, EigenVector3& normal, ScalarType& dist) { return false; }

	inline std::vector<glm::vec3>& GetPositions() { return positionsList; }
	inline std::vector<unsigned short>& GetTriangulation() { return indicesList; }

public:
	glm::vec3 m_pos;

protected:
	std::vector<glm::vec3> positionsList, normalsList;
	std::vector<unsigned short> indicesList;
};

class Plane
{
public:
	Plane() : normal(glm::vec3(0.0, 1.0, 0.0)){ init(); };
	Plane(const glm::vec3 norm, float value) : normal(norm), centerPos(glm::vec3(0.0, value, 0.0)) { init(); };
	Plane(const Plane& other) : normal(other.normal) { init(); };
	~Plane() { 
		positionsList.clear();
		indicesList.clear();
	};

	void init();
	bool intersectionTest(const EigenVector3& p, EigenVector3& normal, ScalarType& dist);

public:
	glm::vec3 centerPos;

protected:	
	std::vector<glm::vec3> positionsList, normalsList;
	glm::vec3 normal; // plane has the same single normal
	std::vector<unsigned short> indicesList;
};

#endif