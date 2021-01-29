
#ifndef _MESH_LOADER_H_
#define _MESH_LOADER_H_

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <glm/glm.hpp>

class MeshLoader {
public:

	struct Tet {
		unsigned int id1, id2, id3, id4;  // four vertices that construct one tet
		Tet() {}
		Tet(int a, int b, int c, int d) : id1(a), id2(b), id3(c), id4(d) {}
		void idMinusMinus() { id1--; id2--; id3--; id4--; }
	};

	MeshLoader();
	MeshLoader(std::string filename);
	virtual ~MeshLoader();

	inline bool info() { return loadSuccess; }

	//Vertices, edges, and faces information
	std::vector<glm::vec3> verticesList;
	std::vector<Tet> tetsList;

	bool loadSuccess;

private:
	std::vector<std::string> dataStringList;  // intermediate vector storing strings from file
	void splitString(const std::string & src, std::vector<std::string>& v, const std::string & split);
	std::vector<std::string> splitString(const std::string & src, const std::string & split);
};

#endif