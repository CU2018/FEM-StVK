#include "MeshLoader.h"

MeshLoader::MeshLoader()
{
	verticesList.clear();
	tetsList.clear();
	loadSuccess = false;
}

MeshLoader::MeshLoader(std::string filename)
{
	std::fstream in;
	in.open(filename, std::ios::in);
	if (!in.is_open()) {
		printf("Error opening the file\n");
		exit(1);
	}

	verticesList.clear();
	tetsList.clear();
	loadSuccess = false;

	std::string buffer;
	int lineNum = 0;
	while (std::getline(in, buffer)) {
		// string to char *
		const std::string firstSplit = ":";
		const std::string secondSplit = ",";
		std::string dataLine = splitString(buffer, firstSplit)[1];
		dataStringList = splitString(dataLine, secondSplit);
		if (lineNum == 0)  // first line of the input file: position
		{
			assert(dataStringList.size() % 3 == 0);
			for (unsigned int i = 0; i < dataStringList.size(); i+=3)
			{
				std::cout << dataStringList[i] << std::endl;
				glm::vec3 vertexPos(std::stof(dataStringList[i]),
									std::stof(dataStringList[i+1]),
									std::stof(dataStringList[i+2]));
				printf("vec %d: %f, %f, %f\n", i/3, vertexPos[0], vertexPos[1], vertexPos[2]);
				verticesList.push_back(vertexPos);
			}
			assert(verticesList.size() == dataStringList.size() / 3);
		}
		else  // second line of the input file: vertices tet
		{
			assert(dataStringList.size() % 4 == 0);
			for (unsigned int i = 0; i < dataStringList.size(); i += 4)
			{
				std::cout << dataStringList[i] << std::endl;
				Tet newTet(std::stoi(dataStringList[i]), std::stoi(dataStringList[i + 1]),
						std::stoi(dataStringList[i + 2]), std::stoi(dataStringList[i + 3]));
				printf("tet %d: %d, %d, %d, %d\n", i/4, newTet.id1, newTet.id2, newTet.id3, newTet.id4);
				tetsList.push_back(newTet);
			}
			assert(tetsList.size() == dataStringList.size() / 4);
		}
		++lineNum;
	}
	in.close();
	loadSuccess = true;
	printf("Read File Done!\n");
}

MeshLoader::~MeshLoader()
{
	verticesList.clear();
	tetsList.clear();
}

void MeshLoader::splitString(const std::string & src, std::vector<std::string>& v, const std::string & split)
{
	std::string::size_type pos1, pos2;
	pos2 = src.find(split);
	pos1 = 0;
	while (std::string::npos != pos2)
	{
		v.push_back(src.substr(pos1, pos2 - pos1));
		pos1 = pos2 + split.size();
		pos2 = src.find(split, pos1);
	}
	if (pos1 != src.length())
		v.push_back(src.substr(pos1));
}

std::vector<std::string> MeshLoader::splitString(const std::string & src, const std::string & split)
{
	std::vector<std::string> _ret = std::vector<std::string>();
	splitString(src, _ret, split);
	return _ret;
}
