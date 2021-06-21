#pragma once
#include <string>
#include <vector>
namespace BaseDataStructure
{
	class QuadMesh;

	class QuadMeshIO
	{
	public:
		QuadMeshIO();
		~QuadMeshIO();

		void ReadQuadMesh(QuadMesh *quadMesh, std::string fileName);
		void WriteQuadMesh(QuadMesh *quadMesh, std::string fileName);

		void SplitString(const std::string& s, std::vector<std::string>& v, const std::string& c);
	};
}

