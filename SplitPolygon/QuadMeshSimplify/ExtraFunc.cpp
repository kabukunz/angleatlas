
#include "ExtraFunc.h"
#include <cstdlib>
#include <ctime>

void SplitPathAndFile(std::string input, std::string &path, std::string &file)
{
	std::string str0("/"), str1("\\");

	auto it0 = std::find_end(input.begin(), input.end(), str0.begin(), str0.end());
	auto it1 = std::find_end(input.begin(), input.end(), str1.begin(), str1.end());

	std::string::iterator finalIt;
	if (it0 == input.end() && it1 == input.end())
	{
		path = std::string("./");
		file = input;
		return;
	}
	else if (it0 == input.end())
		finalIt = it1;
	else if (it1 == input.end())
		finalIt = it0;
	else
	{
		int dis0 = std::distance(input.begin(), it0);
		int dis1 = std::distance(input.begin(), it1);
		if (dis0 < dis1)
			finalIt = it1;
		else
			finalIt = it0;
	}

	path.clear();
	file.clear();

	path.insert(path.end(), input.begin(), finalIt + 1);
	file.insert(file.end(), finalIt + 1, input.end());
}

void SplitFileAndFormat(std::string &input, std::string &fileName, std::string &formatName)
{
	std::string str(".");
	auto it = std::find_end(input.begin(), input.end(), str.begin(), str.end());
	if (it == input.end())
		return;

	fileName.clear();
	formatName.clear();
	fileName.insert(fileName.end(), input.begin(), it + 1);
	fileName.pop_back();
	formatName.insert(formatName.end(), it + 1, input.end());
}

int GetRandomNum(int num1, int num2)
{
	//std::srand((unsigned int)std::time(NULL));
	return (rand() % (num2 - num1)) + num1;
}