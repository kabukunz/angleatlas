#pragma once
#include <sstream>
#include <string>
#include <algorithm>

void SplitPathAndFile(std::string input, std::string &path, std::string &file);
void SplitFileAndFormat(std::string &input, std::string &fileName, std::string &formatName);
int GetRandomNum(int num1, int num2);	//[num1, num2)