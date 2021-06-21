#include "VKString.h"
#include "VKStringList.h"
#include <iostream>

//	// friend io functions
//std::ostream& operator<< (std::ostream& o, VKString const& str)
//{
//	return o << str.c_str();
//}
//
//std::ofstream& operator<< (std::ofstream& o, VKString const& str)
//{
//	o<<str.c_str();
//	return o;
//}
//
//std::ifstream& operator>> (std::ifstream& o, VKString & str)
//{
//	o>>str.m_Data;
//	str.readToken(o);
//	return o;
//}


// member functions

char * VKString::m_pTempBuffer = new char[VKString::TEMP_BUFFER_SIZE];

VKString::VKString()
{
	m_Data = "";
}

VKString::VKString(int number)
{
  std::ostringstream os;
  os << number ;
  m_Data = os.str();
}


VKString::VKString(long int number)
{
	std::ostringstream os;
	os << number ;
	m_Data = os.str();
}


VKString::VKString(double number, int decimals)
{
	assert(decimals==-1);	//TODO
	std::ostringstream os;
	os << number;
	m_Data = os.str();	
}

VKString::VKString(const VKString & otherStr)
{
	m_Data.assign(otherStr.m_Data);	
}

VKString::VKString(const char * strchar)
{
	m_Data.assign(strchar);	
}

VKString::VKString(std::ifstream & fileStream)
{
	m_Data = "";	
	readFromFile(fileStream);	
}

const char * VKString::c_str() const
{
	return m_Data.c_str();
}

void VKString::readFromFile(std::ifstream & fileStream)
{
	m_Data="";
	while (!fileStream.eof())
	{
		char temp = (char)fileStream.get();
		if (fileStream.eof())
			break;
		m_Data.insert(m_Data.size(), &temp, 1);
	}
}

void VKString::readLine(std::ifstream & fileStream)
{
	fileStream.getline(m_pTempBuffer, TEMP_BUFFER_SIZE, '\n' );
	m_Data = m_pTempBuffer;
}

void VKString::readToken(std::ifstream & fileStream)
{
	m_Data = "";
	fileStream>>m_Data;	
}

bool VKString::isEmpty() const
{
	return m_Data.size()==0;
}

bool VKString::operator == (const VKString & other) const
{
	return strcmp(m_Data.c_str(), other.m_Data.c_str())==0;
}

bool VKString::operator != (const VKString & other) const
{
	return strcmp(m_Data.c_str(), other.m_Data.c_str())!=0;	
}

bool VKString::operator == (const char * other) const
{
	return strcmp(m_Data.c_str(), other)==0;
}

bool VKString::operator != (const char * other) const
{
//	std::cout<<"VKSTRING Comparing: |"<<m_Data.c_str()<<"| to |"<<other<<"| val = "<<(strcmp(m_Data.c_str(), other)!=0)<<std::endl;
	return strcmp(m_Data.c_str(), other)!=0;	
}


bool VKString::operator < (const VKString & other) const
{
	return m_Data < other.m_Data;
}

VKString & VKString::operator = (const VKString & other)
{
	m_Data = other.m_Data;
	return *this;
}

VKString VKString::operator + (const VKString & other) const
{
	std::string newData = m_Data;
	newData += other.m_Data;
	return VKString(newData.c_str());
}

VKString VKString::operator + (const char * other) const
{
	return *this + VKString(other);
}

VKString VKString::operator + (const int other) const
{
	return *this + VKString(other);
}

VKString VKString::operator + (const double other) const
{
	return *this + VKString(other);
	
}

VKString VKString::operator + (const float other) const
{
	return *this + VKString(other);
}


VKString & VKString::operator += (const VKString & other)
{
	m_Data += other.m_Data;
	return (*this);
}

VKString & VKString::operator += (const char * other) 
{
	VKString otherStr(other);
	return ((*this)+=otherStr);
}

VKString & VKString::operator += (const int other)
{
	VKString otherStr(other);
	return ((*this)+=otherStr);
}

VKString & VKString::operator += (const double other)
{
	VKString otherStr(other);
	return ((*this)+=otherStr);
}

VKString & VKString::operator += (const float other)
{
	VKString otherStr(other);
	return ((*this)+=otherStr);	
}


bool VKString::startsWith(const char * other) const
{
	VKString tmp(other);
	return startsWith(tmp);
}

bool VKString::endsWith(const char * other) const
{
	VKString tmp(other);
	return endsWith(tmp);
}

bool VKString::startsWith(VKString & other) const
{
	for (int i=0; i<(int)other.m_Data.size(); i++)
	{
		if (other.m_Data[i]!=m_Data[i])
			return false;
	}
	return true;
}

bool VKString::endsWith(VKString & other) const
{
	for (int i=0; i<(int)other.m_Data.size(); i++)
	{
		if (other.m_Data[i]!=m_Data[m_Data.size() - other.m_Data.size() + i])
			return false;
	}
	return true;
}

std::string & VKString::toStdString()
{
	return m_Data;
}

int VKString::toInt(bool * ok) const
{
	std::stringstream ss(m_Data); // Could of course also have done ss("1234") directly.

	int i = 0;

	if (ok!=NULL)
		*ok = !(ss >> i).fail();
	else
		ss >> i;

	return i;
}

double VKString::toDouble(bool * ok) const
{
	std::stringstream ss(m_Data); // Could of course also have done ss("1234") directly.

	double i = 0;

	if (ok!=NULL)
		*ok = !(ss >> i).fail();
	else
		ss >> i;

	return i;

}

VKString & VKString::replace(const char * replaceMe, const VKString & replaceWith)
{
	return replace(VKString(replaceMe), replaceWith);
}

VKString & VKString::replace(const VKString & replaceMe, const char * replaceWith)
{
	return replace(replaceMe, VKString(replaceWith));
}

VKString & VKString::replace(const char * replaceMe, const char * replaceWith)
{
	return replace(VKString(replaceMe), VKString(replaceWith));
}

VKString & VKString::replace(const VKString & replaceMe, const VKString & replaceWith)
{
	for (int i=0; i<(int)m_Data.size() - (int)replaceMe.m_Data.size()+1; i++)
	{
		bool replace = true;
		for (int j=0; j<(int)replaceMe.m_Data.size() && replace; j++)
		{
			if (m_Data[i+j]!=((VKString&)replaceMe).m_Data[j])
				replace = false;
		}
		if (replace)
		{
			m_Data.erase(i, (int)replaceMe.m_Data.size());
			m_Data.insert(i, replaceWith.m_Data);
		}
	}
	return *this;
}

VKStringList VKString::split(const char * separator) const
{
	VKString tmp(separator);
	return split(tmp);
}

VKStringList VKString::split(VKString & separator) const
{
//	std::cout<<"Splitting"<<c_str()<<std::endl;
	VKStringList newStringList;

	int before = 0;
	for (int i=0; i<(int)m_Data.size() - (int)separator.m_Data.size()+1; i++)
	{
		bool same = true;
		for (int j=0; j<(int)separator.m_Data.size() && same; j++)
		{
			if (m_Data[i+j]!=separator.m_Data[j])
				same= false;
		}

		if (same)
		{
			std::string substr = m_Data.substr(before, i-before);
			if (i-before>0)
				newStringList.m_Data.push_back(substr.c_str());

			i += (int)separator.m_Data.size();
			before = i;
		}
	}
	
	std::string substr = m_Data.substr(before, (int)m_Data.size()-before);
	if ((int)m_Data.size()-before>0)
		newStringList.m_Data.push_back(substr.c_str());

	for (int i=0; i<newStringList.count(); i++)
		newStringList[i].replace(separator.c_str(), "");
	
	for (int i=0; i<newStringList.count(); i++)
	{
		if (newStringList[i]=="")
		{
			newStringList.remove(i);
			i--;
		}
	}
	
//	std::cout<<"LAST CHAR = "<<newStringList[newStringList.count()-1].c_str()<<std::endl;
	
	return newStringList;
}

int VKString::lastIndexOf(char symb) const
{
	char symbols[2];
	symbols[0] = symb;
	symbols[1] = 0;
	return lastIndexOf(VKString(symbols));
}

int VKString::lastIndexOf(const VKString & str) const
{
	int lastIndex = -1;
	for (int i=0; i<(int)m_Data.size() - (int)str.m_Data.size()+1; i++)
	{
		bool same = true;
		for (int j=0; j<(int)str.m_Data.size() && same; j++)
		{
			if (m_Data[i+j]!=str.m_Data[j])
			{
				same= false;
			}
		}
		if (same)
			lastIndex = i;
	}
	return lastIndex;
}	

VKString VKString::mid(int first, int numSymb) const
{
	VKString finalStr="";
	for (int i=first; i<first+numSymb; i++)
	{
		assert(i<(int)m_Data.size());
		finalStr.m_Data.append(1, m_Data[i]);
	}
	return finalStr;
}

VKString VKString::left(int numSymb) const
{
	VKString finalStr="";
	for (int i=(int)m_Data.size()-numSymb; i<(int)m_Data.size(); i++)
	{
		assert(i<(int)m_Data.size());		
		finalStr.m_Data.append(1, m_Data[i]);
	}
	return finalStr;	
}

int VKString::length() const
{
	return (int)m_Data.length();
}

VKString VKString::number(int val)
{
	return VKString(val);
}

VKString VKString::number(long int val)
{
	return VKString(val);
}

VKString VKString::number(double val)
{
	return VKString(val);
}

VKString VKString::number(float val)
{
	return VKString(val);
}

