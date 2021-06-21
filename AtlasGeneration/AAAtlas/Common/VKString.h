#include <stdio.h>
#include <vector>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>

#ifndef __VKSTRING_H
#define __VKSTRING_H

class VKStringList;

class VKString
{
	public:
//		friend std::ostream& operator<< (std::ostream& o, VKString const& str);
//		friend std::ofstream& operator<< (std::ofstream& o, VKString const& str);
//		friend std::ifstream& operator>> (std::ifstream& o, VKString & str);
	
		VKString();
		VKString(const int number);
		VKString(const long int number);	
		VKString(const double number, int decimals = -1);
		VKString(const VKString & otherStr);
		VKString(const char * strchar);
		VKString(std::ifstream & fileStream);

	// some things I find useful
		const char * c_str() const;
		void readFromFile(std::ifstream & fileStream);
		void readLine(std::ifstream & fileStream);	
		void readToken(std::ifstream & fileStream);

	// Qt interface
		bool operator < (const VKString & other) const;
		VKString & operator = (const VKString & other);
		bool operator == (const VKString & other) const;
		bool operator != (const VKString & other) const;
		bool operator == (const char * other) const;
		bool operator != (const char * other) const;

		VKString operator + (const VKString & other) const;
		VKString operator + (const char * other) const;	
		VKString operator + (int other) const;
		VKString operator + (double other) const;			
		VKString operator + (float other) const;

		VKString & operator += (const VKString & other);
		VKString & operator += (const char * other);
		VKString & operator += (const int other);
		VKString & operator += (const double other);
		VKString & operator += (const float other);	
	
		bool isEmpty() const;
		bool startsWith(VKString & other) const;
		bool endsWith(VKString & other) const;
		bool startsWith(const char * other) const;
		bool endsWith(const char * other) const;
		VKString & replace(const VKString & replaceMe, const VKString & replaceWith);
		VKString & replace(const char * replaceMe, const char * replaceWith);
		VKString & replace(const char * replaceMe, const VKString & replaceWith);
		VKString & replace(const VKString & replaceMe, const char * replaceWith);
		std::string & toStdString();
		int toInt(bool * ok=NULL) const;
		double toDouble(bool * ok=NULL) const;
		VKStringList split(VKString & separator) const;
		VKStringList split(const char * separator) const;
		int lastIndexOf(char symb) const;	
		int lastIndexOf(const VKString & str) const;
		VKString mid(int first, int numSymb) const;
		VKString left(int numSymb) const;	
		int length() const;

		static VKString number(int val);
		static VKString number(long int val);		
		static VKString number(double val);	
		static VKString number(float val);		
	
	protected:
		std::string m_Data;
		static char * m_pTempBuffer;
		static const int TEMP_BUFFER_SIZE = 10000;
};

#endif

