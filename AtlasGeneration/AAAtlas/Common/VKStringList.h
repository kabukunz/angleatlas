#include "VKString.h"
#ifndef __VKSTRINGLIST_H
#define __VKSTRINGLIST_H

class VKStringList
{
	friend class VKString;
	public:
		VKStringList();
		VKStringList(const VKStringList & other);
	
	// Qt stuff
		void push_back(VKString obj);
		VKString & operator[](int i);
		void clear();
		int count() const;
		bool contains(const VKString & trg) const;
		VKString join(const VKString&  sep) const;
		void remove(int elementID);
		void remove(const VKString & trg);	


	protected:
		std::vector<VKString> m_Data;
};

#endif

