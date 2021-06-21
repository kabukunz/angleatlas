#include "VKStringList.h"

VKStringList::VKStringList()
{
}
VKStringList::VKStringList(const VKStringList & other)
{
	m_Data = other.m_Data;
}


void VKStringList::push_back(VKString obj)
{
	m_Data.push_back(obj);
}

VKString & VKStringList::operator[](int i)
{
	return m_Data[i];
}

void VKStringList::clear()
{
	m_Data.clear();
}

int VKStringList::count() const
{
	return (int)m_Data.size();
}

bool VKStringList::contains(const VKString & trg) const
{
	for (int i=0; i<(int)m_Data.size(); i++)
	{
		if (m_Data[i]==trg)
			return true;
	}
	return false;
}

VKString VKStringList::join(const VKString & sep) const
{
	VKString joinedStr;
	for (int i=0; i<(int)m_Data.size(); i++)
		joinedStr += m_Data[i] + (i<((int)m_Data.size()-1) ? sep : "");

	return joinedStr;
}

void VKStringList::remove(int elementID)
{
	assert(elementID<(int)m_Data.size() && elementID>=0);
	m_Data.erase(m_Data.begin()+elementID);
}

void VKStringList::remove(const VKString & trg)
{
	for (int i=0; i<(int)m_Data.size(); i++)
	{
		if (m_Data[i]==trg)
		{
			remove(i);
			break;
		}
	}
}
