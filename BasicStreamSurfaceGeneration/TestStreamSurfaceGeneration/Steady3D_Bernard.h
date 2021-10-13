#pragma once
#include <Steady3D.h>

/* Bernard dataset */
class Steady3D_Bernard :
	public Steady3D
{
public:
	Steady3D_Bernard(string inputFileName);
	~Steady3D_Bernard();
	void loadBernardDataset(string fileName);
	const char* locateSubString(const char* inputString, const char* searchString);
};

