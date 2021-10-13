#include "Steady3D_Bernard.h"



Steady3D_Bernard::Steady3D_Bernard(string inputFileName)
{
	// Set the configuration for the dataset
	dims[2] = 64;
	dims[1] = 32;
	dims[0] = 128;
	xRange[0] = 0;
	xRange[1] = 4;
	yRange[0] = 0;
	yRange[1] = 1;
	zRange[0] = 0;
	zRange[1] = 2;

	spacing[0] = (xRange[1] - xRange[0]) / (dims[0] - 1);
	spacing[1] = (yRange[1] - yRange[0]) / (dims[1] - 1);
	spacing[2] = (zRange[1] - zRange[0]) / (dims[2] - 1);

	origin[0] = origin[1] = origin[2] = 0;
	vfData = new float[dims[2] * dims[1] * dims[0] * DIMENSION_SIZE];
	loadBernardDataset(inputFileName);
}

/**
* Find a substring in a string
* @return a pointer to the first occurrence of searchString in the inputString
*/
const char* Steady3D_Bernard::locateSubString(const char* inputString, const char* searchString)
{
	const char* foundLoc = strstr(inputString, searchString);
	if (foundLoc) return foundLoc + strlen(searchString);
	return inputString;
}

/**
Read the raw Bernard data file and store in the uniform rectilinear grid array loadData
*/
void Steady3D_Bernard::loadBernardDataset(string fileName) {

	FILE* pFile2 = fopen(fileName.c_str(), "rb");
	if (pFile2 == NULL) { fputs("File error", stderr); exit(1); }
	fseek(pFile2, 0L, SEEK_SET);

	//Read the data
	// - how much to read
	const size_t NumToRead = dims[0] * dims[1] * dims[2] * 3;
	// - prepare memory; use malloc() if you're using pure C
	unsigned char* pData = new unsigned char[NumToRead];
	if (pData)
	{
		// - do it
		const size_t ActRead = fread((void*)pData, sizeof(unsigned char), NumToRead, pFile2);
		// - ok?
		if (NumToRead != ActRead)
		{
			printf("Something went wrong while reading the binary data section.\nPremature end of file?\n");
			delete[] pData;
			fclose(pFile2);
			return;
		}

		//Test: Print all data values
		//Note: Data runs x-fastest, i.e., the loop over the x-axis is the innermost
		//printf("\nPrinting all values in the same order in which they are in memory:\n");
		int Idx(0);
		float tmp[3];
		for (int k = 0; k < dims[2]; k++)
		{
			for (int j = 0; j < dims[1]; j++)
			{
				for (int i = 0; i < dims[0]; i++)
				{
					//Note: Random access to the value (of the first component) of the grid point (i,j,k):
					// pData[((k * yDim + j) * xDim + i) * NumComponents]
					//assert(pData[((k * resolutionY + j) * resolutionX + i) * 3] == pData[Idx * 3]);


					for (int c = 0; c < 3; c++)
					{
						tmp[c] = (float)pData[Idx * 3 + c] / 255. - 0.5;
					}
					float dist = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]);
					for (int c = 0; c < 3; c++)
					{
						vfData[getTupleIdx(i,j,k,c)] = tmp[c] / dist;
					}

					Idx++;
				}
			}
		}

		delete[] pData;

	}

	fclose(pFile2);

}


Steady3D_Bernard::~Steady3D_Bernard()
{
	delete[] vfData;
}
