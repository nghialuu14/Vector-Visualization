#pragma once


#include <iostream>
#include <stdio.h>
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <cstring>
#include <vector>
#include <unordered_set>
#include <queue>
#include <fstream>
#include <assert.h>
#include <string.h>
#include <string>
#include <cstdlib>
#include <omp.h>
#include <algorithm>
#include <time.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

#define PI 3.14159265

using namespace std;
using namespace Eigen;



namespace Util
{
	
	const float BIG_FLOAT = 99999999.0;
	const float SMALL_FLOAT = -99999999.0;
	class Point2D
	{
	public:
		float x, y;
		Point2D() {};
		Point2D(float xx, float yy) { x = xx; y = yy; };
		~Point2D(void) {};
	};


	class Point3D : public Point2D
	{
	public:
		float z;
		Point3D(void) {};
		Point3D(float xx, float yy, float zz) { x = xx; y = yy; z = zz; };
		Point3D(float xx, float yy) { x = xx; y = yy; z = 0.; };
		~Point3D(void) {};
	};

	class Point4D : public Point3D
	{
	public:
		float t;
		Point4D(void) {};
		Point4D(float xx, float yy, float zz) { x = xx; y = yy; z = zz; };
		Point4D(float xx, float yy) { x = xx; y = yy; z = 0.; t = 0.; };
		Point4D(float xx, float yy, float zz, float tt) { x = xx; y = yy; z = zz; t = tt; };
		~Point4D(void) {};
	};

	class Triangle {
	public:
		Point3D p1, p2, p3;
		Triangle(void) {};
		Triangle(Point3D pp1, Point3D pp2, Point3D pp3) { p1 = pp1; p2 = pp2; p3 = pp3; };
		~Triangle(void) {};
	};

	void savePointToVTKFile(vector<Vector3f> results, string corelineFileName);

	void saveLinesToVTKFile(vector<vector<Vector4f>> results, string corelineFileName);
	void saveLinesToVTKFile(vector<vector<Vector3f>> results, string corelineFileName);

	void generateObjSurfaceFromLines(vector<vector<Vector3f>> lines, string outputFileName);

	void outputSurfaceToOBJ(vector<Vector3f> &vertices, vector<Vector3i> &faces, string fileName);

}