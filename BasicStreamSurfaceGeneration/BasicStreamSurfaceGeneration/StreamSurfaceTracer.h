#pragma once

#include <iostream>
#include <vector>
#include "Utilities.h"
#include "Steady3D.h"
#include "StreamlineTracer.h"

using namespace Util;

class StreamSurfaceTracer
{
public:
	StreamSurfaceTracer(Steady3D &steady3D);
	~StreamSurfaceTracer();
	void setStepSize(const float& _stepSize);
	void setLength(const int& _length);
	void setForward(const bool& _forward);

	void trace(vector<Vector3f> seedingCurve, int streamlineLength, vector<Vector3f> &vertices, vector<Vector3i>& faces);
	void addTriangle(vector<Vector3f> &vertices, vector<Vector3i> &faces, vector<Vector3f> points);
	void getSurfaceFromTwoPoints(Vector3f p1, Vector3f p2, int streamlineLength, vector<Vector3f> &vertices, vector<Vector3i> &faces);


private:
	
	StreamlineTracer* mStreamlineTracer;
	int length;
	bool forward;
	float stepSize;
};

