#define _USE_MATH_DEFINES
#include <iostream>
#include <omp.h>
#include <Utilities.h>
#include <Steady3D.h>
#include <StreamlineTracer.h>
#include <StreamSurfaceTracer.h>
#include <cmath>
#include <fstream>
#include <string>

#include "Steady3D_Bernard.h"


using namespace std;

struct point {
	float x;
	float y;
	float z;
};

struct vect {
	float i;
	float j;
	float k;
};

int nSeedPoints;
vector<point> seedPoints;

enum TestCases {
	SingleStreamLine_TestCase = 0,
	MultipleStreamlines_TestCase,
	StreamSurfaceWithSeedingCurve_TestCase
};

void testSingleStreamline();
void testMultipleStreamline();
void testStreamSurfaceWithSeedingCurve();


void extractBernard(vector<Vector3f> &seedPoints) {
	ifstream bernardFile("../../sample_data/bernard3d.vtk");
	string tmp;
	for (int i = 0; i < 9; i++) getline(bernardFile,tmp);
	for (int i = 0; i < 262144; i++) {
		float x, y, z;
		bernardFile >> x >> y >> z;
		seedPoints.push_back(Vector3f(x, y, z));
	}
	bernardFile.close();
};

void printSeedpoints(vector<Vector3f> seedingCurve, string corelineFileName) {
	FILE* pFile;
	pFile = fopen(corelineFileName.c_str(), "w");
	fprintf(pFile, "# vtk DataFile Version 3.0\n");
	fprintf(pFile, "VortexCoreLine\n");
	fprintf(pFile, "ASCII\n");
	fprintf(pFile, "DATASET POLYDATA\n");
	fprintf(pFile, "POINTS %d float\n", seedingCurve.size());
	for (int i = 0; i < seedingCurve.size(); i++) {
		fprintf(pFile, "%f %f %f\n", seedingCurve[i](0), seedingCurve[i](1), seedingCurve[i](2));

	}
};

/* Generate a single streamline from a predefined seeding point*/
void testSingleStreamline() {
	std::cout << "Start testSingleStreamline()..." << std::endl;
	// Step 1: Load a vector field. In this example, we use the Bernard dataset
	Steady3D_Bernard bernardVectorField("../../sample_data/bernard.raw");

	// Step 2: Initialize the streamline tracer
	StreamlineTracer streamlineTracer(bernardVectorField);
	streamlineTracer.setForward(true);
	streamlineTracer.setStepSize(0.1);
	int streamline_length = 50;

	// Step 3: Create a seeding point
	Vector3f sampleSeedingPoint(2.51969, 0.548387, 0.444444);

	// Step 4: Trace the streamline
	vector<Vector3f> streamline = streamlineTracer.trace(sampleSeedingPoint, streamline_length);

	vector<Vector4f> streamline2;
	Vector4f p0 = Vector4f(streamline[0](0), streamline[0](1), streamline[0](2), 0);
	streamline2.push_back(p0);
	for (int i = 1; i < streamline.size() - 1; i++) {
		float x1 = streamline[i - 1](0);
		float y1 = streamline[i - 1](1);
		float z1 = streamline[i - 1](2);
		float x2 = streamline[i](0);
		float y2 = streamline[i](1);
		float z2 = streamline[i](2);
		float x3 = streamline[i + 1](0);
		float y3 = streamline[i + 1](1);
		float z3 = streamline[i + 1](2);
		vect vt1; vect vt2;
		vt1.i = x2 - x1; vt1.j = y2 - y1; vt1.k = z2 - z1;

		vt2.i = x3 - x2; vt2.j = y3 - y2; vt2.k = z3 - z2;
		float tmp1 = vt1.i * vt2.i + vt1.j * vt2.j + vt1.k * vt2.k;
		float tmp2 = sqrt(vt1.i * vt1.i + vt1.j * vt1.j + vt1.k * vt1.k) * sqrt(vt2.i * vt2.i + vt2.j * vt2.j + vt2.k * vt2.k);
		float theta = acos(tmp1 / tmp2)*180/M_PI;
		Vector4f p0 = Vector4f(x2, y2, z2, theta);
		streamline2.push_back(p0);
	}
	int last = streamline.size();
	p0 = Vector4f(streamline[last - 1](0), streamline[last - 1](1), streamline[last - 1](2), 0);
	streamline2.push_back(p0);

	// Step 5: Save the streamline to a vtk file
	vector<vector<Vector4f>> outputLines;
	outputLines.push_back(streamline2);
	Util::saveLinesToVTKFile(outputLines, "../../output/singleStreamline.vtk");

	std::cout << "Finish testSingleStreamline()!" << std::endl;
	
}


void connect(vector<Vector3f>& seedingCurve, float x1, float y1, float z1, float x2, float y2, float z2, int n) {
	float tmpx = abs(x1 - x2) / n;
	float tmpy = abs(y1 - y2) / n;
	float tmpz = abs(z1 - z2) / n;
	float x = x1; float y = y1; float z = z1;
	for (int i = 0; i <= n; i++) {
		seedingCurve.push_back(Vector3f(x, y, z));
		if (x1 < x2) x += tmpx; else x -= tmpx;
		if (y1 < y2) y += tmpy; else y -= tmpy;
		if (z1 < z2) z += tmpz; else z -= tmpz;
	}

};


void fig1(vector<Vector3f>& seedingCurve) {
	seedingCurve.push_back(Vector3f(0.51969, 0.548387, 0.444444));
	seedingCurve.push_back(Vector3f(0.51432, 0.504238, 0.486014));
	seedingCurve.push_back(Vector3f(0.50522, 0.478441, 0.534027));
	seedingCurve.push_back(Vector3f(0.49611, 0.452904, 0.581806));
	seedingCurve.push_back(Vector3f(0.48316, 0.42857, 0.6298));
	seedingCurve.push_back(Vector3f(0.46891, 0.405556, 0.677777));
	seedingCurve.push_back(Vector3f(0.45691, 0.388489, 0.725037));
	seedingCurve.push_back(Vector3f(0.44247, 0.365658, 0.772294));
	seedingCurve.push_back(Vector3f(0.42122, 0.343246, 0.819314));
	seedingCurve.push_back(Vector3f(0.40604, 0.327284, 0.866114));
};

void fig2(vector<Vector3f>& seedingCurve) {
	seedingCurve.push_back(Vector3f(1.51969, 0.548387, 0.444444));
	seedingCurve.push_back(Vector3f(1.51432, 0.504238, 0.486014));
	seedingCurve.push_back(Vector3f(1.50522, 0.478441, 0.534027));
	seedingCurve.push_back(Vector3f(1.49611, 0.452904, 0.581806));
	seedingCurve.push_back(Vector3f(1.48316, 0.42857, 0.6298));
	seedingCurve.push_back(Vector3f(1.46891, 0.405556, 0.677777));
	seedingCurve.push_back(Vector3f(1.45691, 0.388489, 0.725037));
	seedingCurve.push_back(Vector3f(1.44247, 0.365658, 0.772294));
	seedingCurve.push_back(Vector3f(1.42122, 0.343246, 0.819314));
	seedingCurve.push_back(Vector3f(1.40604, 0.327284, 0.866114));
};

void fig3(vector<Vector3f>& seedingCurve) {
	/*connect(seedingCurve, 2.51969, 0.548387, 0.444444, 2.51432, 0.504238, 0.486014, 10);
	connect(seedingCurve, 2.51432, 0.504238, 0.486014, 2.50522, 0.478441, 0.534027, 10);
	connect(seedingCurve, 2.50522, 0.478441, 0.534027, 2.49611, 0.452904, 0.581806, 10);
	connect(seedingCurve, 2.49611, 0.452904, 0.581806, 2.48316, 0.42857, 0.6298, 10);
	connect(seedingCurve, 2.48316, 0.42857, 0.6298, 2.46891, 0.405556, 0.677777, 10);
	connect(seedingCurve, 2.46891, 0.405556, 0.677777, 2.45691, 0.388489, 0.725037, 10);
	connect(seedingCurve, 2.45691, 0.388489, 0.725037, 2.44247, 0.365658, 0.772294, 10);
	connect(seedingCurve, 2.44247, 0.365658, 0.772294, 2.42122, 0.343246, 0.819314, 10);
	connect(seedingCurve, 2.42122, 0.343246, 0.819314, 2.40604, 0.327284, 0.866114, 10);*/

	seedingCurve.push_back(Vector3f(2.51969, 0.548387, 0.444444));
	seedingCurve.push_back(Vector3f(2.51432, 0.504238, 0.486014));
	seedingCurve.push_back(Vector3f(2.50522, 0.478441, 0.534027));
	seedingCurve.push_back(Vector3f(2.49611, 0.452904, 0.581806)); 
	seedingCurve.push_back(Vector3f(2.48316, 0.42857, 0.6298));
	seedingCurve.push_back(Vector3f(2.46891, 0.405556, 0.677777));
	seedingCurve.push_back(Vector3f(2.45691, 0.388489, 0.725037));
	seedingCurve.push_back(Vector3f(2.44247, 0.365658, 0.772294));
	seedingCurve.push_back(Vector3f(2.42122, 0.343246, 0.819314));
	seedingCurve.push_back(Vector3f(2.40604, 0.327284, 0.866114));
};

void fig4(vector<Vector3f> &seedingCurve) {
	seedingCurve.push_back(Vector3f(3.51969, 0.548387, 0.444444));
	seedingCurve.push_back(Vector3f(3.51432, 0.504238, 0.486014));
	seedingCurve.push_back(Vector3f(3.50522, 0.478441, 0.534027));
	seedingCurve.push_back(Vector3f(3.49611, 0.452904, 0.581806));
	seedingCurve.push_back(Vector3f(3.48316, 0.42857, 0.6298));
	seedingCurve.push_back(Vector3f(3.46891, 0.405556, 0.677777));
	seedingCurve.push_back(Vector3f(3.45691, 0.388489, 0.725037));
	seedingCurve.push_back(Vector3f(3.44247, 0.365658, 0.772294));
	seedingCurve.push_back(Vector3f(3.42122, 0.343246, 0.819314));
	seedingCurve.push_back(Vector3f(3.40604, 0.327284, 0.866114));
};

void fig5(vector<Vector3f>& seedingCurve) {
	seedingCurve.push_back(Vector3f(0.51969, 0.548387, 1.444444));
	seedingCurve.push_back(Vector3f(0.51432, 0.504238, 1.494444));
	seedingCurve.push_back(Vector3f(0.50522, 0.478441, 1.544444));
	seedingCurve.push_back(Vector3f(0.49611, 0.452904, 1.594444));
	seedingCurve.push_back(Vector3f(0.48316, 0.42857, 1.644444));
	seedingCurve.push_back(Vector3f(0.46891, 0.405556, 1.694444));
	seedingCurve.push_back(Vector3f(0.45691, 0.388489, 1.744444));
	seedingCurve.push_back(Vector3f(0.44247, 0.365658, 1.794444));
	seedingCurve.push_back(Vector3f(0.42122, 0.343246, 1.844444));
	seedingCurve.push_back(Vector3f(0.42122, 0.343246, 1.894444));
};


void fig6(vector<Vector3f> &seedingCurve) {
	seedingCurve.push_back(Vector3f(1.51969, 0.548387, 1.444444));
	seedingCurve.push_back(Vector3f(1.51432, 0.548387, 1.494444));
	seedingCurve.push_back(Vector3f(1.50522, 0.548387, 1.544444));
	seedingCurve.push_back(Vector3f(1.49611, 0.548387, 1.594444));
	seedingCurve.push_back(Vector3f(1.48316, 0.548387, 1.644444));
	seedingCurve.push_back(Vector3f(1.46891, 0.548387, 1.694444));
	seedingCurve.push_back(Vector3f(1.45691, 0.548387, 1.744444));
	seedingCurve.push_back(Vector3f(1.44247, 0.548387, 1.794444));
	seedingCurve.push_back(Vector3f(1.42122, 0.548387, 1.844444));
	seedingCurve.push_back(Vector3f(1.40604, 0.548387, 1.894444));
};

void fig7(vector<Vector3f>& seedingCurve) {
	seedingCurve.push_back(Vector3f(2.51969, 0.548387, 1.444444));
	seedingCurve.push_back(Vector3f(2.51432, 0.548387, 1.494444));
	seedingCurve.push_back(Vector3f(2.50522, 0.548387, 1.544444));
	seedingCurve.push_back(Vector3f(2.49611, 0.548387, 1.594444));
	seedingCurve.push_back(Vector3f(2.48316, 0.548387, 1.644444));
	seedingCurve.push_back(Vector3f(2.46891, 0.548387, 1.694444));
	seedingCurve.push_back(Vector3f(2.45691, 0.548387, 1.744444));
	seedingCurve.push_back(Vector3f(2.44247, 0.548387, 1.794444));
	seedingCurve.push_back(Vector3f(2.42122, 0.548387, 1.844444));
	seedingCurve.push_back(Vector3f(2.40604, 0.548387, 1.894444));
};

void fig8(vector<Vector3f>& seedingCurve) {
	seedingCurve.push_back(Vector3f(3.51969, 0.548387, 0.444444));
	seedingCurve.push_back(Vector3f(3.51432, 0.504238, 1.494444));
	seedingCurve.push_back(Vector3f(3.50522, 0.478441, 1.544444));
	seedingCurve.push_back(Vector3f(3.49611, 0.452904, 1.594444));
	seedingCurve.push_back(Vector3f(3.48316, 0.42857, 1.644444));
	seedingCurve.push_back(Vector3f(3.46891, 0.405556, 1.694444));
	seedingCurve.push_back(Vector3f(3.45691, 0.388489, 1.744444));
	seedingCurve.push_back(Vector3f(3.44247, 0.365658, 1.794444));
	seedingCurve.push_back(Vector3f(3.42122, 0.343246, 1.844444));
	seedingCurve.push_back(Vector3f(3.40604, 0.327284, 1.894444));
};


/* Generate multiple streamlines from a seeding curve*/
void testMultipleStreamline() {

	std::cout << "Start testMultipleStreamline()..." << std::endl;
	// Step 1: Load a vector field. In this example, we use the Bernard dataset
	Steady3D_Bernard bernardVectorField("../../sample_data/bernard.raw");

	// Step 2: Initialize the streamline tracer
	StreamlineTracer streamlineTracer(bernardVectorField);
	streamlineTracer.setForward(true);
	streamlineTracer.setStepSize(0.1);
	int streamline_length = 50;

	// Step 3: Create a seeding curve with 4 points
	vector<Vector3f> seedingCurve;

	fig3(seedingCurve);
	//fig7(seedingCurve);
	//fig2(seedingCurve);
	//fig6(seedingCurve);
	//fig4(seedingCurve);
	//fig1(seedingCurve);
	//fig5(seedingCurve);
	//fig8(seedingCurve);
	printSeedpoints(seedingCurve, "../../output/seedpoints.vtk");
	
	// Step 4: Generate streamlines 
	vector<vector<Vector3f>> allStreamlines;
	vector<Vector3f> streamline;
	for (int i = 0; i < seedingCurve.size(); i++) {
		streamlineTracer.setForward(true);
		streamline = streamlineTracer.trace(seedingCurve[i], streamline_length);
		allStreamlines.push_back(streamline);


		streamlineTracer.setForward(false);
		streamline = streamlineTracer.trace(seedingCurve[i], streamline_length);
		//reverse
		reverse(streamline.begin(), streamline.end());
		allStreamlines.push_back(streamline);
		//strongly recommend reversing the backward streamline
	}
	vector<vector<Vector4f>> allAngles;
	for (int t = 0; t < allStreamlines.size(); t++) {
		Vector4f p0 = Vector4f(allStreamlines[t][0](0), allStreamlines[t][0](1), allStreamlines[t][0](2), 0);
		vector<Vector4f> angles;
		angles.push_back(p0);
		for (int i = 1; i < allStreamlines[t].size() - 1; i++) {
			float x1 = allStreamlines[t][i - 1](0);
			float y1 = allStreamlines[t][i - 1](1);
			float z1 = allStreamlines[t][i - 1](2);
			float x2 = allStreamlines[t][i](0);
			float y2 = allStreamlines[t][i](1);
			float z2 = allStreamlines[t][i](2);
			float x3 = allStreamlines[t][i + 1](0);
			float y3 = allStreamlines[t][i + 1](1);
			float z3 = allStreamlines[t][i + 1](2);
			vect vt1; vect vt2;
			vt1.i = x2 - x1; vt1.j = y2 - y1; vt1.k = z2 - z1;

			vt2.i = x3 - x2; vt2.j = y3 - y2; vt2.k = z3 - z2;
			float tmp1 = vt1.i * vt2.i + vt1.j * vt2.j + vt1.k * vt2.k;
			float tmp2 = sqrt(vt1.i * vt1.i + vt1.j * vt1.j + vt1.k * vt1.k) * sqrt(vt2.i * vt2.i + vt2.j * vt2.j + vt2.k * vt2.k);
			float theta = acos(tmp1 / tmp2) * 180 / M_PI;
			p0 = Vector4f(x2, y2, z2, theta);
			angles.push_back(p0);
		}
		int last = allStreamlines[t].size();
		p0 = Vector4f(allStreamlines[t][last-1](0), allStreamlines[t][last-1](1), allStreamlines[t][last-1](2), 0);
		angles.push_back(p0);
		allAngles.push_back(angles);
	}
	/*vector<vector<Vector3f>> allConnections;
	for (int i = 0; i < seedingCurve.size() - 1; i++) {
		vector<Vector3f> connection;
		int start = i; int next = i+1;
		if (allStreamlines[i].size() < allStreamlines[i + 1].size()) {
			start = i + 1; next = i;
		}
		for (int j = 0; j < min(allStreamlines[i].size(),allStreamlines[i+1].size()); j++) {
			connection.push_back(allStreamlines[start][j]);
			connection.push_back(allStreamlines[next][j]);
		}
		connection.push_back(allStreamlines[start][allStreamlines[start].size()-1]);
		allConnections.push_back(connection); 
	}*/
	// Step 4: Save the streamlines to a vtk file
	

	Util::saveLinesToVTKFile(allAngles, "../../output/multipleStreamline.vtk");
	/*Util::saveLinesToVTKFile(allConnections, "../../output/connect.vtk");*/
	std::cout << "Finish testMultipleStreamline()!" << std::endl;

}

/* Generate a surface given a seeding curve*/
void testStreamSurfaceWithSeedingCurve() {
	std::cout << "Start testStreamSurfaceWithSeedingCurve()..." << std::endl;
	// Step 1: Load a vector field. In this example, we use the Bernard dataset
	Steady3D_Bernard bernardVectorField("../../sample_data/bernard.raw");

	// Step 2: Initialize the streamline tracer
	StreamSurfaceTracer streamSurfaceTracer(bernardVectorField);
	streamSurfaceTracer.setForward(true);
	streamSurfaceTracer.setStepSize(0.1);
	int streamline_length = 50;

	// Step 3: Create a seeding curve with 4 points
	vector<Vector3f> seedingCurve;
	//fig3(seedingCurve);
	//fig7(seedingCurve);
	//fig2(seedingCurve);
	//fig6(seedingCurve);
	//fig4(seedingCurve);
	//fig1(seedingCurve);
	//fig5(seedingCurve);
	//fig8(seedingCurve);

	// Step 4: Generate a stream surface 
	vector<Vector3f> vertices;
	vector<Vector3i> faces;
	streamSurfaceTracer.trace(seedingCurve, streamline_length, vertices, faces);


	// Step 4: Save the streamlines to a vtk file

	Util::outputSurfaceToOBJ(vertices, faces, "../../output/streamsurface.obj");

	std::cout << "Finish testStreamSurfaceWithSeedingCurve()!" << std::endl;
}

void readData() {
	cout << "Number of seed points: ";
	cin >> nSeedPoints;
};

void init() {
	float x = 0;
	float y = 0;
	float z = 0;
	float tmpx = 4.0 / (std::cbrt(nSeedPoints) + 1);
	float tmpy = 1.0 / (std::cbrt(nSeedPoints) + 1);
	float tmpz = 2.0 / (std::cbrt(nSeedPoints) + 1);
	point tmp;
	for (int i = 0; i < std::cbrt(nSeedPoints); i++) {
		x += tmpx;
		for (int j = 0; j < std::cbrt(nSeedPoints); j++) {
			y += tmpy;
			for (int k = 0; k < std::cbrt(nSeedPoints); k++) {
				z += tmpz;
				tmp.x = x;
				tmp.y = y;
				tmp.z = z;
				seedPoints.push_back(tmp);
			}
			z = 0;
		}
		y = 0;
	}
};

void init2() {
	float x = 0;
	float y = 0;
	float z = 0;
	float tmpx = 4.0 / (128 + 1);
	float tmpy = 1.0 / (32 + 1);
	float tmpz = 2.0 / (64 + 1);
	point tmp;
	int count = 0;
	for (int i = 0; i < 128; i++) {
		x += tmpx;
		for (int j = 0; j < 32; j++) {
			y += tmpy;
			for (int k = 0; k < 64; k++) {
				z += tmpz;
				tmp.x = x;
				tmp.y = y;
				tmp.z = z;
				seedPoints.push_back(tmp);
				count++;
			}
			z = 0;
		}
		y = 0;
	}
	nSeedPoints = seedPoints.size();
}

void generateStreamlines() {
	init();
	Steady3D_Bernard bernardVectorField("../../sample_data/bernard.raw");

	StreamlineTracer streamlineTracer(bernardVectorField);
	streamlineTracer.setForward(true);
	streamlineTracer.setStepSize(0.1);
	int streamline_length = 50;

	vector<Vector3f> seedingCurve;

	
	for (int i = 0; i < nSeedPoints; i++) {
		seedingCurve.push_back(Vector3f(seedPoints[i].x, seedPoints[i].y, seedPoints[i].z));
	}
	printSeedpoints(seedingCurve, "../../output/seedpoints.vtk");
	// Step 4: Generate streamlines 
	vector<vector<Vector3f>> allStreamlines;
	vector<Vector3f> streamline;
	for (int i = 0; i < seedingCurve.size(); i++) {
		streamlineTracer.setForward(true);
		streamline = streamlineTracer.trace(seedingCurve[i], streamline_length);
		if (streamline.size() > 10) allStreamlines.push_back(streamline);

		streamlineTracer.setForward(false);
		streamline = streamlineTracer.trace(seedingCurve[i], streamline_length);
		//reverse
		reverse(streamline.begin(), streamline.end());
		if (streamline.size() > 10) allStreamlines.push_back(streamline);
		
	}


	vector<vector<Vector4f>> allAngles;
	for (int t = 0; t < allStreamlines.size(); t++) {
		Vector4f p0 = Vector4f(allStreamlines[t][0](0), allStreamlines[t][0](1), allStreamlines[t][0](2), 0);
		vector<Vector4f> angles;
		angles.push_back(p0);
		for (int i = 1; i < allStreamlines[t].size() - 1; i++) {
			float x1 = allStreamlines[t][i - 1](0);
			float y1 = allStreamlines[t][i - 1](1);
			float z1 = allStreamlines[t][i - 1](2);
			float x2 = allStreamlines[t][i](0);
			float y2 = allStreamlines[t][i](1);
			float z2 = allStreamlines[t][i](2);
			float x3 = allStreamlines[t][i + 1](0);
			float y3 = allStreamlines[t][i + 1](1);
			float z3 = allStreamlines[t][i + 1](2);
			vect vt1; vect vt2;
			vt1.i = x2 - x1; vt1.j = y2 - y1; vt1.k = z2 - z1;

			vt2.i = x3 - x2; vt2.j = y3 - y2; vt2.k = z3 - z2;
			float tmp1 = vt1.i * vt2.i + vt1.j * vt2.j + vt1.k * vt2.k;
			float tmp2 = sqrt(vt1.i * vt1.i + vt1.j * vt1.j + vt1.k * vt1.k) * sqrt(vt2.i * vt2.i + vt2.j * vt2.j + vt2.k * vt2.k);
			float theta = acos(tmp1 / tmp2) * 180 / M_PI;
			p0 = Vector4f(x2, y2, z2, theta);
			angles.push_back(p0);
		}
		int last = allStreamlines[t].size();
		p0 = Vector4f(allStreamlines[t][last - 1](0), allStreamlines[t][last - 1](1), allStreamlines[t][last - 1](2), 0);
		angles.push_back(p0);
		allAngles.push_back(angles);
	}
	Util::saveLinesToVTKFile(allAngles, "../../output/generate.vtk");
};

void BernardStreamlines() {
	vector<Vector3f> seedingCurve;
	extractBernard(seedingCurve);
	printSeedpoints(seedingCurve, "../../output/seedpoints.vtk");
	Steady3D_Bernard bernardVectorField("../../sample_data/bernard.raw");

	StreamlineTracer streamlineTracer(bernardVectorField);
	streamlineTracer.setForward(true);
	streamlineTracer.setStepSize(0.1);
	int streamline_length = 50;

	// Step 4: Generate streamlines 
	vector<vector<Vector3f>> allStreamlines;
	for (int i = 0; i < seedingCurve.size(); i++) {
		vector<Vector3f> streamline = streamlineTracer.trace(seedingCurve[i], streamline_length);
		if (streamline.size() > 40) allStreamlines.push_back(streamline);
	}

	vector<vector<Vector4f>> allAngles;
	for (int t = 0; t < allStreamlines.size(); t++) {
		Vector4f p0 = Vector4f(allStreamlines[t][0](0), allStreamlines[t][0](1), allStreamlines[t][0](2), 0);
		vector<Vector4f> angles;
		angles.push_back(p0);
		for (int i = 1; i < allStreamlines[t].size() - 1; i++) {
			float x1 = allStreamlines[t][i - 1](0);
			float y1 = allStreamlines[t][i - 1](1);
			float z1 = allStreamlines[t][i - 1](2);
			float x2 = allStreamlines[t][i](0);
			float y2 = allStreamlines[t][i](1);
			float z2 = allStreamlines[t][i](2);
			float x3 = allStreamlines[t][i + 1](0);
			float y3 = allStreamlines[t][i + 1](1);
			float z3 = allStreamlines[t][i + 1](2);
			vect vt1; vect vt2;
			vt1.i = x2 - x1; vt1.j = y2 - y1; vt1.k = z2 - z1;

			vt2.i = x3 - x2; vt2.j = y3 - y2; vt2.k = z3 - z2;
			float tmp1 = vt1.i * vt2.i + vt1.j * vt2.j + vt1.k * vt2.k;
			float tmp2 = sqrt(vt1.i * vt1.i + vt1.j * vt1.j + vt1.k * vt1.k) * sqrt(vt2.i * vt2.i + vt2.j * vt2.j + vt2.k * vt2.k);
			float theta = acos(tmp1 / tmp2) * 180 / M_PI;
			p0 = Vector4f(x2, y2, z2, theta);
			angles.push_back(p0);
		}
		int last = allStreamlines[t].size();
		p0 = Vector4f(allStreamlines[t][last - 1](0), allStreamlines[t][last - 1](1), allStreamlines[t][last - 1](2), 0);
		angles.push_back(p0);
		allAngles.push_back(angles);                                                          
	}
	Util::saveLinesToVTKFile(allAngles, "../../output/Bernard.vtk");

}


void grid_points() {
	init2();
	Steady3D_Bernard bernardVectorField("../../sample_data/bernard.raw");

	StreamlineTracer streamlineTracer(bernardVectorField);
	streamlineTracer.setForward(true);
	streamlineTracer.setStepSize(0.1);
	int streamline_length = 50;

	vector<Vector3f> seedingCurve;

	for (int i = 0; i < nSeedPoints; i++) {
		seedingCurve.push_back(Vector3f(seedPoints[i].x, seedPoints[i].y, seedPoints[i].z));
	}
	//printSeedpoints(seedingCurve, "../../output/seedpoints.vtk");
	// Step 4: Generate streamlines 
	vector<vector<Vector3f>> allStreamlines;
	vector<Vector3f> streamline;
//#pragma omp parallel
//#pragma omp parallel for num_threads(20)
	//#pragma omp parallel for
	for (int i = 0; i < seedingCurve.size(); i++) {
		streamlineTracer.setForward(true);
		//#pragma omp critical
		streamline = streamlineTracer.trace(seedingCurve[i], streamline_length);
		//#pragma omp critical
		allStreamlines.push_back(streamline);


	}

	struct grid_point {
		float x;
		float y;
		float z;
		float total_angle;
	};
	vector<grid_point> allGridPoints;
	vector<vector<Vector4f>> allAngles;
	for (int t = 0; t < allStreamlines.size(); t++) {
		float angle_sum = 0;
		for (int i = 1; i < allStreamlines[t].size() - 1; i++) {
			float x1 = allStreamlines[t][i - 1.0](0);
			float y1 = allStreamlines[t][i - 1.0](1);
			float z1 = allStreamlines[t][i - 1.0](2);
			float x2 = allStreamlines[t][i](0);
			float y2 = allStreamlines[t][i](1);
			float z2 = allStreamlines[t][i](2);
			float x3 = allStreamlines[t][i + 1.0](0);
			float y3 = allStreamlines[t][i + 1.0](1);
			float z3 = allStreamlines[t][i + 1.0](2);
			vect vt1; vect vt2;
			vt1.i = x2 - x1; vt1.j = y2 - y1; vt1.k = z2 - z1;
			vt2.i = x3 - x2; vt2.j = y3 - y2; vt2.k = z3 - z2;
			float tmp1 = vt1.i * vt2.i + vt1.j * vt2.j + vt1.k * vt2.k;
			float tmp2 = sqrt(vt1.i * vt1.i + vt1.j * vt1.j + vt1.k * vt1.k) * sqrt(vt2.i * vt2.i + vt2.j * vt2.j + vt2.k * vt2.k);
			float theta = 0.0;
			if (tmp1 <= tmp2){
				theta = acos(tmp1 / tmp2);
				theta = theta * 180.0 / M_PI;
			}
			angle_sum += theta;
		}
		grid_point tmp;
		tmp.x = allStreamlines[t][0](0);
		tmp.y = allStreamlines[t][0](1);
		tmp.z = allStreamlines[t][0](2);
		tmp.total_angle = angle_sum;
		/*if (angle_sum >= 36.4 && angle_sum <= 36.42) {
			cout << '1 ';
		}*/
		allGridPoints.push_back(tmp);
	}
	//Util::saveLinesToVTKFile(allAngles, "../../output/generate.vtk");

	ofstream pFile;
	string file = "../../output/gridpoints.vtk";
	pFile.open(file);
	
	pFile << "# vtk DataFile Version 3.0" << endl;
	pFile << "bernard3D dataset" << endl;
	pFile << "ASCII" << endl;
	pFile << "DATASET STRUCTURED_POINTS" << endl;
	pFile << "DIMENSIONS 128 32 64" << endl;
	pFile << "ASPECT_RATIO 0.03125 0.03125 0.03125" << endl;
	pFile << "ORIGIN 0.000000 0.000000 0.000000" << endl;
	pFile << "POINT_DATA " << allGridPoints.size() << endl;
	pFile << "SCALARS volume_scalars float 1" << endl;
	pFile << "LOOKUP_TABLE default" << endl;
	
	/*pFile << "# vtk DataFile Version 3.0" << endl;
	pFile << "	VortexCoreLine";
	pFile << "ASCII" << endl;
	pFile << "DATASET POLYDATA" << endl;
	pFile << "POINTS " << allGridPoints.size() << " float" << endl;*/
	/*for (int i = 0; i < allGridPoints.size(); i++)
		pFile << allGridPoints[i].total_angle << endl;*/

	int c = 0;
	int count = 0;
	for (int i = 0; i < 128; i++) {
		for (int j = 0; j < 32; j++) {
			for (int k = 0; k < 64; k++) {
				/*(if (count % 100 == 0) {
					pFile << allGridPoints[c].total_angle << " ";
					c++;
				}
				else pFile << "0 ";
				count++;*/
				/*if (allGridPoints[c].total_angle >= 36.4 && allGridPoints[c].total_angle <= 36.42) {
					cout << '1 ';
				}*/
				
				pFile << allGridPoints[c].total_angle << endl;
				c++;
			}
		}
		pFile << endl;
	}
};
int main(int argc, char **argv)
{
	//BernardStreamlines();
	//testSingleStreamline();
	//readData();
	//generateStreamlines();
	//testMultipleStreamline();
	//testStreamSurfaceWithSeedingCurve();
	grid_points();
	system("PAUSE");
	return 0;
}