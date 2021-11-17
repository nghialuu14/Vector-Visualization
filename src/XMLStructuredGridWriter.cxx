#define _USE_MATH_DEFINES
#include <iostream>
#include <omp.h>
#include "Utilities.h"
#include "Steady3D.h"
#include "StreamlineTracer.h"
#include "StreamSurfaceTracer.h"
#include <cmath>
#include <fstream>
#include <string>
#include <vtkCellArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkFloatArray.h>
#include <vtkStructuredGridReader.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>



#include "Steady3D_Bernard.h"

using namespace std;

struct point {
	double x;
	double y;
	double z;
};

struct vect {
	double i;
	double j;
	double k;
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


void extractBernard(vector<Vector3f>& seedPoints) {
	ifstream bernardFile("../../sample_data/bernard3d.vtk");
	string tmp;
	for (int i = 0; i < 9; i++) getline(bernardFile, tmp);
	for (int i = 0; i < 262144; i++) {
		double x, y, z;
		bernardFile >> x >> y >> z;
		seedPoints.push_back(Vector3f(x, y, z));
	}
	bernardFile.close();
};

void printSeedpoints(vector<Vector3f> seedingCurve, string corelineFileName) {
	ofstream pFile;
	pFile.open(corelineFileName.c_str());
	pFile << "# vtk DataFile Version 3.0" << endl;
	pFile << "VortexCoreLine" << endl;
	pFile << "ASCII" << endl;
	pFile << "DATASET POLYDATA" << endl;
	pFile << "POINTS " << seedingCurve.size() << " float" << endl;
	for (int i = 0; i < seedingCurve.size(); i++)
		pFile << seedingCurve[i](0) << " " << seedingCurve[i](1) << " " << seedingCurve[i](2) << endl;
	pFile.close();
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
		double x1 = streamline[i - 1](0);
		double y1 = streamline[i - 1](1);
		double z1 = streamline[i - 1](2);
		double x2 = streamline[i](0);
		double y2 = streamline[i](1);
		double z2 = streamline[i](2);
		double x3 = streamline[i + 1](0);
		double y3 = streamline[i + 1](1);
		double z3 = streamline[i + 1](2);
		vect vt1; vect vt2;
		vt1.i = x2 - x1; vt1.j = y2 - y1; vt1.k = z2 - z1;

		vt2.i = x3 - x2; vt2.j = y3 - y2; vt2.k = z3 - z2;
		double tmp1 = vt1.i * vt2.i + vt1.j * vt2.j + vt1.k * vt2.k;
		double tmp2 = sqrt(vt1.i * vt1.i + vt1.j * vt1.j + vt1.k * vt1.k) * sqrt(vt2.i * vt2.i + vt2.j * vt2.j + vt2.k * vt2.k);
		double theta = acos(tmp1 / tmp2) * 180 / M_PI;
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


void connect(vector<Vector3f>& seedingCurve, double x1, double y1, double z1, double x2, double y2, double z2, int n) {
	double tmpx = abs(x1 - x2) / n;
	double tmpy = abs(y1 - y2) / n;
	double tmpz = abs(z1 - z2) / n;
	double x = x1; double y = y1; double z = z1;
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

void fig4(vector<Vector3f>& seedingCurve) {
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


void fig6(vector<Vector3f>& seedingCurve) {
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
			double x1 = allStreamlines[t][i - 1](0);
			double y1 = allStreamlines[t][i - 1](1);
			double z1 = allStreamlines[t][i - 1](2);
			double x2 = allStreamlines[t][i](0);
			double y2 = allStreamlines[t][i](1);
			double z2 = allStreamlines[t][i](2);
			double x3 = allStreamlines[t][i + 1](0);
			double y3 = allStreamlines[t][i + 1](1);
			double z3 = allStreamlines[t][i + 1](2);
			vect vt1; vect vt2;
			vt1.i = x2 - x1; vt1.j = y2 - y1; vt1.k = z2 - z1;

			vt2.i = x3 - x2; vt2.j = y3 - y2; vt2.k = z3 - z2;
			double tmp1 = vt1.i * vt2.i + vt1.j * vt2.j + vt1.k * vt2.k;
			double tmp2 = sqrt(vt1.i * vt1.i + vt1.j * vt1.j + vt1.k * vt1.k) * sqrt(vt2.i * vt2.i + vt2.j * vt2.j + vt2.k * vt2.k);
			double theta = acos(tmp1 / tmp2) * 180 / M_PI;
			p0 = Vector4f(x2, y2, z2, theta);
			angles.push_back(p0);
		}
		int last = allStreamlines[t].size();
		p0 = Vector4f(allStreamlines[t][last - 1](0), allStreamlines[t][last - 1](1), allStreamlines[t][last - 1](2), 0);
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
	double x = 0;
	double y = 0;
	double z = 0;
	double tmpx = 4.0 / (std::cbrt(nSeedPoints) + 1);
	double tmpy = 1.0 / (std::cbrt(nSeedPoints) + 1);
	double tmpz = 2.0 / (std::cbrt(nSeedPoints) + 1);
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
	double x = 0;
	double y = 0;
	double z = 0;
	double tmpx = 4.0 / (128 + 1);
	double tmpy = 1.0 / (32 + 1);
	double tmpz = 2.0 / (64 + 1);
	point tmp;
	int count = 0;
	for (int k = 0; k < 64; k++) {
		for (int j = 0; j < 32; j++) {
			for (int i = 0; i < 128; i++) {
				tmp.x = x;
				tmp.y = y;
				tmp.z = z;
				seedPoints.push_back(tmp);
				x += 0.03125;
			}
			x = 0;
			y += 0.03125;
		}
		x = 0;
		y = 0;
		z += 0.03125;
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
			double x1 = allStreamlines[t][i - 1](0);
			double y1 = allStreamlines[t][i - 1](1);
			double z1 = allStreamlines[t][i - 1](2);
			double x2 = allStreamlines[t][i](0);
			double y2 = allStreamlines[t][i](1);
			double z2 = allStreamlines[t][i](2);
			double x3 = allStreamlines[t][i + 1](0);
			double y3 = allStreamlines[t][i + 1](1);
			double z3 = allStreamlines[t][i + 1](2);
			vect vt1; vect vt2;
			vt1.i = x2 - x1; vt1.j = y2 - y1; vt1.k = z2 - z1;

			vt2.i = x3 - x2; vt2.j = y3 - y2; vt2.k = z3 - z2;
			double tmp1 = vt1.i * vt2.i + vt1.j * vt2.j + vt1.k * vt2.k;
			double tmp2 = sqrt(vt1.i * vt1.i + vt1.j * vt1.j + vt1.k * vt1.k) * sqrt(vt2.i * vt2.i + vt2.j * vt2.j + vt2.k * vt2.k);
			double theta = acos(tmp1 / tmp2) * 180 / M_PI;
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
			double x1 = allStreamlines[t][i - 1](0);
			double y1 = allStreamlines[t][i - 1](1);
			double z1 = allStreamlines[t][i - 1](2);
			double x2 = allStreamlines[t][i](0);
			double y2 = allStreamlines[t][i](1);
			double z2 = allStreamlines[t][i](2);
			double x3 = allStreamlines[t][i + 1](0);
			double y3 = allStreamlines[t][i + 1](1);
			double z3 = allStreamlines[t][i + 1](2);
			vect vt1; vect vt2;
			vt1.i = x2 - x1; vt1.j = y2 - y1; vt1.k = z2 - z1;

			vt2.i = x3 - x2; vt2.j = y3 - y2; vt2.k = z3 - z2;
			double tmp1 = vt1.i * vt2.i + vt1.j * vt2.j + vt1.k * vt2.k;
			double tmp2 = sqrt(vt1.i * vt1.i + vt1.j * vt1.j + vt1.k * vt1.k) * sqrt(vt2.i * vt2.i + vt2.j * vt2.j + vt2.k * vt2.k);
			double theta = acos(tmp1 / tmp2) * 180 / M_PI;
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

struct grid_point {
	double x;
	double y;
	double z;
	double scalar;
};

void printVolume(string file_name, vector<grid_point> arr) {
	ofstream pFile;
	string file = file_name;
	pFile.open(file);

	pFile << "# vtk DataFile Version 3.0" << endl;
	pFile << "bernard3D dataset" << endl;
	pFile << "ASCII" << endl;
	pFile << "DATASET STRUCTURED_POINTS" << endl;
	pFile << "DIMENSIONS 128 32 64" << endl;
	pFile << "ASPECT_RATIO 0.03125 0.03125 0.03125" << endl;
	pFile << "ORIGIN 0.000000 0.000000 0.000000" << endl;
	pFile << "POINT_DATA " << arr.size() << endl;
	pFile << "SCALARS volume_scalars double 1" << endl;
	pFile << "LOOKUP_TABLE default" << endl;

	for (int i = 0; i < arr.size(); i++)
		pFile << arr[i].scalar << endl;
	pFile.close();
}

void printVolume2(char* file_name, vector<grid_point> arr, Steady3D &dataset) {
	// Create a grid
	vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	int dims[3];
	dataset.getDomainDimension(dims);
	unsigned int numi = dims[0];
	unsigned int numj = dims[1];
	unsigned int numk = dims[2];

	double x = 0;
	double y = 0;
	double z = 0;
	float xRange[2];
	float yRange[2];
	float zRange[2];
	dataset.getXRange(xRange);
	dataset.getYRange(yRange);
	dataset.getZRange(zRange);
	double tmpx = xRange[1] / numi;
	double tmpy = yRange[1] / numj;
	double tmpz = zRange[1] / numk;
	point tmp;
	for (int k = 0; k < numk; k++) {
		for (int j = 0; j < numj; j++) {
			for (int i = 0; i < numi; i++) {
				points->InsertNextPoint(x, y, z);
				x += tmpx;
			}
			x = 0;
			y += tmpy;
		}
		x = 0;
		y = 0;
		z += tmpz;
	}

	// Specify the dimensions of the grid
	structuredGrid->SetDimensions(numi, numj, numk);
	structuredGrid->SetPoints(points);

	vtkFloatArray* scalars = vtkFloatArray::New();
	for (int i = 0; i < numi * numj * numk; i++) scalars->InsertTuple1(i, arr[i].scalar);
	structuredGrid->GetPointData()->SetScalars(scalars);

	// Write file
	/*vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
	writer->SetFileName(file_name);
	writer->SetInputData(structuredGrid);
	writer->Write();*/
	vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();
	writer->SetFileName(file_name);
	writer->SetInputData(structuredGrid);
	writer->Write();
}

void grid_points(int type) {
	init2();
	Steady3D_Bernard bernardVectorField("../src/sample_data/bernard.raw");

	StreamlineTracer streamlineTracer(bernardVectorField);
	streamlineTracer.setForward(true);
	streamlineTracer.setStepSize(0.5);
	int streamline_length = 50;

	vector<Vector3f> seedingCurve;

	for (int i = 0; i < nSeedPoints; i++) {
		seedingCurve.push_back(Vector3f(seedPoints[i].x, seedPoints[i].y, seedPoints[i].z));
	}

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

	int x = 1;
	int y = 1;
	int z = 1;

	vector<vector<vector<grid_point>>> coordinates;
	for (int i = 0; i <= 129; i++) {
		vector<vector<grid_point>> u;
		for (int j = 0; j <= 33; j++) {
			vector<grid_point> v;
			for (int k = 0; k <= 65; k++) {
				grid_point tmp;
				tmp.x = 0;
				tmp.y = 0;
				tmp.z = 0;
				tmp.scalar = 0;
				v.push_back(tmp);
			}
			u.push_back(v);
		}
		coordinates.push_back(u);
	}
	vector<grid_point> allGridPoints;
	vector<vector<Vector4f>> allAngles;
	vector<Vector3f> seedingCurve2;
	vector<Vector3f> seedingCurve3;
	for (int t = 0; t < allStreamlines.size(); t++) {
		double angle_sum = 0;

		for (int i = 1; i < allStreamlines[t].size() - 1; i++) {
			double x1 = allStreamlines[t][i - 1.0](0);
			double y1 = allStreamlines[t][i - 1.0](1);
			double z1 = allStreamlines[t][i - 1.0](2);
			double x2 = allStreamlines[t][i](0);
			double y2 = allStreamlines[t][i](1);
			double z2 = allStreamlines[t][i](2);
			double x3 = allStreamlines[t][i + 1.0](0);
			double y3 = allStreamlines[t][i + 1.0](1);
			double z3 = allStreamlines[t][i + 1.0](2);
			vect vt1; vect vt2;
			vt1.i = x2 - x1; vt1.j = y2 - y1; vt1.k = z2 - z1;
			vt2.i = x3 - x2; vt2.j = y3 - y2; vt2.k = z3 - z2;
			double tmp1 = vt1.i * vt2.i + vt1.j * vt2.j + vt1.k * vt2.k;
			double tmp2 = sqrt(vt1.i * vt1.i + vt1.j * vt1.j + vt1.k * vt1.k) * sqrt(vt2.i * vt2.i + vt2.j * vt2.j + vt2.k * vt2.k);
			double theta = 0.0;
			if (tmp1 <= tmp2) {
				theta = acos(tmp1 / tmp2);
				theta = theta * 180.0 / M_PI;
			}
			angle_sum += theta;
		}
		grid_point tmp;
		tmp.x = allStreamlines[t][0](0);
		tmp.y = allStreamlines[t][0](1);
		tmp.z = allStreamlines[t][0](2);
		tmp.scalar = angle_sum;
		//override with vx only
		/*float vx, vy, vz;
		bernardVectorField.get_vector_at(tmp.x, tmp.y, tmp.z, vx, vy, vz);
		tmp.scalar = vx;*/
		coordinates[x][y][z] = tmp;
		allGridPoints.push_back(tmp);
		if (type == 1 && tmp.scalar >= 400
			&& tmp.x >= 0.1 && tmp.x <= 3.9
			&& tmp.y >= 0.1 && tmp.y <= 0.9
			&& tmp.z >= 0.1 && tmp.z <= 1.9) {
			seedingCurve2.push_back(Vector3f(tmp.x, tmp.y, tmp.z));
		}
		else if (type == 2 && 700 <= tmp.scalar && tmp.scalar <= 800 
			&& tmp.x >= 0.2 && tmp.x <= 3.8 
			&& tmp.y >= 0.2 && tmp.y <= 0.8 
			&& tmp.z >= 0.2 && tmp.z <= 1.8) {
			seedingCurve2.push_back(Vector3f(tmp.x, tmp.y, tmp.z));
		}
		x++;
		if (x > 128) { x = 1; y++; }
		if (y > 32) { x = 1; y = 1; z++; }
	}
	vector<vector<Vector3f>> allStreamlines2;
	for (int i = 0; i < seedingCurve2.size(); i++) {
		streamlineTracer.setForward(true);
		streamline = streamlineTracer.trace(seedingCurve2[i], streamline_length);
		allStreamlines2.push_back(streamline);
	}
	if (type == 1) {
		Util::saveLinesToVTKFile(allStreamlines2, "../src/output/streamlineLARGEANGLE.vtk");
		printSeedpoints(seedingCurve2, "../src/output/seedpointsLARGE.vtk");
	}
	else if (type == 2) {
		Util::saveLinesToVTKFile(allStreamlines2, "../src/output/streamlineSMALLANGLE.vtk");
		printSeedpoints(seedingCurve2, "../src/output/seedpointsSMALL.vtk");
	}
	int c = 0;

	vector<grid_point> allGradient;
	for (int k = 1; k <= 64; k++)
		for (int j = 1; j <= 32; j++)
			for (int i = 1; i <= 128; i++) {
				double vtx = abs(coordinates[i + 1][j][k].scalar - coordinates[i - 1][j][k].scalar) / abs(coordinates[i + 1][j][k].x - coordinates[i - 1][j][k].x);
				double vty = abs(coordinates[i][j + 1][k].scalar - coordinates[i][j - 1][k].scalar) / abs(coordinates[i][j + 1][k].y - coordinates[i][j - 1][k].y);
				double vtz = abs(coordinates[i][j][k + 1].scalar - coordinates[i][j][k - 1].scalar) / abs(coordinates[i][j][k + 1].z - coordinates[i][j][k - 1].z);
				double gradient = sqrt(vtx * vtx + vty * vty + vtz * vtz);
				grid_point tmp;
				tmp.x = coordinates[i][j][k].x;
				tmp.y = coordinates[i][j][k].y;
				tmp.z = coordinates[i][j][k].z;
				tmp.scalar = gradient;
				allGradient.push_back(tmp);
				if (type == 1 && tmp.scalar >= 10000
					&& tmp.x >= 0.1 && tmp.x <= 3.9
					&& tmp.y >= 0.1 && tmp.y <= 0.9
					&& tmp.z >= 0.1 && tmp.z <= 1.9) {
					seedingCurve3.push_back(Vector3f(tmp.x, tmp.y, tmp.z));
				}
				else if (type == 2 && 1500 <= tmp.scalar && tmp.scalar <= 2000
					&& tmp.x >= 0.2 && tmp.x <= 3.8
					&& tmp.y >= 0.2 && tmp.y <= 0.8
					&& tmp.z >= 0.2 && tmp.z <= 1.8) {
					seedingCurve3.push_back(Vector3f(tmp.x, tmp.y, tmp.z));
				}
			}
	vector<vector<Vector3f>> allStreamlines3;
	for (int i = 0; i < seedingCurve3.size(); i++) {
		streamlineTracer.setForward(true);
		streamline = streamlineTracer.trace(seedingCurve3[i], streamline_length);
		allStreamlines3.push_back(streamline);
	}
	if (type == 1) {
		Util::saveLinesToVTKFile(allStreamlines3, "../src/output/largeGradient.vtk");
		printSeedpoints(seedingCurve3, "../src/output/seedpointsLargeGradient.vtk");
	}
	else if (type == 2) {
		Util::saveLinesToVTKFile(allStreamlines3, "../src/output/smallGradient.vtk");
		printSeedpoints(seedingCurve3, "../src/output/seedpointsSmallGradient.vtk");
	}

	printVolume2("../src/output/gridpoints.vtk", allGridPoints, bernardVectorField);
	printVolume2("../src/output/gradients.vtk", allGradient, bernardVectorField);

};

void surfaceGenerate() {
	int type;
	cout << "Please choose twisted (1) or smooth (2) surface: ";
	cin >> type;
	grid_points(type);
}

int main(int argc, char** argv)
{
	//BernardStreamlines();
	//testSingleStreamline();
	//readData();
	//generateStreamlines();
	//testMultipleStreamline();
	//testStreamSurfaceWithSeedingCurve();
	surfaceGenerate();
	system("PAUSE");
	return 0;
}