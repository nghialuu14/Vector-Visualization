#include "utilities.h"

namespace Util
{

	// Save the 3D points to a VTK file
	void savePointToVTKFile(vector<Vector3f> results, string corelineFileName) {
		FILE * pFile;

		int totalPoints = results.size();

		pFile = fopen(corelineFileName.c_str(), "w");
		fprintf(pFile, "# vtk DataFile Version 3.0\n");
		fprintf(pFile, "VortexCoreLine\n");
		fprintf(pFile, "ASCII\n");
		fprintf(pFile, "DATASET POLYDATA\n");
		fprintf(pFile, "POINTS %d float\n", totalPoints);
		for (int i = 0; i < results.size(); i++) {

			fprintf(pFile, "%f %f %f\n", results[i](0), results[i](1), results[i](2));

		}

		// Save points
		//int numb_vertices = totalPoints + 1;
		//fprintf(pFile, "VERTICES 1 %d\n %d", numb_vertices, totalPoints);
		//for (int i = 0; i < results.size(); i++) {
		//	fprintf(pFile, " %d", i);
		//}

		int numb_vertices = totalPoints + 1;
		fprintf(pFile, "LINES 1 %d\n %d", numb_vertices, totalPoints);
		for (int i = 0; i < results.size(); i++) {
			fprintf(pFile, " %d", i);
		}
		//
		/*int lines = results.size();
		int points = lines + totalPoints;
		fprintf(pFile, "LINES %d %d \n", lines, points);

		int pointIndex = 0;
		for (int i = 0; i < pathlineList.size(); i++) {
		fprintf(pFile, "%d ", pathlineList[i].size());
		for (int j = 0; j < pathlineList[i].size(); j++) {
		fprintf(pFile, "%d ", pointIndex);
		pointIndex++;
		}

		fprintf(pFile, "\n");
		}*/



		fclose(pFile);
	}


	// Save a line containing multiple 3D points to a VTK file
	void saveLinesToVTKFile(vector<vector<Vector4f>> results, string corelineFileName) {
		FILE * pFile;

		int totalPoints = 0;
		int line_size = 0;
		int line_numb = 0;
		for (int i = 0; i < results.size(); i++) {


			if (results[i].size() > 1) {
				for (int j = 0; j < results[i].size(); j++)
					totalPoints++;

				line_numb += 1;
			}
		}

		pFile = fopen(corelineFileName.c_str(), "w");
		fprintf(pFile, "# vtk DataFile Version 3.0\n");
		fprintf(pFile, "VortexCoreLine\n");
		fprintf(pFile, "ASCII\n");
		fprintf(pFile, "DATASET POLYDATA\n");
		fprintf(pFile, "POINTS %d float\n", totalPoints);
		for (int i = 0; i < results.size(); i++) {
			if (results[i].size() > 1) {
				for (int j = 0; j < results[i].size(); j++)
					fprintf(pFile, "%f %f %f\n", results[i][j](0), results[i][j](1), results[i][j](2));
			}

		}



		int numb_vertices = totalPoints + 1;
		line_size = line_numb + totalPoints;

		fprintf(pFile, "LINES %d %d\n", line_numb, line_size);
		int index = 0;
		for (int i = 0; i < results.size(); i++) {
			if (results[i].size() > 1) {
				fprintf(pFile, "%d ", results[i].size());
				for (int j = 0; j < results[i].size(); j++) {
					fprintf(pFile, "%d ", index);
					index++;
				}
				fprintf(pFile, "\n");
			}
		}
		fprintf(pFile, "POINT_DATA %d\n", totalPoints);
		fprintf(pFile, "SCALARS sample_scalars float 1\n");
		fprintf(pFile, "LOOKUP_TABLE lookup_table\n");
		for (int i = 0; i < results.size(); i++) {
			if (results[i].size() > 1) {
				for (int j = 0; j < results[i].size(); j++)
					fprintf(pFile, "%f\n", results[i][j](3));
			}
		}
		//
		/*int lines = results.size();
		int points = lines + totalPoints;
		fprintf(pFile, "LINES %d %d \n", lines, points);

		int pointIndex = 0;
		for (int i = 0; i < pathlineList.size(); i++) {
		fprintf(pFile, "%d ", pathlineList[i].size());
		for (int j = 0; j < pathlineList[i].size(); j++) {
		fprintf(pFile, "%d ", pointIndex);
		pointIndex++;
		}

		fprintf(pFile, "\n");
		}*/



		fclose(pFile);
	}
	
	void saveLinesToVTKFile(vector<vector<Vector3f>> results, string corelineFileName) {
		FILE* pFile;

		int totalPoints = 0;
		int line_size = 0;
		int line_numb = 0;
		for (int i = 0; i < results.size(); i++) {


			if (results[i].size() > 1) {
				for (int j = 0; j < results[i].size(); j++)
					totalPoints++;

				line_numb += 1;
			}
		}

		pFile = fopen(corelineFileName.c_str(), "w");
		fprintf(pFile, "# vtk DataFile Version 3.0\n");
		fprintf(pFile, "VortexCoreLine\n");
		fprintf(pFile, "ASCII\n");
		fprintf(pFile, "DATASET POLYDATA\n");
		fprintf(pFile, "POINTS %d float\n", totalPoints);
		for (int i = 0; i < results.size(); i++) {
			if (results[i].size() > 1) {
				for (int j = 0; j < results[i].size(); j++)
					fprintf(pFile, "%f %f %f\n", results[i][j](0), results[i][j](1), results[i][j](2));
			}

		}



		int numb_vertices = totalPoints + 1;
		line_size = line_numb + totalPoints;

		fprintf(pFile, "LINES %d %d\n", line_numb, line_size);
		int index = 0;
		for (int i = 0; i < results.size(); i++) {
			if (results[i].size() > 1) {
				fprintf(pFile, "%d ", results[i].size());
				for (int j = 0; j < results[i].size(); j++) {
					fprintf(pFile, "%d ", index);
					index++;
				}
				fprintf(pFile, "\n");
			}
		}
		fclose(pFile);
	}

	void generateObjSurfaceFromLines(vector<vector<Vector3f>> lines, string outputFileName)
	{
		FILE * pFile;
		pFile = fopen(outputFileName.c_str(), "w");
		fprintf(pFile, "g object\n");
		vector<int> streamline_start_pos;

		// Output vertices
		int start_pos = 1;
		for (int i = 0; i < lines.size(); i++) {
			for (int j = 0; j < lines[i].size(); j++) {
				fprintf(pFile, "v %f %f %f 1.0 1.0 0.0\n", lines[i][j](0), lines[i][j](1), lines[i][j](2));
			}
			streamline_start_pos.push_back(start_pos);
			start_pos += lines[i].size();

		}
		// Output faces
		for (int streamlineID = 0; streamlineID < lines.size() - 1; streamlineID++) {
			int minLength = min(lines[streamlineID].size(), lines[streamlineID + 1].size());
			for (int j = 0; j < minLength - 1; j++) {
				/*
					streamline1            streamline2
						idx0 -----------------idx2
						 |						|
						idx1 -----------------idx3
				*/
				int idx0 = streamline_start_pos[streamlineID] + j;
				int idx1 = streamline_start_pos[streamlineID] + j + 1;
				int idx2 = streamline_start_pos[streamlineID + 1] + j;
				int idx3 = streamline_start_pos[streamlineID + 1] + j + 1;
				fprintf(pFile, "f %d//%d %d//%d %d//%d \n", idx0, idx0, idx2, idx2, idx1, idx1);
				fprintf(pFile, "f %d//%d %d//%d %d//%d \n", idx2, idx2, idx3, idx3, idx1, idx1);
			}
		}
		fclose(pFile);
	}


	/* Save a surface to an OBJ file. The surface is represented in the triangle mesh structure which has vertices and triangle connections*/
	void outputSurfaceToOBJ(vector<Vector3f> &vertices, vector<Vector3i> &faces, string fileName) {
		FILE * pFile;
		pFile = fopen(fileName.c_str(), "w");
		fprintf(pFile, "g object\n");

		// Output vertices
		for (int i = 0; i < vertices.size(); i++) {
			fprintf(pFile, "v %f %f %f 1.0 1.0 0.0\n", vertices[i](0), vertices[i](1), vertices[i](2));

		}

		// output faces
		for (int i = 0; i < faces.size(); i++) {
			fprintf(pFile, "f %d//%d %d//%d %d//%d \n", faces[i](0), faces[i](0), faces[i](1), faces[i](1), faces[i](2), faces[i](2));

		}

		fclose(pFile);
	}

}