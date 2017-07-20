// MomSolver.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "Mesh.h"
#include <iostream>
#include <Eigen\Dense>
#include "MOM.h"
#include "PO.h"
#include "PlaneWave.h"
#include <time.h>
using namespace std;
using namespace Eigen;
int main()
{
	string FileName;
	/*cout << "Please input the mesh file name:" << endl;
	cin >> FileName;*/
	clock_t start, end;
	double duration;
	FileName = "C:/Users/sszz/Documents/Visual Studio 2015/Projects/MomSolver/rect.nas";
	Mesh mesh;
	mesh.open(FileName);
	mesh.PrintMeshInfo();
	if (mesh.WriteMesh())
	{
		cout << "Mesh output successfully!" << endl;
	}
	//this is the procedure for MoM
	mesh.EdgeInit();
	mesh.WriteEdge();
	mesh.PrintEdgeInfo();
	PlaneWave Einc(1.0, 3e8, PI / 4.0, 0.0, H);
	//MOM momsolver(&mesh, Einc);
	//start = clock();
	//cout << "Filling Matrix..." << endl;
	//momsolver.FillMatrix();
	//momsolver.FillRhs();
	//end = clock();
	//duration = (double)(end - start) / CLOCKS_PER_SEC;
	//cout << "Filling Matrix use " << duration << " second." << endl;
	//start = clock();
	//cout << "Solving system equation..." << endl;
	//momsolver.Solver();
	//end = clock();
	//duration = (double)(end - start) / CLOCKS_PER_SEC;
	//cout << "Solving system equation use " << duration << " second." << endl;
	//momsolver.RCS(PI / 4.0);
	//this is the procedure for PO
	PO posolver(&mesh, Einc);
	posolver.judgeLitPatch(Einc.kDir);
	posolver.WriteLitPatch();
	posolver.Solver();
	posolver.RCS(PI / 4.0);
	cout << "done." << endl;
    return 0;
}

