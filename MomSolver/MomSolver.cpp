// MomSolver.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "Mesh.h"
#include <iostream>
#include <Eigen\Dense>
#include "MOM.h"
#include "PlaneWave.h"
using namespace std;
using namespace Eigen;
int main()
{
	string FileName;
	/*cout << "Please input the mesh file name:" << endl;
	cin >> FileName;*/
	FileName = "C:/Users/sszz/Documents/Visual Studio 2015/Projects/MomSolver/sphere.nas";
	Mesh mesh;
	mesh.open(FileName);
	mesh.PrintMeshInfo();
	if (mesh.WriteMesh())
	{
		cout << "Mesh output successfully!" << endl;
	}
	mesh.EdgeInit();
	mesh.WriteEdge();
	mesh.PrintEdgeInfo();
	PlaneWave Einc(1.0, 3e8, PI / 2.0, 0.0, H);
	MOM momsolver(&mesh, Einc);
	momsolver.FillMatrix();
	momsolver.FillRhs();
	momsolver.Solver();
	momsolver.RCS();
    return 0;
}

