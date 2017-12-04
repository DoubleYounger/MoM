#pragma once
#include "Source.h"
#include "constant.h"
#include "PlaneWave.h"
#include "Mesh.h"
#include "PO.h"
#include "MOM.h"
#include <fstream>
#include <Eigen\Dense>
using namespace Eigen;
using namespace std;


class Post :
	public Source
{
private:
	Mesh * mesh;
	PO posolver;
	MOM momsolver;
	PlaneWave planewave;
	int POedgeN;
	int MOMedgeN;
	vector<Triangle> Triangles;
	vector<Edge*> Edges;
	vector<int> MOMEdges;
	vector<int> POEdges;
	vector<Vector3d> Vertexes;
	double k, eta, mu, epsilon, omega;
	VectorXcd POCurrent;
	VectorXcd MOMCurrent;
public:
	complex<double> green(double R)
	{
		return exp(-im*k*R) / R;
	}
	Post(Mesh *mesh, PO &posolver, MOM &momsolver, PlaneWave &planewave);
	void RCS(double theta);
	Vector3cd ScatteredField(double phi, double theta);
	Vector3cd IntegralEfieldPO(int m, double phi, double theta);
	Vector3cd IntegralEfieldMOM(int m, double phi, double theta);
	~Post();
};

