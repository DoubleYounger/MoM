#pragma once
#include "Source.h"
#include "Mesh.h"
#include "PlaneWave.h"
#include "constant.h"
#include <fstream>
#include <Eigen\Dense>
using namespace Eigen;
using namespace std;
class PO :
	public Source
{
private:
	VectorXcd KaJ;
	int *delta;
	int *RWGdelta;
	Mesh *mesh;
	int facetN, edgeN, verN, litPatchN;
	vector<Triangle> Triangles;
	vector<Edge*> Edges;
	vector<Vector3d> Vertexes;
	PlaneWave planewave;
	double k, eta, mu, epsilon, omega;
public:
	PO(Mesh *mesh, PlaneWave &planewave);

	//give the direction of incident wave, find the lit patch
	void judgeLitPatch(Vector3d kDir);
	void WriteLitPatch();
	void LitPatchInfo();
	complex<double> green(double R)
	{
		return exp(-im*k*R) / R;
	}
	void Solver();
	Vector3cd IntegralEfield(int m, double phi, double theta);
	Vector3cd ScatteredField(double phi, double theta);
	void RCS(double theta);
	~PO();
};

