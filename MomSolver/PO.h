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
	VectorXcd deltaKaJ;
	int *delta;
	int *RWGdelta;
	int **RWGPOdelta;
	Mesh *mesh;
	int facetN, edgeN, verN, litPatchN;
	vector<Triangle> Triangles;
	vector<Edge*> Edges;
	vector<int> MOMEdges;
	vector<int> POEdges;
	vector<Vector3d> Vertexes;
	PlaneWave planewave;
	double k, eta, mu, epsilon, omega;
public:
	PO(Mesh *mesh, PlaneWave &planewave);
	void judgeLitEdge();
	//give the direction of incident wave, find the lit patch
	void judgeLitPatch(Vector3d kDir);
	void WriteLitPatch();
	void LitPatchInfo();
	complex<double> green(double R)
	{
		return exp(-im*k*R) / R;
	}
	complex<double> gradientG(double R);
	VectorXcd Solver();
	complex<double> nearHField(Vector3d r, VectorXcd MOMJ, int n, Vector3d, Vector3d, Vector3d , Vector3d);
	VectorXcd updateCurrent(VectorXcd MOMJ);
	VectorXcd getCurrent();
	Vector3cd IntegralEfield(int m, double phi, double theta);
	Vector3cd ScatteredField(double phi, double theta);
	void RCS(double theta);
	~PO();
};

