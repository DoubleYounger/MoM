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
	double *delta;
	Mesh *mesh;
	int facetN, edgeN, verN;
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
	~PO();
};

