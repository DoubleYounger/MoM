#pragma once
#include "Source.h"
#include "Mesh.h"
#include "PlaneWave.h"
#include "constant.h"
#include <Eigen\Dense>
using namespace Eigen;
class MOM :
	public Source
{
private:
	MatrixXcd ZMatrix;
	VectorXcd rhs;
	VectorXcd X;
	Mesh *mesh;
	int edgeN;
	vector<Triangle> Triangles;
	vector<Edge*> Edges;
	vector<Vector3d> Vertexes;
	PlaneWave planewave;
	double k, eta, mu, epsilon, omega;
public:
	MOM(Mesh *mesh, PlaneWave &planewave);
	void FillMatrix();
	complex<double> Integral(int m, int n, int trim, int trin, int vm3, int vn3);
	complex<double> SingularIntegral(int m, int n, int trim, int trin, int vm3, int vn3);
	complex<double> green(double R);
	complex<double> Singulargreen(double R);
	void FillRhs();
	complex<double> IntegralRhs(int m);
	Vector3cd IntegralEfield(int m, double phi, double theta);
	void Solver();
	Vector3cd ScatteredField(double phi, double theta);
	void RCS();
	~MOM();
};

