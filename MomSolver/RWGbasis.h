#pragma once
#include "Mesh.h"
#include "constant.h"
#include <Eigen\Dense>
using namespace Eigen;
class RWGbasis
{
private:
	Mesh *mesh;
public:
	Vector3d r[4];
	double l;
	double S_plus, S_minus;
	Vector3d rho_c_plus, rho_c_minus;
	Vector3d r_center_plus, r_center_minus;
public:
	RWGbasis(Mesh *mesh,int n);
	~RWGbasis();
};

