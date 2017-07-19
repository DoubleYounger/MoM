#pragma once
#include "Source.h"
#include "constant.h"
#include <Eigen\Dense>
using namespace Eigen;
enum Polarization
{
	H, V
};
class PlaneWave :
	public Source
{
public:
	double magnitude, theta_i, phi_i;
	double fre;
	double omega, wavenumber, lambda;
	Polarization polar;
	Vector3d polarDir, kDir;
public:
	PlaneWave(double magnitude, double fre, double theta, double phi, Polarization p);
	~PlaneWave();
	Vector3cd Eincident(Vector3d r);
	Vector3cd Hincident(Vector3d r);
};

