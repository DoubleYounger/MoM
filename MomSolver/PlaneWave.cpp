#include "stdafx.h"
#include "PlaneWave.h"


PlaneWave::PlaneWave(double magnitude, double fre, double theta, double phi, Polarization p = H)
{
	this->magnitude = magnitude;
	this->fre = fre;
	this->theta_i = theta;
	this->phi_i = phi;
	this->omega = 2 * PI*fre;
	this->wavenumber = omega*sqrt(epsilon*mu);
	this->lambda = C0 / fre;
	this->polar = p;
	kDir << -sin(theta_i)*cos(phi_i), -sin(theta_i)*sin(phi_i), -cos(theta_i);
	if (polar == H)
	{
		polarDir = kDir.cross(Vector3d(0, 0, 1)).normalized();
	}
	else if (polar == V)
	{
		Vector3d t = kDir.cross(Vector3d(0, 0, 1)).normalized();
		polarDir = t.cross(kDir).normalized();
	}
}


PlaneWave::~PlaneWave()
{
}

Vector3cd PlaneWave::Eincident(Vector3d r)
{
	return magnitude*polarDir*exp(-im*wavenumber*kDir.dot(r));
}

Vector3cd PlaneWave::Hincident(Vector3d r)
{

	return kDir.cross(Eincident(r)) / eta;
}
