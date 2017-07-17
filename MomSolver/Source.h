#pragma once
#define PI 3.14159265358979323846264338327950288419716939937510
#define EIGEN_USE_MKL_ALL
#include <Eigen\Dense>
#include <complex>
using namespace std;
using namespace Eigen;
typedef struct ConstantStruct
{
	double epsilon = 8.85418781761e-12;
	double mu = 4 * PI*1e-7;
	double C0 = 1.0 / sqrt(epsilon*mu);
	double eta = sqrt(mu / epsilon);
}ElecConst;
class Source:
	public ElecConst
{
public:
	Source();
	~Source();
};

