#include "stdafx.h"
#include "MOM.h"


MOM::MOM(Mesh *mesh, PlaneWave &planewave):mesh(mesh), planewave(planewave)
{
	edgeN = mesh->getEdgeCount();
	ZMatrix = MatrixXcd::Zero(edgeN, edgeN);
	rhs = VectorXcd::Zero(edgeN);
	X = VectorXcd::Zero(edgeN);
	Triangles = mesh->getTriangles();
	Edges = mesh->getEdges();
	Vertexes = mesh->getVertexes();
	k = planewave.wavenumber;
	eta = planewave.eta;
	mu = planewave.mu;
	epsilon = planewave.epsilon;
	omega = planewave.omega;
}

void MOM::FillMatrix()
{
	complex<double> zll = 0, zlr = 0, zrl = 0, zrr = 0;
	for (int m = 0;m < edgeN;m++)
	{
		//cout << m << endl;
		for (int n = 0;n < edgeN;n++)
		{
			int tml, tmr, tnl, tnr, vm3, vn3;
			tml = Edges[m]->Triangle1;
			tmr = Edges[m]->Triangle2;
			tnl = Edges[n]->Triangle1;
			tnr = Edges[n]->Triangle2;
			vm3 = Edges[m]->EdgeNode3;
			vn3 = Edges[n]->EdgeNode3;
			if (tml != tnl)
			{
				zll= Integral(m, n, tml, tnl, vm3, vn3);
			}
			else
			{
				zll = SingularIntegral(m, n, tml, tnl, vm3, vn3);
			}
			vm3 = Edges[m]->EdgeNode3;
			vn3 = Edges[n]->EdgeNode4;
			if (tml != tnr)
			{
				zlr = Integral(m, n, tml, tnr, vm3, vn3);
			}
			else
			{
				zlr = SingularIntegral(m, n, tml, tnr, vm3, vn3);
			}
			vm3 = Edges[m]->EdgeNode4;
			vn3 = Edges[n]->EdgeNode3;
			if (tmr != tnl)
			{
				zrl = Integral(m, n, tmr, tnl, vm3, vn3);
			}
			else
			{
				zrl =  SingularIntegral(m, n, tmr, tnl, vm3, vn3);
			}
			vm3 = Edges[m]->EdgeNode4;
			vn3 = Edges[n]->EdgeNode4;
			if (tmr != tnr)
			{
				zrr = Integral(m, n, tmr, tnr, vm3, vn3);
			}
			else
			{
				zrr = SingularIntegral(m, n, tmr, tnr, vm3, vn3);
			}
			ZMatrix(m, n) = zll - zlr - zrl + zrr;
		}
	}
}

complex<double> MOM::Integral(int m, int n, int trim, int trin, int vm3, int vn3)
{
	Vector3d rm[N], rn[N];
	complex<double> IA = 0;
	Vector3d vm[3], vn[3], vm3c, vn3c;
	double R = 0, lm, ln;
	vm[0] = Vertexes[Triangles[trim].Node1];
	vm[1] = Vertexes[Triangles[trim].Node2];
	vm[2] = Vertexes[Triangles[trim].Node3];
	vn[0] = Vertexes[Triangles[trin].Node1];
	vn[1] = Vertexes[Triangles[trin].Node2];
	vn[2] = Vertexes[Triangles[trin].Node3];
	vm3c = Vertexes[vm3];
	vn3c = Vertexes[vn3];
	lm = Edges[m]->len;
	ln = Edges[n]->len;
	for (int i = 0;i < N;i++)
	{
		rm[i] = vm[0] * xx4[i] + vm[1] * yy4[i] + vm[2] * zz4[i];
		rn[i] = vn[0] * xx4[i] + vn[1] * yy4[i] + vn[2] * zz4[i];
	}
	for (int i = 0;i < N;i++)
	{
		for (int j = 0;j < N;j++)
		{
			R = (rm[i] - rn[j]).norm();
			IA = IA + ww4[i] * ww4[j] * ((rm[i] - vm3c).dot(rn[j] - vn3c) - 4.0 / k / k)*green(R);
		}
	}
	return im*k*eta*lm*ln*IA / (4.0*PI);
}

complex<double> MOM::SingularIntegral(int m, int n, int trim, int trin, int vm3, int vn3)
{
	Vector3d rm[M], rn[M];
	complex<double> IA = 0;
	Vector3d vm[3], vn[3], vm3c, vn3c;
	double R = 0, lm, ln;
	vm[0] = Vertexes[Triangles[trim].Node1];
	vm[1] = Vertexes[Triangles[trim].Node2];
	vm[2] = Vertexes[Triangles[trim].Node3];
	vn[0] = Vertexes[Triangles[trin].Node1];
	vn[1] = Vertexes[Triangles[trin].Node2];
	vn[2] = Vertexes[Triangles[trin].Node3];
	vm3c = Vertexes[vm3];
	vn3c = Vertexes[vn3];
	lm = Edges[m]->len;
	ln = Edges[n]->len;
	for (int i = 0;i < M;i++)
	{
		rm[i] = vm[0] * xx7[i] + vm[1] * yy7[i] + vm[2] * zz7[i];
		rn[i] = vn[0] * xx7[i] + vn[1] * yy7[i] + vn[2] * zz7[i];
	}
	for (int i = 0;i < M;i++)
	{
		for (int j = 0;j < M;j++)
		{
			R = (rm[i] - rn[j]).norm();
			IA = IA + ww7[i] * ww7[j] * ((rm[i] - vm3c).dot(rn[j] - vn3c) - 4.0 / k / k)*Singulargreen(R);
		}
	}
	assert(trim == trin);
	//Singular part of the singular integral
	Vector3d rho[3], r12, r23, r31;
	Vector3d normal[3];
	double area = Triangles[trim].Area;
	double l_plus[3], l_minus[3], R_plus[3], R_minus[3], f[3], P[3];
	r12 = (vm[1] - vm[0]).normalized();
	r23 = (vm[2] - vm[1]).normalized();
	r31 = (vm[0] - vm[2]).normalized();
	double Is = 0;
	for (int i = 0;i < M;i++)
	{
		rho[0] = vm[0] - rm[i];
		rho[1] = vm[1] - rm[i];
		rho[2] = vm[2] - rm[i];
		l_minus[0] = rho[0].dot(r12);
		l_plus[0] = rho[1].dot(r12);
		l_minus[1] = rho[1].dot(r23);
		l_plus[1] = rho[2].dot(r23);
		l_minus[2] = rho[2].dot(r31);
		l_plus[2] = rho[0].dot(r31);
		normal[0] = rho[0] - l_minus[0] * r12;
		normal[1] = rho[1] - l_minus[1] * r23;
		normal[2] = rho[2] - l_minus[2] * r31;
		for (int j = 0;j < 3;j++)
		{
			P[j] = normal[j].norm();
			normal[j].normalize();
			R_plus[j] = sqrt(l_plus[j] * l_plus[j] + P[j] * P[j]);
			R_minus[j] = sqrt(l_minus[j] * l_minus[j] + P[j] * P[j]);
			f[j] = log((l_plus[j] + R_plus[j]) / (l_minus[j] + R_minus[j]));
		}
		Vector3d Ip = Vector3d::Zero();
		double IR = 0;
		for (int j = 0;j < 3;j++)
		{
			Ip += normal[j] * (P[j] * P[j] * f[j] + l_plus[j] * R_plus[j] - l_minus[j] * R_minus[j]) / 2.0;
			if (P[j]>1e-10)
				IR += P[j] * f[j];
		}
		Is = Is + ww7[i] * (((rm[i] - vm3c).dot(rm[i] - vn3c) - 4.0 / k / k)*IR + (rm[i] - vm3c).dot(Ip));
	}
	return im*k*eta*lm*ln*Is / (4.0*PI) / 2.0 / area + im*k*eta*lm*ln*IA / (4.0*PI);
	//return im*k*eta*lm*ln*IA / (4.0*PI);
}

complex<double> MOM::green(double R)
{
	return exp(-im*k*R) / R;
}

complex<double> MOM::Singulargreen(double R)
{
	//Ì©ÀÕÕ¹¿ª
	return -im*k - 0.5*k*k*R + im*pow(k, 3)*pow(R, 2) / 6.0;
}

void MOM::FillRhs()
{
	for (int m = 0;m < edgeN;m++)
	{
		rhs(m) = IntegralRhs(m);
	}
}

complex<double> MOM::IntegralRhs(int m)
{
	Vector3d rm[M], rn[M], vm[3], vn[3];
	complex<double> I = 0;
	int trim = Edges[m]->Triangle1;
	int trin = Edges[m]->Triangle2;
	vm[0] = Vertexes[Triangles[trim].Node1];
	vm[1] = Vertexes[Triangles[trim].Node2];
	vm[2] = Vertexes[Triangles[trim].Node3];
	vn[0] = Vertexes[Triangles[trin].Node1];
	vn[1] = Vertexes[Triangles[trin].Node2];
	vn[2] = Vertexes[Triangles[trin].Node3];
	Vector3d Node3 = Vertexes[Edges[m]->EdgeNode3];
	Vector3d Node4 = Vertexes[Edges[m]->EdgeNode4];
	for (int i = 0;i < M;i++)
	{
		rm[i] = vm[0] * xx7[i] + vm[1] * yy7[i] + vm[2] * zz7[i];
		rn[i] = vn[0] * xx7[i] + vn[1] * yy7[i] + vn[2] * zz7[i];
		I = I + ww7[i] * ((rm[i] - Node3).dot(planewave.Eincident(rm[i])) + (Node4 - rn[i]).dot(planewave.Eincident(rn[i])));
	}
	return I*Edges[m]->len;
}

Vector3cd MOM::IntegralEfield(int m, double phi, double theta)
{
	Vector3d rm[M], rn[M], vm[3], vn[3];
	Vector3cd I = Vector3cd::Zero();
	int trim = Edges[m]->Triangle1;
	int trin = Edges[m]->Triangle2;
	vm[0] = Vertexes[Triangles[trim].Node1];
	vm[1] = Vertexes[Triangles[trim].Node2];
	vm[2] = Vertexes[Triangles[trim].Node3];
	vn[0] = Vertexes[Triangles[trin].Node1];
	vn[1] = Vertexes[Triangles[trin].Node2];
	vn[2] = Vertexes[Triangles[trin].Node3];
	Vector3d Node3 = Vertexes[Edges[m]->EdgeNode3];
	Vector3d Node4 = Vertexes[Edges[m]->EdgeNode4];
	double R = 1000;
	Vector3d rDir = Vector3d::Zero();
	rDir << sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta);
	for (int i = 0;i < M;i++)
	{
		rm[i] = vm[0] * xx7[i] + vm[1] * yy7[i] + vm[2] * zz7[i];
		rn[i] = vn[0] * xx7[i] + vn[1] * yy7[i] + vn[2] * zz7[i];
		//I = I + ww7[i] * (rDir.cross(rDir.cross((rm[i] - Node3)))*exp(im*k*rDir.dot(rm[i])) + rDir.cross(rDir.cross((Node4 - rn[i])))*exp(im*k*rDir.dot(rn[i])));
		I = I + ww7[i] * ((rm[i] - Node3)*exp(im*k*rDir.dot(rm[i])) + (Node4 - rn[i])*exp(im*k*rDir.dot(rn[i])));
	}
	return -im*omega*mu*I*Edges[m]->len*X[m] * green(R) / (4 * PI);
}

void MOM::Solver()
{
	this->X = ZMatrix.householderQr().solve(rhs);
	double relative_error = (ZMatrix*X - rhs).norm() / rhs.norm(); // norm() is L2 norm
	cout << "The relative error is:\n" << relative_error << endl;
}

Vector3cd MOM::ScatteredField(double phi, double theta)
{
	Vector3cd Efield = Vector3cd::Zero();
	for (int m = 0;m < edgeN;m++)
	{
		Efield += IntegralEfield(m, phi, theta);
	}
	return Efield;
}

void MOM::RCS()
{
	double phi_start = 0.0, phi_end = 2 * PI, phi;
	double theta = PI / 2.0;
	double step_phi = (phi_end - phi_start) / 360.0;
	double R = 1000;
	Vector3cd Efield = Vector3cd::Zero();
	double RCS[361] = { 0 };
	ofstream RCSFile("RCS.txt", ios::out);
	Vector3d theta_dir = Vector3d(0, 0, 1);
	Vector3d phi_dir;
	for (int i = 0;i < 361;i++)
	{
		phi = i*step_phi;
		phi_dir << -sin(phi), cos(phi), 0;
		Efield = ScatteredField(phi, theta);
		RCS[i] = 4.0 * PI*R*R*std::norm(Efield.dot(phi_dir));
		RCSFile << 10 * log10(RCS[i]) << endl;
	}
	return;
}


MOM::~MOM()
{
}
