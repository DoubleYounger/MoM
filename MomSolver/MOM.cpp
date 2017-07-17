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
	complex<double> IA = 0;
	for (int n = 0;n < edgeN;n++)
	{
		cout << n << endl;
		for (int m = 0;m < edgeN;m++)
		{
			complex<double> IA = 0;
			int tml, tmr, tnl, tnr, vm3, vn3;
			tml = Edges[m]->Triangle1;
			tmr = Edges[m]->Triangle2;
			tnl = Edges[n]->Triangle1;
			tnr = Edges[n]->Triangle2;
			vm3 = Edges[m]->EdgeNode3;
			vn3 = Edges[n]->EdgeNode3;
			if (tml != tnl)
			{
				IA = IA + Integral(m, n, tml, tnl, vm3, vn3);
			}
			else
			{
				IA = IA + SingularIntegral(m, n, tml, tnl, vm3, vn3);
			}
			vm3 = Edges[m]->EdgeNode3;
			vn3 = Edges[n]->EdgeNode4;
			if (tml != tnr)
			{
				IA = IA - Integral(m, n, tml, tnr, vm3, vn3);
			}
			else
			{
				IA = IA - SingularIntegral(m, n, tml, tnr, vm3, vn3);
			}
			vm3 = Edges[m]->EdgeNode4;
			vn3 = Edges[n]->EdgeNode3;
			if (tmr != tnl)
			{
				IA = IA - Integral(m, n, tmr, tnl, vm3, vn3);
			}
			else
			{
				IA = IA - SingularIntegral(m, n, tmr, tnl, vm3, vn3);
			}
			vm3 = Edges[m]->EdgeNode4;
			vn3 = Edges[n]->EdgeNode4;
			if (tmr != tnr)
			{
				IA = IA + Integral(m, n, tmr, tnr, vm3, vn3);
			}
			else
			{
				IA = IA + SingularIntegral(m, n, tmr, tnr, vm3, vn3);
			}
			ZMatrix(m, n) = IA;
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
	complex<double> nonsingularI = im*k*eta*lm*ln*IA / (4.0*PI);
	assert(trim == trin);
	Vector3d center = Triangles[trim].Center;

	//Singular part of the singular integral
	Vector3d rho[3], r12, r23, r31;
	Vector3d normal[3];
	double area = Triangles[trim].Area;
	double l_plus[3], l_minus[3], R_plus[3], R_minus[3], f[3], P[3];
	r12 = (vm[1] - vm[0]).normalized();
	r23 = (vm[2] - vm[1]).normalized();
	r31 = (vm[0] - vm[2]).normalized();
	Vector3d Ip = Vector3d::Zero();
	double IR = 0, Is = 0;
	for (int i = 0;i < M;i++)
	{
		rho[0] = vm[0] - rm[0];
		rho[1] = vm[1] - rm[1];
		rho[2] = vm[2] - rm[2];
		l_minus[0] = rho[0].dot(r12);
		l_plus[0] = rho[1].dot(r12);
		l_minus[1] = rho[1].dot(r23);
		l_plus[1] = rho[2].dot(r23);
		l_minus[2] = rho[2].dot(r31);
		l_plus[2] = rho[0].dot(r31);
		assert(l_minus[0] < 1e-15&&l_minus[1] < 1e-15&&l_minus[2] < 1e-15);
		normal[0] = rho[0] - l_minus[0] * r12;
		normal[1] = rho[1] - l_minus[1] * r23;
		normal[2] = rho[2] - l_minus[2] * r31;
		for (int j = 0;j < 3;j++)
		{
			P[i] = normal[i].norm();
			normal[i].normalize();
		}
		for (int j = 0;j < 3;j++)
		{
			Ip += normal[j] * (P[j] * P[j] * f[j] + l_plus[j] * R_plus[j] - l_minus[j] * R_minus[j]) / 2.0;
			IR += P[j] * f[j];
		}
		Is = Is + ww7[i] * (((rm[i] - vm3c).dot(rm[i] - vn3c) - 4.0 / k / k)*IR + (rm[i] - vm3c).dot(Ip));
	}
	return im*k*eta*lm*ln*Is / (4.0*PI) / 2.0 / area + im*k*eta*lm*ln*IA / (4.0*PI);
}

//complex<double> MOM::SingularIntegral(int m, int n, int trim, int trin, int vm3, int vn3)
//{
//	assert(trim.Node1 == trin.Node1);
//	assert(trim.Node2 == trin.Node2);
//	assert(trim.Node3 == trin.Node3);
//	assert(RWGm->l - RWGn->l < 1e-40);
//	Vector3d r[7];
//	complex<double> IA = 0, IF = 0;
//	complex<double> im(0, 1);
//	double k = planewave.wavenumber;
//	double eta = planewave.eta;
//	for (int i = 0;i < 7;i++)
//	{
//		r[i] = trin.Center + Vertexes[trin.Node1] * (-u[i] + sqrt(3.0)*v[i]) / 3.0 + Vertexes[trin.Node2] * (-u[i] - sqrt(3.0)*v[i]) / 3.0 + 2.0*u[i] * Vertexes[trin.Node3] / 3.0;
//		IA = IA + w[i] * ((trim.Center - vm).dot(r[i] - vn) - 4.0 / k / k)*Singulargreen(trim.Center, r[i]);
//	}
//	Vector3d rho_m_c = trim.Center - vm,
//		rho_n_c = trin.Center - vn;
//	Vector3d rho[3];
//	rho[0] = Vertexes[trin.Node1] - trin.Center;
//	rho[1] = Vertexes[trin.Node2] - trin.Center;
//	rho[2] = Vertexes[trin.Node3] - trin.Center;
//	Vector3d r12, r23, r31;
//	r12 = Vertexes[trin.Node2] - Vertexes[trin.Node1];
//	r23 = Vertexes[trin.Node3] - Vertexes[trin.Node2];
//	r31 = Vertexes[trin.Node1] - Vertexes[trin.Node3];
//	double sin_beta[3], cos_beta[3], sin_beta_i[3], cos_beta_i[3];
//	sin_beta[0] = -rho[0].dot(r12) / (rho[0].norm()*r12.norm());
//	sin_beta[1] = -rho[1].dot(r23) / (rho[1].norm()*r23.norm());
//	sin_beta[2] = -rho[2].dot(r31) / (rho[2].norm()*r31.norm());
//	cos_beta[0] = sqrt(1 - pow(sin_beta[0], 2));
//	cos_beta[1] = sqrt(1 - pow(sin_beta[1], 2));
//	cos_beta[2] = sqrt(1 - pow(sin_beta[2], 2));
//	sin_beta_i[0] = rho[1].dot(r12) / (rho[1].norm()*r12.norm());
//	sin_beta_i[1] = rho[2].dot(r23) / (rho[2].norm()*r23.norm());
//	sin_beta_i[2] = rho[0].dot(r31) / (rho[0].norm()*r31.norm());
//	cos_beta_i[0] = sqrt(1 - pow(sin_beta_i[0], 2));
//	cos_beta_i[1] = sqrt(1 - pow(sin_beta_i[1], 2));
//	cos_beta_i[2] = sqrt(1 - pow(sin_beta_i[2], 2));
//	double A[3], B[3], p[3];
//	p[0] = rho[0].norm()*cos_beta[0];
//	p[1] = rho[1].norm()*cos_beta[1];
//	p[2] = rho[2].norm()*cos_beta[2];
//	A[0] = rho_m_c.dot(rho[0].normalized() + r12.normalized()*rho[0].norm()*sin_beta[0]);
//	A[1] = rho_m_c.dot(rho[1].normalized() + r23.normalized()*rho[1].norm()*sin_beta[1]);
//	A[2] = rho_m_c.dot(rho[2].normalized() + r31.normalized()*rho[2].norm()*sin_beta[2]);
//	B[0] = p[0] * rho_m_c.dot(r12.normalized());
//	B[1] = p[1] * rho_m_c.dot(r23.normalized());
//	B[2] = p[2] * rho_m_c.dot(r31.normalized());
//	double SF = 0, SA = 0;
//	for (int i = 0;i < 3;i++)
//	{
//		SF += p[i] * log(cos_beta[i] * (1 + sin_beta_i[i]) / (cos_beta_i[i] * (1 - sin_beta[i])));
//		SA += p[i] * (A[i] * log(cos_beta[i] * (1 + sin_beta_i[i]) / (cos_beta_i[i] * (1 - sin_beta[i]))) + B[i] * (1.0 / cos_beta_i[i] - 1.0 / cos_beta[i]));
//	}
//	SA += rho_m_c.dot(rho_n_c)*SF;
//	IA = (SA - 4.0*SF / k / k) / trin.Area;
//	return im*k*eta*IA*RWGm->l*RWGn->l / 4.0 / (4.0*PI);
//}

complex<double> MOM::green(double R)
{
	return exp(-im*planewave.wavenumber*R) / R;
}

complex<double> MOM::Singulargreen(double R)
{
	//Ì©ÀÕÕ¹¿ª
	return -im*k - 0.5*k*k*R + im*pow(k, 3)*pow(R, 2) / 6.0 + pow(k, 4)*pow(R, 3) / 24.0;
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
	Vector3d rm[7], rn[7], vm[3], vn[3];
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

void MOM::Solver()
{
	this->X = ZMatrix.fullPivHouseholderQr().solve(rhs);
	double relative_error = (ZMatrix*X - rhs).norm() / rhs.norm(); // norm() is L2 norm
	cout << X << endl;
	cout << "The relative error is:\n" << relative_error << endl;
}

Vector3cd MOM::ScatteredField(double phi, double theta)
{
	double R = 1000;
	Vector3d rDir = Vector3d::Zero();
	Vector3d r_plus = Vector3d::Zero();
	Vector3d r_minus = Vector3d::Zero();
	rDir << sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta);
	Vector3cd Efield = Vector3cd::Zero();
	double l = 0;
	Vector3d rho_plus = Vector3d::Zero();
	Vector3d rho_minus = Vector3d::Zero();
	complex<double> im(0, 1);
	double k = planewave.wavenumber;
	for (int i = 0;i < edgeN;i++)
	{
		l = (Vertexes[Edges[i]->EdgeNode1] - Vertexes[Edges[i]->EdgeNode2]).norm();
		r_plus = Triangles[Edges[i]->Triangle1].Center;
		r_minus = Triangles[Edges[i]->Triangle2].Center;
		rho_plus = Triangles[Edges[i]->Triangle1].Center-Vertexes[Edges[i]->EdgeNode3];
		rho_minus = Vertexes[Edges[i]->EdgeNode4] - Triangles[Edges[i]->Triangle2].Center;
		Efield = Efield + X[i] * (l/2.0*rho_plus*exp(im*planewave.wavenumber*rDir.dot(r_plus)) + l/2.0*rho_minus*exp(im*k*rDir.dot(r_minus)))*exp(-im*k*R);
	}
	double omega = planewave.omega;
	double mu = planewave.mu;
	return Efield = -im*omega*mu*Efield / (4.0*PI*R);
}

void MOM::RCS()
{
	double phi_start = 0.0, phi_end = 2 * PI, phi;
	double theta = PI / 2.0;
	double step_phi = (phi_start - phi_end) / 360.0;
	double R = 1000;
	Vector3d rDir = Vector3d::Zero();
	Vector3d r = Vector3d::Zero();
	Vector3cd Efield = Vector3cd::Zero();
	double RCS[180] = { 0 };
	ofstream RCSFile("RCS.txt", ios::out);
	for (int i = 0;i < 360;i++)
	{
		phi = i*step_phi;
		Efield = ScatteredField(phi, theta);
		RCS[i] = 4.0 * PI*R*R*pow(Efield.norm(), 2);
		RCSFile << 10 * log10(RCS[i]) << endl;
	}
	return;
}


MOM::~MOM()
{
}
