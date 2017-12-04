#include "stdafx.h"
#include "PO.h"

// declaration of determine whether point P in triangle ABC
bool PointinTriangle(Vector3d A, Vector3d B, Vector3d C, Vector3d P);

PO::PO(Mesh * mesh, PlaneWave & planewave)
	:mesh(mesh), planewave(planewave)
{
	//edgeN = mesh->getEdgeCount();
	edgeN = mesh->getPOEdgesCount();
	facetN = mesh->getTriangleCount();
	verN = mesh->getVertexCount();
	KaJ = VectorXcd::Zero(edgeN);
	this->delta = new int[facetN];
	for (int i = 0;i < facetN;i++)
	{
		delta[i] = 1; //defaulted, the whole PO region is lit
	}
	this->RWGdelta = new int[edgeN];
	for (int i = 0;i < edgeN;i++)
	{
		RWGdelta[i] = 1;
	}
	Triangles = mesh->getTriangles();
	Edges = mesh->getEdges();
	MOMEdges = mesh->getMOMEdges();
	POEdges = mesh->getPOEdges();
	Vertexes = mesh->getVertexes();
	deltaKaJ = VectorXcd::Zero(edgeN);
	judgeLitEdge();
	k = planewave.wavenumber;
	eta = planewave.eta;
	mu = planewave.mu;
	epsilon = planewave.epsilon;
	omega = planewave.omega;
}

void PO::judgeLitEdge()
{
	RWGPOdelta = new int*[edgeN];
	Vector3d kdir, momedgecenter;
	Triangle patch2;
	Vector3d rOM, rOV, rOP, n2;
	double alpha;
	for (int i = 0;i < edgeN;i++)
	{
		RWGPOdelta[i] = new int[mesh->getMOMEdgesCount()];
	}
	for (int i = 0;i < edgeN;i++)
	{
		for (int j = 0;j < mesh->getMOMEdgesCount();j++)
		{
			RWGPOdelta[i][j] = 1;
		}
	}
	for (int i = 0;i < edgeN;i++)
	{
		rOM = (Vertexes[Edges[POEdges[i]]->EdgeNode1] + Vertexes[Edges[POEdges[i]]->EdgeNode2]) / 2.0;
		for (int j = 0;j < mesh->getMOMEdgesCount();j++)
		{
			assert(Edges[POEdges[i]] != Edges[MOMEdges[j]]);
			momedgecenter= (Vertexes[Edges[MOMEdges[j]]->EdgeNode1] + Vertexes[Edges[MOMEdges[j]]->EdgeNode2]) / 2.0;
			kdir = (rOM - momedgecenter).normalized();
			for (int k = 0;k < facetN;k++)
			{
				if (k == Edges[POEdges[i]]->Triangle1 || k == Edges[POEdges[i]]->Triangle2 || k == Edges[MOMEdges[j]]->Triangle1 || k == Edges[MOMEdges[j]]->Triangle2)
					continue;
					
				patch2 = Triangles[k];
				n2 = patch2.Normal;
				rOV = Vertexes[patch2.Node1];
				alpha = (rOM - rOV).dot(n2) / kdir.dot(n2);

				rOP = rOM - alpha*kdir;
				if ((rOM - rOP).dot(kdir) > 1e-15 && (rOP - momedgecenter).dot(kdir) > 1e-15)
				{
					if (PointinTriangle(Vertexes[patch2.Node1], Vertexes[patch2.Node2], Vertexes[patch2.Node3], rOP))
					{
						RWGPOdelta[i][j] = 0;
						break;
					}
				}
			}
		}
	}
}

void PO::judgeLitPatch(Vector3d kDir)
//A facet if considered to be shadowed only if its central point is invisible to the source
//test passed!!!!!!!
{
	Triangle patch1, patch2;
	Vector3d rOM, rOV, rOP, n2;
	double alpha;
	for (int i = 0;i < facetN;i++)
	{
		patch1 = mesh->getTriangles()[i];
		rOM = patch1.Center;
		for (int j = 0;j < facetN;j++)
		{
			if (j == i)
				continue;
			patch2 = mesh->getTriangles()[j];
			n2 = patch2.Normal;
			rOV = Vertexes[patch2.Node1];
			alpha = (rOM - rOV).dot(n2) / kDir.dot(n2);
			if (alpha > 1e-10)
			//duration 1e-10~1e-15 is a reasonable choice for FEKO mesh
			{
				rOP = rOM - alpha*kDir;
				if (PointinTriangle(Vertexes[patch2.Node1], Vertexes[patch2.Node2], Vertexes[patch2.Node3], rOP))
				{
					delta[i] = 0;
					break;
				}	
			}
		}
	}
}

void PO::WriteLitPatch()
{
	fstream LitFileName;
	LitFileName.open(mesh->getFileName().substr(0, mesh->getFileName().find(".")) + ".lit", ios::out);
	if (LitFileName.fail())
	{
		throw "File open failed";
	}
	LitFileName << verN << endl;
	for (int i = 0;i < verN;i++)
	{
		LitFileName << Vertexes[i].transpose() << endl;
	}
	int i = 0;
	litPatchN = 0;
	for (int i = 0;i < facetN;i++)
	{
		if (delta[i] == 1)
		{
			LitFileName << Triangles[i].Node1 << " " << Triangles[i].Node2 << " " << Triangles[i].Node3 << endl;
			litPatchN++;
		}

	}
	LitFileName.close();
}

void PO::LitPatchInfo()
{
	cout << "Number of lit patches: " << litPatchN << endl;
	cout << "----------------------------------------------------" << endl;
	return;
}

complex<double> PO::gradientG(double R)
{
	return (1.0 + im*k*R)*exp(-im*k*R) / R / R / R;
}

VectorXcd PO::Solver()
//compute the current over the lit patch based on RWG basis
{
	Vector3d edgecenter, n_plus, n_minus, t_plus, t_minus, l, rho_plus, rho_minus;
	for (int i = 0;i < edgeN; i++)
	{
		RWGdelta[i] = 0;
		KaJ[i] = 0;
		if (delta[Edges[POEdges[i]]->Triangle1]==1 && delta[Edges[POEdges[i]]->Triangle2]==1)
		{
			RWGdelta[i] = 1;
			edgecenter = (Vertexes[Edges[POEdges[i]]->EdgeNode1] + Vertexes[Edges[POEdges[i]]->EdgeNode2]) / 2.0;
			n_plus = Triangles[Edges[POEdges[i]]->Triangle1].Normal;
			n_minus = Triangles[Edges[POEdges[i]]->Triangle2].Normal;
			l = (Vertexes[Edges[POEdges[i]]->EdgeNode1] - Vertexes[Edges[POEdges[i]]->EdgeNode2]).normalized();
			rho_plus = Triangles[Edges[POEdges[i]]->Triangle1].Center - Vertexes[Edges[POEdges[i]]->EdgeNode3];
			rho_minus = Vertexes[Edges[POEdges[i]]->EdgeNode4] - Triangles[Edges[POEdges[i]]->Triangle2].Center;
			t_plus = (rho_plus - rho_plus.dot(l)*l).normalized();
			t_minus = (rho_minus - rho_minus.dot(l)*l).normalized();
			KaJ[i] = (t_plus.dot(n_minus.cross(planewave.Hincident(edgecenter))) + t_minus.dot(n_minus.cross(planewave.Hincident(edgecenter))));
		}
	}
	return KaJ;
}

complex<double> PO::nearHField(Vector3d r, VectorXcd MOMJ, int n, Vector3d n_plus, Vector3d n_minus, Vector3d t_plus, Vector3d t_minus)
{
	Vector3d rm[M], rn[M], vm[3], vn[3];
	complex<double> I, result;
	int trim = 0, trin = 0;
	Vector3d Node3 = Vector3d::Zero();
	Vector3d Node4 = Vector3d::Zero();
	double R1 = 0, R2 = 0;
	for (int m = 0;m < mesh->getMOMEdgesCount();m++)
	{
		if (RWGPOdelta[n][m] == 1)
		{
			trim = Edges[MOMEdges[m]]->Triangle1;
			trin = Edges[MOMEdges[m]]->Triangle2;
			vm[0] = Vertexes[Triangles[trim].Node1];
			vm[1] = Vertexes[Triangles[trim].Node2];
			vm[2] = Vertexes[Triangles[trim].Node3];
			vn[0] = Vertexes[Triangles[trin].Node1];
			vn[1] = Vertexes[Triangles[trin].Node2];
			vn[2] = Vertexes[Triangles[trin].Node3];
			Node3 = Vertexes[Edges[MOMEdges[m]]->EdgeNode3];
			Node4 = Vertexes[Edges[MOMEdges[m]]->EdgeNode4];
			I = 0;
			for (int i = 0;i < M;i++)
			{
				rm[i] = vm[0] * xx7[i] + vm[1] * yy7[i] + vm[2] * zz7[i];
				rn[i] = vn[0] * xx7[i] + vn[1] * yy7[i] + vn[2] * zz7[i];
				R1 = (rm[i] - r).norm();
				R2 = (rn[i] - r).norm();
				I = I + ww7[i] * (t_plus.dot(n_plus.cross((rm[i] - Node3).cross(r - rm[i])*gradientG(R1))) + t_minus.dot(n_minus.cross((Node4 - rn[i]).cross(r - rn[i])*gradientG(R2))))*Edges[MOMEdges[m]]->len;
			}
			result -= I*MOMJ[m] / (4.0*PI);
		}
	}
	return result;
}

VectorXcd PO::updateCurrent(VectorXcd MOMJ)
{
	//cout << MOMJ.size() << endl;
	deltaKaJ = VectorXcd::Zero(edgeN);
	VectorXcd J = VectorXcd::Zero(edgeN);
	Vector3d edgecenter, n_plus, n_minus, t_plus, t_minus, l, rho_plus, rho_minus;
	for (int i = 0;i < edgeN; i++)
	{
		edgecenter = (Vertexes[Edges[POEdges[i]]->EdgeNode1] + Vertexes[Edges[POEdges[i]]->EdgeNode2]) / 2.0;
		n_plus = Triangles[Edges[POEdges[i]]->Triangle1].Normal;
		n_minus = Triangles[Edges[POEdges[i]]->Triangle2].Normal;
		l = (Vertexes[Edges[POEdges[i]]->EdgeNode1] - Vertexes[Edges[POEdges[i]]->EdgeNode2]).normalized();
		rho_plus = Triangles[Edges[POEdges[i]]->Triangle1].Center - Vertexes[Edges[POEdges[i]]->EdgeNode3];
		rho_minus = Vertexes[Edges[POEdges[i]]->EdgeNode4] - Triangles[Edges[POEdges[i]]->Triangle2].Center;
		t_plus = (rho_plus - rho_plus.dot(l)*l).normalized();
		t_minus = (rho_minus - rho_minus.dot(l)*l).normalized();
		deltaKaJ[i] = this->nearHField(edgecenter, MOMJ, i, n_plus, n_minus, t_plus, t_minus);//????????
		//cout << KaJ[i] << " " << deltaKaJ[i] << endl;
	}
	return deltaKaJ;
}

VectorXcd PO::getCurrent()
{
	return KaJ + deltaKaJ;
}

Vector3cd PO::IntegralEfield(int m, double phi, double theta)
{
	Vector3d rm[M], rn[M], vm[3], vn[3];
	Vector3cd I = Vector3cd::Zero();
	int trim = Edges[POEdges[m]]->Triangle1;
	int trin = Edges[POEdges[m]]->Triangle2;
	vm[0] = Vertexes[Triangles[trim].Node1];
	vm[1] = Vertexes[Triangles[trim].Node2];
	vm[2] = Vertexes[Triangles[trim].Node3];
	vn[0] = Vertexes[Triangles[trin].Node1];
	vn[1] = Vertexes[Triangles[trin].Node2];
	vn[2] = Vertexes[Triangles[trin].Node3];
	Vector3d Node3 = Vertexes[Edges[POEdges[m]]->EdgeNode3];
	Vector3d Node4 = Vertexes[Edges[POEdges[m]]->EdgeNode4];
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
	return -im*omega*mu*I*Edges[POEdges[m]]->len*KaJ[m] * green(R) / (4.0 * PI);
}

Vector3cd PO::ScatteredField(double phi, double theta)
{
	Vector3cd Efield = Vector3cd::Zero();
	for (int m = 0;m < edgeN;m++)
	{
		Efield += IntegralEfield(m, phi, theta);
	}
	return Efield;
}

void PO::RCS(double theta)
{
	double phi_start = 0.0, phi_end = 2 * PI, phi;
	double R = 1000.0;
	double step_phi = (phi_end - phi_start) / 360.0;
	Vector3cd Efield = Vector3cd::Zero();
	double RCS[361] = { 0 };
	ofstream RCSFile("RCS_PO.txt", ios::out);
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

// Determine whether point P in triangle ABC
bool PointinTriangle(Vector3d A, Vector3d B, Vector3d C, Vector3d P)
{
	Vector3d v0 = C - A;
	Vector3d v1 = B - A;
	Vector3d v2 = P - A;
	double dot00 = v0.dot(v0);
	double dot01 = v0.dot(v1);
	double dot02 = v0.dot(v2);
	double dot11 = v1.dot(v1);
	double dot12 = v1.dot(v2);
	double inverDeno = 1 / (dot00 * dot11 - dot01 * dot01);
	double u = (dot11 * dot02 - dot01 * dot12) * inverDeno;
	if (u < 0.0 || u > 1.0) // if u out of range, return directly
	{
		return false;
	}
	double v = (dot00 * dot12 - dot01 * dot02) * inverDeno;
	if (v < 0.0 || v > 1.0) // if v out of range, return directly
	{
		return false;
	}
	return u + v <= 1.0;
}

PO::~PO()
{
}
