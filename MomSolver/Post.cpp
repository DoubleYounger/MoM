#include "stdafx.h"
#include "Post.h"


Post::Post(Mesh * mesh, PO & posolver, MOM & momsolver, PlaneWave &planewave):
	mesh(mesh),posolver(posolver),momsolver(momsolver), planewave(planewave)
{
	POedgeN = mesh->getPOEdgesCount();
	MOMedgeN = mesh->getMOMEdgesCount();
	Triangles = mesh->getTriangles();
	Edges = mesh->getEdges();
	MOMEdges = mesh->getMOMEdges();
	POEdges = mesh->getPOEdges();
	Vertexes = mesh->getVertexes();
	k = planewave.wavenumber;
	eta = planewave.eta;
	mu = planewave.mu;
	epsilon = planewave.epsilon;
	omega = planewave.omega;
	POCurrent = posolver.getCurrent();
	MOMCurrent = momsolver.getCurrent();
}

void Post::RCS(double theta)
{
	double phi_start = 0.0, phi_end = 2 * PI, phi;
	double R = 1000.0;
	double step_phi = (phi_end - phi_start) / 360.0;
	Vector3cd Efield = Vector3cd::Zero();
	double RCS[361] = { 0 };
	ofstream RCSFile("RCS_PO_MOM.txt", ios::out);
	//Vector3d theta_dir = Vector3d(0, 0, 1);
	//Vector3d phi_dir;
	for (int i = 0;i < 361;i++)
	{
		phi = i*step_phi;
		//phi_dir << -sin(phi), cos(phi), 0;
		Efield = ScatteredField(phi, theta);
		RCS[i] = 2.0 * sqrt(PI)*R*Efield.norm();
		//RCS[i] = 4.0 * PI*R*R*std::norm(Efield.dot(phi_dir));
		RCSFile << 20 * log10(RCS[i]) << endl;
	}
	return;
}

Vector3cd Post::ScatteredField(double phi, double theta)
{
	Vector3cd Efield = Vector3cd::Zero();
	for (int m = 0;m < POedgeN;m++)
	{
		Efield += IntegralEfieldPO(m, phi, theta);
	}
	for (int m = 0;m < MOMedgeN;m++)
	{
		Efield += IntegralEfieldMOM(m, phi, theta);
	}
	return Efield;
}

Vector3cd Post::IntegralEfieldPO(int m, double phi, double theta)
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
	return -im*omega*mu*I*Edges[POEdges[m]]->len*POCurrent[m] * green(R) / (4.0 * PI);
}

Vector3cd Post::IntegralEfieldMOM(int m, double phi, double theta)
{
	Vector3d rm[M], rn[M], vm[3], vn[3];
	Vector3cd I = Vector3cd::Zero();
	int trim = Edges[MOMEdges[m]]->Triangle1;
	int trin = Edges[MOMEdges[m]]->Triangle2;
	vm[0] = Vertexes[Triangles[trim].Node1];
	vm[1] = Vertexes[Triangles[trim].Node2];
	vm[2] = Vertexes[Triangles[trim].Node3];
	vn[0] = Vertexes[Triangles[trin].Node1];
	vn[1] = Vertexes[Triangles[trin].Node2];
	vn[2] = Vertexes[Triangles[trin].Node3];
	Vector3d Node3 = Vertexes[Edges[MOMEdges[m]]->EdgeNode3];
	Vector3d Node4 = Vertexes[Edges[MOMEdges[m]]->EdgeNode4];
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
	return -im*omega*mu*I*Edges[MOMEdges[m]]->len*MOMCurrent[m] * green(R) / (4 * PI);
}

Post::~Post()
{
}
