#include "stdafx.h"
#include "PO.h"

// declaration of determine whether point P in triangle ABC
bool PointinTriangle(Vector3d A, Vector3d B, Vector3d C, Vector3d P);

PO::PO(Mesh * mesh, PlaneWave & planewave)
	:mesh(mesh), planewave(planewave)
{
	edgeN = mesh->getEdgeCount();
	facetN = mesh->getTriangleCount();
	verN = mesh->getVertexCount();
	KaJ = VectorXcd::Zero(edgeN);
	this->delta = new double[facetN];
	for (int i = 0;i < facetN;i++)
	{
		delta[i] = 1.0; //defaulted, the whole PO region is lit
	}
	Triangles = mesh->getTriangles();
	Edges = mesh->getEdges();
	Vertexes = mesh->getVertexes();
	k = planewave.wavenumber;
	eta = planewave.eta;
	mu = planewave.mu;
	epsilon = planewave.epsilon;
	omega = planewave.omega;
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
			patch2 = mesh->getTriangles()[j];
			n2 = patch2.Normal;
			rOV = Vertexes[patch2.Node1];
			alpha = (rOM - rOV).dot(n2) / kDir.dot(n2);
			if (alpha > 1e-10)
			{
				rOP = rOM - alpha*kDir;
				if (PointinTriangle(Vertexes[patch2.Node1], Vertexes[patch2.Node2], Vertexes[patch2.Node3], rOP))
				{
					delta[i] = 0.0;
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
	for (int i = 0;i < facetN;i++)
	{
		if(delta[i]==1.0)
			LitFileName << Triangles[i].Node1 << " " << Triangles[i].Node2 << " " << Triangles[i].Node3 << endl;
	}
	LitFileName.close();
}

void PO::LitPatchInfo()
{
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
