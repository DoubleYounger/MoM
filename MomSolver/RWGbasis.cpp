#include "stdafx.h"
#include "RWGbasis.h"

RWGbasis::RWGbasis(Mesh * mesh, int n)
{
	this->mesh = mesh;
	r[0] = mesh->getVertexes()[mesh->getEdges()[n]->EdgeNode1];
	r[1] = mesh->getVertexes()[mesh->getEdges()[n]->EdgeNode2];
	r[2] = mesh->getVertexes()[mesh->getEdges()[n]->EdgeNode3];
	r[3] = mesh->getVertexes()[mesh->getEdges()[n]->EdgeNode4];
	l= (r[0] - r[1]).norm();
	S_plus = mesh->getTriangles()[mesh->getEdges()[n]->Triangle1].Area;
	S_minus = mesh->getTriangles()[mesh->getEdges()[n]->Triangle2].Area;
	rho_c_plus = mesh->getTriangles()[mesh->getEdges()[n]->Triangle1].Center - r[2];
	rho_c_minus = r[3] - mesh->getTriangles()[mesh->getEdges()[n]->Triangle2].Center;
	r_center_plus = mesh->getTriangles()[mesh->getEdges()[n]->Triangle1].Center;
	r_center_minus = mesh->getTriangles()[mesh->getEdges()[n]->Triangle2].Center;
}

RWGbasis::~RWGbasis()
{
}
