#pragma once
#include<string>
#include<vector>
#include <fstream>
#include <iostream>
#include <math.h>
#define EIGEN_USE_MKL_ALL
#include <Eigen\Dense>
using namespace Eigen;
using namespace std;

typedef struct TriangleStruct
{
	//No of the triangle in the coordinate list
	int Node1;
	int Node2;
	int Node3;
	//normal direct to the out space
	Vector3d Normal;
	//area of the triangle
	double Area;
	Vector3d Center;
}Triangle;

typedef struct EdgeStruct
{
	int EdgeNode1;
	int EdgeNode2;
	int EdgeNode3; //left node for left triangle
	int EdgeNode4; //right node for right triangle
	int Triangle1;
	int Triangle2;
	double len;
}Edge;

class Mesh
{
private:
	vector<Triangle> Triangles;
	vector<Vector3d> Vertexes;
	int VertexCount, TriangleCount, EdgeCount;
	//file name of the mesh
	string FileName;
	vector<Edge*> Edges;
public:
	Mesh();
	bool WriteMesh();
	bool open(string FileName);
	void PrintMeshInfo();
	int getVertexCount() 
	{ 
		return VertexCount; 
	}
	int getTriangleCount() 
	{ 
		return TriangleCount; 
	}
	int getEdgeCount()
	{
		return EdgeCount;
	}
	vector<Triangle>& getTriangles() 
	{
		return Triangles;
	}
	vector<Edge*>& getEdges()
	{
		return Edges;
	}
	vector<Vector3d>& getVertexes()
	{
		return Vertexes;
	}
	string getFileName()
	{
		return FileName;
	}
	void EdgeInit();
	void InsertEdge(Edge *edge, int i);
	void WriteEdge();
	void PrintEdgeInfo();
	~Mesh();
};

