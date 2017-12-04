#include "stdafx.h"
#include "Mesh.h"

double getArea(Vector3d v1, Vector3d v2, Vector3d v3);
Mesh::Mesh()
{
	VertexCount = 0;
	TriangleCount = 0;
	EdgeCount = 0;
}


Mesh::~Mesh()
{
}

bool Mesh::open(string FileName_)
{
	FileName = FileName_;
	string str;
	char *chs = new char[200];
	fstream MeshFile;
	MeshFile.open(FileName, ios::in);
	if (MeshFile.fail())
	{
		throw "File open failed";
	}
	VertexCount = 0;
	do
	{
		MeshFile >> str;
		if (str == "GRID*")
			VertexCount++;
	} while (str != "CTRIA3");
	MeshFile.seekg(0, ios::beg);
	for (size_t i = 0; i < 6; i++)
	{
		MeshFile.getline(chs, 200);
	}
	int temp;
	MeshFile >> str >> str >> TriangleCount;
	MeshFile.get();
	MeshFile.getline(chs, 200);
	MeshFile.getline(chs, 200);
	//read map file first
	if (!mapFileReader())
		cout << "map file read error!" << endl;

	Triangle Triangle;
	double x, y, z;
	Vector3d Vertex;
	for (int i = 0; i < VertexCount; i++)
	{
		MeshFile >> str >> temp >> x >> y >> temp >> str >> temp >> z;
		Vertex << x, y, z;
		Vertexes.push_back(Vertex);
	}
	for (int i = 0; i < TriangleCount; i++)
	{
		MeshFile >> str >> temp >> Triangle.patch >> Triangle.Node1 >> Triangle.Node2 >> Triangle.Node3;
		Triangle.RegionType = mapPatch[Triangle.patch];
		Triangle.Node1--;
		Triangle.Node2--;
		Triangle.Node3--;
		//the normal is pointing to the out space
		//test passed!!!!!!
		Triangle.Normal = (Vertexes[Triangle.Node1] - Vertexes[Triangle.Node2]).cross(Vertexes[Triangle.Node2] - Vertexes[Triangle.Node3]);
		Triangle.Normal.normalize();
		Triangle.Area = getArea(Vertexes[Triangle.Node1], Vertexes[Triangle.Node2], Vertexes[Triangle.Node3]);
		Triangle.Center = (Vertexes[Triangle.Node1] + Vertexes[Triangle.Node2] + Vertexes[Triangle.Node3]) / 3.0;
		Triangles.push_back(Triangle);
	}

	return false;
}

bool Mesh::mapFileReader()
//to determine the region type of the whole patches must read the .map file
{
	string mapFileName = FileName.substr(0, FileName.find(".")) + ".map";
	fstream mapFile(mapFileName, ios::in);
	if (mapFile.fail())
	{
		cout << "map File open failed... (please make sure if it doesn't matter!)" << endl;
	}
	int patch;
	string type;
	string patchName;
	EnumParser<Region> fieldTypeParser;
	while (!mapFile.eof())
	{
		mapFile >> patch >> patchName >> type;
		mapPatch[patch] = fieldTypeParser.ParseSomeEnum(type);
	}
	return true;
}

double getArea(Vector3d v1, Vector3d v2, Vector3d v3)
//°Ã[p(p-a)(p-b)(p-c) ]∆‰÷–p=1/2(a+b+c)
{
	double l[3], p;
	l[0] = (v1 - v2).norm();
	l[1] = (v1 - v3).norm();
	l[2] = (v2 - v3).norm();
	p = (l[0] + l[1] + l[2]) / 2.0;
	return sqrt(p*(p - l[0])*(p - l[1])*(p - l[2]));
}
bool Mesh::WriteMesh()
{
	fstream MeshFile;
	MeshFile.open(FileName.substr(0, FileName.find(".")) + ".msh", ios::out);
	if (MeshFile.fail())
	{
		throw "File open failed";
	}
	MeshFile << VertexCount << "\t" << TriangleCount << endl;
	for (int i=0;i<VertexCount;i++)
	{
		MeshFile << Vertexes[i][0] << " " << Vertexes[i][1] << " " << Vertexes[i][2] << endl;
	}
	for (int i = 0;i<TriangleCount;i++)
	{
		MeshFile << Triangles[i].Node1 << " " << Triangles[i].Node2 << " " << Triangles[i].Node3 << " " << Triangles[i].Normal[0] << " " << Triangles[i].Normal[1] << " " << Triangles[i].Normal[2] << endl;
	}
	MeshFile.close();
	return true;
}

void Mesh::PrintMeshInfo()
{
	cout << "----------------------------------------------------" << endl;
	cout << "Number of the vertices: " << VertexCount << endl;
	cout << "Number of the triangles: " << TriangleCount << endl;
	cout << "----------------------------------------------------" << endl;
	return;
}

void Mesh::EdgeInit()
{
	for (int i = 0;i < this->getTriangleCount();i++)
	{
		for (int j = 0;j < 3;j++)
		{
			Edge *edge1 = new Edge();
			edge1->EdgeNode1 = Triangles[i].Node1;
			edge1->EdgeNode2 = Triangles[i].Node2;
			InsertEdge(edge1, i);
			Edge *edge2 = new Edge();
			edge2->EdgeNode1 = Triangles[i].Node2;
			edge2->EdgeNode2 = Triangles[i].Node3;
			InsertEdge(edge2, i);
			Edge *edge3 = new Edge();
			edge3->EdgeNode1 = Triangles[i].Node1;
			edge3->EdgeNode2 = Triangles[i].Node3;
			InsertEdge(edge3, i);
		}
	}
	EdgeCount = (int)Edges.size();
}

class Comp
{
private:
	Edge *edge;
public:
	Comp(Edge *edge_) :edge(edge_) {};
	bool operator()(Edge* lhs)
	{
		return ((lhs->EdgeNode1 == edge->EdgeNode1) && (lhs->EdgeNode2 == edge->EdgeNode2)) || ((lhs->EdgeNode2 == edge->EdgeNode1) && (lhs->EdgeNode1 == edge->EdgeNode2));
	}
};

void Mesh::InsertEdge(Edge *edge, int i)
{
	vector<Edge*>::iterator pos;
	pos = find_if(Edges.begin(), Edges.end(), Comp(edge));
	if (pos == Edges.end())
	{
		for (int j = 0;j < this->getTriangleCount();j++)
		{
			if (i == j)
				continue;
			if ((edge->EdgeNode1 == Triangles[j].Node1 || edge->EdgeNode1 == Triangles[j].Node2 || edge->EdgeNode1 == Triangles[j].Node3) && (edge->EdgeNode2 == Triangles[j].Node1 || edge->EdgeNode2 == Triangles[j].Node2 || edge->EdgeNode2 == Triangles[j].Node3))
			{
				edge->Triangle1 = i < j ? i : j;
				edge->Triangle2 = i > j ? i : j;
				edge->EdgeNode3 = Triangles[edge->Triangle1].Node1 + Triangles[edge->Triangle1].Node2 + Triangles[edge->Triangle1].Node3 - edge->EdgeNode1 - edge->EdgeNode2;
				edge->EdgeNode4 = Triangles[edge->Triangle2].Node1 + Triangles[edge->Triangle2].Node2 + Triangles[edge->Triangle2].Node3 - edge->EdgeNode1 - edge->EdgeNode2;
				edge->len = (Vertexes[edge->EdgeNode1] - Vertexes[edge->EdgeNode2]).norm();
				if (Triangles[edge->Triangle1].RegionType == Triangles[edge->Triangle2].RegionType)
				{
					edge->RegionType = Triangles[edge->Triangle1].RegionType;
				}
				else
				{
					edge->RegionType = po;
				}
				Edges.push_back(edge);
				break;
			}
		}
	}
}

void Mesh::WriteEdge()
{
	fstream EdgeFile;
	EdgeFile.open(FileName.substr(0, FileName.find(".")) + ".edg", ios::out);
	if (EdgeFile.fail())
	{
		throw "File open failed";
	}
	EdgeFile << Edges.size() << endl;
	for (int i=0;i<Edges.size();i++)
	{
		EdgeFile << Edges[i]->EdgeNode1 << "\t" << Edges[i]->EdgeNode2 << "\t" << Edges[i]->Triangle1 << "\t" << Edges[i]->Triangle2 << endl;
	}
	EdgeFile.close();
}

void Mesh::PrintEdgeInfo()
{
	cout << "Number of the edges: " << EdgeCount << endl;
	cout << "----------------------------------------------------" << endl;
	return;
}

void Mesh::divideEdges()
{
	for (int i = 0;i < EdgeCount; i++)
	{
		if (Edges[i]->RegionType == po)
			POEdges.push_back(i);
		else if (Edges[i]->RegionType == mom)
			MOMEdges.push_back(i);
		else
			MOMEdges[i] = i;
	}
}
