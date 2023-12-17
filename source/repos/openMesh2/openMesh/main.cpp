

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <OpenMesh/Core/Utils/Predicates.hh>
#include <OpenMesh/Core/Utils/Property.hh>


bool is_divisible_by_3(OpenMesh::FaceHandle vh) { return vh.idx() % 3 == 0; }
using namespace OpenMesh::Predicates;
using namespace std;

#include <GL/glut.h>

typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;


GLfloat MyVertices[8][3] = {
	{ -0.25, -0.25, 0.25 },
	{ -0.25, 0.25, 0.25   },
	{ 0.25, 0.25, 0.25    },
	{ 0.25, -0.25, 0.25   },
	{ -0.25, -0.25, -0.25 },
	{ -0.25, 0.25, -0.25   },
	{ 0.25, 0.25, -0.25    },
	{ 0.25, -0.25, -0.25    }
};
GLfloat MyColors[8][3] = {
	{ 0.2, 0.2, 0.2 },
	{ 1.0, 0.0, 0.0 },
	{ 1.0, 1.0, 0.0 },
	{ 0.0, 1.0, 0.0 },
	{ 0.0, 0.0, 1.0 },
	{ 1.0, 0.0, 1.0 },
	{ 1.0, 1.0, 1.0 },
	{ 0.0, 1.0, 1.0 }
};
unsigned short MyVertexList[36] = {
	0,3,2,
	2,1,0,
	2,5,1,
	2,6,5,
	2,3,6,
	3,7,6,
	0,1,5,
	0,5,4,
	3,0,4,
	7,3,4,
	4,5,7,
	5,6,7
};
class float3 {
public:
	float x, y, z;
	float3(int a, int b, int c) { x = a; y = b; z = c; }
	float3(float a, float b, float c) { x = a; y = b; z = c; }
	float3(double a, double b, double c) { x = a; y = b; z = c; }

	float3 operator-(const float3& a) const {
		return { x - a.x, y - a.y, z - a.z };
	}
	float3 operator+(const float3& a) const {
		return { x + a.x, y + a.y, z + a.z };
	}
	float3 operator/(const float3& a)const {
		return { x / a.x, y / a.y, z / a.z };
	}
	float3 operator/(const float& a) const {
		return { x / a, y / a, z / a };
	}
	float3 operator*(const float3& a) const {
		return { x * a.x, y * a.y, z * a.z };
	}
	float3 operator*(const float& a)const {
		return { x * a, y * a, z * a };
	}
	bool operator>(const float& a) const {
		if (x > a && y > a && z > a)return true;
		else return false;
	}
	bool operator<(const float& a) const {
		if (x < a && y < a && z < a)return true;
		else return false;
	}
	float3 cross(const float3& a) const {
		return { y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x }; //T.T 
	}
	void normalize() {
		float slen = x * x + y * y + z * z; // ���� ����
		if (slen == 0)
			return;
		float ret = 1.0 / sqrtf(slen);
		x = x * ret;
		y = y * ret;
		z = z * ret;
	}
	float vectorlength() {
		return sqrt(x * x + y * y + z * z);
	}

	float dotProduct(const float3& dir) const {
		return{ x * dir.x + y * dir.y + z * dir.z };
	}
};

/*
int checkOppositeEdge(int iEdge) {
		if (getOppositeEdgeId(iEdge) == -1) return false; // oppositeEdge�� ���� �����
		else return getOppositeEdgeId(iEdge); // �ִ� ���

	}
	int getOppositeEdgeId(int iEdge) {
		if (edges[iEdge].deleted)
			return -1; // ������

		int vert1 = edges[iEdge].v0;
		int vert2 = edges[iEdge].v1;

		//for (auto it = EdgeArrForVertex[v1].begin(); it != EdgeArrForVertex[v1].end(); it++) {
		//	if (edges[*it].first == v2 && edges[*it].second == v1)
		//		return *it;
		//}
		for (int i = 0; i < EdgeArrForVertex[vert1].size(); i++) {
			int edge = EdgeArrForVertex[vert1][i];

			if (edges[edge].v0 == vert2 && edges[edge].v1 == vert1 && edges[edge].deleted == 0)
				return edge;
		}
		return -1; // ��ã����
	}
	void getOneRingNeighbor(int vindex, vector<int>& arr) { // arr�� ä���� ��
		arr.clear();
		int edgeIndex = 0;


		// vindex���� ������ edge�� ã���� �̿��� ��� ã�� �� ����
		for (int i = 0; i < EdgeArrForVertex[vindex].size(); i++) {
			edgeIndex = EdgeArrForVertex[vindex][i];
			if (edges[edgeIndex].v0 != vindex && !edges[edgeIndex].deleted)
				arr.push_back(edges[edgeIndex].v0);
			//if (edges[edgeIndex].second != vindex)
			//	arr.push_back(edges[edgeIndex].second);

		}
	}
	void getOneRingTriNeighborFromVertex(int vertNum, vector<int>& arr) {
		for (int i = 0; i < EdgeArrForVertex[vertNum].size(); i++) {
			int newE = EdgeArrForVertex[vertNum][i];
			if (edges[newE].v0 == vertNum && !edges[newE].deleted)
				arr.push_back(edges[newE].tri);
		}
	}
	float TriArea(const float3& point0, const float3& point1, const float3& point2) {
		float3 e1 = point1 - point0; // point1.operator-(point0)
		float3 e2 = point2 - point0;
		float3 n = e1.cross(e2);
		float area = n.vectorlength() * 0.5; //n ���̿� 0.5
		return area;
	}
	float3 getVertexNormal(int vertNum) {
		vector<int> arr; // �ﰢ���� ��ȣ�� ������ �迭
		vector<int> Earr;
		getOneRingTriNeighborFromVertex(vertNum, arr);
		// �� �ﰢ���� ���� normal, �׵��� ���, return
		float3 sum{ 0,0,0 };
		float total = 0.0;
		for (int i = 0; i < arr.size(); i++) {
			// arr[i] = �ﰢ�� ��ȣ
			int vertNum0 = vlistArr[(arr[i] * 3) + 0];
			float3 point0 = vertArr[vertNum0];
			int vertNum1 = vlistArr[(arr[i] * 3) + 1];
			float3 point1 = vertArr[vertNum1];				  //getOneRingTriNeighbor(arr[i], Earr);
			int vertNum2 = vlistArr[(arr[i] * 3) + 2];
			float3 point2 = vertArr[vertNum2];

			//float area = TriArea(point0, point1, point2);
			float3 e1 = point1 - point0;
			float3 e2 = point2 - point0;
			float3 n = e1.cross(e2);
			float area = n.vectorlength() * 0.5; //n ���̿� 0.5

			n.normalize();
			sum = sum + (n * area); // �ﰢ�� ���̿� ����ϴ� ����������� ��ġ��.
			total = total + area;
		}

		//float3 ret{ sum.x / arr.size(), sum.y / arr.size(), sum.z / arr.size() };
		if (total > 0.0) {
			sum = sum / total;
			sum.normalize();
		}
		//ret.normalize();
		return sum;
	}

	// TriNum: �ﰢ�� ��ȣ
	void getOneRingTriNeighbor(int TriNum, vector<int>& arr) {
		if (edges[TriNum * 3].deleted)
			return;
		// trinum���� edge��ȣ�� ����
		//0, 1, 2 -> 0
		//	3, 4, 5->1
			// edge->tri
			// tri->���� , 3��
			// ����->edge, 3*... = 20
			// for... 20 tri �Ҽ�����?
		//int e0 = TriNum * 3; // �ﰢ���� �����ϴ� edge��ȣ
		//int e1 = TriNum * 3 + 1;
		//int e2 = TriNum * 3 + 2;

		for (int i = 0; i < 3; i++) {
			int e0 = TriNum * 3 + i;
			int oe0 = getOppositeEdgeId(e0);
			if (oe0 != -1)
				arr.push_back(edges[oe0].tri);
		}
	}
	void Draw() { // ��� �������� �̿��ؼ� �׸� �׸�
		//;
		glVertexPointer(3, GL_FLOAT, 0, vertArr.data()); // &vertArr[0]);
		for (GLint i = 0; i < 12; i++) {
			if (!edges[3 * i].deleted) {
				glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_SHORT, &vlistArr[3 * i]);
			}
		}
	}
	// �ﰢ���� ���� ������ �������� ����� �˷���
	void DeleteTriangle(int triNum) {
		edges[triNum * 3].deleted = edges[triNum * 3 + 1].deleted = edges[triNum * 3 + 2].deleted = 1;

	}
	// �޽ÿ��� �Ƿ翧(2d)�� �ش��ϴ� edge��ȣ���� �����Ѵ�.
	void GetSilhouetteEdges2D(vector<int>& arr) {
		for (int i = 0; i < edges.size(); i++) {
			if (edges[i].deleted == 0 && checkOppositeEdge(i) == -1) //�ٵ� oppositeEdge ���� �� delete Ȯ���ؼ� �� ���Ǹ� �־ �� ��
				arr.push_back(i);

		}

	}



	void GetSilhouetteEdges3D(vector<int>& arr, const float3& dir) {
		for (int i = 0; i < edges.size(); i++) {
			float3 compareVector{ 0,0,0 };
			int oppositeEdge = checkOppositeEdge(i);
			if (oppositeEdge != -1){ // �ݴ� ������ ������ ���
				if (dir.dotProduct(getVertexNormal(edges[oppositeEdge].tri))* dir.dotProduct(getVertexNormal(edges[i].tri)) <0)
					arr.push_back(i);
			}
		}

	}


*/




//void printOppositeEdges(const MyMesh& mesh) {
//	for (MyMesh::EdgeIter edge_it = mesh.edges_begin(); edge_it != mesh.edges_end(); ++edge_it) {
//		MyMesh::HalfedgeHandle heh0 = mesh.halfedge_handle(*edge_it, 0);
//		MyMesh::HalfedgeHandle oppositeHeh0 = mesh.opposite_halfedge_handle(heh0);
//		MyMesh::EdgeHandle oppositeEdgeHandle0 = mesh.edge_handle(oppositeHeh0);
//
//		cout << "Edge " << (*edge_it).idx() << " -> Opposite Edge: ";
//		if (oppositeEdgeHandle0.is_valid())
//			cout << oppositeEdgeHandle0.idx() << " ";
//		cout << endl;
//	}
//}

float3 dir(1.0f, 0.0f, 0.0f);

struct EdgeInfo {
	int idx; // ���� �ε���
	int oppositeEdgeIdx; // �ݴ� ���� �ε��� (���� ��� -1)
};

// ������ �ݴ� ������ Ȯ���ϴ� �Լ�
EdgeInfo checkOppositeEdge(const MyMesh& mesh, int edgeIndex) {
	const MyMesh::EdgeHandle edgeHandle = mesh.edge_handle(edgeIndex);

	MyMesh::HalfedgeHandle heh0 = mesh.halfedge_handle(edgeHandle, 0);
	MyMesh::HalfedgeHandle heh1 = mesh.halfedge_handle(edgeHandle, 1);

	// �ݴ� ������ �ε����� ã��
	int oppositeEdgeIdx = -1;
	if (heh0.is_valid()) {
		oppositeEdgeIdx = mesh.edge_handle(heh0).idx();
	}
	else if (heh1.is_valid()) {
		oppositeEdgeIdx = mesh.edge_handle(heh1).idx();
	}

	return { edgeHandle.idx(), oppositeEdgeIdx };
}

vector<int> GetSilhouetteEdges3D(const MyMesh& mesh, const float3 dir) {
	std::vector<int> arr;

	for (MyMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
		const MyMesh::EdgeHandle edgeHandle = *e_it;

		EdgeInfo oppositeEdge = checkOppositeEdge(mesh, edgeHandle.idx());
		if (oppositeEdge.oppositeEdgeIdx != -1) {
			// Assuming you have functions getVertexNormal and dotProduct defined

			// if (dir.dotProduct(getVertexNormal(edges[oppositeEdge].tri))* dir.dotProduct(getVertexNormal(edges[i].tri)) <0)

			const float3 normal1 = getVertexNormal(mesh.property(points_, mesh.from_vertex_handle(edgeHandle)));
			const float3 normal2 = getVertexNormal(mesh.property(points_, mesh.from_vertex_handle(mesh.edge_handle(oppositeEdge))));

			if (dir.dotProduct(normal1) * dir.dotProduct(normal2) < 0)
				arr.push_back(edgeHandle.idx());
		}
	}

	return arr;
}

void printOppositeEdges(const MyMesh& mesh) {
	// iter -> handle.... *iter==handle
	for (MyMesh::HalfedgeIter edge_it = mesh.halfedges_begin(); edge_it != mesh.halfedges_end(); ++edge_it) {

		MyMesh::HalfedgeHandle oppositeHeh0 = mesh.opposite_halfedge_handle(*edge_it);
		cout << "Edge " << (*edge_it).idx() << " -> Opposite Edge: ";
		if (oppositeHeh0.is_valid())
			cout << oppositeHeh0.idx() << " ";
		cout << endl;
	}

}

void getAllNeighborVertices( MyMesh& mesh, MyMesh::VertexHandle *vhandle)
{
	//vector<MyMesh::VertexHandle> neighbors;
	cout << "= getNeighborVertices = \n" ;
	for(int i=0; i<8; i++){
		cout << i << " :";
		for (MyMesh::VertexVertexIter vv_it = mesh.vv_begin(vhandle[i]); vv_it != mesh.vv_end(vhandle[i]); ++vv_it)
		{
			MyMesh::VertexHandle neighbor_vertex = *vv_it;
			cout << neighbor_vertex << ",";
			//neighbors.push_back(neighbor_vertex);
		}
		cout << "\n";
	}
	
}

void getNeighborVertices(MyMesh& mesh, MyMesh::VertexHandle target_vertex)
{
	//vector<MyMesh::VertexHandle> neighbors;
	cout << "= getNeighborVertices = \n";
	for (MyMesh::VertexVertexIter vv_it = mesh.vv_begin(target_vertex); vv_it != mesh.vv_end(target_vertex); ++vv_it)
	{
		cout << target_vertex << " :" << endl;
		MyMesh::VertexHandle neighbor_vertex = *vv_it;
		cout << neighbor_vertex << ",";
		//neighbors.push_back(neighbor_vertex);
	}

}

void OppositeEdges(const MyMesh& mesh) {
	for (MyMesh::HalfedgeIter edge_it = mesh.halfedges_begin(); edge_it != mesh.halfedges_end(); edge_it++) {
		MyMesh::HalfedgeHandle oppositeHandle = mesh.opposite_halfedge_handle(*edge_it);
		cout << "edge" << (*edge_it).idx() << "-> opposite ";
		if (oppositeHandle.is_valid()) cout << oppositeHandle.idx() << "";
		cout << endl;
	}

	//for (MyMesh::HalfedgeIter edge_it = mesh.halfedges_begin(); edge_it != mesh.halfedges_end(); ++edge_it) {
	//
	//	MyMesh::HalfedgeHandle oppositeHeh0 = mesh.opposite_halfedge_handle(*edge_it);
	//	cout << "Edge " << (*edge_it).idx() << " -> Opposite Edge: ";
	//	if (oppositeHeh0.is_valid())
	//		cout << oppositeHeh0.idx() << " ";
	//	cout << endl;
	//}

}

//void findOneRingNeighbors( MyMesh& mesh, MyMesh::VertexHandle vh, vector<int>& neighbors) {
//	neighbors.clear();
//	for (MyMesh::VertexVertexIter vv_it = mesh.vv_begin(vh); vv_it != mesh.vv_end(vh); ++vv_it) {
//		neighbors.push_back((*vv_it).idx());
//	}
//}

// �ٵ� ���ʿ� one ring neighbor�� ���� �ʿ䰡 ����..
// vertex_vertex_iter �� ��ü: �־��� ���� �������� ���� ���� �迭�� �ִµ�, �� iter
// �־��� ������ ����? (1) VertexIter �̴�. (2) VertexHandle �̴�.
void findOneRingNeighbors23(MyMesh& mesh, MyMesh::VertexHandle vh, vector<int> neighbors) {
	neighbors.clear();
	// 1 VertexIter ���� vv_iter ���ϸ� ��
	//MyMesh::VertexIter v_it = mesh.vertices_begin();
	//MyMesh::VertexVertexIter    vv_it;
	//for (vv_it = mesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it) {
	//
	//}

	// 2 VertexHandle���� vv_begin�� ���ϸ� ��
	for (MyMesh::VertexVertexIter vv_it = mesh.vv_begin(vh); vv_it != mesh.vv_end(vh); ++vv_it) {
		//neighbors.push_back((*vv_it));
	}
}

// geodesic distance
/*
float getDistance(MyMesh& mesh, MyMesh::VertexHandle v0, MyMesh::VertexHandle v1) {

	// top-down : ��ü ������ ������ ���鼭 �κ��� �ϼ� - ������ ���� ����
	// bottom - up : �˰� �ִ� �ɰ����� �����ؼ� ��ü�� ���� - ���� ������ �� �ȸ³�...?
	//float arr[8] = { 0, };

	//float arr[8] = { 0, }; //
	const float INF = 1000000;
	float distance[8] = { INF,INF,INF,INF,INF,INF,INF,INF };// ���� �������� ���� ������ ������ �ִܰŸ� ����ġ
	bool visit[8] = { false };//
	MyMesh::VertexHandle nowVertex = v0;
	MyMesh::VertexHandle destVertex = v1;
	//int <vector>NeighborsArr;

	//MyMesh::VertexHandle vh(v0);
	vector<MyMesh::VertexHandle> NeighborsArr;
	distance[nowVertex.idx()] = 0;
	while (nowVertex != v1) {
		visit[nowVertex.idx()] = true;
		findOneRingNeighbors(mesh, nowVertex, NeighborsArr);
		MyMesh::Point p0 = mesh.point(nowVertex);
		for (int i = 0; i < NeighborsArr.size(); i++) {
			//mesh.vertex(NeighborsArr[i])
			MyMesh::VertexHandle vhandle = NeighborsArr[i];
			MyMesh::Point p1 = mesh.point(vhandle);

			float distance1 = (p1 - p0).length();
			distance[vhandle.idx()] = distance1; // ���� ������ ����. ���� �߰ߵ� ������� ������Ʈ
		}
		float minV = 100000; // distance[0];
		int nxtid = 0;
		for (it �� ��� ������ ���� handle) {
			id = handle.id;
			if (visited[id])
				continue;
			distance[id]; �̿� �ּڰ�

		}
		�ּҶ�� �˷��� nexthandle = v1 �̸� ����. return;
		//for (int i = 1; i < 8; i++) {
		//	if (minV < distance[i]) nxtid = i;
		//}
	}

	//�迭 �ʱ�ȭ; �Ÿ��� �������;
	//�湮���� �ʱ�ȭ; �湮������ bool �迭.false�ʱ�ȭ;
	//���� ����V = v0���� ����;
	//while
	//	���� ������ �湮ǥ��;
	//	���� �������� ���� ���� list Ȯ��; findOneRingNeighbors
	//	�� list ����� ���� �Ÿ����; mesh ������ǥ
	//	�Ÿ� ����� ���� �Ÿ� �迭�� ������Ʈ;
	//
	//	���� ������ ã�Ƽ� ���� ���� V = �迭���� �ּҰŸ��� ������.�� �湮������ �װ� ����.;





	// ����...
	// ���� �̰��ϰ�
	// ����ϰ�
	// �ݺ��ϸ鼭
	neighArr
		findOneRingNeighbors(mesh, [v0], NeighborsArr)
}

*/

void addProperty(MyMesh& mesh, MyMesh::VertexIter v_it, MyMesh::VertexVertexIter vv_it, MyMesh::VertexIter v_end, OpenMesh::VPropHandleT<MyMesh::Point> cogs, MyMesh::Scalar valence) {
	
	// �޽��� ��� ������ �ݺ��ϸ鼭, �� ������ ���� �߽��� �Ӽ��� �ʱ�ȭ�ϰ� �̿� �������� ��ġ�� �����Ͽ� �߽����� ���
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) {
		mesh.property(cogs, *v_it).vectorize(0.0f); //�Ӽ� �κ� �ʱ�ȭ
		valence = 0.0;
		for (vv_it = mesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it) {
			mesh.property(cogs, *v_it) += mesh.point(*vv_it); //�߽��� �����ϰ�
			++valence;
		}
		mesh.property(cogs, *v_it) /= valence;
	}
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) {
		if (!mesh.is_boundary(*v_it)) //��� �κ� �����ϰ� ���鼭
			mesh.set_point(*v_it, mesh.property(cogs, *v_it));
	}
	//�Ӽ� ��� �κ�
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) {
		MyMesh::Point cog = mesh.property(cogs, *v_it);
		cout << "Vertex " << *v_it << " COG: (" << cog[0] << ", " << cog[1] << ", " << cog[2] << ")" << std::endl;
	}
}

//chap4

template <class Mesh> class SmootherT
{
public:
	typedef typename Mesh::Point  cog_t;
	typedef OpenMesh::VPropHandleT< cog_t > Property_cog;

public:
	explicit SmootherT(Mesh& _mesh) : mesh_(_mesh)
	{
		mesh_.add_property(cog_);
	}

	~SmootherT()
	{
		mesh_.remove_property(cog_);
	}

	void smooth(unsigned int _iterations)
	{
		for (unsigned int i = 0; i < _iterations; ++i)
		{
			for_each(mesh_.vertices_begin(), mesh_.vertices_end(), ComputeCOG(mesh_, cog_));
			for_each(mesh_.vertices_begin(), mesh_.vertices_end(), SetCOG(mesh_, cog_));
		}
	}

private:
	class ComputeCOG
	{
	public:
		ComputeCOG(Mesh& _mesh, Property_cog& _cog) : mesh_(_mesh), cog_(_cog)
		{}

		void operator()(const typename Mesh::VertexHandle& _vh)
		{
			typename Mesh::VertexVertexIter vv_it;
			typename Mesh::Scalar valence(0.0);

			cog_[_vh] = typename Mesh::Point(0.0, 0.0, 0.0);
			for (vv_it = mesh_.vv_iter(_vh); vv_it.is_valid(); ++vv_it)
			{
				cog_[_vh] += mesh_.point(*vv_it);
				++valence;
			}
			cog_[_vh] /= valence;
		}

	private:
		Mesh& mesh_;
		Property_cog& cog_;
	};

	class SetCOG
	{
	public:
		SetCOG(Mesh& _mesh, Property_cog& _cog)	: mesh_(_mesh), cog_(_cog) {}

		void operator()(const typename Mesh::VertexHandle& _vh)
		{
			if (!mesh_.is_boundary(_vh))
				mesh_.set_point(_vh, cog_[_vh]);
		}

	private:
		Mesh& mesh_;
		Property_cog& cog_;
	};

private:
	Mesh& mesh_;
	Property_cog cog_;
};

// chap 5
void opposite_vertices_arr(MyMesh& mesh) {
	
	vector<OpenMesh::VertexHandle> opposite_vertices;
	
	for (auto vh : mesh.vertices())
	{
		// iterate over all outgoing halfedges
		for (auto heh : vh.outgoing_halfedges())
		{
			// navigate to the opposite vertex and store it in the vector
			opposite_vertices.push_back(heh.next().opp().next().to());
		}
	}

	//for (auto h_it = opposite_vertices.begin(); h_it != opposite_vertices.end(); ++h_it) {
	//	cout << (*h_it)<< endl;
	//}
	
	//�־��� �ڵ�� �޽��� �� ������ ���� �ݺ��ϸ鼭, 
	// �� ������ ���� �ݴ� ������ �������� ã�� ���� �ڵ�
	// iterate over vertices of the mesh
	for (auto vh : mesh.vertices())
	{

		// create lambda that returns opposite vertex
		auto opposite_vertex = [](OpenMesh::SmartHalfedgeHandle heh) { 
			cout << "heh : " << heh << endl;
			cout << "heh.next() : " << heh.next() << endl;
			cout << "heh.next().opp() : " << heh.next().opp() << endl;
			cout << "heh.next().opp().next() : " << heh.next().opp().next() << endl;
			cout << "heh.next().opp().next().to() : " << heh.next().opp().next().to() << endl;


			return heh.next().opp().next().to(); 
		};
		// create vector containing all opposite vertices
		auto opposite_vertices = vh.outgoing_halfedges().to_vector(
			opposite_vertex);
	}




}

void calc_laplace(MyMesh& mesh){
	// Add a vertex property storing the laplace vector
	auto laplace = OpenMesh::VProp<MyMesh::Point>(mesh);

	// Add a vertex property storing the laplace of the laplace
	auto bi_laplace = OpenMesh::VProp<MyMesh::Point>(mesh);

	// Get a propertymanager of the points property of the mesh to use as functor
	auto points = OpenMesh::getPointsProperty(mesh);


	//for (int i = 0; i < iterations; ++i) { // ���� iterations ���� �𸣰���...
	for (auto i : mesh.vertices()){
		// Iterate over all vertices to compute laplace vector
		for (const auto& vh : mesh.vertices())
			laplace(vh) = vh.vertices().avg(points) - points(vh);

		// Iterate over all vertices to compute the laplace vector of the laplace vectors
		for (const auto& vh : mesh.vertices())
			bi_laplace(vh) = (vh.vertices().avg(laplace) - laplace(vh));

		// update points by substracting the bi-laplacian damped by a factor of 0.5
		for (const auto& vh : mesh.vertices())
			points(vh) += -0.5 * bi_laplace(vh);
	}
} // The laplace and update properties are removed from the mesh at the end of this scope.

void calc_barycenters(MyMesh mesh, MyMesh::VertexIter v_it, MyMesh::VertexIter v_end, MyMesh::VertexVertexIter vv_it) {
	
	MyMesh::Point         cog;
	MyMesh::Scalar        valence;
	vector <MyMesh::Point> cogs;
	vector <MyMesh::Point>::iterator cog_it;

	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) {
		cog[0] = cog[1] = cog[2] = valence = 0.0;
		for (vv_it = mesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it) {
			cog += mesh.point(*vv_it);
			++valence;
		}
		cogs.push_back(cog / valence);
	}
	for (v_it = mesh.vertices_begin(), cog_it = cogs.begin(); v_it != v_end; ++v_it, ++cog_it) {
		if (!mesh.is_boundary(*v_it)) mesh.set_point(*v_it, *cog_it);
	}

}

void OneRingNeighbor(MyMesh& mesh, MyMesh::VertexIter v_it, MyMesh::VertexIter v_end, MyMesh::VertexVertexIter vv_it, vector<int>& neighbors) {
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) {
		for (vv_it = mesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it) {
			neighbors.push_back(vv_it->idx());
		}
	}
}

void NeighborsOfVertex(MyMesh& mesh, MyMesh::VertexHandle vh) {
	MyMesh::VertexVertexIter vv_it;
	std::cout << "Vertex " << vh.idx() << "�� �̿� ����: ";
	for (vv_it = mesh.vv_iter(vh); vv_it.is_valid(); ++vv_it) {
		std::cout << vv_it->idx() << " ";
	}
	std::cout << std::endl;
}

void MyDisplay() {
	glClear(GL_COLOR_BUFFER_BIT);
	glFrontFace(GL_CCW); // ���̰˻�(8��)�� ���� �ʰ�, front face���
	//glEnable(GL_CULL_FACE);
	glEnableClientState(GL_COLOR_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);
	glColorPointer(3, GL_FLOAT, 0, MyColors);
	//glVertexPointer(3, GL_FLOAT, 0, MyVertices);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRotatef(30.0, 1.0, 1.0, 1.0);

	//mesh.Draw();

	glFlush();
}

//int main() {
//	MyMesh mesh;
//	MyMesh::VertexHandle vhandle[8];
//	MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
//	MyMesh::VertexVertexIter vv_it;
//	vector<MyMesh::VertexHandle> face_vhandles;
//
//	for (int i = 0; i < 8; i++) {
//		vhandle[i] = mesh.add_vertex(MyMesh::Point(MyVertices[i][0], MyVertices[i][1], MyVertices[i][2]));
//	}
//	for (int i = 0; i < 10; i++) {
//		face_vhandles.clear();
//		face_vhandles.push_back(vhandle[MyVertexList[i*3]]);
//		face_vhandles.push_back(vhandle[MyVertexList[i*3+1]]);
//		face_vhandles.push_back(vhandle[MyVertexList[i*3+2]]);
//		mesh.add_face(face_vhandles);
//	}
//	
//	vector<int>neighbors;
//	neighbors.clear();
//
//	OneRingNeighbor(mesh, v_it, v_end, vv_it, neighbors);
//	for (size_t i = 0; i < neighbors.size(); ++i) {
//		cout << "Vertex " << i << "�� �̿�: " << neighbors[i] << endl;
//	}
//	//calc_barycenters(mesh, v_it, v_end, vv_it);
//
//}

	

int main()
{
	MyMesh mesh;
	// generate vertices

	MyMesh::VertexHandle vhandle[8];
	for (int i = 0; i < 8; i++) {
		vhandle[i] = mesh.add_vertex(MyMesh::Point(MyVertices[i][0], MyVertices[i][1], MyVertices[i][2]));
	}


	std::vector<MyMesh::VertexHandle> face_vhandles;

	for (int i = 0; i < 10; i++) {
		face_vhandles.clear();
		face_vhandles.push_back(vhandle[MyVertexList[i * 3 + 0]]);
		face_vhandles.push_back(vhandle[MyVertexList[i * 3 + 1]]);
		face_vhandles.push_back(vhandle[MyVertexList[i * 3 + 2]]);
		mesh.add_face(face_vhandles);
	}

	vector<int> neighbors;
	MyMesh::VertexHandle vertexHandle = vhandle[4];
	//findOneRingNeighbors(mesh, vertexHandle, neighbors);

	//getOppositeEdgeId(mesh, 4);
	
	cout << " == vertex ==" << endl;

	for (int i = 0; i < 8; i++) {
		cout << i << " : " << mesh.point(vhandle[i]) << endl;

	}

	// OppositeEdges(mesh);
	getAllNeighborVertices(mesh, vhandle);
//	getNeighborVertices(mesh, vhandle[4]);

	// chap3
	OpenMesh::VPropHandleT<MyMesh::Point> cogs;
	mesh.add_property(cogs);

	MyMesh::VertexIter          v_it, v_end(mesh.vertices_end());
	MyMesh::VertexVertexIter    vv_it;
	MyMesh::Point               cog;
	MyMesh::Scalar              valence;

	valence = 0.0;
	
	neighbors.clear();


	cout << "\n\n == OneRingNeighbor ==" << endl;
	OneRingNeighbor(mesh, v_it, v_end, vv_it, neighbors);
	for (size_t i = 0; i < neighbors.size(); ++i) {
		cout << "Vertex " << i << "�� �̿�: " << neighbors[i] << endl;
	}

	MyMesh::VertexHandle targetVertex = vhandle[4];
	// ������ ������ �̿� ���� ���
	//NeighborsOfVertex(mesh, targetVertex);
	//targetVertex = vhandle[0];
	//NeighborsOfVertex(mesh, targetVertex);

	//for (int i = 0; i < 8; i++) {
	//	targetVertex = vhandle[i];
	//	NeighborsOfVertex(mesh, targetVertex);
	//}
	


	//addProperty(mesh, v_it, vv_it, v_end, cogs, valence);
	//opposite_vertices_arr(mesh);

	//printOppositeEdges(mesh);

	//printOppositeEdges(mesh);
/*
  //Chap Filtering ranges with predicates
  // Count boundary vertices
	cout << "Mesh contains " << mesh.vertices().count_if(Boundary()) << " boundary vertices";
	// Selected inner vertices
	cout << "These are the selected inner vertices: " << std::endl;
	for (auto vh : mesh.vertices().filtered(!Boundary() && Selected()))
		cout << vh.idx() << ", ";
	cout << std::endl;
	

	// Faces whose id is divisible by 3
	auto vec = mesh.faces().filtered(is_divisible_by_3).to_vector();
	cout << "There are " << vec.size() << " faces whose id is divisible by 3" << std::endl;

	// Faces which are tagged or whose id is not divisible by 3
	auto vec2 = mesh.faces().filtered(Tagged() || !make_predicate(is_divisible_by_3)).to_vector();
	cout << "There are " << vec2.size() << " faces which are tagged or whose id is not divisible by 3" << std::endl;

	// Edges that are longer than 10 or shorter than 2
	OpenMesh::EProp<bool> longer_than_10(mesh);
	for (auto eh : mesh.edges())
		longer_than_10[eh] = mesh.calc_edge_length(eh) > 10;
	cout << "There are " <<
		mesh.edges().count_if(make_predicate(longer_than_10) || make_predicate([&](OpenMesh::EdgeHandle eh) { return mesh.calc_edge_length(eh) < 2; })) <<
		" edges which are shorter than 2 or longer than 10" << std::endl;
	*/
	
	try
	{
		if (!OpenMesh::IO::write_mesh(mesh, "output.off"))
		{
			cerr << "Cannot write mesh to file 'output.off'" << std::endl;
			return 1;
		}
	}
	catch (std::exception& x)
	{
		cerr << x.what() << endl;
		return 1;
	}

	return 0;
}


