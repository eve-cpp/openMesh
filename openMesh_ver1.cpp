#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
using namespace std;

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <GL/glut.h>

GLfloat MyVertices[8][3] = {
	{ -0.25, -0.25, 0.25 },
	{ -0.25, 0.25, 0.25 },
	{ 0.25, 0.25, 0.25 },
	{ 0.25, -0.25, 0.25 },
	{ -0.25, -0.25, -0.25 },
	{ -0.25, 0.25, -0.25 },
	{ 0.25, 0.25, -0.25 },
	{ 0.25, -0.25, -0.25 }
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
	0,4,5,
	3,0,4,
	7,3,4,
	4,5,7,
	5,6,7
};
class float3 {
public:
	float x, y, z;
	float3 operator-(const float3&a) {
		return { x - a.x, y - a.y, z - a.z };
	}
	float3 operator+(const float3& a) {
		return { x + a.x, y + a.y, z + a.z };
	}
	float3 cross(const float3& a) {
		return { y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x }; //T.T 
	}
	void normalize() {
		float slen = x * x + y * y + z * z; // 길이 제곱
		if (slen == 0)
			return;
		float ret = 1.0 / sqrtf(slen);
		x= x  * ret;
		y = y * ret;
		z = z * ret;
	}


};



class Mesh {
public:
	vector<float3> vertArr; // 정점 번호를 알면 정점 좌표가 있음
	vector<vector<int>> EdgeArrForVertex;  //EdgeArrForVertex; // 각 정점에 연결되는 edge들 번호를 모은 집합
	vector<short> vlistArr; // 삼각형 번호 -> 정점 번호가 나옴
	//vlistArr[3*1+0] vlistArr[3*1+1] vlistArr[3*1+2]
	//vector< pair<short, short >> edges; // 6,1

	class Edge {
	public:
		int v0, v1;
		int tri;
	};
	vector<Edge> edges; // edge번호를 알면 edge정보가 있음
	//vector< pair< pair<short, short>, int >> edges; // 6,1

	Mesh(float* Vert, int vnum, short *index, int tnum) { 
		EdgeArrForVertex.clear();

		for (int i = 0; i < vnum; i++) {
			//float3 vtx { Vert[i * 3], Vert[i * 3+1] , Vert[i * 3+2] }; // (&Vert[i * 3]);
			vertArr.push_back({ Vert[i * 3], Vert[i * 3 + 1] , Vert[i * 3 + 2] });
		}
		for (int i = 0; i < tnum * 3; i++) {
			vlistArr.push_back(index[i]);
		}
		EdgeArrForVertex.resize(tnum*3);
		for (int tid = 0; tid < tnum; tid++) { //12

			int edgeid = tid * 3; // 30번 edge
			// 삼각형을 구성하는 정점 번호들 (감겨있으므로, 첫점이 중복됨)
			// 삼각형 edge를 구성하는 각 정점들
			int vIds[4] = { index[tid * 3] , index[tid * 3 + 1], index[tid * 3 + 2], index[tid * 3] };
			for(int e=0; e<3; e++) { // 삼각형의 에지 3개
				int vid0 = vIds[e], vid1 = vIds[e + 1];
				edges.push_back({ vid0, vid1, tid}); // 정점번호 16, 72
				EdgeArrForVertex[vid0].push_back(edgeid+e);
				EdgeArrForVertex[vid1].push_back(edgeid+e);
			}
			//int edgeid = i * 3; // 30번 edge
			//edges.push_back({ vlist[i * 3], vlist[i * 3 + 1] }); // 정점번호 16, 72
			//EdgeArrForVertex[vlist[i * 3]].push_back(edgeid);
			//EdgeArrForVertex[vlist[i * 3 + 1]].push_back(edgeid);
			//edgeid++;
			//
			//edges.push_back(make_pair(vlist[i * 3 + 1], vlist[i * 3 + 2]));
			//EdgeArrForVertex[vlist[i * 3 + 1]].push_back(edgeid);
			//EdgeArrForVertex[vlist[i * 3 + 2]].push_back(edgeid);
			//edgeid++;
			//
			//edges.push_back(make_pair(vlist[i * 3 + 2], vlist[i * 3]));
			//EdgeArrForVertex[vlist[i * 3 + 2]].push_back(edgeid);
			//EdgeArrForVertex[vlist[i * 3]].push_back(edgeid);
		}
	}
	int getOppositeEdgeId(int iEdge) {
		int vert1 = edges[iEdge].v0;
		int vert2 = edges[iEdge].v1;

		//for (auto it = EdgeArrForVertex[v1].begin(); it != EdgeArrForVertex[v1].end(); it++) {
		//	if (edges[*it].first == v2 && edges[*it].second == v1)
		//		return *it;
		//}
		for (int i = 0; i < EdgeArrForVertex[vert1].size(); i++) {
			int edge =EdgeArrForVertex[vert1][i];
			if (edges[edge].v0 == vert2 && edges[edge].v1 == vert1)
				return edge;
		}
		return -1; // 못찾았음
	}
	void getOneRingNeighbor(int vindex, vector<int>& arr) { // arr에 채워야 함
		arr.clear();
		int edgeIndex = 0;

		for (int i = 0; i < EdgeArrForVertex[vindex].size(); i++) {
			edgeIndex = EdgeArrForVertex[vindex][i];
			if(edges[edgeIndex].v0 != vindex)
				arr.push_back(edges[edgeIndex].v0);
			//if (edges[edgeIndex].second != vindex)
			//	arr.push_back(edges[edgeIndex].second);

		}
	}
	void getOneRingTriNeighborFromVertex(int vertNum, vector<int>& arr) {
		for (int i = 0; i < EdgeArrForVertex[vertNum].size(); i++) {
			int newE = EdgeArrForVertex[vertNum][i];
			if(edges[newE].v0==vertNum)
				arr.push_back(edges[newE].tri);
		}
	}
	float3 getVertexNormal(int vertNum) {
		vector<int> arr; // 삼각형의 번호를 저장한 배열
		vector<int> Earr;
		getOneRingTriNeighborFromVertex(vertNum, arr);
		// 각 삼각형에 대해 normal, 그들을 평균, return
		float3 sum{ 0,0,0 };
		for (int i = 0; i < arr.size(); i++) {
			// arr[i] = 삼각형 번호
			int vertNum0 = vlistArr[(arr[i] * 3) + 0];
			float3 point0 = vertArr[vertNum0];
			int vertNum1 = vlistArr[(arr[i] * 3) + 1];
			float3 point1 = vertArr[vertNum1];				  //getOneRingTriNeighbor(arr[i], Earr);
			int vertNum2 = vlistArr[(arr[i] * 3) + 2];
			float3 point2 = vertArr[vertNum2];
			
			float3 e1 = point1 - point0;
			float3 e2 = point2 - point0;
			float3 n = e1.cross(e2);
			n.normalize();
			sum =sum + n; // 삼각형 넓이에 비례하는 가중평균으로 고치기.
		}
		float3 ret{ sum.x / arr.size(), sum.y / arr.size(), sum.z / arr.size() };
		ret.normalize();
		return ret;
	}

	// TriNum: 삼각형 번호
	void getOneRingTriNeighbor(int TriNum, vector<int>& arr) {
		// trinum에서 edge번호를 추출
		//0, 1, 2 -> 0
		//	3, 4, 5->1
			// edge->tri
			// tri->정점 , 3개
			// 정점->edge, 3*... = 20
			// for... 20 tri 소속인지?
		//int e0 = TriNum * 3; // 삼각형을 구성하는 edge번호
		//int e1 = TriNum * 3 + 1;
		//int e2 = TriNum * 3 + 2;

		for (int i = 0; i < 3; i++) {
			int e0 = TriNum * 3 + i;
			int oe0 = getOppositeEdgeId(e0);
			arr.push_back(edges[oe0].tri);
		}
	}
};
void MyDisplay() {
	glClear(GL_COLOR_BUFFER_BIT);
	glFrontFace(GL_CCW); // 깊이검사(8장)를 하지 않고, front face사용
	//glEnable(GL_CULL_FACE);
	glEnableClientState(GL_COLOR_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);
	glColorPointer(3, GL_FLOAT, 0, MyColors);
	glVertexPointer(3, GL_FLOAT, 0, MyVertices);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRotatef(30.0, 1.0, 1.0, 1.0);
	//mesh.Draw();
	//for (GLint i = 0; i < 12; i++)
	//	glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_SHORT, &MyVertexList[3 * i]);
	glFlush();
}

int main(int argc, char** argv) {
	Mesh mesh((float*)MyVertices, 8, (short*)MyVertexList, 12);
	//mesh.init(MyVertices, 8, MyVertexList, 12); // 점은 8개 삼각형은 12개
	vector<int> arr;
	mesh.getOneRingNeighbor(3, arr);
	for (vector<int>::iterator it = arr.begin(); it != arr.end(); it++) {
		cout << *it << ' ';
	}
	cout << endl << mesh.getOppositeEdgeId(3) << endl;

	arr.clear();
	mesh.getOneRingTriNeighbor(0, arr); // 0번 삼각형의 이웃 삼각형(아마도 3개)
	for (vector<int>::iterator it = arr.begin(); it != arr.end(); it++) {
		cout << *it << ' ';
	}

	cout << endl;
	arr.clear();
	mesh.getOneRingTriNeighborFromVertex(0, arr); // 0번 정점의 이웃 삼각형(아마도 많음 )
	for (vector<int>::iterator it = arr.begin(); it != arr.end(); it++) {
		cout << *it << ' ';
	}

	float3 n = mesh.getVertexNormal(0);
	cout << endl << n.x << ' ' << n.y << ' ' << n.z << endl;
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB);
	glutInitWindowSize(300, 300);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("OpenGL Drawing Example");
	glPolygonMode(GL_FRONT, GL_LINE);
	glPolygonMode(GL_BACK, GL_LINE);
	glClearColor(0.5, 0.5, 0.5, 1.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
	glutDisplayFunc(MyDisplay);
	glutMainLoop();
	return 0;
}
