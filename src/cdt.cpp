#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <boost/chrono.hpp>
#include <boost/graph/iteration_macros.hpp>

#include "GraphDefs.h"
#include "timer.h"
#include "GraphUtil.h"
#include "GraphParse.h"
#include "GraphGen.h"

#include <fstream>
#include <string>

#define SHOW_DEBUG true

const int vertexIdx = 0; // Arbitrary index id

void printGraph(const char* title, VertexVector* vertices, EdgeVector* edges, bool printVertices = false) {
	std::cout << std::endl << "=== " << title << std::endl;

	if (printVertices) {
		// Iterate through the vertices and print them out
		for (int i = 0; i < vertices->size(); i++) {
			CGALPoint* v = (*vertices)[i];
			std::cout << "index: " << i << " (" << v->x() << ", " << v->y() << ")" << std::endl;
		}
	}

	// Iterate through the edges and print them out
	std::cout << std::endl << "Edges: " << std::endl;
	for (int i = 0; i < edges->size(); i++) {
		SimpleEdge* e = (*edges)[i];
		CGALPoint* src = (*vertices)[e->u];
		CGALPoint* tar = (*vertices)[e->v];
		//std::cout << edgeWeightMap->at(e).weight << " (" << (*g)[src].pt << ") (" << (*g)[tar].pt << ")" << std::endl;
		std::cout << e->u << " " << e->v << " (" << (*src) << ") (" << (*tar) << ")" << std::endl;
	}

	std::cout << std::endl;
}

TriVertexHandle OppositeOfEdge(TriVertexHandle ev0, TriVertexHandle ev1, TriFaceHandle f) {
	TriVertexHandle v0 = f->vertex(vertexIdx);
	if (v0 != ev0 && v0 != ev1) {
		return v0;
	}

	TriVertexHandle v1 = f->vertex(f->cw(vertexIdx));
	if (v1 != ev0 && v1 != ev1) {
		return v1;
	}

	TriVertexHandle v2 = f->vertex(f->ccw(vertexIdx));
	if (v2 != ev0 && v2 != ev1) {
		return v2;

	}

	assert(false); // Unreachable
	return NULL;
}

bool IsLocallyDelaunay(CGALPoint* p, CGALPoint* q, CGALPoint* r, CGALPoint* t) {
	if (SHOW_DEBUG) {
		//std::cout << "(" << (*p) << ") (" << (*q) << ") (" << (*r) << ") (" << (*t) << ")" << std::endl;
	}
	CGALCircle c(*p, *q, *r);
	CGAL::Bounded_side side = c.bounded_side(*t);

	// v is outside or on the circumcircle of f
	return side == CGAL::Bounded_side::ON_UNBOUNDED_SIDE
		|| side == CGAL::Bounded_side::ON_BOUNDARY;
}

CDT* computeCdt(VertexVector* vertices, EdgeVector* edges, boost::unordered_map<TriVertexHandle, VertexIndex>** handlesToIndex) {
	CDT* cdt = new CDT();

	(*handlesToIndex) = new boost::unordered_map<TriVertexHandle, VertexIndex>();
	boost::unordered_map<VertexIndex, TriVertexHandle> vertexHandles;
	for (int i = 0; i < vertices->size(); i++) {
		CGALPoint* pt = (*vertices)[i];
		TriVertexHandle vHandle = cdt->insert(*pt);
		vertexHandles.emplace(i, vHandle);
		(*handlesToIndex)->emplace(vHandle, i);
	}

	// Insert constraint edges
	for (int i = 0; i < edges->size(); i++) {
		SimpleEdge* edge = (*edges)[i];
		VertexIndex u = edge->u;
		VertexIndex v = edge->v;
		TriVertexHandle uH = vertexHandles[u];
		TriVertexHandle vH = vertexHandles[v];
		cdt->insert_constraint(uH, vH);
	}

	assert(cdt->is_valid());

	return cdt;
}

// Given 3 colinear points, if the 2 furthest points are constrained, how is this delaunay triangulated?
// It seems the outcome is implementation dependent. CGAL will split the constraint into 2 edges.
// Other implementations might allow for overlapping, collinear edges. Either way, the current plan
// is to assume that an input forest, F, with collinear edges will be replaced by smaller constraint edges.
EdgeVector* newConstraintSetFromCt(CDT* cdt, VertexVector* originalVertices) {
	EdgeVector* newEdgeVector = new EdgeVector();

	boost::unordered_map<CGALPoint, VertexIndex> vertexIndex;
	for (int i = 0; i < originalVertices->size(); i++) {
		vertexIndex[*(*originalVertices)[i]] = i;
	}

	// Add edges to graph
	for (CDT::Edge_iterator eit = cdt->edges_begin(); eit != cdt->edges_end(); ++eit) {
		CDT::Edge cgal_e = *eit;

		if (cdt->is_constrained(cgal_e)) {
			// Assumes point coord are unique
			CGALSegment segement = cdt->segment(cgal_e);
			CGALPoint cgal_u = segement.point(0);
			CGALPoint cgal_v = segement.point(1);
			VertexIndex u = vertexIndex[cgal_u];
			VertexIndex v = vertexIndex[cgal_v];

			SimpleEdge* edge = new SimpleEdge(u, v, 0);
			newEdgeVector->push_back(edge);
		}
	}
	return newEdgeVector;
}

void computeNonLocallyDelaunay(
	VertexVector* vertices,
	EdgeVector* edges,
	EdgeVector** NewEdges,
	EdgeVector** S_Edges) {

	boost::unordered_map<TriVertexHandle, VertexIndex>* handlesToIndex;
	CDT* cdt = computeCdt(vertices, edges, &handlesToIndex);
	(*NewEdges) = newConstraintSetFromCt(cdt, vertices);

	boost::unordered_set<TriEdge>* S = new boost::unordered_set<TriEdge>();

	TriVertexHandle infiniteVertex = cdt->infinite_vertex();
	TriFaceHandle infiniteFace = cdt->infinite_face();
	int edgeCount = 0;

	boost::chrono::high_resolution_clock::time_point startTotal = boost::chrono::high_resolution_clock::now();

	for (FiniteEdgeIter iter = cdt->finite_edges_begin(); iter != cdt->finite_edges_end(); ++iter) {
		edgeCount++;

		// typedef std::pair<Face_handle, int> Edge;
		TriEdge e = *iter;
		int eIndex = e.second;

		// Edge shared by faces f0 and f1
		TriFaceHandle f0 = e.first;
		TriFaceHandle f1 = e.first->neighbor(eIndex);

		// Vertex opposite of edge e in f0
		TriVertexHandle opp0 = f0->vertex(eIndex);

		if (SHOW_DEBUG) {
			// Vertex endpoint v0, v1 of edge e for f0 (doesn't seem to be possible to identify the same edge for f1 in a similar way)
			TriVertexHandle v0f0 = f0->vertex(f0->cw(eIndex));
			TriVertexHandle v1f0 = f0->vertex(f0->ccw(eIndex));

			//std::cout << eIndex << std::endl;
			std::cout << "(" << *v0f0 << ") (" << *v1f0 << ")" << "(" << *opp0 << ")" << std::endl;
		}

		if (opp0 != infiniteVertex && f1 != infiniteFace) {
			TriVertexHandle pH = f1->vertex(0);
			TriVertexHandle qH = f1->vertex(1);
			TriVertexHandle rH = f1->vertex(2);

			if (pH != infiniteVertex && qH != infiniteVertex && rH != infiniteVertex) {
				CGALPoint p = pH->point();
				CGALPoint q = qH->point();
				CGALPoint r = rH->point();
				CGALPoint t = opp0->point();
				if (!IsLocallyDelaunay(&p, &q, &r, &t)) {
					S->emplace(e);
				}
			}
		}

		//// Vertex of edge e
		//TriVertexHandle e0 = f0->vertex(f0->cw(eIndex));
		//TriVertexHandle e1 = f0->vertex(f0->ccw(eIndex));
		//
		//// Vertex opposite of edge e in f1
		//TriVertexHandle opp1 = OppositeOfEdge(e0, e1, f1);
		//
		//if (opp1 != infiniteVertex && f1 != infiniteFace) {
		//	TriVertexHandle pH = f1->vertex(0);
		//	TriVertexHandle qH = f1->vertex(1);
		//	TriVertexHandle rH = f1->vertex(2);

		//	if (pH != infiniteVertex && qH != infiniteVertex && rH != infiniteVertex) {
		//		CGALPoint p = pH->point();
		//		CGALPoint q = qH->point();
		//		CGALPoint r = rH->point();
		//		CGALPoint t = opp1->point();
		//		if (!IsLocallyDelaunay(&p, &q, &r, &t)) {
		//			S->emplace(e);
		//		}
		//	}
		//}
	}

	(*S_Edges) = new EdgeVector();
	for (boost::unordered_set<TriEdge>::iterator iter = S->begin(); iter != S->end(); ++iter) {
		TriEdge e = *iter;
		int eIndex = e.second;
		TriFaceHandle f0 = e.first;
		TriVertexHandle e0 = f0->vertex(f0->cw(eIndex));
		TriVertexHandle e1 = f0->vertex(f0->ccw(eIndex));

		VertexIndex u = (*handlesToIndex)[e0];
		VertexIndex v = (*handlesToIndex)[e1];

		(*S_Edges)->push_back(new SimpleEdge(u, v, 0));
	}

	boost::chrono::high_resolution_clock::time_point endTotal = boost::chrono::high_resolution_clock::now();
	boost::chrono::milliseconds total = (boost::chrono::duration_cast<boost::chrono::milliseconds>(endTotal - startTotal));
	printDuration("computeCDT duration", total);

	delete S;
	delete cdt;
	delete handlesToIndex;
}

EdgeVector* intersectInputSetWithConstraintSet(EdgeVector* inputSet, EdgeVector* constraintSet) {
	boost::unordered_set<SimpleEdge>* NewEdges_Hashset = new boost::unordered_set<SimpleEdge>();
	for (int i = 0; i < inputSet->size(); i++) {
		SimpleEdge* e = (*inputSet)[i];
		NewEdges_Hashset->emplace((*e));
	}

	EdgeVector* intersect = new EdgeVector();
	for (int i = 0; i < constraintSet->size(); i++) {
		SimpleEdge* se = (*constraintSet)[i];
		if (NewEdges_Hashset->count(*se) > 0) {
			intersect->push_back(new SimpleEdge(se->u, se->v, se->weight));
		}
	}

	return intersect;
}

bool containsEdge(boost::unordered_set<SimpleEdge>* edgeSet, SimpleEdge* edge) {
	SimpleEdge se(edge->u, edge->v, 0);
	return edgeSet->count(se) > 0;
}

// True if A a subgraph of B
bool isSubgraph(VertexVector* vertices, EdgeVector* a, EdgeVector* b) {
	boost::unordered_set<SimpleEdge>* bEdgeSet = createSimpleEdgeSet(b);

	// Iterate through the edges
	for (int i = 0; i < a->size(); i++) {
		SimpleEdge* edge = (*a)[i];
		if (!containsEdge(bEdgeSet, edge)) {
			CGALPoint* u = (*vertices)[edge->u];
			CGALPoint* v = (*vertices)[edge->v];
			EdgeWeight weight = CGAL::squared_distance(*u, *v);
			std::cout << "b contains edge not in a: " << weight << " (" << *u << ") (" << *v << ")" << std::endl;
			return false;
		}
	}

	return true;
}

EdgeVector* convertCdtToGraph(VertexVector* vertices, CDT* cdt) {
	EdgeVector* edgeVec = new EdgeVector();
	// edgeVec->reserve(cdt->number_of_faces() * 3); // Upper bound on number of edges (typically too much)

	// Map CGALPoint -> VertexIndex
	boost::unordered_map<CGALPoint, VertexIndex> vertexIndex;
	for (int i = 0; i < vertices->size(); i++) {
		vertexIndex[*(*vertices)[i]] = i;
	}

	// Add edges to graph
	for (CDT::Edge_iterator eit = cdt->edges_begin(); eit != cdt->edges_end(); ++eit) {
		CDT::Edge cgal_e = *eit;
		CGALSegment segement = cdt->segment(cgal_e);
		CGALPoint cgal_u = segement.point(0);
		CGALPoint cgal_v = segement.point(1);
		VertexIndex u = vertexIndex[cgal_u];
		VertexIndex v = vertexIndex[cgal_v];

		SimpleEdge* edge = new SimpleEdge(u, v, 0);
		edgeVec->push_back(edge);
	}

	return edgeVec;
}

bool isCdtSubgraph(VertexVector* vertices, EdgeVector* edgesF, EdgeVector* edgesS) {
	boost::unordered_map<TriVertexHandle, VertexIndex>* handlesToIndex;
	CDT* cdtS = computeCdt(vertices, edgesS, &handlesToIndex); // CDT of mimimum edge constraint
	EdgeVector* ev_cdtS = convertCdtToGraph(vertices, cdtS);
	bool res = isSubgraph(vertices, edgesF, ev_cdtS);

	if (!res) {
		printGraph("Input", vertices, edgesF, true);
		std::cout << "foo";
	}

	delete ev_cdtS;
	delete cdtS;
	delete handlesToIndex;

	return res;
}

int main(int argc, char* argv[]) {
	const char* vertFile = (argc > 2) ? argv[1] : NULL;
	const char* edgeFile = (argc > 2) ? argv[2] : NULL;
	for (;;) {
		VertexVector* vertices = NULL;
		EdgeVector* edges = NULL;

		if (vertFile == NULL || edgeFile == NULL) {
			// Random graph
			createRandomPlaneForest(10, 10, 10, &vertices, &edges);
		}
		else {
			// Load graph from file
			parseGraph(vertFile, edgeFile, &vertices, &edges);
		}

		/*vertices = new VertexVector();
		vertices->push_back(new CGALPoint(1.44669, 4.81919));
		vertices->push_back(new CGALPoint(2.91838, 17.0605));
		vertices->push_back(new CGALPoint(11.9717, 0.196299));
		vertices->push_back(new CGALPoint(0.264731, 7.71428));
		vertices->push_back(new CGALPoint(18.7201, 14.895));
		vertices->push_back(new CGALPoint(0.102537, 8.57164));
		vertices->push_back(new CGALPoint(1.12243, 5.39687));
		vertices->push_back(new CGALPoint(9.96529, 6.02341e-005));
		vertices->push_back(new CGALPoint(19.8294, 11.8395));
		vertices->push_back(new CGALPoint(6.86348, 19.4954));

		edges = new EdgeVector();
		edges->push_back(new SimpleEdge(0, 8, 0));
		edges->push_back(new SimpleEdge(6, 9, 0));
		edges->push_back(new SimpleEdge(8, 1, 0));*/

		if (SHOW_DEBUG) {
			printGraph("Input", vertices, edges, true);
		}

		// Compute CT(V, E)
		EdgeVector* NewEdges;
		EdgeVector* cdtS;
		computeNonLocallyDelaunay(vertices, edges, &NewEdges, &cdtS);
		//EdgeVector* S = intersectInputSetWithConstraintSet(NewEdges, cdtS);
		//std::cout << "Edges in E: " << NewEdges->size() << " Edges in S: " << S->size() << " Ratio: " << (double)((double)S->size() / (double)NewEdges->size()) << std::endl;
		std::cout << "Edges in E: " << NewEdges->size() << " Edges in S: " << cdtS->size() << " Ratio: " << (double)((double)cdtS->size() / (double)NewEdges->size()) << std::endl;

		if (SHOW_DEBUG) {
			printGraph("computeNonLocallyDelaunay", vertices, cdtS, false);
			//printGraph("S", vertices, S, false);
		}

		// Validatation:
		// S ⊆ E
		//if (!isSubgraph(vertices, S, NewEdges)) {
			//std::cout << "Error: isSubgraph is false" << std::endl;
		//}

		// F ⊆ CDT(V, S)
		//if (!isCdtSubgraph(vertices, NewEdges, S)) {
		if (!isCdtSubgraph(vertices, NewEdges, cdtS)) {
			std::cout << "Error: isCdtSubgraph is false" << std::endl;
		}

		//deleteEdgeVector(S);
		deleteEdgeVector(cdtS);
		deleteEdgeVector(NewEdges);
		deleteEdgeVector(edges);
		deleteVerticesVector(vertices);
	}
	return 0;
}