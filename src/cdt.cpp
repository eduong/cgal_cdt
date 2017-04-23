#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <boost/chrono.hpp>
#include <boost/graph/iteration_macros.hpp>

#include "GraphDefs.h"
#include "timer.h"
#include "verify.h"

#include <fstream>
#include <string>

#define SHOW_DEBUG false

const int vertexIdx = 0; // Arbitrary index id

BoostGraph* CreateRandomVertices(int numVertices, int bounds) {
	BoostGraph* g = new BoostGraph();

	// Define bounds random gen
	boost::uniform_int<> boundsRange(0, bounds);
	boost::variate_generator<boost::minstd_rand, boost::uniform_int<>> boundsDice(gen, boundsRange);
	boundsDice.engine().seed(static_cast<unsigned int>(std::time(0)));

	// To maintain unique x, y
	boost::unordered_map<int, boost::unordered_map<int, int>> xBounds;

	// Generate vertices with random coordinated within bounds
	for (int i = 0; i < numVertices; i++) {
		int x = boundsDice();
		int y = boundsDice();
		if (xBounds.count(x)) { // x exists
			while (xBounds[x].count(y)) { // y exists, reroll (TODO: range check to prevent infinite loop)
				y = boundsDice();
			}
			xBounds[x].emplace(y, NULL);
		}
		else {
			boost::unordered_map<int, int> yBounds;
			yBounds.emplace(y, NULL);
			xBounds.emplace(x, yBounds);
		}
		Vertex v = add_vertex(*g);
		(*g)[v].pt = CGALPoint(x, y);
	}

	return g;
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

bool IsLocallyDelaunay(Triangulation t, TriFaceHandle f, TriVertexHandle v) {
	CGAL::Oriented_side side = t.side_of_oriented_circle(f, v->point());
	// v is outside or on the circumcircle of f, i.e. ON_NEGATIVE_SIDE or ON_ORIENTED_BOUNDARY, respectively
	return side == CGAL::Oriented_side::ON_NEGATIVE_SIDE
		|| side == CGAL::Oriented_side::ON_ORIENTED_BOUNDARY;
}

int main(int argc, char* argv[]) {
	Triangulation t; // Could use Triangulation_hierarchy_2 for better performance?

	BoostGraph* g = CreateRandomVertices(1000, 1000);

	std::pair<VertexIter, VertexIter> vp;
	for (vp = boost::vertices(*g); vp.first != vp.second; ++vp.first) {
		Vertex v = *vp.first;
		t.insert((*g)[v].pt);
	}

	// Simple example
	/*CGALPoint p0(1, 1);
	CGALPoint p1(5, 5);
	CGALPoint p2(4, 1);
	CGALPoint p3(2, 4);

	t.insert(p0);
	t.insert(p1);
	t.insert(p2);
	t.insert(p3);*/

	boost::unordered_set<TriEdge>* s = new boost::unordered_set<TriEdge>();

	int edgeCount = 0;
	for (FiniteEdgeIter iter = t.finite_edges_begin(); iter != t.finite_edges_end(); ++iter) {
		edgeCount++;

		// typedef std::pair<Face_handle, int> Edge;
		TriEdge e = *iter;
		int eIndex = e.second;

		// Edge shared by faces f0 and f1
		TriFaceHandle f0 = e.first;
		TriFaceHandle f1 = e.first->neighbor(eIndex);

		// Vertex opposite of edge e in f0, f1
		TriVertexHandle opp0 = f0->vertex(eIndex);

		if (SHOW_DEBUG) {
			// Vertex endpoint v0, v1 of edge e for f0 (doesn't seem to be possible to identify the same edge for f1 in a similar way)
			TriVertexHandle v0f0 = f0->vertex(f0->cw(eIndex));
			TriVertexHandle v1f0 = f0->vertex(f0->ccw(eIndex));

			std::cout << eIndex << std::endl;
			std::cout << "(" << *v0f0 << ") (" << *v1f0 << ")" << std::endl;
			std::cout << "(" << *opp0 << ")" << std::endl;
		}

		if (opp0 != t.infinite_vertex() && !IsLocallyDelaunay(t, f1, opp0)) {
			s->emplace(e);
		}
	}

	std::cout << "Edges in T: " << edgeCount << " Edges in S: " << s->size() << " Ratio: " << (double)((double)s->size() / (double)edgeCount) << std::endl;

	delete s;

	return 0;
}