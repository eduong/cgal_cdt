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
		std::cout << "(" << (*p) << ") (" << (*q) << ") (" << (*r) << ") (" << (*t) << ")" << std::endl;
	}
	CGALCircle c(*p, *q, *r);
	CGAL::Bounded_side side = c.bounded_side(*t);

	// v is outside or on the circumcircle of f
	return side == CGAL::Bounded_side::ON_UNBOUNDED_SIDE
		|| side == CGAL::Bounded_side::ON_BOUNDARY;
}

int main(int argc, char* argv[]) {
	std::cout << "insertion of random points" << std::endl;
	Triangulation* t = new Triangulation();
	CGAL::Random_points_on_circle_2<TriPoint, Creator> graph(100.);
	CGAL::cpp11::copy_n(graph, 100000, std::back_inserter(*t));
	//verbose mode of is_valid ; shows the number of vertices at each  level
	std::cout << "The number of vertices at successive levels" << std::endl;
	assert(t->is_valid(true));

	boost::unordered_set<TriEdge>* s = new boost::unordered_set<TriEdge>();

	TriVertexHandle infiniteVertex = t->infinite_vertex();
	TriFaceHandle infiniteFace = t->infinite_face();
	int edgeCount = 0;

	boost::chrono::high_resolution_clock::time_point startTotal = boost::chrono::high_resolution_clock::now();

	for (FiniteEdgeIter iter = t->finite_edges_begin(); iter != t->finite_edges_end(); ++iter) {
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
					s->emplace(e);
				}
			}
		}
	}

	boost::chrono::high_resolution_clock::time_point endTotal = boost::chrono::high_resolution_clock::now();
	boost::chrono::milliseconds total = (boost::chrono::duration_cast<boost::chrono::milliseconds>(endTotal - startTotal));
	printDuration("Total duration", total);

	std::cout << "Edges in T: " << edgeCount << " Edges in S: " << s->size() << " Ratio: " << (double)((double)s->size() / (double)edgeCount) << std::endl;

	delete s;
	delete t;

	return 0;
}