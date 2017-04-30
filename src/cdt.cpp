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

TriVertexHandle OppositeOfEdge(TriVertexHandle ev0, TriVertexHandle ev1, TriFaceHandle f) {
	const int vertexIdx = 0; // Arbitrary index id
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
	// v is outside (ON_NEGATIVE_SIDE) or on the circumcircle (ON_ORIENTED_BOUNDARY) of f
	return side == CGAL::Oriented_side::ON_NEGATIVE_SIDE
		|| side == CGAL::Oriented_side::ON_ORIENTED_BOUNDARY;
}

int main(int argc, char* argv[]) {
	std::cout << "insertion of 1000 random points" << std::endl;
	Triangulation t;
	CGAL::Random_points_on_circle_2<TriPoint, Creator> graph(1.);
	CGAL::cpp11::copy_n(graph, 1000, std::back_inserter(t));
	//verbose mode of is_valid ; shows the number of vertices at each  level
	std::cout << "The number of vertices at successive levels" << std::endl;
	assert(t.is_valid(true));

	boost::unordered_set<TriEdge>* s = new boost::unordered_set<TriEdge>();

	TriVertexHandle infiniteVertex = t.infinite_vertex();
	int edgeCount = 0;

	boost::chrono::high_resolution_clock::time_point startTotal = boost::chrono::high_resolution_clock::now();

	boost::chrono::high_resolution_clock::time_point start;
	boost::chrono::high_resolution_clock::time_point end;
	boost::chrono::milliseconds duration(0);

	for (FiniteEdgeIter iter = t.finite_edges_begin(); iter != t.finite_edges_end(); ++iter) {
		edgeCount++;

		start = boost::chrono::high_resolution_clock::now();

		// typedef std::pair<Face_handle, int> Edge;
		TriEdge e = *iter;
		int eIndex = e.second;

		// Edge shared by faces f0 and f1
		TriFaceHandle f0 = e.first;
		TriFaceHandle f1 = e.first->neighbor(eIndex);

		end = boost::chrono::high_resolution_clock::now();
		duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
		printDuration("Find 2 shared face of edge e", duration);

		start = boost::chrono::high_resolution_clock::now();

		// Vertex opposite of edge e in f0, f1
		TriVertexHandle opp0 = f0->vertex(eIndex);

		end = boost::chrono::high_resolution_clock::now();
		duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
		printDuration("Get vertex that is not e in f0", duration);

		if (SHOW_DEBUG) {
			start = boost::chrono::high_resolution_clock::now();

			// Vertex endpoint v0, v1 of edge e for f0 (doesn't seem to be possible to identify the same edge for f1 in a similar way)
			TriVertexHandle v0f0 = f0->vertex(f0->cw(eIndex));
			TriVertexHandle v1f0 = f0->vertex(f0->ccw(eIndex));

			std::cout << eIndex << std::endl;
			std::cout << "(" << *v0f0 << ") (" << *v1f0 << ")" << std::endl;
			std::cout << "(" << *opp0 << ")" << std::endl;

			end = boost::chrono::high_resolution_clock::now();
			duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
			printDuration("Debug data", duration);
		}

		start = boost::chrono::high_resolution_clock::now();
		if (opp0 != infiniteVertex && !IsLocallyDelaunay(t, f1, opp0)) {
			s->emplace(e);
		}
		end = boost::chrono::high_resolution_clock::now();
		duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
		printDuration("Is not locally delaunay check", duration);
	}

	boost::chrono::high_resolution_clock::time_point endTotal = boost::chrono::high_resolution_clock::now();
	boost::chrono::milliseconds total = (boost::chrono::duration_cast<boost::chrono::milliseconds>(endTotal - startTotal));
	printDuration("Total duration", total);

	std::cout << "Edges in T: " << edgeCount << " Edges in S: " << s->size() << " Ratio: " << (double)((double)s->size() / (double)edgeCount) << std::endl;

	delete s;

	return 0;
}