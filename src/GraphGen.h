#ifndef GRAPH_GEN_H_
#define GRAPH_GEN_H_

#include "GraphDefs.h"

#include <boost/pending/disjoint_sets.hpp>

/**
* Naive linear time intersection
* Returns true if edge (u, v) intersects an edge in g, otherwise false
**/
bool doesIntersect(VertexVector* vertices, EdgeVector* edges, VertexIndex u, VertexIndex v) {
	CGALPoint* uPt = (*vertices)[u];
	CGALPoint* vPt = (*vertices)[v];
	CGALSegment segUV((*uPt), (*vPt));

	for (int i = 0; i < edges->size(); i++) {
		SimpleEdge* e = (*edges)[i];
		CGALPoint* src = (*vertices)[e->u];
		CGALPoint* tar = (*vertices)[e->v];
		CGALSegment seg((*src), (*tar));

		CGAL::cpp11::result_of<CGALIntersect(CGALSegment, CGALSegment)>::type result = intersection(seg, segUV);
		if (result) {
			if (const CGALSegment* s = boost::get<CGALSegment>(&*result)) {
				//std::cout << *s << std::endl;
				return true;
			}
			else if (const CGALPoint* p = boost::get<CGALPoint >(&*result)) {
				//std::cout << " i " << *p;
				// Ignore intersection at segment endpoints
				if (*p != *uPt && *p != *vPt) {
					return true;
				}
			}
		}
	}
	return false;
}

void createRandomPlaneForest(int numVertices, int radius, int upToNumEdges, VertexVector** vertices, EdgeVector** edges) {
	(*vertices) = new VertexVector();
	(*edges) = new EdgeVector();

	//CGAL::Random_points_in_disc_2<CGALPoint, Creator> randPts(radius);
	CGAL::Random_points_on_circle_2<CGALPoint, Creator> randPts(radius);

	// Generate vertices with random coordinated within bounds
	for (int i = 0; i < numVertices; i++) {
		CGALPoint* pt = new CGALPoint((*randPts++));
		CGALPoint* ptTranslated = new CGALPoint(pt->x() + radius, pt->y() + radius);
		(*vertices)->push_back(ptTranslated);
		delete pt;
	}

	// Define edge random gen
	boost::uniform_int<> vertexRange(0, numVertices - 1);
	boost::variate_generator<boost::minstd_rand, boost::uniform_int<>> vertexDice(gen, vertexRange);
	vertexDice.engine().seed(static_cast<unsigned int>(std::time(0)));

	std::vector<int> rank(numVertices);
	std::vector<int> parent(numVertices);
	boost::disjoint_sets<int*, int*, boost::find_with_full_path_compression> ds(&rank[0], &parent[0]);
	for (int i = 0; i < rank.size(); i++) {
		rank[i] = i;
		parent[i] = i;
	}

	// Select random vertices u, v for edgeRolls number of times
	// An edge connects u, v:
	//		1. u != v
	//		3. adding edge(u, v) does not create a cycle
	//		4. edge(u, v) does not intersect any existing edge
	for (int i = 0; i < upToNumEdges; i++) {
		VertexIndex u = vertexDice();
		VertexIndex v = vertexDice();
		if (u != v
			&& ds.find_set(u) != ds.find_set(v)
			&& !doesIntersect((*vertices), (*edges), u, v)) {

			// Add edge(u, v)
			//std::pair<Edge, bool> result = add_edge(u, v, *g);
			//assert(result.second);
			(*edges)->push_back(new SimpleEdge(u, v, 0));
			ds.link(u, v);
			//std::cout << " - added";
		}
	}
}

#endif  // GRAPH_GEN_H_