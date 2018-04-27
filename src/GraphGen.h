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

void createRandomNearTriangulation(int numVertices, int radius, VertexVector** vertices, EdgeVector** edges) {
	(*vertices) = new VertexVector();
	(*edges) = new EdgeVector();

	CGAL::Random_points_in_disc_2<CGALPoint, Creator> randPts(radius);
	//CGAL::Random_points_on_circle_2<CGALPoint, Creator> randPts(radius);

	// Generate vertices with random coordinated within bounds
	for (int i = 0; i < numVertices; i++) {
		CGALPoint* pt = new CGALPoint((*randPts++));
		(*vertices)->push_back(pt);
	}

	CDT* cdt = new CDT();
	boost::unordered_map<VertexIndex, Vertex_handle> vertexHandles;
	for (int i = 0; i < (*vertices)->size(); i++) {
		CGALPoint* pt = (**vertices)[i];
		Vertex_handle vHandle = cdt->insert(*pt);
		vertexHandles.emplace(i, vHandle);
	}

	// Map CGALPoint -> VertexIndex
	boost::unordered_map<CGALPoint, VertexIndex> vertexIndex;
	for (int i = 0; i < (*vertices)->size(); i++) {
		vertexIndex[*(**vertices)[i]] = i;
	}

	// Disjoint set for forest property
	std::vector<int> rank(numVertices);
	std::vector<int> parent(numVertices);
	boost::disjoint_sets<int*, int*, boost::find_with_full_path_compression> ds(&rank[0], &parent[0]);
	for (int i = 0; i < rank.size(); i++) {
		rank[i] = i;
		parent[i] = i;
	}

	// Define edge random gen
	boost::uniform_int<> randomRange(0, 9);
	boost::variate_generator<boost::minstd_rand, boost::uniform_int<>> edgeDice(gen, randomRange);
	edgeDice.engine().seed(static_cast<unsigned int>(std::time(0)));

	// Add edges to graph
	for (CDT::Edge_iterator eit = cdt->edges_begin(); eit != cdt->edges_end(); ++eit) {
		CDT::Edge cgal_e = *eit;
		CGALSegment segement = cdt->segment(cgal_e);
		CGALPoint cgal_u = segement.point(0);
		CGALPoint cgal_v = segement.point(1);
		VertexIndex u = vertexIndex[cgal_u];
		VertexIndex v = vertexIndex[cgal_v];
		int rand = edgeDice();

		unsigned long long setU = ds.find_set(u);
		unsigned long long setV = ds.find_set(v);

		if (ds.find_set(u) != ds.find_set(v)
			&& rand == 9) {
			ds.link(u, v);
			SimpleEdge* edge = new SimpleEdge(u, v, 0);
			(*edges)->push_back(edge);
		}
	}
}

void createRandomPlaneForest(int numVertices, int radius, int numEdges, VertexVector** vertices, EdgeVector** edges) {
	(*vertices) = new VertexVector();
	(*edges) = new EdgeVector();

	//CGAL::Random_points_in_disc_2<CGALPoint, Creator> randPts(radius);
	CGAL::Random_points_on_circle_2<CGALPoint, Creator> randPts(radius);

	// Generate vertices with random coordinated within bounds
	for (int i = 0; i < numVertices; i++) {
		CGALPoint* pt = new CGALPoint((*randPts++));
		(*vertices)->push_back(pt);
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
	int edgeCount = 0;
	while (edgeCount < numEdges) {
		VertexIndex u = vertexDice();
		VertexIndex v = vertexDice();
		// An edge connects u, v:
		//		1. u != v
		//		3. adding edge(u, v) does not create a cycle
		//		4. edge(u, v) does not intersect any existing edge
		if (u != v
			&& ds.find_set(u) != ds.find_set(v)
			&& !doesIntersect((*vertices), (*edges), u, v)) {

			// Add edge(u, v)
			(*edges)->push_back(new SimpleEdge(u, v, 0));
			ds.link(u, v);
			edgeCount++;
			//std::cout << " - added";
		}
	}
}

#endif  // GRAPH_GEN_H_