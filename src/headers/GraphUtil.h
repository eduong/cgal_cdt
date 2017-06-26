#ifndef GRAPH_UTIL_H_
#define GRAPH_UTIL_H_

#include "GraphDefs.h"

boost::unordered_set<SimpleEdge>* createSimpleEdgeSet(EdgeVector* edges) {
	boost::unordered_set<SimpleEdge>* contraintEdgeSet = new boost::unordered_set<SimpleEdge>();

	// Map each edge into a hash-able edge
	for (int i = 0; i < edges->size(); i++) {
		SimpleEdge* edge = (*edges)[i];
		contraintEdgeSet->insert(*edge);
	}
	return contraintEdgeSet;
}

#endif  // GRAPH_UTIL_H_