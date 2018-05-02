#ifndef GRAPH_DEFS_H_
#define GRAPH_DEFS_H_

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Lazy_exact_nt.h>
#include <cassert>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/random/linear_congruential.hpp>

// CGAL typedefs (2 space)
//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2 CGALPoint;
typedef K::Segment_2 CGALSegment;
typedef K::Intersect_2 CGALIntersect;
typedef K::Circle_2 CGALCircle;

typedef CGAL::Triangulation_vertex_base_2<K>				Vbb;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb>	Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>		Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>		TDS;
typedef CGAL::No_intersection_tag							Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
typedef CGAL::Triangulation_hierarchy_2<CDT> Triangulation;

typedef CDT::Vertex_handle Vertex_handle;

//typedef CGAL::Triangulation_2<K> Triangulation;
typedef Triangulation::Finite_faces_iterator FiniteFaceIter;
typedef Triangulation::Finite_edges_iterator FiniteEdgeIter;
typedef Triangulation::Finite_vertices_iterator FiniteVertexIter;

typedef Triangulation::Face_handle TriFaceHandle;
typedef Triangulation::Vertex_handle TriVertexHandle;

typedef Triangulation::Face TriFace;
typedef Triangulation::Edge TriEdge;
typedef Triangulation::Vertex TriVertex;
typedef Triangulation::Point TriPoint;

typedef CGAL::Creator_uniform_2<double, TriPoint> Creator;

bool operator==(TriEdge const& e1, TriEdge const& e2)
{
	return e1.first == e2.first && e1.second == e2.second;
}

namespace boost {
	template <> struct hash < TriEdge > {
		size_t operator()(TriEdge const& e) const {
			std::size_t seed = 31;
			boost::hash_combine(seed, e.first);
			boost::hash_combine(seed, e.second);
			return seed;
		}
	};
}

typedef double EdgeWeight;

bool operator==(CGALPoint const& p1, CGALPoint const& p2)
{
	return p1.x() == p2.x() && p1.y() == p2.y();
}

namespace boost {
	template <> struct hash < CGALPoint > {
		size_t operator()(CGALPoint const& p) const {
			std::size_t seed = 31;
			// Warning: Possibly area of precision loss. The hash is dependent on to_double returning the same
			// values that we initialized our CGALPoint with
			double x = CGAL::to_double(p.x());
			double y = CGAL::to_double(p.y());
			boost::hash_combine(seed, x);
			boost::hash_combine(seed, y);
			return seed;
		}
	};
}

// Custom vertex
typedef size_t VertexIndex;
typedef std::vector<CGALPoint*> VertexVector;

// A hashable wrapper for an Edge
struct SimpleEdge {
	VertexIndex u;
	VertexIndex v;
	EdgeWeight weight;

	SimpleEdge(VertexIndex aU_idx, VertexIndex aV_idx, int aWeight) {
		u = aU_idx;
		v = aV_idx;
		weight = aWeight;
	}
};

bool operator==(SimpleEdge const& e1, SimpleEdge const& e2)
{
	return (e1.u == e2.u && e1.v == e2.v) || (e1.u == e2.v && e1.v == e2.u);
}

namespace boost {
	template <> struct hash < SimpleEdge > {
		size_t operator()(SimpleEdge const& e) const {
			std::size_t seed = 31;
			boost::hash<VertexIndex> index_hasher;

			// Add hash so edge endpoints {(x1, y1) (x2, y2)}
			// have the same hash as {(x2, y2) (x1, y1)}
			return index_hasher(e.u) + index_hasher(e.v);
		}
	};
}

// Custom edge
typedef std::vector<SimpleEdge*> EdgeVector;

void deleteVerticesVector(VertexVector* vv) {
	for (VertexVector::iterator it = vv->begin(); it != vv->end(); ++it) {
		delete (*it);
	}
	delete vv;
}

void deleteEdgeVector(EdgeVector* ev) {
	for (EdgeVector::iterator it = ev->begin(); it != ev->end(); ++it) {
		delete (*it);
	}
	delete ev;
}

// Random gen
boost::minstd_rand gen;

inline bool SegmentIntersect(CGALPoint existingU, CGALPoint existingV, CGALPoint inputU, CGALPoint inputY) {
	CGALSegment seg1(existingU, existingV);
	CGALSegment seg2(inputU, inputY);

	CGAL::cpp11::result_of<CGALIntersect(CGALSegment, CGALSegment)>::type result = intersection(seg1, seg2);
	if (result) {
		if (const CGALSegment* s = boost::get<CGALSegment>(&*result)) {
			//std::cout << *s << std::endl;
			return true;
		}
		else if (const CGALPoint* p = boost::get<CGALPoint >(&*result)) {
			//std::cout << " i " << *p;
			// Ignore intersection at segment endpoints
			if (*p != inputU && *p != inputY) {
				return true;
			}
		}
	}
	return false;
}

#endif  // GRAPH_DEFS_H_