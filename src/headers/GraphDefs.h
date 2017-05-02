#ifndef GRAPH_DEFS_H_
#define GRAPH_DEFS_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <cassert>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/random/linear_congruential.hpp>

// CGAL typedefs (2 space)
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 CGALPoint;
typedef K::Segment_2 CGALSegment;
typedef K::Intersect_2 CGALIntersect;
typedef K::Circle_2 CGALCircle;

typedef CGAL::Triangulation_vertex_base_2<K>				Vbb;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb>	Vb;
typedef CGAL::Triangulation_face_base_2<K>					Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>		Tds;
typedef CGAL::Triangulation_2<K, Tds>						T2;
typedef CGAL::Triangulation_hierarchy_2<T2>					Triangulation;

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

typedef long long EdgeWeight;

struct VertexProperties
{
	CGALPoint pt;
};

struct EdgeProperties
{
	EdgeWeight weight;
};

bool operator==(CGALPoint const& p1, CGALPoint const& p2)
{
	return p1.x() == p2.x() && p1.y() == p2.y();
}

namespace boost {
	template <> struct hash < CGALPoint > {
		size_t operator()(CGALPoint const& p) const {
			std::size_t seed = 31;
			boost::hash_combine(seed, p.x());
			boost::hash_combine(seed, p.y());
			return seed;
		}
	};
}

struct SimpleEdge {
	CGALPoint u;
	CGALPoint v;
	int u_idx;
	int v_idx;
	EdgeWeight weight;

	SimpleEdge(CGALPoint aU, CGALPoint aV, int aU_idx, int aV_idx, int aWeight) {
		u = aU;
		v = aV;
		u_idx = aU_idx;
		v_idx = aV_idx;
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
			boost::hash<CGALPoint> point_hasher;

			// Add hash so edge endpoints {(x1, y1) (x2, y2)}
			// have the same hash as {(x2, y2) (x1, y1)}
			return point_hasher(e.u) + point_hasher(e.v);
		}
	};
}

typedef boost::adjacency_list <
	boost::hash_setS, // OutEdgeList
	boost::vecS, // VertexList
	boost::undirectedS, // Undirected edges
	VertexProperties, // Vertex obj representation
	EdgeProperties, // Edge obj representation
	boost::no_property, // Graph obj representation
	boost::listS > // EdgeList (there is limited documentation on how to change this)
	BoostGraph;

// Boost typdefs
typedef boost::graph_traits<BoostGraph>::vertex_descriptor Vertex;
typedef boost::graph_traits<BoostGraph>::edge_descriptor Edge;

typedef boost::graph_traits<BoostGraph>::vertex_iterator VertexIter;
typedef boost::graph_traits<BoostGraph>::edge_iterator EdgeIter;

// Random gen
boost::minstd_rand gen;

#endif  // GRAPH_DEFS_H_