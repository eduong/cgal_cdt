#ifndef GRAPH_PARSE_H_
#define GRAPH_PARSE_H_

#include <fstream>

#include <boost/pending/disjoint_sets.hpp>
#include <boost/chrono.hpp>
#include <boost/tokenizer.hpp>

#include "GraphDefs.h"
#include "Timer.h"

#define PERFORM_RESTRICTION_CHECKS false

typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

void parseGraph(const char* vertFile, const char* edgeFile, VertexVector** vertices, EdgeVector** edges) {
	std::ifstream vStream(vertFile);
	std::ifstream eStream(edgeFile);

	(*vertices) = new VertexVector();
	(*edges) = new EdgeVector();

	boost::chrono::high_resolution_clock::time_point start;
	boost::chrono::high_resolution_clock::time_point end;
	boost::chrono::milliseconds duration(0);

	start = boost::chrono::high_resolution_clock::now();

	std::string line;
	boost::char_separator<char> sep(" ");
	
	// Parse vertices count
	// Note: insertion of vertices/edges directly into BoostGraph is very slow
	std::getline(vStream, line);
	int verticesCount = std::stoi(line);
	(*vertices) = new VertexVector();
	(*vertices)->reserve(verticesCount);

	// Parse vertices
	// Note: insertion of vertices/edges directly into BoostGraph is very slow
	while (std::getline(vStream, line))
	{
		std::istringstream iss(line);
		std::string line = iss.str();
		// Vertex lines, e.g.: "v 1 -73530767 41085396"
		tokenizer tokens(line, sep);
		tokenizer::iterator beg = tokens.begin();
		beg++; // Skip the index (first number)
		int x = std::stoi(*beg);
		beg++;
		int y = std::stoi(*beg);
		(*vertices)->push_back(new CGALPoint(x, y));
	}

	vStream.close();

	end = boost::chrono::high_resolution_clock::now();
	duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
	printDuration("Read vertices from file", duration);

	start = boost::chrono::high_resolution_clock::now();

	// Parse edges count
	std::getline(eStream, line);
	int edgeCount = std::stoi(line);
	(*edges) = new EdgeVector();
	(*edges)->reserve(edgeCount);

	// Parse edges
	while (std::getline(eStream, line))
	{
		std::istringstream iss(line);
		std::string line = iss.str();
		// Edge (arc) lines, e.g.: "1 2"
		// Where 1 2 are the vertex index, starting from 0 
		tokenizer tokens(line, sep);
		tokenizer::iterator beg = tokens.begin();
		int u = std::stoi(*beg);
		beg++;
		int v = std::stoi(*beg);
		SimpleEdge* edge = new SimpleEdge(u, v, 0);

		(*edges)->push_back(edge);
	}

	eStream.close();

	end = boost::chrono::high_resolution_clock::now();
	duration = (boost::chrono::duration_cast<boost::chrono::milliseconds>(end - start));
	printDuration("Read edges from file", duration);
}

#endif  // GRAPH_PARSE_H_