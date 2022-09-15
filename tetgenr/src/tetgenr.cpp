#include <cstddef> // for std::size_t
#include <cassert> // for assert
#include <Rcpp.h>
#include "tetgen.h"
using namespace Rcpp;

// [[Rcpp::export]]
List mesh_internal(NumericMatrix input_geometry, IntegerMatrix input_hull, std::string options) {
	if (input_geometry.nrow() < 4) {
		stop("input_geometry.nrow() < 4 (is %i)", input_geometry.nrow());
	}
	if (input_geometry.ncol() != 3) {
		stop("input_geometry.ncol() != 3 (is %i)", input_geometry.ncol());
	}
	if (input_hull.nrow() < 4) {
		stop("input_hull.nrow() < 4 (is %i)", input_hull.nrow());
	}
	if (input_hull.ncol() != 3) {
		stop("input_hull.ncol() != 3 (is %i)", input_hull.ncol());
	}
	std::size_t n_nodes = input_geometry.nrow();
	
	tetgenio input;
	tetgenio output;
	
	// construct input
	input.firstnumber = 1;
	input.numberofpoints = n_nodes;
	input.pointlist = new REAL[3 * n_nodes];
	
	for (int p = 0; p < n_nodes; ++p) {
		for (int d = 0; d < 3; ++d) {
			input.pointlist[3 * p + d] = input_geometry(p, d);
		}
	}
	
	input.numberoffacets = input_hull.nrow();
	input.facetlist = new tetgenio::facet[input.numberoffacets];
	input.facetmarkerlist = new int[input.numberoffacets];
	
	for (int s = 0; s < input.numberoffacets; ++s) {
		input.facetmarkerlist[s] = s + 1;
		tetgenio::facet & f = input.facetlist[s];
		f.numberofpolygons = 1;
		f.polygonlist = new tetgenio::polygon[1];
		f.numberofholes = 0;
		f.holelist = nullptr;
		tetgenio::polygon & p = f.polygonlist[0];
		p.numberofvertices = 3;
		p.vertexlist = new int[3];
		for (int i = 0; i < 3; ++i) {
			p.vertexlist[i] = input_hull(s, i);
		}
	}
	
	// call
	char * input_options = new char[options.size() + 1];
	options.copy(input_options, options.size());
	input_options[options.size()] = '\0';
	
	tetrahedralize(input_options, &input, &output);
	
	delete[] input_options;
	input_options = nullptr;
	
	delete[] input.pointlist;
	input.pointlist = nullptr;
	
	// parse output
	n_nodes = output.numberofpoints;
	NumericMatrix output_geometry(n_nodes, 3);
	for (int p = 0; p < n_nodes; ++p) {
		for (int d = 0; d < 3; ++d) {
			output_geometry(p, d) = output.pointlist[3 * p + d];
		}
	}
	
	assert(output.numberofcorners == 4);
	std::size_t n_elements = output.numberoftetrahedra;
	IntegerMatrix output_topology(n_elements, 4);
	for (int e = 0; e < n_elements; ++e) {
		for (int n = 0; n < 4; ++n) {
			output_topology(e, n) = output.tetrahedronlist[4 * e + n];
		}
	}
	
	List result;
	result["geometry"] = output_geometry;
	result["topology"] = output_topology;
	
	return result;
}