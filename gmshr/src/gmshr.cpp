#include <cstddef> // for std::size_t
#include <vector> // for std::vector
#include <unordered_map> // for std::unordered_map
#include <stdexcept> // for std::length_error, std::domain_error
#include <string> // for std::string
#include <boost/container_hash/hash.hpp>
#include <boost/dynamic_bitset.hpp>
#include <Rcpp.h>
#include <gmsh.h>
using namespace Rcpp;

// [[Rcpp::export]]
List mesh_internal(NumericMatrix input_geometry, IntegerMatrix input_hull, List options) {
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
	
	if (gmsh::isInitialized()) {
		gmsh::finalize();
	}
	gmsh::initialize();
	
	try {
		std::vector<std::string> option_keys = options.names();
		for (int i = 0; i < option_keys.size(); ++i) {
			auto key = option_keys[i];
			RObject value_any = options[key];
			if (is<NumericVector>(value_any)) {
				NumericVector value = as<NumericVector>(value_any);
				if (value.length() != 1) {
					throw std::length_error(key);
				}
				gmsh::option::setNumber(
					key,
					value[0]
				);
			} else if (is<IntegerVector>(value_any)) {
				IntegerVector value = as<IntegerVector>(value_any);
				if (value.size() != 1) {
					throw std::length_error(key);
				}
				gmsh::option::setNumber(
					key,
					value[0]
				);
			} else if (is<LogicalVector>(value_any)) {
				LogicalVector value = as<LogicalVector>(value_any);
				if (value.size() != 1) {
					throw std::length_error(key);
				}
				gmsh::option::setNumber(
					key,
					value[0] ? 1 : 0
				);
			} else if (is<CharacterVector>(value_any)) {
				CharacterVector value = as<CharacterVector>(value_any);
				if (value.size() != 1) {
					throw std::length_error(key);
				}
				gmsh::option::setString(
					key,
					static_cast<std::string>(value[0])
				);
			} else {
				throw std::domain_error(key);
			}
		}
	
		for (int p = 0; p < n_nodes; ++p) {
			gmsh::model::geo::addPoint(
				input_geometry(p, 0),
				input_geometry(p, 1),
				input_geometry(p, 2),
				0,
				p + 1
			);
		}
		
		std::vector<int> surface_tags;
		std::unordered_map<
			std::pair<std::size_t, std::size_t>,
			int,
			boost::hash<std::pair<std::size_t, std::size_t>>
		> line_tags;
		boost::dynamic_bitset<> boundary_node_tags(n_nodes + 1);
		for (int s = 0; s < input_hull.nrow(); ++s) {
			std::vector<int> l = {0, 0, 0};
			
			for (int i = 0; i < 3; ++i) {
				int start = input_hull(s, i);
				int end = input_hull(s, (i + 1) % 3);
				boundary_node_tags.set(start);
				boundary_node_tags.set(end);
				
				if (line_tags.count({start, end}) > 0) {
					l[i] = line_tags[{start, end}];
				} else if (line_tags.count({end, start}) > 0) {
					l[i] = -1 * line_tags[{end, start}];
				} else {
					l[i] = line_tags[{start, end}] = gmsh::model::geo::addLine(start, end);
				}
			}
			
			gmsh::model::geo::addCurveLoop(l, s + 1, true);
			gmsh::model::geo::addPlaneSurface({s + 1}, s + 1);
			
			surface_tags.push_back(s + 1);
		}
		gmsh::model::geo::addSurfaceLoop(surface_tags, 1);
		gmsh::model::geo::addVolume({1}, 1);
		gmsh::model::geo::synchronize();
		
		std::vector<int> internal_node_tags;
		for (int p = 1; p < n_nodes + 1; ++p) {
			if (!boundary_node_tags.test(p)) {
				internal_node_tags.push_back(p);
			}
		}
		gmsh::model::mesh::embed(0, internal_node_tags, 3, 1);
		
		// gmsh::option::setNumber("Mesh.MeshSizeMin", 2.0);
		gmsh::model::mesh::generate(3);
		gmsh::model::mesh::renumberNodes();
		gmsh::model::mesh::renumberElements();
		
		std::vector<std::size_t> node_tags;
		std::vector<double> node_coordinates;
		std::vector<double> node_parametric_coordinates;
		gmsh::model::mesh::getNodes(node_tags, node_coordinates, node_parametric_coordinates);
		
		n_nodes = node_tags.size();
		
		NumericMatrix output_geometry(n_nodes, 3);
		for (int i = 0; i < n_nodes; ++i) {
			int p = node_tags[i] - 1;
			for (int d = 0; d < 3; ++d) {
				output_geometry(p, d) = node_coordinates[3 * p + d];
			}
		}
		
		int element_type = gmsh::model::mesh::getElementType("tetrahedron", 1);
		std::string element_name;
		int element_dimension;
		int element_order;
		int element_n_nodes;
		std::vector<double> element_local_coordinates;
		int element_n_nodes_primary;
		gmsh::model::mesh::getElementProperties(
			element_type,
			element_name,
			element_dimension,
			element_order,
			element_n_nodes,
			element_local_coordinates,
			element_n_nodes_primary
		);
		
		if (element_dimension != 3) {
			stop("element_dimension != 3 (is %i)", element_dimension);
		}
		if (element_order != 1) {
			stop("element_order != 1 (is %i)", element_order);
		}
		if (element_n_nodes != 4) {
			stop("element_n_nodes != 4 (is %i)", element_n_nodes);
		}
		
		std::vector<int> element_types;
		std::vector<std::vector<std::size_t>> element_tags;
		std::vector<std::vector<std::size_t>> element_node_tags;
		gmsh::model::mesh::getElements(
			element_types,
			element_tags,
			element_node_tags,
			-1
		);
		
		int element_index = -1;
		for (int i = 0; i < element_types.size(); ++i) {
			if (element_types[i] == element_type) {
				element_index = i;
				break;
			}
		}
		
		if (element_index == -1) {
			stop("element_index == -1");
		}
		
		if (element_tags[element_index].size() == 0) {
			stop("element_tags[element_index].size() == 0");
		}
		if (element_tags[element_index].size() * element_n_nodes != element_node_tags[element_index].size()) {
			stop("element_tags[element_index].size() * element_n_nodes != element_node_tags[element_index].size()");
		}
		
		std::size_t n_elements = element_tags[element_index].size();
		
		IntegerMatrix output_topology(n_elements, element_n_nodes);
		for (int i = 0; i < n_elements; ++i) {
			for (int n = 0; n < element_n_nodes; ++n) {
				output_topology(i, n) = element_node_tags[element_index][element_n_nodes * i + n];
			}
		}
		
		gmsh::clear();
		gmsh::finalize();
		
		List result;
		result["geometry"] = output_geometry;
		result["topology"] = output_topology;
		
		return result;
	} catch (...) {
		gmsh::vectorpair dim_tags;
		std::vector<std::size_t> node_tags;
		
		gmsh::model::mesh::getLastEntityError(dim_tags);
		gmsh::model::mesh::getLastNodeError(node_tags);
		
		for (auto & e : dim_tags) {
			Rcerr << "{\"dim\": " << e.first << ", \"tag\"" << e.second << "}\n";
		}
		
		gmsh::clear();
		gmsh::finalize();
		throw;
	}
}