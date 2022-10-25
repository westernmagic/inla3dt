#' Write a mesh to an XDMF file
#' 
#' @param mesh_file
#'   Path to the output HDF5 file. Required.
#' @param xdmf_file
#'   Path to the output XDMF file. Required.
#' @param geometry
#'   The mesh geometry. Required.
#' @param topology
#'   The mesh topology. Required.
#' @param ...
#'   Further nodal or elemental attributes. Optional.
#' 
#' @importFrom assertthat assert_that
#' @importFrom hdf5r      H5File
#' @importFrom glue       glue
#' @importFrom stringr    str_c
#' @importFrom dplyr      case_when
#' @importFrom xml2       as_xml_document
#' @importFrom xml2       write_xml
#' @export
write_xdmf <- function(
	mesh_file,
	xdmf_file,
	geometry,
	topology,
	times = NULL,
	node_attributes = list(),
	cell_attributes = list()
) {
	if (is.null(node_attributes)) {
		node_attributes <- list()
	}
	if (is.null(cell_attributes)) {
		cell_attributes <- list()
	}
	
	assert_that(is.matrix(geometry))
	assert_that(is.matrix(topology))
	assert_that(is.vector(times) || is.null(times))
	assert_that(is.list(node_attributes))
	assert_that(is.list(cell_attributes))
	
	# number of nodes
	n_nodes <- dim(geometry)[1]
	# number of dimensions
	d <- dim(geometry)[2]
	# number of cells
	n_cells <- dim(topology)[1]
	# number of time steps
	n_t <- length(times)
	# nodes_per_cell
	nodes_per_cell <- dim(topology)[2]
	
	# check number of dimensions
	assert_that(d %in% c(2, 3))
	# check nodes per cell
	assert_that(nodes_per_cell %in% c(3, 4))
	# check dimension matches nodes per cell
	assert_that(d + 1 == nodes_per_cell)
	
	# check data types
	assert_that(is.double(geometry) || is.integer(geometry))
	assert_that(is.integer(topology))
	assert_that(is.double(times) || is.integer(times) || is.null(times))
	
	# check that topology only references valid nodes
	assert_that(all(1 <= topology & topology <= n_nodes))
	# check that there are no useless orphaned nodes
	# assert_that(all(1:n_nodes %in% topology))
	
	# check that no nodes overlap
	assert_that(!anyDuplicated(geometry, MARGIN = 1))
	# check that all cells are unique
	assert_that(!anyDuplicated(apply(topology, 1, sort), MARGIN = 2))
	
	if (n_t > 0) {
		# check that times are sorted
		assert_that(!is.unsorted(times))
	}
	
	if (length(node_attributes) > 0) {
		# check that all attributes are named
		assert_that(!is.null(names(node_attributes)) && !any(names(node_attributes) == ""))
		# check that all attribute names are unique
		assert_that(!anyDuplicated(names(node_attributes)))
		
		for (key in names(node_attributes)) {
			value <- node_attributes[[key]]
			
			assert_that(is.vector(value) || is.matrix(value) || is.array(value))
			assert_that(is.double(value) || is.integer(value))
			
			n <- if (is.vector(value)) {
				length(value)
			} else if (is.matrix(value)) {
				dim(value)
			} else if (is.array(value)) {
				dim(value)
			}
			
			assert_that(n[1] == n_nodes)
			
			if (n_t == 0) {
				# assert_that(is.vector(value) || is.matrix(value) || is.array(value))
				assert_that(length(n) %in% c(1, 2))
				
				if (is.array(value)) {
					assert_that(length(n) == 2)
					value <- as.matrix(value)
				}
				if (is.matrix(value)) {
					assert_that(n[2] == d)
				}
			} else {
				assert_that(is.matrix(value) || is.array(value))
				assert_that(length(n) %in% c(2, 3))
				
				if (length(n) == 2) {
					value <- as.matrix(value)
				}
				
				if (is.matrix(value)) {
					assert_that(n[2] == n_t)
				} else if (is.array(value)) {
					assert_that(n[2] == d)
					assert_that(n[3] == n_t)
				}
			}
			node_attributes[[key]] <- value
		}
	}
	
	if (length(cell_attributes) > 0) {
		# check that all attributes are named
		assert_that(!is.null(names(cell_attributes)) && !any(names(cell_attributes) == ""))
		# check that all attribute names are unique
		assert_that(!anyDuplicated(names(cell_attributes)))
		
		for (key in names(cell_attributes)) {
			value <- cell_attributes[[key]]
			
			assert_that(is.vector(value) || is.matrix(value) || is.array(value))
			assert_that(is.double(value) || is.integer(value))
			
			n <- if (is.vector(value)) {
				length(value)
			} else if (is.matrix(value)) {
				dim(value)
			} else if (is.array(value)) {
				dim(value)
			}
			
			assert_that(n[1] == n_cells)
			
			if (n_t == 0) {
				# assert_that(is.vector(value) || is.matrix(value) || is.array(value))
				assert_that(length(n) %in% c(1, 2))
				
				if (is.array(value)) {
					assert_that(length(n) == 2)
					value <- as.matrix(value)
				}
				if (is.matrix(value)) {
					assert_that(n[2] == d)
				}
			} else {
				assert_that(is.matrix(value) || is.array(value))
				assert_that(length(n) %in% c(2, 3))
				
				if (length(n) == 2) {
					value <- as.matrix(value)
				}
				
				if (is.matrix(value)) {
					assert_that(n[2] == n_t)
				} else if (is.array(value)) {
					assert_that(n[2] == d)
					assert_that(n[3] == n_t)
				}
			}
			cell_attributes[[key]] <- value
		}
	}
	
	# write mesh file
	hdf_file <- H5File$new(mesh_file, mode = "w")
	hdf_file$create_group("mesh")
	
	hdf_file[["mesh/geometry"]] <- t(geometry)
	hdf_file[["mesh/topology"]] <- t(topology - 1)
	
	if (n_t == 0) {
		if (length(node_attributes) > 0) {
			hdf_file$create_group("mesh/node_attributes")
			for (key in names(node_attributes)) {
				value <- node_attributes[[key]]
				
				if (is.vector(value)) {
					hdf_file[[glue("mesh/node_attributes/{key}")]] <- value
				} else if (is.matrix(value)) {
					hdf_file[[glue("mesh/node_attributes/{key}")]] <- t(value)
				} else if (is.array(value)) {
					hdf_file[[glue("mesh/node_attributes/{key}")]] <- aperm(value, length(dim(value)):1)
				}
			}
		}
		
		if (length(cell_attributes) > 0) {
			hdf_file$create_group("mesh/cell_attributes")
			for (key in names(cell_attributes)) {
				value <- cell_attributes[[key]]
				
				if (is.vector(value)) {
					hdf_file[[glue("mesh/cell_attributes/{key}")]] <- value
				} else if (is.matrix(value)) {
					hdf_file[[glue("mesh/cell_attributes/{key}")]] <- t(value)
				} else if (is.array(value)) {
					hdf_file[[glue("mesh/cell_attributes/{key}")]] <- aperm(value, length(dim(value)):1)
				}
			}
		}
	} else {
		for (i in 1:n_t) {
			hdf_file$create_group(glue("mesh/{times[i]}"))
		}
		
		if (length(node_attributes) > 0) {
			for (i in 1:n_t) {
				hdf_file$create_group(glue("mesh/{times[i]}/node_attributes"))
				for (key in names(node_attributes)) {
					value <- node_attributes[[key]]
					
					if (is.matrix(value)) {
						hdf_file[[glue("mesh/{times[i]}/node_attributes/{key}")]] <- as.vector(value[, i])
					} else if (is.array(value)) {
						hdf_file[[glue("mesh/{times[i]}/node_attributes/{key}")]] <- aperm(value[,, i])
					}
				}
			}
		}
		
		if (length(cell_attributes) > 0) {
			for (i in 1:n_t) {
				hdf_file$create_group(glue("mesh/{times[i]}/cell_attributes"))
				for (key in names(cell_attributes)) {
					value <- cell_attributes[[key]]
					
					if (is.matrix(value)) {
						hdf_file[[glue("mesh/{times[i]}/cell_attributes/{key}")]] <- as.vector(value[, i])
					} else if (is.array(value)) {
						hdf_file[[glue("mesh/{times[i]}/cell_attributes/{key}")]] <- aperm(value[,, i])
					}
				}
			}
		}
	}
	
	hdf_file$close_all()
	
	# write xdmf file
	xdmf_geometry <- structure(
		Type = case_when(
			d == 2 ~ "XY",
			d == 3 ~ "XYZ"
		),
		list(
			DataItem = structure(
				Format = "HDF",
				Dimensions = glue("{n_nodes} {d}"),
				NumberType = case_when(
					is.double(geometry)  ~ "Float",
					is.integer(geometry) ~ "Integer"
				),
				list(
					glue("{mesh_file}:/mesh/geometry")
				)
			)
		)
	)
	
	xdmf_topology <- structure(
		TopologyType = case_when(
			d == 2 ~ "Triangle",
			d == 3 ~ "Tetrahedron"
		),
		list(
			DataItem = structure(
				Format = "HDF",
				Dimensions = glue("{n_cells} {nodes_per_cell}"),
				NumberType = "Int",
				list(
					glue("{mesh_file}:/mesh/topology")
				)
			)
		)
	)
	
	if (n_t == 0) {
		xdmf_grid <- list(
			Geometry = xdmf_geometry,
			Topology = xdmf_topology
		)
		
		for (key in names(node_attributes)) {
			value <- node_attributes[[key]]
			
			n <- if (is.vector(value)) {
				length(value)
			} else if (is.matrix(value)) {
				dim(value)
			} else if (is.array(value)) {
				dim(value)
			}
			
			xdmf_grid <- append(
				xdmf_grid,
				list(Attribute = structure(
					Name = key,
					AttributeType = case_when(
						is.vector(value) ~ "Scalar",
						is.matrix(value) ~ "Vector"
					),
					Center = "Node",
					list(
						DataItem = structure(
							Format = "HDF",
							Dimensions = str_c(n, collapse = " "),
							NumberType = case_when(
								is.double(value)  ~ "Float",
								is.integer(value) ~ "Int"
							),
							list(
								glue("{mesh_file}:/mesh/node_attributes/{key}")
							)
						)
					)
				))
			)
		}
		
		for (key in names(cell_attributes)) {
			value <- cell_attributes[[key]]
			
			n <- if (is.vector(value)) {
				length(value)
			} else if (is.matrix(value)) {
				dim(value)
			} else if (is.array(value)) {
				dim(value)
			}
			
			xdmf_grid <- append(
				xdmf_grid,
				list(Attribute = structure(
					Name = key,
					AttributeType = case_when(
						is.vector(value) ~ "Scalar",
						is.matrix(value) ~ "Vector"
					),
					Center = "Cell",
					list(
						DataItem = structure(
							Format = "HDF",
							Dimensions = str_c(n, collapse = " "),
							NumberType = case_when(
								is.double(value)  ~ "Float",
								is.integer(value) ~ "Int"
							),
							list(
								glue("{mesh_file}:/mesh/cell_attributes/{key}")
							)
						)
					)
				))
			)
		}
		
		xdmf_domain <- list(
			Grid = structure(
				Name = "mesh",
				GridType = "Uniform",
				xdmf_grid
			)
		)
	} else {
		xdmf_grids <- list()
		for (i in 1:n_t) {
			xdmf_time <- structure(
				TimeType = "single",
				Value    = times[i],
				list()
			)
			
			xdmf_grid <- list(
				Time     = xdmf_time,
				Geometry = xdmf_geometry,
				Topology = xdmf_topology
			)
			
			for (key in names(node_attributes)) {
				value <- node_attributes[[key]]
				
				n <- if (is.vector(value)) {
					length(value)
				} else if (is.matrix(value)) {
					dim(value)
				} else if (is.array(value)) {
					dim(value)
				}
				
				xdmf_grid <- append(
					xdmf_grid,
					list(Attribute = structure(
						Name = key,
						AttributeType = case_when(
							is.matrix(value) ~ "Scalar",
							is.array(value)  ~ "Vector"
						),
						Center = "Node",
						list(
							DataItem = structure(
								Format     = "HDF",
								Dimensions = str_c(n[1:(length(n) - 1)], collapse = " "),
								NumberType = case_when(
									is.double(value)  ~ "Float",
									is.integer(value) ~ "Int"
								),
								list(
									glue("{mesh_file}:/mesh/{times[i]}/node_attributes/{key}")
								)
							)
						)
					))
				)
			}
			
			for (key in names(cell_attributes)) {
				value <- cell_attributes[[key]]
				
				n <- if (is.vector(value)) {
					length(value)
				} else if (is.matrix(value)) {
					dim(value)
				} else if (is.array(value)) {
					dim(value)
				}
				
				xdmf_grid <- append(
					xdmf_grid,
					list(Attribute = structure(
						Name = key,
						AttributeType = case_when(
							is.matrix(value) ~ "Scalar",
							is.array(value)  ~ "Vector"
						),
						Center = "Cell",
						list(
							DataItem = structure(
								Format     = "HDF",
								Dimensions = str_c(n[1:(length(n) - 1)], collapse = " "),
								NumberType = case_when(
									is.double(value)  ~ "Float",
									is.integer(value) ~ "Int"
								),
								list(
									glue("{mesh_file}:/mesh/{times[i]}/cell_attributes/{key}")
								)
							)
						)
					))
				)
			}
			
			xdmf_grids <- append(
				xdmf_grids,
				list(
					Grid = structure(
						Name = glue("mesh {times[i]}"),
						GridType = "Uniform",
						xdmf_grid
					)
				)
			)
		}
		
		xdmf_domain <- list(
			Grid = structure(
				Name = "mesh",
				GridType = "Collection",
				CollectionType = "Temporal",
				xdmf_grids
			)
		)
	}
	
	xdmf <- as_xml_document(
		list(
			Xdmf = structure(
				Version = "2.0",
				list(
					Domain = xdmf_domain
				)
			)
		)
	)
	
	write_xml(xdmf, xdmf_file, options = "format")
	
	invisible()
}
