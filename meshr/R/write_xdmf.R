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
#' @param open
#'   Whether to open the result in ParaView, if available.
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
	open = TRUE,
	...
) {
	assert_that(is.matrix(geometry))
	assert_that(is.matrix(topology))
	
	# number of nodes
	n_nodes <- dim(geometry)[1]
	# number of dimensions
	d <- dim(geometry)[2]
	# number of cells
	n_cells <- dim(topology)[1]
	# nodes_per_cell
	nodes_per_cell <- dim(topology)[2]
	# capture attributes
	values <- list(...)
	
	# check number of dimensions
	assert_that(d %in% c(2, 3))
	# check nodes per cell
	assert_that(nodes_per_cell %in% c(3, 4))
	# check dimension matches nodes per cell
	assert_that(d + 1 == nodes_per_cell)
	
	# check data types
	assert_that(is.double(geometry) || is.integer(geometry))
	assert_that(is.integer(topology))
	
	# check that topology only references valid nodes
	assert_that(all(1 <= topology & topology <= n_nodes))
	# check that there are no useless orphaned nodes
	# assert_that(all(1:n_nodes %in% topology))
	
	# check that no nodes overlap
	assert_that(!anyDuplicated(geometry, MARGIN = 1))
	# check that all cells are unique
	assert_that(!anyDuplicated(apply(topology, 1, sort), MARGIN = 2))
	
	if (length(values) > 0) {
		# check that we can differentiate between node and cell attributes
		assert_that(n_nodes != n_cells)
		# check that all attributes are named
		assert_that(!is.null(names(values)) && !any(names(values) == ""))
		
		for (key in names(values)) {
			value <- values[[key]]
			
			assert_that(is.vector(value) || is.matrix(value))
			assert_that(is.double(value) || is.integer(value))
			
			if (is.matrix(value)) {
				assert_that(length(dim(values)) == 2)
				assert_that(dim(value)[1] %in% c(n_nodes, n_cells))
				
				# normalize 1-d matrices to vectors
				if (dim(value)[2] == 1) {
					value <- as.vector(value)
					values[key] <- value
				} else {
					assert_that(dim(value)[2] == d)
				}
			} else if (is.vector(value)) {
				assert_that(length(value) %in% c(n_nodes, n_cells))
			}
		}
	}
	
	# write mesh file
	hdf_file <- H5File$new(mesh_file, mode = "w")
	hdf_file$create_group("mesh")
	hdf_file[["mesh/geometry"]] <- t(geometry)
	hdf_file[["mesh/topology"]] <- t(topology - 1)
	hdf_file$create_group("mesh/values")
	for (key in names(values)) {
		value <- values[[key]]
		if (is.matrix(value)) {
			hdf_file[[glue("mesh/values/{key}")]] <- t(value)
		} else if (is.vector(value)) {
			hdf_file[[glue("mesh/values/{key}")]] <- value
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
	
	xdmf_grid <- list(
		Geometry = xdmf_geometry,
		Topology = xdmf_topology
	)
	
	for (key in names(values)) {
		value <- values[[key]]
		assert_that(is.vector(value) || is.matrix(value))
		
		if (is.vector(value)) {
			n <- length(value)
		} else if (is.matrix(value)) {
			n <- dim(value)
		}
		
		assert_that(n[1] == n_nodes || n[1] == n_cells)
		
		if (is.matrix(value)) {
			assert_that(length(n) == 2)
			assert_that(n[2] == d)
		}
		
		v <- structure(
			Name = key,
			AttributeType = case_when(
				is.vector(value) ~ "Scalar",
				is.matrix(value) ~ "Vector"
			),
			Center = case_when(
				n[1] == n_nodes ~ "Node",
				n[1] == n_cells ~ "Cell"
			),
			list(
				DataItem = structure(
					Format = "HDF",
					Dimensions = str_c(n, collapse = " "),
					NumberType = case_when(
						is.double(value)  ~ "Float",
						is.integer(value) ~ "Int"
					),
					list(
						glue("{mesh_file}:/mesh/values/{key}")
					)
				)
			)
		)
		
		xdmf_grid <- append(
			xdmf_grid,
			list(Attribute = v)
		)
	}
	
	xdmf <- as_xml_document(
		list(
			Xdmf = structure(
				Version = "2.0",
				list(
					Domain = list(
						Grid = structure(
							Name = "mesh",
							GridType = "Uniform",
							xdmf_grid
						)
					)
				)
			)
		)
	)
	
	write_xml(xdmf, xdmf_file, options = "format")
	
	paraview <- Sys.which("paraview")
	if (open && paraview != "") {
		system2(paraview, xdmf_file, wait = FALSE)
	}
	
	invisible()
}
