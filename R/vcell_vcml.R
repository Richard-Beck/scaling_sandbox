read_vcml_model <- function(path = "biomodel.vcml", simulation = "simulation2") {
  if (!requireNamespace("xml2", quietly = TRUE)) stop("Package 'xml2' is required")
  document <- xml2::read_xml(path)
  xml2::xml_ns_strip(document)

  value <- function(xpath) xml2::xml_text(xml2::xml_find_first(document, xpath))
  attr_num <- function(xpath, attribute) {
    as.numeric(xml2::xml_attr(xml2::xml_find_first(document, xpath), attribute))
  }

  geometry_node <- xml2::xml_find_first(document, ".//SimulationSpec/Geometry")
  image_node <- xml2::xml_find_first(geometry_node, ".//ImageData")
  image_nx <- as.integer(xml2::xml_attr(image_node, "X"))
  image_ny <- as.integer(xml2::xml_attr(image_node, "Y"))
  compressed_hex <- gsub("[[:space:]]", "", xml2::xml_text(image_node))
  byte_starts <- seq.int(1L, nchar(compressed_hex), by = 2L)
  compressed <- as.raw(strtoi(substring(compressed_hex, byte_starts,
                                        byte_starts + 1L), 16L))
  pixels <- memDecompress(compressed, type = "gzip")
  if (length(pixels) != image_nx * image_ny) stop("Decoded image has unexpected size")
  image <- matrix(as.integer(pixels), nrow = image_nx, ncol = image_ny)

  simulation_xpath <- sprintf(".//Simulation[@Name='%s']", simulation)
  simulation_node <- xml2::xml_find_first(document, simulation_xpath)
  if (inherits(simulation_node, "xml_missing")) stop("Simulation not found: ", simulation)
  mesh_node <- xml2::xml_find_first(simulation_node, ".//MeshSpecification/Size")
  time_node <- xml2::xml_find_first(simulation_node, ".//TimeBound")
  output_node <- xml2::xml_find_first(simulation_node, ".//OutputOptions")

  constants <- xml2::xml_find_all(document, ".//MathDescription/Constant")
  parameter_values <- as.numeric(xml2::xml_text(constants))
  names(parameter_values) <- xml2::xml_attr(constants, "Name")

  pde_nodes <- xml2::xml_find_all(document,
    ".//MathDescription/CompartmentSubDomain[@Name='Cell']/PdeEquation")
  species <- xml2::xml_attr(pde_nodes, "Name")
  initial_names <- xml2::xml_text(xml2::xml_find_all(pde_nodes, "./Initial"))
  diffusion_names <- xml2::xml_text(xml2::xml_find_all(pde_nodes, "./Diffusion"))
  initial <- parameter_values[initial_names]
  diffusion <- parameter_values[diffusion_names]
  names(initial) <- names(diffusion) <- species

  extent_node <- xml2::xml_find_first(geometry_node, "./Extent")
  origin_node <- xml2::xml_find_first(geometry_node, "./Origin")
  list(
    name = xml2::xml_attr(xml2::xml_find_first(document, ".//BioModel"), "Name"),
    vcml_version = xml2::xml_attr(xml2::xml_find_first(document, "/vcml"), "Version"),
    geometry = list(
      extent = c(x = as.numeric(xml2::xml_attr(extent_node, "X")),
                 y = as.numeric(xml2::xml_attr(extent_node, "Y"))),
      origin = c(x = as.numeric(xml2::xml_attr(origin_node, "X")),
                 y = as.numeric(xml2::xml_attr(origin_node, "Y"))),
      image = image,
      cell_pixel_value = 1L
    ),
    mesh = c(x = as.integer(xml2::xml_attr(mesh_node, "X")),
             y = as.integer(xml2::xml_attr(mesh_node, "Y"))),
    time = c(start = as.numeric(xml2::xml_attr(time_node, "StartTime")),
             end = as.numeric(xml2::xml_attr(time_node, "EndTime")),
             output = as.numeric(xml2::xml_attr(output_node, "OutputTimeStep"))),
    species = species,
    initial = initial,
    diffusion = diffusion,
    constants = parameter_values,
    equations = c(
      J_r0 = value(".//MathDescription/Function[@Name='J_r0']"),
      J_r1 = value(".//MathDescription/Function[@Name='J_r1']"),
      J_r2 = value(".//MathDescription/Function[@Name='J_r2']"),
      J_r3 = value(".//MathDescription/Function[@Name='J_r3']"),
      J_r4 = value(".//MathDescription/Function[@Name='J_r4']"),
      J_r5 = value(".//MathDescription/Function[@Name='J_r5']"),
      Ic = value(".//MathDescription/Function[@Name='Ic']")
    )
  )
}

resample_vcell_mask <- function(model) {
  nx <- model$mesh[["x"]]
  ny <- model$mesh[["y"]]
  image <- model$geometry$image
  # VCell image samples are x-fastest. Mesh voxel centres select the nearest
  # containing source-image pixel.
  ix <- pmin(nrow(image), floor(((seq_len(nx) - 0.5) / nx) * nrow(image)) + 1L)
  iy <- pmin(ncol(image), floor(((seq_len(ny) - 0.5) / ny) * ncol(image)) + 1L)
  image[ix, iy, drop = FALSE] == model$geometry$cell_pixel_value
}
