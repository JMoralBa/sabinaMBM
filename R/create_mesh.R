#' @name create_mesh
#'
#' @title Create INLA mesh for spatial modeling
#'
#' @description Generates a 2D mesh object using the convex hull of the global species presences
#' from a processed \code{nsdm.vinput} object (output of \code{\link{NSDM.SelectCovariates}}).
#'
#' @param nsdm_obj An object of class `nsdm.vinput` as returned by `sabinaNSDM::NSDM.SelectCovariates()`.
#' @param edge A numeric vector of length 2. Maximum triangle edge lengths for the inner and outer parts of the mesh (default: `c(0.5, 1)`). #@@@JMB los valores por defecto son para crs en grados, habrá que ajustar cosas aquí para cuando los rasters sean metros
#' @param offset A numeric vector of length 2. Offsets to expand the domain inward and outward from the boundary (default: `c(0.25, 0.5)`).
#' @param buffer A numeric scalar. Amount to buffer (expand) the convex hull before mesh creation (default: `0.01`).
#' @param plot Logical. If TRUE, plots the mesh with context (presence points and subtitle). Default is FALSE.
#'
#' @return An INLA mesh object (`inla.mesh`). If `plot = TRUE`, the plot is also printed for inspection.
#' @export
create_mesh <- function(nsdm_obj, edge = c(0.5, 1), offset = c(0.25, 0.5), buffer = 0.01, plot = FALSE) {
  if (!inherits(nsdm_obj, "nsdm.vinput")) {
    stop("The 'nsdm_obj' must be an object of class 'nsdm.vinput'.\nUse NSDM.SelectCovariates() to obtain it.")
  }

  points_global <- sf::st_as_sf(nsdm_obj$SpeciesData.XY.Global, coords = c("x", "y"))
  raster_unwrapped <- terra::unwrap(nsdm_obj$IndVar.Global.Selected)
  crs <- sf::st_crs(raster_unwrapped)
  sf::st_crs(points_global) <- crs

  # Convex hull + buffer
  boundary <- sf::st_convex_hull(sf::st_union(points_global))
  boundary <- sf::st_buffer(boundary, buffer)

  mesh <- fmesher::fm_mesh_2d(
    boundary = boundary,
    max.edge = 4 * edge,
    offset = offset
  )
  fmesher::fm_crs(mesh) <- crs

  if (plot) {
    bbox_mesh <- sf::st_bbox(fm_as_sfc(mesh))
    world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    world <- sf::st_transform(world, crs)
    sf::st_crs(points_global) <- crs

    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(ggtext))
    suppressPackageStartupMessages(library(inlabru))

    print(
      ggplot() +
        geom_sf(data = world, fill = "grey92", color = "white", linewidth = 0.2) +
        gg(mesh) +
        geom_sf(data = points_global, color = "darkgreen", alpha = 0.5, size = 1.2) +
        coord_sf(
          xlim = c(bbox_mesh["xmin"], bbox_mesh["xmax"]),
          ylim = c(bbox_mesh["ymin"], bbox_mesh["ymax"])
        ) +
        ggtitle(
          "INLA MESH",
          subtitle = paste0(
            "<span style='color:darkgreen;'>Presence sp. points</span> │ ",
            "<span style='color:grey40;'>Mesh triangles</span> │ ",
            "<span style='color:blue;'>Inner mesh boundary</span> │ ",
            "<span style='color:black;'>Outer mesh boundary</span>"
          )
        ) +
        labs(x = "Longitude", y = "Latitude") +
        theme_minimal(base_size = 11) +
        theme(
          legend.position = "none",
          plot.subtitle = ggtext::element_markdown(size = 10)
        )
    )
  }

  return(mesh)
}
