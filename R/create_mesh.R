#' @name create_mesh
#'
#' @title Create INLA mesh for spatial modeling
#'
#' @description Generates a 2D mesh object using the convex hull of the global species presences
#' from a processed \code{nsdm.vinput} object (output of \code{\link{NSDM.SelectCovariates}}).
#'
#' @param nsdm_obj An object of class `nsdm.vinput` as returned by `sabinaNSDM::NSDM.SelectCovariates()`.
#' @param edge A numeric vector of length 2. Maximum triangle edge lengths for the inner and outer parts of the mesh (default: c(0.5, 1)). #@@@JMB los valores por defecto son para crs en grados, habrá que ajustar cosas aquí para cuando los rasters sean metros
#' @param offset A numeric vector of length 2. Offsets to expand the domain inward and outward from the boundary (default: c(0.25, 0.5)).
#' @param buffer A numeric scalar. Amount to buffer (expand) the convex hull before mesh creation (default: 0.01).
#' @param boundary.method Character. Either `"convex_hull"` (default), "concave_hull" or `"raster_mask"`. #@@@JMB in process
#' @param plot Logical. If TRUE, plots the mesh for inspection (default: FALSE).
#'
#' @return An INLA mesh object (`inla.mesh`).
#'
#' @export
create_mesh <- function(nsdm_obj, 
                        edge = c(0.5, 1),
                        offset = c(0.25, 0.5),
                        buffer = 0.01,
                        boundary.method = "convex_hull", #"convex_hull", "concave_hull" o "raster_mask"
                        concavity =3, # ajusta envolvente cuando concave_hull
                       plot = FALSE) {

  if(!inherits(nsdm_obj, "nsdm.vinput")) {
    stop("The 'nsdm_obj' must be an object of class 'nsdm.vinput'.\nUse sabinaNSDM::NSDM.SelectCovariates() to obtain it.")
  }
  if(!boundary.method %in% c("convex_hull", "raster_mask", "concave_hull")) {
    stop("Invalid 'boundary.method'. Please, select 'convex_hull', 'raster_mask' or 'concave_hull' .")
  }

  # Prepare data
  pp <- do.call(rbind, list(
    nsdm_obj$SpeciesData.XY.Global,
    nsdm_obj$SpeciesData.XY.Regional
  ))
  pp_sf <- sf::st_as_sf(pp, coords = c("x", "y"))
  pp_sf$type <- "Presence"

  ap <- list(
    nsdm_obj$Absences.XY.Global,
    nsdm_obj$Absences.XY.Regional,
    nsdm_obj$Background.XY.Global,
    nsdm_obj$Background.XY.Regional
  )
  ap <- ap[!sapply(ap, is.null)]

  ap_sf <- if (length(ap) > 0) {
    absences <- do.call(rbind, ap)
    temp_sf <- sf::st_as_sf(absences, coords = c("x", "y"))
    temp_sf$type <- "Absence"
    temp_sf
  } else {
    NULL
  }

  # CRS
  raster_unwrapped <- terra::unwrap(nsdm_obj$IndVar.Global.Selected)
  crs <- sf::st_crs(raster_unwrapped)
  sf::st_crs(pp_sf) <- crs
  if(!is.null(ap_sf)) sf::st_crs(ap_sf) <- crs

  points_sf <- dplyr::bind_rows(pp_sf, ap_sf)

  # boudary method
  boundary <- switch(boundary.method,
    convex_hull = boundary_convex_hull(points_sf, buffer),
    raster_mask = boundary_raster_mask(nsdm_obj, buffer),
    concave_hull  = boundary_concave_hull(nsdm_obj, concavity = 2, buffer = buffer)
  )

  # # ANTIGUO enfoque (solo presencias globales):
  # points_sf <- sf::st_as_sf(nsdm_obj$SpeciesData.XY.Global, coords = c("x", "y"))
  # sf::st_crs(points_sf) <- crs
  # boundary <- sf::st_convex_hull(sf::st_union(points_sf))

  # create mesh
  mesh <- fmesher::fm_mesh_2d(
    boundary = boundary,
    max.edge = 4 * edge,
    offset = offset
  )
  fmesher::fm_crs(mesh) <- crs

  # plot mesh
  if (plot) {
    bbox_mesh <- sf::st_bbox(fm_as_sfc(mesh))
    world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    world <- sf::st_transform(world, crs)
    sf::st_crs(points_sf) <- crs

    print(
      ggplot() +
        geom_sf(data = world, fill = "grey92", color = "white", linewidth = 0.2) +
        gg(mesh) +
        geom_sf(data = points_sf, aes(color = type), alpha = 0.5, size = 1.2) +
        scale_color_manual(values = c("Presence" = "darkgreen", "Absence" = "red")) +
        coord_sf(
          xlim = c(bbox_mesh["xmin"], bbox_mesh["xmax"]),
          ylim = c(bbox_mesh["ymin"], bbox_mesh["ymax"])
        ) +
        ggtitle(
          "INLA MESH",
          subtitle = paste0(
            "<span style='color:darkgreen;'>Presences</span> │ ",
            "<span style='color:red;'>Absences/Background</span> │ ",
            "<span style='color:grey40;'>Mesh triangles</span> │ ",
            "<span style='color:blue;'>Inner mesh boundary</span> │ ",
            "<span style='color:black;'>Outer mesh boundary</span>"
          )
        ) +
        labs(x = "Longitude", y = "Latitude") +
        theme_minimal(base_size = 11) +
        theme(
          legend.position = "none",
          plot.title = element_text(size = 10),
          plot.subtitle = ggtext::element_markdown(size = 9)
        )
    )
  }

  return(mesh)
}


## Different mesh
# mesh convex hull
boundary_convex_hull <- function(points_sf, buffer) {
  boundary <- sf::st_convex_hull(sf::st_union(points_sf))
  sf::st_buffer(boundary, buffer)
}


# mesh concave_hull
boundary_concave_hull <- function(nsdm_obj, concavity, buffer) {
  # Raster global
  r_glo <- terra::unwrap(nsdm_obj$IndVar.Global.Selected)
  crs <- sf::st_crs(r_glo)

  # puntos global + regional
  all_points <- do.call(rbind, list(
    nsdm_obj$SpeciesData.XY.Global,
    nsdm_obj$SpeciesData.XY.Regional,
    nsdm_obj$Absences.XY.Global,
    nsdm_obj$Absences.XY.Regional,
    nsdm_obj$Background.XY.Global,
    nsdm_obj$Background.XY.Regional
  ))
  all_points <- all_points[!sapply(all_points, is.null), , drop = FALSE]
  if (nrow(all_points) == 0) stop("No points available.")

  # to sf
  pts_sf <- sf::st_as_sf(all_points, coords = c("x", "y"), crs = crs)

  # Filtrar puntos no-NA
  vals <- terra::extract(r_glo, all_points[, c("x", "y")], ID = FALSE)
  valid <- apply(vals, 1, function(row) any(!is.na(row)))
  pts_valid <- pts_sf[valid, ]
  if (nrow(pts_valid) == 0) stop("No valid points found over raster.")

  # Concave hull
  boundary <- concaveman::concaveman(pts_valid, concavity = concavity)

  # Buffer
  if (buffer > 0) {
    boundary <- sf::st_buffer(boundary, dist = buffer)
  }

  boundary <- sf::st_make_valid(boundary)
  boundary <- sf::st_cast(boundary, "MULTIPOLYGON")

  return(boundary)
}


