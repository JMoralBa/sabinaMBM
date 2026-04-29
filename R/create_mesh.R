#' @name create_mesh
#'
#' @title Create INLA mesh for spatial modeling
#'
#' @description Generates a 2D mesh object using the convex hull, concave_hull, or raster_mask of the species presences/absences (plus optional new scenario rasters), from a processed \code{nsdm.vinput} object (output of \code{\link{NSDM.SelectCovariates}}).
#'
#' @param nsdm_obj An object of class `nsdm.vinput` as returned by `sabinaNSDM::NSDM.SelectCovariates()`.
#' @param edge A numeric vector of length 2. Maximum triangle edge lengths for the inner and outer mesh domains (default: \code{c(0.5, 1)}). Units must match the CRS of the input rasters (degrees for geographic CRS, metres for projected CRS). #@@@JMB los valores por defecto son para crs en grados, habrá que ajustar cosas aquí para cuando los rasters sean metros
#' @param offset A numeric vector of length 2. Offsets to expand the domain inward and outward from the boundary (default: c(0.25, 0.5)).
#' @param buffer A numeric scalar. Amount to buffer (expand) the convex hull before mesh creation (default: 0.01).
#' @param boundary.method Character. One of \code{"convex_hull"} (default), \code{"concave_hull"}, or \code{"raster_mask"}. \code{"raster_mask"} is recommended for fragmented distributions or when prediction extends to coastal/island areas. #@@@JMB raster mask funciona bien para fragmentos separados y para hacer predicciones en oceano. Para esto último preguntar a Virgilio si triangulos pequeños en borde del raster es necesario...
#' @param concavity Numeric. Required when `boundary.method = "concave_hull"`. Controls how tightly the concave hull wraps the points; lower values (e.g., 2–3) produce tighter boundaries.
#' @param remove_holes Logical. Only used when `boundary.method = "raster_mask"`; if TRUE, internal holes are removed before mesh creation (default: FALSE).
#' @param proj.new.env Logical. If TRUE, includes the extent of new scenarios in the mesh domain if they exist in `nsdm_obj$Scenarios` (default: TRUE). Only applies if \code{boundary.method = "raster_mask"}. 
#' @param plot Logical. If TRUE, plots the mesh for inspection (default: FALSE).
#'
#' @return An \code{fm_mesh_2d} object (class \code{inla.mesh}) suitable for use with \code{INLA} and \code{inlabru} SPDE components. Pass directly to \code{NSBM.pure(spde.mesh = ...)}. #@@@JMB pensar si ponemos clase propia a este objeto (inla.mesh??)
#'
#' @seealso \code{\link{NSBM.pure}}
#'
#' @references
#' Lindgren, F., Rue, H. & Lindström, J. (2011). An explicit link between Gaussian fields
#' and Gaussian Markov random fields: the stochastic partial differential equation approach.
#' \emph{Journal of the Royal Statistical Society B}, 73, 423–498.
#'
#' @export
create_mesh <- function(nsdm_obj, 
                        edge = c(0.5, 1),
                        offset = c(0.25, 0.5),
                        buffer = 0.01,
                        boundary.method = "convex_hull", #"convex_hull", "concave_hull" o "raster_mask"
                        concavity = NULL,
                        remove_holes = FALSE,
                        proj.new.env = TRUE,
                        plot = FALSE) {

  if(!inherits(nsdm_obj, "nsdm.vinput")) {
    stop("The 'nsdm_obj' must be an object of class 'nsdm.vinput'.\nPlease, use sabinaNSDM::NSDM.SelectCovariates() to obtain it.")
  }
  if(!boundary.method %in% c("convex_hull", "raster_mask", "concave_hull")) {
    stop("Invalid 'boundary.method'. Please, select 'convex_hull', 'raster_mask' or 'concave_hull'.")
  }
  if(boundary.method == "concave_hull") {
    if(is.null(concavity)) {
      stop("boundary.method = 'concave_hull' requires a 'concavity' value.\n",
        "  Example: concavity = 3  # lower values produce tighter boundaries.")
    }
  } else if(!is.null(concavity)) {
    warning("concavity is ignored unless boundary.method = 'concave_hull'.")
  }
  if (boundary.method != "raster_mask" && remove_holes) {
    warning("remove_holes only applies when boundary.method = 'raster_mask'.")
  }
  if (proj.new.env && boundary.method != "raster_mask") {
    warning(
      "proj.new.env = TRUE is ignored unless boundary.method = 'raster_mask'; ",
      "'convex_hull' and 'concave_hull' derive their boundary from presence/absence/backgorund points only."
    )
  }

  # Prepare points
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

  points_sf <- dplyr::bind_rows(ap_sf, pp_sf)

  # boundary method
  boundary <- switch(boundary.method,
    convex_hull = { boundary_convex_hull(points_sf, buffer = buffer) },
    concave_hull = { boundary_concave_hull(nsdm_obj, concavity = concavity, buffer = buffer) },
    raster_mask = { 
    b <- boundary_raster_mask(nsdm_obj, buffer = buffer, remove_holes = remove_holes)
      if(proj.new.env && !is.null(nsdm_obj$Scenarios)) {
        b <- add_new_scenario_boundary(b, nsdm_obj$Scenarios, remove_holes = remove_holes)
      }
      b
    }
  )

  # create mesh
  mesh <- fmesher::fm_mesh_2d(
    boundary = boundary,
    max.edge = edge, #4 * edge,  #@@@JMB este 4* hay que quitarlo para la versión final. (Virgilio?). Aumenta arificialmente el tamaño los triangulos... 
    offset = offset
  )
  fmesher::fm_crs(mesh) <- crs

  # plot mesh
  if(plot) {
    bbox_mesh <- sf::st_bbox(fmesher::fm_as_sfc(mesh))
    world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    world <- sf::st_transform(world, crs)
    sf::st_crs(points_sf) <- crs

    print(
      ggplot2::ggplot() +
        ggplot2::geom_sf(data = world, fill = "grey92", color = "white",
                         linewidth = 0.2) +
        inlabru::gg(mesh) +
        ggplot2::geom_sf(data = points_sf,
                         ggplot2::aes(color = type),
                         alpha = 0.5, size = 1.2) +
        ggplot2::scale_color_manual(
          values = c("Presence" = "darkgreen", "Absence" = "red")) +
        ggplot2::coord_sf(
          xlim = c(bbox_mesh["xmin"], bbox_mesh["xmax"]),
          ylim = c(bbox_mesh["ymin"], bbox_mesh["ymax"])
        ) +
        ggplot2::ggtitle(
          paste0("MESH SPDE (", boundary.method, ")"),
          subtitle = paste0(
            "<span style='color:darkgreen;'>Presences</span> │ ",
            "<span style='color:red;'>Absences/Background</span> │ ",
            "<span style='color:grey40;'>Mesh triangles</span> │ ",
            "<span style='color:blue;'>Inner mesh boundary</span> │ ",
            "<span style='color:black;'>Outer mesh boundary</span>"
          )
        ) +
        ggplot2::labs(x = "Longitude", y = "Latitude") +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          legend.position = "none",
          plot.title = ggplot2::element_text(size = 10),
          plot.subtitle = ggtext::element_markdown(size = 9)
        )
    )
  }

  return(mesh)
}



## Different boundary_mesh
#-------------------------
## convex hull
boundary_convex_hull <- function(points_sf, buffer) {
  old_s2 <- suppressMessages(sf::sf_use_s2())
  suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)
  boundary <- sf::st_convex_hull(sf::st_union(points_sf))
  sf::st_buffer(boundary, buffer)
}


## concave_hull
boundary_concave_hull <- function(nsdm_obj, concavity = concavity, buffer = buffer) {
  # Raster global
  r_glo <- terra::unwrap(nsdm_obj$IndVar.Global.Selected)
  crs <- sf::st_crs(r_glo)

  # puntos global + regional
  all_points <- do.call(rbind, Filter(Negate(is.null), list(
    nsdm_obj$SpeciesData.XY.Global,
    nsdm_obj$SpeciesData.XY.Regional,
    nsdm_obj$Absences.XY.Global,
    nsdm_obj$Absences.XY.Regional,
    nsdm_obj$Background.XY.Global,
    nsdm_obj$Background.XY.Regional
  )))
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


## raster mask
boundary_raster_mask <- function(nsdm_obj, buffer = buffer, remove_holes = remove_holes) {
  old_s2 <- suppressMessages(sf::sf_use_s2())
  suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)

  # Raster global y crs
  r_glo <- terra::unwrap(nsdm_obj$IndVar.Global.Selected)
  proj4 <- terra::crs(r_glo, proj = TRUE)
  crs_sf <- sf::st_crs(proj4)

  # Points
  pts_all <- do.call(rbind, Filter(Negate(is.null), list(
    nsdm_obj$SpeciesData.XY.Global,
    nsdm_obj$SpeciesData.XY.Regional,
    nsdm_obj$Absences.XY.Global,
    nsdm_obj$Absences.XY.Regional,
    nsdm_obj$Background.XY.Global,
    nsdm_obj$Background.XY.Regional
  )))
  pts_sf  <- sf::st_as_sf(pts_all, coords = c("x","y"), crs = crs_sf)

  # Mask raster + points
  r_pts <- terra::rasterize(terra::vect(pts_sf), r_glo, field = 1)
  mask_raster <- !is.na(terra::app(r_glo, sum, na.rm = TRUE))
  mask_combined <- mask_raster | !is.na(r_pts)
  mask_combined[!mask_combined] <- NA 

  # Mask a sf
  poly_vec <- terra::as.polygons(mask_combined, dissolve = TRUE)
  terra::crs(poly_vec) <- terra::crs(r_glo)
  boundary_all <- sf::st_as_sf(poly_vec)
  #boundary_all <- sf::st_make_valid(boundary_all)

  boundary_all <- sf::st_buffer(boundary_all, 0)

  # Holes internos
  extract_holes <- function(sf_poly) {
    holes <- list()
    for (g in sf::st_geometry(sf_poly)) {
      tp <- sf::st_geometry_type(g)
      if (tp == "POLYGON" && length(g) > 1) {
        for (i in 2:length(g)) {
          holes <- c(holes, list(sf::st_polygon(list(g[[i]]))))
        }
      } else if (tp == "MULTIPOLYGON") {
        for (poly in g) {
          if (length(poly) > 1) {
            for (i in 2:length(poly)) {
              holes <- c(holes, list(sf::st_polygon(list(poly[[i]]))))
            }
          }
        }
      }
    }
    if (length(holes) == 0) return(NULL)
    st_sfc(holes, crs = sf::st_crs(sf_poly))
  }
  holes_sfc <- extract_holes(boundary_all)

  # Valid polys (those with points)
  keep <- lengths(sf::st_intersects(boundary_all, pts_sf)) > 0
  boundary_useful <- boundary_all[keep, ]

  # boundary para mesh
  if (remove_holes) {
    boundary_mesh <- remove_inner_holes(boundary_useful)
  } else {
    boundary_mesh <- boundary_useful
  }

  # Clean boundary
  #boundary_mesh <- sf::st_make_valid(boundary_mesh)
  #if(!sf::st_is_longlat(boundary_mesh)) {
  #  boundary_mesh <- sf::st_buffer(boundary_mesh, 0)
  #}
  boundary_mesh <- sf::st_union(sf::st_make_valid(boundary_mesh))
  pixel_res <- min(terra::res(r_glo)) / 2
  boundary_mesh <- sf::st_simplify(boundary_mesh, dTolerance = pixel_res)
  boundary_mesh <- sf::st_make_valid(boundary_mesh)
  if(!sf::st_is_longlat(boundary_mesh)) {
    boundary_mesh <- sf::st_buffer(boundary_mesh, 0)
  }

  # buffer
  if (buffer > 0) {
    boundary_mesh <- sf::st_buffer(boundary_mesh, dist = buffer)
  }

  boundary <- boundary_mesh
  return(boundary)
}


# remove inner holes if remove_holes = TRUE
remove_inner_holes <- function(multi_poly_sf) {
  exteriors <- lapply(
    sf::st_geometry(multi_poly_sf),
    function(g) {
      if (sf::st_geometry_type(g) == "POLYGON") {
        sf::st_polygon(list(g[[1]]))
      } else {
        rings <- lapply(g, function(poly) list(poly[[1]]))
        sf::st_multipolygon(rings)
      }
    }
  )
  result <- sf::st_union(
    sf::st_sf(geometry = sf::st_sfc(exteriors, crs = sf::st_crs(multi_poly_sf)))
  ) |>
    sf::st_cast("MULTIPOLYGON") |>
    sf::st_make_valid()
  return(result)
}


# Add boundary of new scenariaos (pixels no-NA) if peoj.new.env = TRUE to extend the mesh if necessary
add_new_scenario_boundary <- function(boundary, scenario_list, remove_holes = FALSE) {
  old_s2 <- suppressMessages(sf::sf_use_s2())
  suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)
  for(r_scen in scenario_list) {
    r_scen <- terra::unwrap(r_scen)

    # Mask raster scenario
    m_scen <- !is.na(r_scen[[1]])
    m_scen[!m_scen] <- NA # only TRUE create polygon

    # Dissolve boundary
    poly_scen <- terra::as.polygons(m_scen, dissolve = TRUE, values = FALSE)
    terra::crs(poly_scen) <- terra::crs(r_scen)
    poly_scen_sf <- sf::st_as_sf(poly_scen) |>
      sf::st_make_valid() |>
      sf::st_cast("MULTIPOLYGON")

    # rm holes
    if (remove_holes) {
      poly_scen_sf <- remove_inner_holes(poly_scen_sf)
    }

    # CRS
    crs_bound   <- sf::st_crs(boundary)
    poly_scen_sf <- sf::st_transform(poly_scen_sf, crs_bound)

    # when scenario is within current boundary, next
    if(all(sf::st_within(poly_scen_sf, boundary, sparse = FALSE))) {
      next
    }

    # Extract only new pixels
    new_frag <- sf::st_difference(poly_scen_sf, boundary) |>
      sf::st_make_valid()

    # merge new pixels + boundary
    if (!all(sf::st_is_empty(new_frag))) {
      boundary <- sf::st_union(boundary, new_frag) |>
                  sf::st_make_valid()
    }
  }
  return(boundary)
}

