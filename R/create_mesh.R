#' @name create_mesh
#'
#' @title Create INLA mesh for spatial modeling
#'
#' @description Generates a 2D mesh object using the convex hull of the global species presences from a processed \code{nsdm.vinput} object (output of \code{\link{NSDM.SelectCovariates}}).
#'
#' @param nsdm_obj An object of class `nsdm.vinput` as returned by `sabinaNSDM::NSDM.SelectCovariates()`.
#' @param edge A numeric vector of length 2. Maximum triangle edge lengths for the inner and outer parts of the mesh (default: c(0.5, 1)). #@@@JMB los valores por defecto son para crs en grados, habrá que ajustar cosas aquí para cuando los rasters sean metros
#' @param offset A numeric vector of length 2. Offsets to expand the domain inward and outward from the boundary (default: c(0.25, 0.5)).
#' @param buffer A numeric scalar. Amount to buffer (expand) the convex hull before mesh creation (default: 0.01).
#' @param boundary.method Character. One of `"convex_hull"` (default), `"concave_hull"` or `"raster_mask"`. #@@@JMB raster mask funciona bien para fragmentos separados y para hacer predicciones en oceano. Para esto último preguntar a Virgilio si triangulos pequeños en borde del raster es apropiado...
#' @param concavity Numeric. Required when `boundary.method = "concave_hull"`. Controls how tightly the concave hull wraps the points; lower values (e.g., 2–3) produce tighter boundaries.
#' @param remove_holes Logical. Only used when `boundary.method = "raster_mask"`; if TRUE, internal holes are removed before mesh creation (default: FALSE).
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
                        concavity = NULL,
                        remove_holes = FALSE,
                        plot = FALSE) {

  if(!inherits(nsdm_obj, "nsdm.vinput")) {
    stop("The 'nsdm_obj' must be an object of class 'nsdm.vinput'.\nUse sabinaNSDM::NSDM.SelectCovariates() to obtain it.")
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
    warning("Param concavity is ignored unless boundary.method = 'concave_hull'.")
  }
  if (boundary.method != "raster_mask" && remove_holes) {
    warning("remove_holes only applies when boundary.method = 'raster_mask'.")
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
    convex_hull = boundary_convex_hull(points_sf, buffer = buffer),
    raster_mask   = boundary_raster_mask(nsdm_obj, buffer = buffer, remove_holes = remove_holes),
    concave_hull  = boundary_concave_hull(nsdm_obj, concavity = concavity, buffer = buffer)
  )

  # # ANTIGUO enfoque Virgilio (solo presencias globales) #@@@JMB debería ser presencias/ausencias, no?
  # points_sf <- sf::st_as_sf(nsdm_obj$SpeciesData.XY.Global, coords = c("x", "y"))
  # sf::st_crs(points_sf) <- crs
  # boundary <- sf::st_convex_hull(sf::st_union(points_sf))

  # create mesh
  mesh <- fmesher::fm_mesh_2d(
    boundary = boundary,
    max.edge = 4 * edge,  #@@@JMB este 4* hay que quitarlo para la versión final. (Virgilio?). Aumenta arificialmente el tamaño los triangulos... 
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
        ggtitle(paste0("MESH SPDE (",boundary.method,")"),
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



## Different boundary_mesh

# convex hull
boundary_convex_hull <- function(points_sf, buffer) {
  boundary <- sf::st_convex_hull(sf::st_union(points_sf))
  sf::st_buffer(boundary, buffer)
}


# concave_hull
boundary_concave_hull <- function(nsdm_obj, concavity = concavity, buffer = buffer) {
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


# raster mask
boundary_raster_mask <- function(nsdm_obj, buffer = buffer, remove_holes = remove_holes) {
  old_s2 <- suppressMessages(sf::sf_use_s2())
  suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)

  # Raster global y crs
  #r_glo <- terra::rast(nsdm_obj$IndVar.Global.Selected)
  r_glo <- terra::unwrap(nsdm_obj$IndVar.Global.Selected)
  proj4  <- terra::crs(r_glo, proj = TRUE)
  crs_sf <- sf::st_crs(proj4)

  # Puntos
  pts_all <- do.call(rbind, list(
    nsdm_obj$SpeciesData.XY.Global,
    nsdm_obj$SpeciesData.XY.Regional,
    nsdm_obj$Absences.XY.Global,
    nsdm_obj$Absences.XY.Regional,
    nsdm_obj$Background.XY.Global,
    nsdm_obj$Background.XY.Regional
  ))
  pts_all <- pts_all[!sapply(pts_all, is.null), , drop = FALSE]
  pts_sf  <- st_as_sf(pts_all, coords = c("x","y"), crs = crs_sf)

  # Mask raster + puntos
  r_pts         <- terra::rasterize(terra::vect(pts_sf), r_glo, field = 1)
  mask_raster   <- !is.na(terra::app(r_glo, sum, na.rm = TRUE))
  mask_combined <- mask_raster | !is.na(r_pts)

  # Mask a sf
  poly_vec     <- terra::as.polygons(mask_combined, dissolve = TRUE)
  terra::crs(poly_vec) <- terra::crs(r_glo)
  boundary_all <- st_as_sf(poly_vec)
  boundary_all <- sf::st_make_valid(boundary_all)
  boundary_all <- sf::st_buffer(boundary_all, 0)

  # Holes internos
  extract_holes <- function(sf_poly) {
    holes <- list()
    for (g in st_geometry(sf_poly)) {
      tp <- st_geometry_type(g)
      if (tp == "POLYGON" && length(g) > 1) {
        for (i in 2:length(g)) {
          holes <- c(holes, list(st_polygon(list(g[[i]]))))
        }
      } else if (tp == "MULTIPOLYGON") {
        for (poly in g) {
          if (length(poly) > 1) {
            for (i in 2:length(poly)) {
              holes <- c(holes, list(st_polygon(list(poly[[i]]))))
            }
          }
        }
      }
    }
    if (length(holes) == 0) return(NULL)
    st_sfc(holes, crs = st_crs(sf_poly))
  }
  holes_sfc <- extract_holes(boundary_all)

  # Fragments válidos (los que tienen puntos)
  keep            <- st_intersects(boundary_all, pts_sf, sparse = FALSE)[,1]
  boundary_useful <- boundary_all[keep, ]

  # boundary para mesh
  if (remove_holes) {
    exteriors <- lapply(st_geometry(boundary_useful), function(g) {
      if (st_geometry_type(g) == "POLYGON") {
        st_polygon(list(g[[1]]))
      } else {
        rings <- lapply(g, function(poly) list(poly[[1]]))
        st_multipolygon(rings)
      }
    })
    boundary_mesh <- st_union(st_sf(geometry = st_sfc(exteriors, crs = crs_sf)))
    boundary_mesh <- st_cast(boundary_mesh, "MULTIPOLYGON")
  } else {
    boundary_mesh <- boundary_useful
  }

  # limpiar boudary
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

