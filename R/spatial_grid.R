
#' Make a spatial grid
#'
#' @param shape A shape to cover. SpatialPolygons* object. 
#' @param structure Sampling and connectivity structure. 
#' One of: quadratic, hexagonal, triangular, and random. 
#' Quadratic puts vertices on a regular grid and connects each vertex its direct 4 neighbours. 
#' Hexagonal samples vertices as the centroids of hexagons and connects each vertex with its direct 6 neighbours. 
#' Triangular samples vertices as the corners of hexagons and connects each vertex with its direct 3 neighbours. 
#' Random samples vertices randomly within the shape, connecting each with its neighbours through a Delaunay triangulation. 
#' @param N_pts Number of points (approximate in structure not random)
#' @param intermediate_proj Intermediate projection to use for sampling. 
#' Can be used to sample from, e.g., an equal area Albers projection that creates spatially more regular grids.
#'
#' @return An igraph object with x and y coordinates encoded as vertex attributes. 
#' @export
spatial_grid <- function(shape,
                        structure,
                        N_pts,
                        intermediate_proj = NULL){
  if(structure == "hexagonal"){
    spatial_hexgrid(shape = shape,
                    N_pts = N_pts,
                    intermediate_proj = NULL)
  } else if(structure == "quadratic"){
    spatial_quadgrid(shape = shape,
                     N_pts = N_pts,
                     intermediate_proj = NULL)
  } else if(structure == "triangular"){
    spatial_trigrid(shape = shape,
                    N_pts = N_pts,
                    intermediate_proj = NULL)
  } else if(structure == "random"){
    spatial_randgrid(shape = shape,
                     N_pts = N_pts,
                     intermediate_proj = NULL)
  } else {
    stop("Invalid structure.")
  }
}



#' Reproject shape for sampling
#'
#' @param shape Shape
#' @param proj4string arget projection
#'
#' @return Reprojected shape
reproject_shape <- function(shape, proj4string){
  # Set projection of shape if necessary
  if(is.null(proj4string(shape)) | is.na(proj4string(shape))){
    warning("Setting projection to WGS84")
    proj4string(shape) <- "+proj=longlat +ellps=WGS84 +no_defs"
  }
  old.proj4string <- proj4string(shape)
  
  # Reproject
  if(old.proj4string != proj4string){
    shape <- spTransform(shape, CRS(proj4string))
  }
  
  # Return
  return(shape)
}


#' Backproject graph
#'
#' @param graph Igraph object with x and y coords
#' @param from source projection
#' @param to target projection
#'
#' @return Reprojected graph
reproject_graph <- function(graph, from, 
                            to = "+proj=longlat +ellps=WGS84 +no_defs",
                            coords = c("x","y")){
  # Make points
  pts <- SpatialPoints(cbind(vertex_attr(graph, coords[1]),
                             vertex_attr(graph, coords[2])),
                       proj4string = CRS(from))
  
  # Set projection of shape if necessary
  if(from != to){
    pts <- spTransform(pts, CRS(to))
  }
  
  # Reset Vertex coordinates
  vertex_attr(graph, coords[1]) <- pts@coords[,1]
  vertex_attr(graph, coords[2]) <- pts@coords[,2]
  
  # Return
  return(graph)
}

#' Delete edges that are too long on regular grid 
#'
#' @param graph graph
#' @param prec precision -- outliers that are longer and still allowed. 
#' @param max.degree Set maximum for degree. Not used currently. 
#'
#' @return A reduced graph
delete_outlier_edges <- function(graph, prec = .01, max.degree = NULL){
  # Delete long edges
  new.pts <- SpatialPoints(cbind(V(graph)$x, V(graph)$y))
  E(graph)$length <- apply(ends(graph, E(graph), names = F), 1,
                           function(e){
                             dist(new.pts@coords[e,])
                           })
  graph <- delete.edges(graph, 
                        which(E(graph)$length > min(E(graph)$length) + prec*min(E(graph)$length)))
  v.degree <- degree(graph)
  
  # Delete high degree if so wished
  if(!is.null(max.degree)){
    graph <- delete.vertices(graph, which(v.degree > max.degree))
  }
  
  
  # Return
  return(graph)
}


#' Make grid with triangular sampling and connectivity structure
#'
#' @param shape A shape to cover. SpatialPolygons* object. 
#' @param N_pts Number of points
#' @param intermediate_proj Intermediate projection to use for sampling.
#'
#' @return An igraph object with x and y coordinates encoded as vertex attributes. 
#' @export
spatial_trigrid <- function(shape, N_pts, intermediate_proj = NULL){
  
  # Reproject
  if(!is.null(intermediate_proj)){
    # Save old projection for later
    old.proj4string <- proj4string(shape)
    
    # Reproject
    shape <- reproject_shape(shape = shape, 
                             proj4string = intermediate_proj)
  }
  
  
  # Sample
  pts <- spsample(shape, 
                  n = N_pts,
                  type = "hexagonal")
  
  # Draw Network
  graph <- graph_from_pts(pts)
  
  # Delete outlier edges
  graph <- delete_outlier_edges(graph)
  
  # Triangles
  triangles <- triangles(graph)
  
  # Average point per triangle
  pts <- SpatialPoints(cbind(V(graph)$x, V(graph)$y))
  new.pts <- do.call(rbind, lapply(seq(1, length(triangles)-2, by = 3), 
                                   function(t){
                                     apply(pts@coords[triangles[t:(t+2)],], 2, mean)
                                   }))
  new.pts <- SpatialPoints(new.pts)
  
  # Make new graph
  graph <- graph_from_pts(new.pts)
  
  # Delete outlier edges
  graph <- delete_outlier_edges(graph)
  
  # Back-project graph
  if(!is.null(intermediate_proj)){
    graph <- reproject_graph(graph, from = intermediate_proj,
                             to = old.proj4string,
                             coords = c("x","y"))
  }
  
  # Check
  if(length(V(graph)) == 0){
    warning("No vertices in graph")
    return(NULL)
  }
  
  # Return
  return(graph)
}




#' Make grid with quadratic sampling and connectivity structure
#'
#' @param shape A shape to cover. SpatialPolygons* object. 
#' @param N_pts Number of points
#' @param intermediate_proj Intermediate projection to use for sampling.
#' @rawNamespace import(igraph, except = c(knn, union))
#' @return An igraph object with x and y coordinates encoded as vertex attributes. 
#' @export
spatial_quadgrid <- function(shape, N_pts, intermediate_proj = NULL){
  
  # Reproject
  if(!is.null(intermediate_proj)){
    # Save old projection for later
    old.proj4string <- proj4string(shape)
    
    # Reproject
    shape <- reproject_shape(shape = shape, 
                             proj4string = intermediate_proj)
  }
  
  
  # Sample
  pts <- spsample(shape, 
                  n = N_pts,
                  type = "regular")
  
  # Draw Network
  graph <- graph_from_pts(pts)
 
  # Delete outlier edges
  graph <- delete_outlier_edges(graph)
  
  
  # Back-project graph
  if(!is.null(intermediate_proj)){
    graph <- reproject_graph(graph, from = intermediate_proj,
                             to = old.proj4string,
                             coords = c("x","y"))
  }
  
  # Check
  if(length(V(graph)) == 0){
    warning("No vertices in graph")
    return(NULL)
  }

  # Return
  return(graph)
}


#' Make grid with hexagonal sampling and connectivity structure
#'
#' @param shape A shape to cover. SpatialPolygons* object. 
#' @param N_pts Number of points
#' @param intermediate_proj Intermediate projection to use for sampling.
#' @rawNamespace import(igraph, except = c(knn, union))
#' @return An igraph object with x and y coordinates encoded as vertex attributes. 
#' @export
spatial_hexgrid <- function(shape, N_pts, intermediate_proj = NULL){
  
  # Reproject
  if(!is.null(intermediate_proj)){
    # Save old projection for later
    old.proj4string <- proj4string(shape)
    
    # Reproject
    shape <- reproject_shape(shape = shape, 
                             proj4string = intermediate_proj)
  }
  
  
  # Sample
  pts <- spsample(shape, 
                  n = N_pts,
                  type = "hexagonal")
  
  # Draw Network
  graph <- graph_from_pts(pts)
  
  # Delete outlier edges
  graph <- delete_outlier_edges(graph)
  
  
  # Back-project graph
  if(!is.null(intermediate_proj)){
    graph <- reproject_graph(graph, from = intermediate_proj,
                             to = old.proj4string,
                             coords = c("x","y"))
  }
  
  # Check
  if(length(V(graph)) == 0){
    warning("No vertices in graph")
    return(NULL)
  }
  
  # Return
  return(graph)
}




#' Make grid with random sampling and Delaunay connectivity. 
#'
#' @param shape A shape to cover. SpatialPolygons* object. 
#' @param N_pts Number of points
#' @param intermediate_proj Intermediate projection to use for sampling.
#' @rawNamespace import(igraph, except = c(knn, union))
#' @return An igraph object with x and y coordinates encoded as vertex attributes. 
#' @export
spatial_randgrid <- function(shape, N_pts, intermediate_proj = NULL){
  if(!is.null(intermediate_proj)){
    # Save old projection for later
    old.proj4string <- proj4string(shape)
    
    # Reproject
    shape <- reproject_shape(shape = shape, 
                             proj4string = intermediate_proj)
  }
  
  
  # Sample
  pts <- spsample(shape, 
                  n = N_pts,
                  type = "random")
  
  # Draw Network
  graph <- graph_from_pts(pts)

  # Back-project graph
  if(!is.null(intermediate_proj)){
    graph <- reproject_graph(graph, from = intermediate_proj,
                             to = old.proj4string,
                             coords = c("x","y"))
  }
  
  # Check
  if(length(V(graph)) == 0){
    warning("No vertices in graph")
    return(NULL)
  }
  
  # Return
  return(graph)
}
