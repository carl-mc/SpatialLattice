#' Polygons to spatial graph
#' 
#' @description Transforms a set of polygons to a spatial graph. 
#' Vertices are encodes as polygons centroids (if pts are not supplied) and 
#' edges are drawn on the basis of the adjacency matrix of polygons (based on \code{sp::gTouches})
#'
#' @param pts Optional SpatialPoints* object to set different x and y coordinates for graph vertices. 
#'
#' @return An igraph object with an "x" and "y" coordinate attribute and all 
#' polygon-level data transferred to vertex attributes. 
#' 
#' @export
#' @import
#' rgeos
#' Matrix
#' @rawNamespace import(igraph, except = c(knn, union))
graph_from_poly <- function(poly, pts = NULL, plot.result = F) {

  # Use centroids if no points given
  if(is.null(pts)){
    pts <- gCentroid(poly, byid = T)
  }
  
  # Neighborhood matrix
  # Get neighbour list
  ngb.ls <- gTouches(poly, byid = T, returnDense = F)
  
  # Make sparse matrix
  ngb.mat <- sparseMatrix(i = rep(1:length(ngb.ls), unlist(lapply(ngb.ls, length))),
                          j = unlist(ngb.ls),
                          x = T)
  
  # Make Graph From ngb matrix
  graph <- graph_from_adjacency_matrix(ngb.mat, mode = "undirected", weighted = NULL, 
                                       diag = F)
  
  # Add vertex attributes
  vertex_attr(graph, name = "x") <- coordinates(pts)[,1]
  vertex_attr(graph, name = "y") <- coordinates(pts)[,2]
  
  # Add other data
  if(class(poly) == "SpatialPolygonsDataFrame"){
    for(v in colnames(poly@data)){
      vertex_attr(graph, name = v) <- poly@data[,v]
    }
  }
  
  # Plot
  if(plot.result){
    asp <- (max( coordinates(pts)[, 2]) - min( coordinates(pts)[, 
                                                                2]))/(max( coordinates(pts)[, 1]) - min( coordinates(pts)[, 
                                                                                                                          1]))
    plot(graph, vertex.size = 0.25, vertex.label = NA, 
         layout = cbind(vertex_attr(graph, name = "x"), vertex_attr(graph, name = "y")), 
         asp = asp)
  }
  
  # Return
  return(graph)
}



#' Points to spatial graph
#' 
#' @description Transforms a set of points to a spatial graph via a Delaunay triangulation. 
#'
#' @param points SpatialPoints* object
#'
#' @return An igraph object with an "x" and "y" coordinate attribute and all 
#' point-level data transferred to vertex attributes. 
#' @export
#' @import BoostLines
#' rgeos
#'
#' @rawNamespace import(igraph, except = c(knn, union))
#' @importFrom rtree knn RTree
graph_from_pts <- function(points){

  # Remove duplicate points
  dupl <- duplicated(points@coords)
  if(any(dupl)){
    print(paste("Removing" , sum(dupl), "duplicate points. "))
    points <- points[!dupl,]
  }
  
  # Delauney Triangulation
  del.tri = rgeos::gDelaunayTriangulation(points, onlyEdges = T)
  
  # remove multilines
  del.tri <- BoostLines::remove_multilines(del.tri)
  
  # To BoostLines
  bl <- BoostLines(del.tri)
  
  # Lines to igraph
  g <- boost2graph(bl, df = data.frame(lid = 1:length(del.tri)),
                   lonlat = T, plot = F)
  
  # Check
  stopifnot(length(V(g)) == length(points))
  
  # Add vertex info
  lat.df <- cbind(x = V(g)$x, y = V(g)$y)
  rt <- rtree::RTree(lat.df)
  ngb <- unlist(rtree::knn(rt, y = points@coords, k = as.integer(1)))
  stopifnot(length(ngb) == length(unique(ngb)))
  
  if(.hasSlot(points, "data")){
    for(v in colnames(points@data)[!colnames(points@data) %in% c("x","y")]){
      vertex_attr(g, name = v, index = ngb) <- points@data[,v]
    }
  }
  
  
  # Return
  return(g)
}
