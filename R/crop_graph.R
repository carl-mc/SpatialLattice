##############################
# CROP GRAPH WITH POLYGON
# Carl MC, 20.9.2019
##############################


#' Checks whether vertices are located in (exntent of) Spatial* object
#'
#' @description Function checks for each whether it is located within 
#' x -- if it is a polygon, its area, otherwise its spatial extent. 
#' 
#' @param g igraph object with vertex attributes x and y, encoding coordinates. 
#' @param x Spatial* Object of the same projection as the graph. 
#'
#' @return A reduced igraph object. 
#' 
#' @import sp
#' @rawNamespace import(igraph, except = c(knn, union))
#' 
#' @importFrom raster extent
#' 
#' @export
vertex_within <- function(g, x){
  # Convert x to polygon if not already a polygon
  if(!grepl("SpatialPolygon",class(x))){
    x <- as(raster::extent(x), "SpatialPolygons")
  }
  # Verteces to Spatial Points
  vertex.pts <- SpatialPoints(cbind(V(g)$x, V(g)$y))
  proj4string(vertex.pts) <- proj4string(x)
  
  # Within area
  return(apply(gContains(x,vertex.pts,  byid = T), 1, any))
}

#' Crop a Spatial Graph to (extent of) Spatial* Object
#'
#' @description Function deletes all vertices that are not located within 
#' x -- if it is a polygon, its area, otherwise its spatial extent. 
#' 
#' @param g igraph object with vertex attributes x and y, encoding coordinates. 
#' @param x Spatial* Object of the same projection as the graph. 
#' @param delete.islands Delete all disconnected components except the largest one. 
#'
#' @return A reduced igraph object. 
#' 
#' @import sp
#' raster
#' @rawNamespace import(igraph, except = c(knn, union))
#' 
#' @export
crop_graph <- function(g, x, delete.islands = F){

  
  # Subset graph
  v.within <- which(!vertex_within(g, x))
  g.sub = delete.vertices(g, v.within)
  
  # Delete islands
  if(delete.islands){
    g.components <- components(g.sub)$membership
    g.sub <- induced_subgraph(g.sub,g.components == which.max(table(g.components)))
  }
  
  # Return
  return(g.sub)
}


#' Delete graph islands
#'
#' @description Deletes all disconnected components except the largest one. 
#' @param g igraph
#'
#' @return A reduced igraph
#' @export
#' @rawNamespace import(igraph, except = c(knn, union))
delete_graph_islands <- function(g){
  g.components <- components(g)$membership
  induced_subgraph(g,g.components == which.max(table(g.components)))
}

