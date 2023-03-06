

#' Encode Polygon ID on graph
#'
#' @param graph An igraph object with an "x" and "y" coordinate column.
#' @param poly A SpatialPolygons* object. Function assumes non-overlapping polygons.
#'  Vertices which intersect with two or more polygons, are assigned an NA value. 
#' @param multi Value to set for intersections with multiple polygons. "NA" = NA, "first" = first, "last" = last
#' @return Vector of length V(graph) which encodes the numberic position of the polygon in poly vertices intersect with. 
#' @export
#' @import sp
#' @rawNamespace import(igraph, except = c(knn, union))
poly2graph <- function(graph, poly, multi = "NA"){
  int.mat <- gContains(poly, 
                       SpatialPoints(cbind(V(graph)$x, V(graph)$y)), byid = T)
  polyid <- apply(int.mat, 1, function(i){
    i <- which(i)
    if(length(i) == 1){
      i
    } else if(length(i) == 0){
      NA
    } else {
      if(multi == "NA"){
        NA
      } else if(multi == "first"){
        i[1]
      } else if(multi == "last"){
        i[length(i)]
      } else {
        stop("Provide valid multi value: NA, first, or last")
      }
    }
  })
  return(polyid)
}


# 
#' Binary edge encoding of vertex differences
#'
#' @param g Graph
#' @param attr_name vertex attribute to compute edge-level differences on
#'
#' @return Binary vector of length E(g)
#' @export
#' @rawNamespace import(igraph, except = c(knn, union))
#'
vid_2_edge_diff <- function(g, attr_name){
  e.ends <- ends(g, E(g), names = F)
  attr <- vertex_attr(g, attr_name)
  as.numeric(attr[e.ends[,1]] != attr[e.ends[,2]])
}



# 
#' Aggregate SpatialLines* Data to Spatial Graph
#'
#' @param sl SpatialLinesDataFrame to aggregate onto spatial graph
#' @param graph An igraph object with an "x" and "y" coordinate column.
#' @param aggvar Variable in sl to aggregate
#' @param aggfun Function to aggregate values if an edge crosses multiple lines. 
#' @param na.val Value to set if an edge crosses no line
#' @param only.uneven Count only intersections where an edge crosses a 
#' line an uneven number of times. 
#' This avoids ``fake intersections, where the starts and ends of the
#' edge are on the same side of the line. 
#' 
#' @rawNamespace import(igraph, except = c(knn, union))
#'
#' @return Vector of length E(g)
#' @export
lines2graph <- function(sl, graph, aggvar, aggfun, 
                       na.val = 0, only.uneven = T){
  ### Spatial Lines
  g.sl <- edges2sl(graph)
  
  ### Intersect with lines
  int <- intersect_lines(sl = g.sl, sl2 = sl, aggvar = aggvar, 
                         aggfun = aggfun, na.val = 0, only.uneven = T)
  
  ### Return
  return(int)
}


#' Extract raster values to spatial graph
#'
#' @param raster A raster* object
#' @param graph An igraph object with an "x" and "y" coordinate column.
#' @param aggfun Function to aggregate values, when edges = TRUE. Defaults to max.
#' @param edges Flag on whether to aggregate to edges or vertices. 
#' @param ... Additional parameters passed to \code{raster::extract()}
#' @rawNamespace import(igraph, except = c(knn, union))
#' @return The results of raster::extract().
#' @export 
#' @importFrom raster extract
raster2graph <- function(raster, graph, edges = TRUE, aggfun = max, ...){
  if(edges){
    ### Spatial Lines
    g.sl <- edges2sl(graph)
    
    ### Intersect with lines
    int <- raster::extract(raster, g.sl, fun = aggfun, ...)
    
  } else {
    ### Spatial Lines
    g.sp <- vertices2sp(graph)
    
    ### Intersect with lines
    int <- raster::extract(raster, g.sp, ...)
    
  }
  
  ### Return
  return(int)
}


#' Edges to SpatialLines
#'
#' @param g An igraph object with an "x" and "y" coordinate column.
#' @rawNamespace import(igraph, except = c(knn, union))
#' @return A SpatialLines object of same length and order as E(g)
#' @export
#' @import sp
edges2sl <- function(g){
  end.mat <- ends(g, E(g))
  SpatialLines(lapply(1:nrow(end.mat), function(e){
    Lines(list(Line(cbind(vertex_attr(g, name = "x", index = end.mat[e,]),
                          vertex_attr(g, name = "y", index = end.mat[e,])))), ID = as.character(e))
  }))
}

#' Vertices to SpatialPoints
#'
#' @param g An igraph object with an "x" and "y" coordinate column.
#'
#' @return A SpatialPoints object of same length and order as V(g)
#' @export
#' @import sp
vertices2sp <- function(g){
  SpatialPoints(cbind(V(g)$x, V(g)$y))
}


#' Intersect Lines with Lines
#'
#' @param sl Line to aggregate
#' @param sl2 Line that queries
#' @param aggvar Variable to aggrgate
#' @param aggfun Function for aggregation
#' @param na.val NA value
#' @param only.uneven only count uneven intersections
intersect_lines <- function(sl, sl2, aggvar, aggfun = function(x){max(x)}, na.val = 0, only.uneven = T){
  # Get decent row.names to be able to work with intersection points
  row.names(sl) <- as.character(1:length(sl))
  row.names(sl2) <- as.character(1:length(sl2))
  
  ### Intersect
  sl.int <- gIntersects(sl, sl2, byid = T, returnDense =  F)
  
  ### Count number of intersections by line x river combination
  sl.int.pts <- gIntersection(sl2, sl, byid = T)
  sl.int.pts.id <- apply(do.call(rbind, strsplit(row.names(sl.int.pts), " ")), 2, as.numeric)
  
  ### Generate aggregate of variable of sl2
  agg.vec <- unlist(lapply(1:length(sl), function(i){
    r.i <- sl.int[[i]]
    if(length(r.i) == 0){ ## No single crossing
      na.val
    } else { ## Keep only those with an *uneven* number of crossings
      r.ids <- unlist(lapply(r.i, function(j){
        if(only.uneven & sum(sl.int.pts.id[,1] == j & sl.int.pts.id[,2] == i) %% 2 == 0){
          NULL
        } else {
          j
        }
      }))
      if(length(r.ids) == 0){
        na.val
      } else {
        aggfun(sl2@data[r.ids, aggvar])
      }
    }
  }))
  
  # Return
  return(agg.vec)
}


