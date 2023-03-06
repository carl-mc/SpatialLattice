#' Polygons to spatial graph
#' 
#' @description Transforms a set of polygons to a spatial graph. 
#' Vertices are encodes as polygons centroids (if pts are not supplied) and 
#' edges are drawn on the basis of the adjacency matrix of polygons (based on \code{sp::gTouches})
#'
#' @param poly SpatialPolygons* object to transform into spatial igraph.  
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
graph_from_poly <- function(poly, pts = NULL) {

  # Use centroids if no points given
  if(is.null(pts)){
    pts <- gCentroid(poly, byid = TRUE)
  }
  
  # Neighborhood matrix
  # Get neighbour list
  ngb.ls <- gTouches(poly, byid = TRUE, returnDense = FALSE)
  
  # Make sparse matrix
  ngb.mat <- sparseMatrix(i = rep(1:length(ngb.ls), unlist(lapply(ngb.ls, length))),
                          j = unlist(ngb.ls),
                          x = TRUE)
  
  # Make Graph From ngb matrix
  graph <- graph_from_adjacency_matrix(ngb.mat, mode = "undirected", weighted = NULL, 
                                       diag = F)
  
  # Add vertex attributes
  vertex_attr(graph, name = "x") <- coordinates(pts)[,1]
  vertex_attr(graph, name = "y") <- coordinates(pts)[,2]
  
  # Add other data
  if(inherits(poly, "SpatialPolygonsDataFrame")){
    for(v in colnames(poly@data)){
      vertex_attr(graph, name = v) <- poly@data[,v]
    }
  }
  
  # Return
  return(graph)
}



#' Eliminate Multilines
#'
#' Removes all multilines from a SpatialLinesDataFrame by replacing them with simple lines.
#' Attribute data for multilines is duplicated accordingly. From BoostLines package by Philipp Hunziker. 
#'
#' @param s A SpatialLinesDataFrame or SpatialLines object.
#'
#' @return A SpatialLinesDataFrame or SpatialLines object.
#'
#' @import sp
#'
remove_multilines <- function(s) {
  
  if (inherits(s, 'SpatialLines') | inherits(s, 'SpatialLinesDataFrame')) {
    lines.ls <- s@lines
    multiline.count <- unlist(lapply(lines.ls, function(x) length(x@Lines)))
    exploded_lines.id <- rep(1:length(lines.ls), times=multiline.count)
    exploded_lines.ls <- lapply(lines.ls, function(x) {
      this.Lines.ls <- x@Lines
      this.id <- x@ID
      cnt <- 0
      this.lines.ls <- lapply(this.Lines.ls, function(x) {
        cnt <<- cnt + 1;
        Lines(list(x), paste0(this.id, '.', cnt));
      })
      return(this.lines.ls);
    })
    exploded_lines.ls <- unlist(exploded_lines.ls)
    exploded.sl <- SpatialLines(exploded_lines.ls, proj4string = s@proj4string)
    
    if (inherits(s, 'SpatialLinesDataFrame')) {
      exploded.sldf <- SpatialLinesDataFrame(exploded.sl, s@data[exploded_lines.id,], FALSE)
      out <- exploded.sldf
    } else {
      out <- exploded.sl
    }
    
  } else {
    stop('s is not a supported object.')
  }
  
  return(out)
}

#' Generates an igraph graph from a SpatialLines* object.
#'
#' @details
#' Generates an igraph graph from a SpatialLines* object.
#' Line attribute data is transformed into edge attribute data
#'
#'
#' @param x A SpatialLines* object.
#' @return A \code{igraph::graph} object.
#'
#' @rawNamespace import(igraph, except = c(knn, union))
#'
sl2graph <- function(x) {
  
  if (!inherits(x, 'SpatialLines')) {
    stop('x is not a SpatialLines object.')
  }
  
  ## Make vertex DF
  coord.ls <- unlist(lapply(x@lines, function(l) lapply(l@Lines, function(ll) ll@coords)),
                     recursive = F)
  endcoords.ls <- lapply(coord.ls, function(x) x[c(1,nrow(x)),])
  endcoord.mat <- unique(do.call("rbind", endcoords.ls))
  vertex.df <- data.frame(vid=1:nrow(endcoord.mat), x=endcoord.mat[,1], y=endcoord.mat[,2])
  
  ## Make edgelist DF
  endcoords.el.ls <- lapply(endcoords.ls, function(x) matrix(c(x[1,], x[2,]), ncol=4))
  el.mat <- do.call("rbind", endcoords.el.ls)
  el.df <- as.data.frame(el.mat)
  names(el.df) <- c('from.x', 'from.y', 'to.x', 'to.y')
  el.df$order <- 1:nrow(el.df)
  fromvertex.df <- vertex.df
  names(fromvertex.df) <- c('from.vid', 'from.x', 'from.y')
  el.df <- merge(el.df, fromvertex.df, all.x=TRUE, all.y=FALSE)
  tovertex <- vertex.df
  names(tovertex) <- c('to.vid', 'to.x', 'to.y')
  el.df <- merge(el.df, tovertex, all.x=TRUE, all.y=FALSE)
  el.df <- el.df[order(el.df$order),]
  el.df <- el.df[,c('from.vid', 'to.vid', 'order')]
  
  ## Add attribute data, euclidean length
  if(inherits(x, "SpatialLinesDataFrame")){
    el.df <- cbind(el.df, x@data)
  }

  ## Make graph
  graph <- graph_from_data_frame(el.df, FALSE, vertex.df)
  
  ## Return
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
#' @import rgeos
#'
#' @rawNamespace import(igraph, except = c(knn, union))
#' @importFrom class knn
#' @importFrom methods .hasSlot
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
  del.tri <- remove_multilines(del.tri)
  
  # # To BoostLines
  # bl <- BoostLines(del.tri)
  # 
  # # Lines to igraph // "snaps" points to grid.
  # g <- boost2graph(bl, df = data.frame(lid = 1:length(del.tri)),
  #                  lonlat = T, plot = F, snap = 1e-12)
  g <- sl2graph(del.tri)
  
  # Check
  stopifnot(length(V(g)) == length(points))
  
  # Add vertex info
  lat.df <- cbind(x = V(g)$x, y = V(g)$y)
  ngb <- class::knn(lat.df, points@coords, 
                    cl = 1:nrow(lat.df),
                    k = as.integer(1))
  stopifnot(length(ngb) == length(unique(ngb)))
  
  if(.hasSlot(points, "data")){
    for(v in colnames(points@data)[!colnames(points@data) %in% c("x","y")]){
      vertex_attr(g, name = v, index = ngb) <- points@data[,v]
    }
  }
  
  # Original x and y coords
  vertex_attr(g, name = "x", index = ngb) <- points@coords[,1]
  vertex_attr(g, name = "y", index = ngb) <- points@coords[,2]
  
  
  # Return
  return(g)
}
