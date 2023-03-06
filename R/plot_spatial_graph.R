#################################
# PLOT IGRAPH AS SPATIAL OBJECT
#################################


#' Transform edges to lines
#'
#' @description Transforms the graph's edges to plottable lines
#' 
#' @param g igraph object with vertex attributes x and y, encoding coordinates. 
#'
#' @return A four-column matrix of point pairs in the same order as g's edges. 
#' @export
#' @rawNamespace import(igraph, except = c(knn, union))
edges2lines <- function(g){
  e.vertices <- ends(g, E(g), names = F)
  pts <- cbind(vertex_attr(g, "x"),vertex_attr(g, "y"))
  cbind(pts[e.vertices[,1, drop = F],, drop = F], pts[e.vertices[,2, drop = F], , drop = F] )
}

#' Coloring Vertices' Partitions
#'
#' @description Implements Philipp Hunzikers \code{MapColoring} approach for partitions of a spatial graphs such that 
#' neighboring partitions receive colors of maximal distance. Returns result at the vertex level.
#' 
#' @param graph igraph object with vertex attributes x and y, encoding coordinates. 
#' @param groupvar vertex attribute name that indicates the grouping / partitioning of vertices
#' @param cand.colors Vector of candidate colors. Defaults to \code{viridis(16)}
#'
#' @return A vector of colors to color vertices with one color per partitioning
#' @export
#' @import MapColoring
#' viridis
#' @importFrom stats na.omit
#' @rawNamespace import(igraph, except = c(knn, union))
graph_coloring <- function(graph, groupvar, cand.colors = viridis(16)){
  # Edge ends
  e.ends <- ends(graph, E(graph), names = F)
  
  # Groups
  unit.vec <- vertex_attr(graph, groupvar)
  all.units <- na.omit(unique(unit.vec))
  
  # NGB Matrix
  pairs <- na.omit(e.ends[unit.vec[e.ends[,2]] != 
                            unit.vec[e.ends[,1]],])
  ngb.mat <- matrix(FALSE, nrow = length(all.units), ncol = length(all.units))
  for(p in seq_len(nrow(pairs))){
    u1 <- which(all.units == unit.vec[pairs[p,1]])
    u2 <- which(all.units == unit.vec[pairs[p,2]])
    ngb.mat[u1, u2] <- TRUE
    ngb.mat[u2, u1] <- TRUE
  }
  
  # Colors
  opt.colors <- MapColoring::getOptimalContrast(x=ngb.mat, 
                                                col=cand.colors)
  
  # To verteces
  vert.col <- unlist(sapply(unit.vec, function(i){
    ifelse(is.na(i), NA, opt.colors[which(all.units == i)])
  }))
  
  # Return
  return(vert.col)
}

#' Plot a spatial graph
#'
#' @description Spatial plotting of a graph using R base plot. 
#' Plots vertices as SpatialPoints and edges as SpatialLines.
#'
#' @param g igraph object with vertex attributes x and y, encoding coordinates. 
#' @param vertex.size Size of vertices
#' @param vertex.color Color of vertices
#' @param vertex.pch Pch of vertices
#' @param edge.color Color of edges
#' @param edge.width Width of edges
#' @param axes Show axes? 
#' @param bty Show box?
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param add Add to previous plot?
#' @param edge.lty Lty of edges
#' @param ... Additional parameters passed to \code{plot(NULL)}
#' @rawNamespace import(igraph, except = c(knn, union))
#' @return Nothing.
#' @export
#' @importFrom graphics segments points
plot_spatial_graph <- function(g, 
                               vertex.size = 0.25, 
                               vertex.color = "black", vertex.pch = 19,
                               edge.color = "grey", edge.width = .25,
                               axes=F, bty = "n", xlab = "", ylab = "",
                               add = F, edge.lty = 1,
                               ...){
  
  # Vertices to SpatialPoints
  pts <- SpatialPoints(cbind(vertex_attr(g, "x"),vertex_attr(g, "y")))
  
  # Edges to SpatialLines
  lines <- edges2lines(g)
  
  # Plot extent
  ext <- extent(pts)
  
  # Plot
  if(!add){
    plot(NULL, xlim = ext[1:2], ylim = ext[3:4], axes = axes, bty = bty, xlab = xlab, ylab = ylab, ...) 
  }
  if(!is.null(lines)){
    segments(lines[,1],lines[,2],lines[,3],lines[,4], lwd = edge.width, col = edge.color, lty = edge.lty)
  }
  points(pts, pch = vertex.pch, cex = vertex.size, col = vertex.color)
}



#' Coloring Graph Partitions
#'
#' @description Implements Philipp Hunzikers \code{MapColoring} approach for partitions of a spatial graphs such that 
#' neighboring partitions receive colors of maximal distance. Returns result at the partition level.
#' 
#' @param graph igraph object with vertex attributes x and y, encoding coordinates. 
#' @param var vertex attribute name that indicates the grouping / partitioning of vertices
#' @param colors Vector of candidate colors. Defaults to \code{viridis(16)}
#'
#' @return A vector of colors one color per partitioning. Ordered in the same order as \code{unique(vertex_attr(g, var))}
#' @export
#' @import MapColoring
#' viridis
#' @rawNamespace import(igraph, except = c(knn, union))
#'
make_opt_color <- function(graph, var, colors = viridis(n = 10, option = "D")){

  
  Y <- as.numeric(factor(vertex_attr(graph, var), 
                            levels = unique(vertex_attr(graph, var)), 
                            ordered = T))
  Y.ngb <- unique(rbind(cbind(Y[ends(graph, E(graph), names = F)[,1]], 
                              Y[ends(graph, E(graph), names = F)[,2]]),
                        cbind(Y[ends(graph, E(graph), names = F)[,2]], 
                              Y[ends(graph, E(graph), names = F)[,1]])))
  adj.mat <- sapply(1:max(Y),function(i){
    sapply(1:max(Y), function(j){
      any(Y.ngb[,1] == i & Y.ngb[, 2] == j)
    })
  })
  diag(adj.mat) <- T
  
  ## Get Optimal contrast colors
  return(getOptimalContrast(x=adj.mat, colors=colors)[Y])
}
