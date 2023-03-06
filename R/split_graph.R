#####################################
# Split Graph function
#####################################

#' Split graph into smaller pieces
#'
#' @param g  igraph object with vertex attributes x and y, encoding coordinates. 
#' @param type Type of split applied. 
#' sprandom, spregular, sphexagonal applies the respective spsample technique to vertices. 
#'  connected_component splits into connected components. 
#'  kmeans applies \code{stats::kmeans} clustering to vertices spatial coordinates.
#'  byvar splits by grouping variable. 
#' @param size Average size of clusteres required. 
#' @param drop.smaller.than Drop clusters of size x.
#' @param variable If \code{type = byvar}, this is the name of the grouping vertex attribute to split by. 
#'
#' @return List of graphs
#' 
#' @import sp
#' @importFrom rtree knn
#' @importFrom stats na.omit
#' @rawNamespace import(igraph, except = c(knn, union))
#' @export
split_graph <- function(g, type = c("sprandom", "spregular", "sphexagonal",
                                    "connected_component", "kmeans", "byvar"), 
                        size = NULL, drop.smaller.than = 1, variable = NULL){
  # Check size
  if(!is.null(size)){
    if(size >= length(V(g))){
      return(list(g))
    }
  }

  if(type == "connected_component"){
    # Find clusters
    clusters <- clusters(g)
    
    # Split
    cl.sum <- table(clusters$membership)
    g_ls <- lapply(which(cl.sum > drop.smaller.than), function(c){
      induced_subgraph(g, vids = which(clusters$membership == c))
    })
  } else if(type %in% c("sprandom", "spregular", "sphexagonal")){
    stopifnot(!is.null(size))
    
    ## Make Spatial Points
    v.pts <- SpatialPoints(cbind(V(g)$x, V(g)$y))
    
    ## Sample
    centers <- spsample(gBuffer(v.pts, width = .1, byid = T), ceiling(length(V(g))/size), 
                        type = gsub("sp","",type))
    
    ## Cluster
    centr.tree <- RTree(centers@coords)
    membership <- unlist(rtree::knn(centr.tree, v.pts@coords, k = 1L))
    
    ## Split graph
    cl.sum <- table(membership)
    g_ls <- lapply(which(cl.sum > drop.smaller.than), function(c){
      induced_subgraph(g, vids = which(membership == c))
    })
  } else if(type == "kmeans"){
    ## membership
    cluster <- stats::kmeans(cbind(V(g)$x, V(g)$y), centers = ceiling(length(V(g))/size))
    membership <- cluster$cluster
    
    ## Split graph
    cl.sum <- table(membership)
    g_ls <- lapply(which(cl.sum > drop.smaller.than), function(c){
      induced_subgraph(g, vids = which(membership == c))
    })
  } else if(type == "byvar"){
    g_ls <- lapply(na.omit(unique(vertex_attr(g, variable))), function(c){
      induced_subgraph(g, vids = which(vertex_attr(g, variable) == c & 
                                         !is.na(vertex_attr(g, variable))))
    })
    names(g_ls) <- paste0(variable, ".", na.omit(unique(vertex_attr(g, variable))))
  }
  
  # Return
  return(g_ls)
}
