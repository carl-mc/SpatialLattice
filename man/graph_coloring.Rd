% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_spatial_graph.R
\name{graph_coloring}
\alias{graph_coloring}
\title{Coloring Vertices' Partitions}
\usage{
graph_coloring(graph, groupvar, cand.colors = viridis(16))
}
\arguments{
\item{graph}{igraph object with vertex attributes x and y, encoding coordinates.}

\item{groupvar}{vertex attribute name that indicates the grouping / partitioning of vertices}

\item{cand.colors}{Vector of candidate colors. Defaults to \code{viridis(16)}}
}
\value{
A vector of colors to color vertices with one color per partitioning
}
\description{
Implements Philipp Hunzikers \code{MapColoring} approach for partitions of a spatial graphs such that 
neighboring partitions receive colors of maximal distance. Returns result at the vertex level.
}
