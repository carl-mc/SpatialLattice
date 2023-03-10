% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_spatial_graph.R
\name{make_opt_color}
\alias{make_opt_color}
\title{Coloring Graph Partitions}
\usage{
make_opt_color(graph, var, colors = viridis(n = 10, option = "D"))
}
\arguments{
\item{graph}{igraph object with vertex attributes x and y, encoding coordinates.}

\item{var}{vertex attribute name that indicates the grouping / partitioning of vertices}

\item{colors}{Vector of candidate colors. Defaults to \code{viridis(16)}}
}
\value{
A vector of colors one color per partitioning. Ordered in the same order as \code{unique(vertex_attr(g, var))}
}
\description{
Implements Philipp Hunzikers \code{MapColoring} approach for partitions of a spatial graphs such that 
neighboring partitions receive colors of maximal distance. Returns result at the partition level.
}
