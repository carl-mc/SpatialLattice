% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crop_graph.R
\name{delete_graph_islands}
\alias{delete_graph_islands}
\title{Delete graph islands}
\usage{
delete_graph_islands(g)
}
\arguments{
\item{g}{igraph}
}
\value{
A reduced igraph
}
\description{
Deletes all disconnected components except the largest one.
}
