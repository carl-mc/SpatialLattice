% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_grid.R
\name{reproject_graph}
\alias{reproject_graph}
\title{Backproject graph}
\usage{
reproject_graph(
  graph,
  from,
  to = "+proj=longlat +ellps=WGS84 +no_defs",
  coords = c("x", "y")
)
}
\arguments{
\item{graph}{Igraph object with x and y coords}

\item{from}{source projection}

\item{to}{target projection}

\item{coords}{character vector with the names of vertex attributes that encode the xand y coordinates. 
Defaults to \code{c("x","y")}.}
}
\value{
Reprojected graph
}
\description{
Backproject graph
}
