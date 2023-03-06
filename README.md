# SpatialLattice: Spatial Graph Manipulation in R

## Introduction

## Installation

``` r
# Download SpatialLattice package
library(devtools)
install_github(repo = "carl-mc/SpatialLattice")
```

## Creating spatial graph data

``` r
# Load some libraries needed below
library(SpatialLattice)
library(igraph)
library(sp)
library(cshapes)
library(sf)

# Load country borders in Africa
africa <- cshp(date = as.Date("2000-01-01"))
africa <- africa[africa$gwcode %in% c(400:626,651),]
africa <- as_Spatial(africa)
```

### Regular grid structures

``` r
# Sample grid with hexagonal structure
grid <- spatial_grid(shape = africa, 
                     structure = "hexagonal",
                     N_pts = 500)

# Print vertex attributes
print(vertex_attr_names(grid))
```

    ## [1] "name" "x"    "y"

``` r
# Plot spatial graph
plot_spatial_graph(grid)
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

### Graphs from polygons

``` r
# Graph from polygons
poly.graph <- graph_from_poly(africa)

# Print vertex attributes
print(vertex_attr_names(poly.graph))
```

    ##  [1] "x"            "y"            "gwcode"       "country_name" "start"       
    ##  [6] "end"          "status"       "owner"        "capname"      "caplong"     
    ## [11] "caplat"       "b_def"        "fid"

``` r
# Plot spatial graph
plot(africa)
plot_spatial_graph(poly.graph, add = T)
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### Graphs from points

``` r
# Graph from points

## Make points as Spatial* object
pts <- SpatialPoints(cbind(africa$caplong,
                           africa$caplat))

## Make graph based on Delaunay Triangulation
pts.graph <- graph_from_pts(points = pts)

## Plot
plot(africa)
points(pts, col = "red")
plot_spatial_graph(pts.graph, add = T)
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Basic manipulation of spatial graphs

``` r
# Simple transformations

## Vertices to SpatialPoints 
pts <- vertices2sp(grid)

## Edges to SpatialLines
lines <- edges2sl(grid)

## Length
E(grid)$length <- geosphere::lengthLine(lines) / 1000
```

## Adding data to spatial graphs

### Adding data from Spatial Polygons

``` r
# Adding poly data to spatial graphs

## From polygons
V(grid)$poly.id <- poly2graph(grid, africa)

## Use IDs to encode data
V(grid)$gwcode <- africa$gwcode[V(grid)$poly.id]

## Encode country-borders on edges
E(grid)$gwcode.diff <- vid_2_edge_diff(grid,
                                       "gwcode")

## Plot countries on grid
gwcode.col <- graph_coloring(graph = grid, groupvar = "gwcode")
plot_spatial_graph(grid, 
                   vertex.color = gwcode.col, vertex.size = .5,
                   edge.color = ifelse(E(grid)$gwcode.diff == 0, "black", "lightgrey"),
                   edge.width = 1)
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### Adding data from spatial lines

``` r
# Adding line data to spatial graphs

## Make lines
gw.borders <- as(africa, "SpatialLines")
gw.borders <- SpatialLinesDataFrame(gw.borders,
                                    data.frame(is.border = rep(1, length(gw.borders))),
                                    match.ID = F)

## From polygons
E(grid)$crosses.border <- lines2graph(gw.borders, grid, aggvar = "is.border",
                               aggfun = max, 
                               na.val = 0,
                               only.uneven = TRUE)

## Plot borders on grid
plot_spatial_graph(grid, 
                   vertex.color = gwcode.col, vertex.size = .5,
                   edge.color = ifelse(E(grid)$crosses.border == 1, "black", "lightgrey"),
                   edge.width = 1)
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### Adding data from a spatial raster

``` r
# Adding raster data to spatial graphs

## Raster library
library(raster)

## Make a raster
rand.r <- raster(ext = extent(africa),
                 res = 1)
rand.r <- raster::setValues(rand.r, 
                            runif(n = ncell(rand.r)))

## Aggregate to edge level 
E(grid)$rand.r <- raster2graph(raster = rand.r, graph = grid, 
                               edges = T, 
                               aggfun = mean)

## Plot
plot_spatial_graph(grid, 
                   vertex.size = .5,
                   edge.color = grey(E(grid)$rand.r),
                   edge.width = 2)
```

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## Analyzing Spatial Partitionings with pspm

``` r
# Analyze raster with PSPM
library(pspm)
model <- fit_pspm_model(gwcode ~ length + rand.r , 
                        g_ls = list(grid))
summary(model)
```

    ## --------------------------------------------
    ## Maximum Likelihood estimation
    ## BFGS maximization, 31 iterations
    ## Return code 0: successful convergence 
    ## Log-Likelihood: -198.0226 
    ## 3  free parameters
    ## Estimates:
    ##           Estimate Std. error t value Pr(> t)
    ## Constant -2.353742   1.575533  -1.494   0.135
    ## length    0.004134   0.005570   0.742   0.458
    ## rand.r    0.237885   0.443821   0.536   0.592
    ## --------------------------------------------

## Feedback, comments, questions

We are very grateful for any bug reports, feedback, questions, or
contributions to this package. Please report any issues here or write to
c.a.muller-crepon \[at\] lse.ac.uk .
