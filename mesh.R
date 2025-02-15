
# description -------------------------------------------------------------

# This script creates the lattices for the mesh. 

# nepal_lattice_sfc/2: hexagon mesh lattices with different edge lengths
# bnd: boundary
# bnd_out: outer boundary
# bnd_buff: smoothed boundary
# mesh_fm: mesh with finer resolution
# mesh_fmc: mesh with coarser resolution

# construct lattice -------------------------------------------------------


system.time({
  nepal_lattice_sfc <- hexagon_lattice(
    bnd = bnd, x_bin = x_bin,
    edge_len_n = edge_len_n
  )
})

system.time({
  nepal_lattice_sfc2 <- hexagon_lattice(
    bnd = bnd, x_bin = x_bin / 2,
    edge_len_n = edge_len_n
  )
})

boundary <- fm_as_segm_list(list(
  bnd,
  bnd_out
))

bnd1 <- fm_nonconvex_hull(bnd,
  nepal_lattice_sfc$edge_len,
  crs = fm_crs(bnd)
)

bnd2 <- fm_nonconvex_hull(bnd,
  10 * nepal_lattice_sfc$edge_len,
  crs = fm_crs(bnd)
)

## fm_rcdt_2d --------------------------------------------------------------
system.time({
  rcdt_fm <- fm_rcdt_2d_inla(
    loc = nepal_lattice_sfc$lattice,
    boundary = bnd1,
    refine = list(
      max.edge = c(1.1 * nepal_lattice_sfc$edge_len),
      min.angle = 1
    ),
    cutoff = 0.5 * nepal_lattice_sfc$edge_len,
    crs = fm_crs(bnd)
  ) 
})


## fm_mesh_2d --------------------------------------------------------------

system.time({
  mesh_fm <- fm_mesh_2d_inla(
    loc = nepal_lattice_sfc$lattice,
    boundary = fm_extensions(bnd, c(
      nepal_lattice_sfc$edge_len,
      10 * nepal_lattice_sfc$edge_len
    )), # a shortcut for bnd1 and bnd2
    max.edge = c(
      1.1 * nepal_lattice_sfc$edge_len, # avoid numerical overflow
      5 * nepal_lattice_sfc$edge_len
    ),
    max.n.strict = c(
      as.integer(0.9 * rcdt_fm$n),
      2000
    ), # arbitrary
    cutoff = 0.9 * nepal_lattice_sfc$edge_len, # Filter away adjacent points to avoid Q matrix not positive definite because the boundary resolution is too high
    crs = fm_crs(bnd)
  )
})

system.time({
  mesh_fmc <- fm_mesh_2d_inla(
    loc = nepal_lattice_sfc2$lattice,
    boundary = fm_extensions(bnd, c(
      nepal_lattice_sfc2$edge_len,
      10 * nepal_lattice_sfc2$edge_len
    )), # a shortcut for bnd1 and bnd2
    max.edge = c(
      1.1 * nepal_lattice_sfc2$edge_len, # avoid numerical overflow
      5 * nepal_lattice_sfc2$edge_len
    ),
    max.n.strict = c(
      as.integer(0.9 * rcdt_fm$n),
      2000
    ), # arbitrary
    cutoff = 0.9 * nepal_lattice_sfc2$edge_len, # Filter away adjacent points to avoid Q matrix not positive definite because the boundary resolution is too high
    crs = fm_crs(bnd)
  )
})


# Plot mesh
if (to_plot) {
  pmesh <- ggplot() +
    geom_fm(data = mesh_fm) +
    geom_sf(data = bnd, col = "red", fill = NA) +
    geom_sf(data = bnd1, col = "purple", fill = NA) +
    geom_sf(data = bnd2, col = "pink", fill = NA) +
    geom_sf(data = bnd_buff, col = "green", fill = NA)

  prcdt <- ggplot() +
    geom_fm(data = rcdt_fm) +
    geom_sf(data = bnd, fill = NA)
}

rm(rcdt_fm, nepal_lattice_sfc, nepal_lattice_sfc2)

