library(fmesher)

vec <- expand.grid(
  max.edge = rev(exp(seq(log(0.1), log(1.5), length.out = 201))),
  min.angle = 21
)
vec <- expand.grid(
  max.edge = rev(exp(seq(log(0.25 - 1e-8), log(0.25 + 1e-8), length.out = 201))),
  min.angle = 34
)
vec <- expand.grid(
  max.edge = rev(exp(seq(log(0.1), log(1.5), length.out = 201))),
  min.angle = 26
)
vec <- expand.grid(
  max.edge = rev(exp(seq(log(0.25 - 1e-8), log(0.25 + 1e-8), length.out = 201))),
  min.angle = 1e-9
)
vec <- expand.grid(
  max.edge = 0.25,
  min.angle = seq(30, 40, by = 1)
)
vec <- expand.grid(
  max.edge = 10,
  min.angle = seq(1, 34, by = 0.1)
)
vec <- expand.grid(
  max.edge = 10,
  min.angle = seq(20, 26, by = 0.01)
)

 loc <- rbind(c(0.1, 0.2), c(0.2, 0.1), c(0.15, -0.1))
#loc <- rbind(c(0.0, 0.0))

res <- data.frame(
  dof = integer(nrow(vec)),
  min.angle0 = vec$min.angle,
  min.angle = numeric(nrow(vec)),
  max.edge0 = vec$max.edge,
  max.edge = numeric(nrow(vec))
)
mesh <- list()
for (k in seq_len(nrow(vec))) {
  mesh[[k]] <- fm_rcdt_2d_inla(loc,
    extend = list(offset = 1, n = 16),
    refine = list(
      max.edge = vec$max.edge[k],
      min.angle = vec$min.angle[k]
    )
  )
  # mesh[[k]] <- fm_rcdt_2d_inla(
  #   globe = 1,
  #   refine = list(
  #     max.edge = vec$max.edge[k],
  #     min.angle = vec$min.angle[k]
  #   )
  # )
  res$dof[k] <- fm_dof(mesh[[k]])
  edges <- list(
    mesh[[k]]$loc[mesh[[k]]$graph$tv[, 2], ] - mesh[[k]]$loc[mesh[[k]]$graph$tv[, 1], ],
    mesh[[k]]$loc[mesh[[k]]$graph$tv[, 3], ] - mesh[[k]]$loc[mesh[[k]]$graph$tv[, 2], ],
    mesh[[k]]$loc[mesh[[k]]$graph$tv[, 1], ] - mesh[[k]]$loc[mesh[[k]]$graph$tv[, 3], ]
  )
  res$min.angle[k] <-
    180 / pi * min(acos(pmin(1, pmax(
      -1, c(
        -rowSums(edges[[1]] * edges[[2]]) / rowSums(edges[[1]]^2)^0.5 / rowSums(edges[[2]]^2)^0.5,
        -rowSums(edges[[2]] * edges[[3]]) / rowSums(edges[[2]]^2)^0.5 / rowSums(edges[[3]]^2)^0.5,
        -rowSums(edges[[3]] * edges[[1]]) / rowSums(edges[[3]]^2)^0.5 / rowSums(edges[[1]]^2)^0.5
      )
    ))))
  res$max.edge[k] <- max(rowSums(do.call(rbind, edges)^2)^0.5)
}


plot(res$dof, res$min.angle, pch = 20, log = "x")
plot(res$dof, res$max.edge, pch = 20, log = "xy")
plot(res$dof, res$min.angle0 - res$min.angle, pch = 20, log = "x")
plot(res$dof, res$max.edge / vec$max.edge, pch = 20, log = "xy")

plot(res$min.angle0, res$dof, pch = 20, log = "xy")
plot(res$min.angle0, res$min.angle, pch = 20, log = "x", cex=0.5)
lines(res$min.angle0, res$min.angle0, col = 2)
plot(res$min.angle0, res$min.angle0 - res$min.angle, pch = 20, log = "x")
plot(res$min.angle0, res$max.edge, pch = 20, log = "xy")
plot(res$min.angle0, res$max.edge / vec$max.edge, pch = 20, log = "xy")
plot(res$min.angle, res$max.edge, pch = 20, log = "xy")

plot(vec$max.edge, res$dof, pch = 20, log = "xy")
plot(vec$max.edge, res$min.angle, pch = 20, log = "xy")
plot(vec$max.edge, res$max.edge, pch = 20, log = "xy")
plot(vec$max.edge, res$max.edge / vec$max.edge, pch = 20, log = "xy")

for (k in seq_len(nrow(vec) - 1) + 1) {
  plot(mesh[[k - 1]])
  points(mesh[[k]]$loc, pch = 20, col = 2)
  points(mesh[[k - 1]]$loc, pch = 20, col = 4)
  title(paste0(
    "min.angle = ", round(vec$min.angle[k], 3), " < ", round(res$min.angle[k], 3),
    ", ",
    "max.edge = ", round(vec$max.edge[k], 3), " > ", round(res$max.edge[k], 3)
  ))
  Sys.sleep(0.1)
}
