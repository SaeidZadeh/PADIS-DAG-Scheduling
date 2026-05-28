################################################################################
# PADIS: Precedence-Aware Deadline-Impact Scheduling
# Alirezazadeh et al., 2025
#
# Requirements: plot3D, vioplot, matrixStats, pcalg (via Bioconductor)
# Install with: source("R/padis.R") after installing dependencies
################################################################################

# ── 0. Dependencies ───────────────────────────────────────────────────────────
if (!require(plot3D))      install.packages("plot3D")
if (!require(vioplot))     install.packages("vioplot")
if (!require(matrixStats)) install.packages("matrixStats")
if (!require(pcalg)) {
  install.packages("BiocManager")
  BiocManager::install("Rgraphviz")
  install.packages("pcalg")
}
library(plot3D); library(vioplot); library(matrixStats); library(pcalg)

set.seed(1)

# ── 1. Graph utilities ────────────────────────────────────────────────────────

#' Convert a pcalg graph object to an adjacency matrix
toadj <- function(g) {
  a <- matrix(0, length(g@nodes), length(g@nodes))
  for (i in seq_along(g@nodes)) {
    k <- as.numeric(g@edgeL[[i]][[1]])
    if (length(k) > 1)              a[i, k] <- 1
    if (length(k) == 1 && k != i)  a[i, k] <- 1
  }
  a
}

#' Precompute DAG structure: adjacency, out-degrees, reachability, and Ov
prep_dag <- function(a) {
  n     <- ncol(a)
  Ab    <- (a > 0) * 1
  outdeg <- rowSums(Ab)
  reach  <- Ab > 0
  if (n >= 1)
    for (k in seq_len(n))
      reach <- reach | outer(reach[, k], reach[k, ], "&")
  Ov <- as.numeric((reach * 1) %*% outdeg)
  list(A = Ab, Ov = Ov, outdeg = outdeg, n = n)
}

# ── 2. PADIS-ZC algorithm ─────────────────────────────────────────────────────

#' Run PADIS-ZC (zero communication) on a single instance.
#'
#' @param pp     Output of prep_dag()
#' @param exec   m x n matrix of average execution times
#' @param deadline  Length-n vector of task deadlines
#' @return Number of accepted tasks
simulate_padis <- function(pp, exec, deadline) {
  m       <- nrow(exec)
  n       <- pp$n
  A       <- pp$A
  Ov      <- pp$Ov
  indeg   <- colSums(A)
  load    <- rep(0, m)
  removed <- rep(FALSE, n)
  acc     <- 0

  while (!all(removed)) {
    cand <- which(!removed & indeg == 0)
    if (length(cand) == 0) break

    # Lexicographic selection: earliest deadline, then highest Ov
    dmin <- min(deadline[cand])
    tie  <- cand[deadline[cand] == dmin]
    tmin <- tie[which.max(Ov[tie])]

    # Assign to feasible processor with minimum completion time
    H    <- load + exec[, tmin]
    feas <- which(H < deadline[tmin])
    if (length(feas) > 0) {
      p       <- feas[which.min(H[feas])]
      load[p] <- H[p]
      acc     <- acc + 1
    }

    # Update residual graph
    removed[tmin] <- TRUE
    succ <- which(A[tmin, ] > 0)
    if (length(succ) > 0) indeg[succ] <- indeg[succ] - 1
  }
  acc
}

run_padis <- function(pp, exec, deadline) simulate_padis(pp, exec, deadline)

# ── 3. Benchmark DAG generators ───────────────────────────────────────────────

make_ge_dag <- function(N) {
  id  <- matrix(0, N, N); cnt <- 0
  for (k in 1:(N - 1)) for (j in k:N) { cnt <- cnt + 1; id[k, j] <- cnt }
  a   <- matrix(0, cnt, cnt)
  for (k in 1:(N - 1)) for (j in (k + 1):N) {
    a[id[k, k], id[k, j]] <- 1
    if (id[k + 1, j] != 0) a[id[k, j], id[k + 1, j]] <- 1
  }
  a
}

make_stencil_dag <- function(W, L) {
  idx <- function(t, x) (t - 1) * W + x
  a   <- matrix(0, W * L, W * L)
  for (t in 2:L) for (x in 1:W) for (dx in c(-1, 0, 1)) {
    xp <- x + dx
    if (xp >= 1 && xp <= W) a[idx(t - 1, xp), idx(t, x)] <- 1
  }
  a
}

make_laplace_dag <- function(G, It) {
  idx <- function(it, x, y) ((it - 1) * G + (x - 1)) * G + y
  a   <- matrix(0, G * G * It, G * G * It)
  nb  <- list(c(0, 0), c(-1, 0), c(1, 0), c(0, -1), c(0, 1))
  for (it in 2:It) for (x in 1:G) for (y in 1:G) for (d in nb) {
    xp <- x + d[1]; yp <- y + d[2]
    if (xp >= 1 && xp <= G && yp >= 1 && yp <= G)
      a[idx(it - 1, xp, yp), idx(it, x, y)] <- 1
  }
  a
}

# ── 4. Simulation grids ───────────────────────────────────────────────────────

run_grid_fixed <- function(R = 50) {
  our <- matrix(0, 20, 50); sdm <- matrix(0, 20, 50)
  for (m in 2:21) {
    for (n in 2:51) {
      deadlines <- runif(n, 1, 40)
      exec_time <- matrix(0, m, n)
      for (i in 1:m) exec_time[i, ] <- runif(n, 1, 20)
      o <- replicate(R, {
        pp <- prep_dag(toadj(randDAG(n, n / 3)))
        run_padis(pp, exec_time, deadlines)
      })
      our[m - 1, n - 1] <- mean(o)
      sdm[m - 1, n - 1] <- sd(o)
    }
    cat("Fixed grid: m =", m, "\n")
  }
  list(mean = our, sd = sdm)
}

run_grid_all_random <- function(R = 50) {
  our <- matrix(0, 20, 50); sdm <- matrix(0, 20, 50)
  for (m in 2:21) {
    for (n in 2:51) {
      o <- replicate(R, {
        deadlines <- runif(n, 1, 40)
        exec_time <- matrix(0, m, n)
        for (i in 1:m) exec_time[i, ] <- runif(n, 1, 20)
        pp <- prep_dag(toadj(randDAG(n, n / 3)))
        run_padis(pp, exec_time, deadlines)
      })
      our[m - 1, n - 1] <- mean(o)
      sdm[m - 1, n - 1] <- sd(o)
    }
    cat("All-random grid: m =", m, "\n")
  }
  list(mean = our, sd = sdm)
}

# ── 5. Benchmark evaluation ───────────────────────────────────────────────────

run_benchmark <- function(a, label, R = 100, mset = 2:21) {
  pp  <- prep_dag(a); n <- pp$n; out <- data.frame()
  for (m in mset) {
    o <- replicate(R, {
      deadlines <- runif(n, 1, 40)
      exec      <- matrix(0, m, n)
      for (i in 1:m) exec[i, ] <- runif(n, 1, 20)
      run_padis(pp, exec, deadlines)
    })
    out <- rbind(out, data.frame(
      benchmark = label, n = n, m = m,
      PADIS = mean(o), sdP = sd(o)
    ))
  }
  out
}

run_all_benchmarks <- function() {
  ge_a      <- make_ge_dag(10)
  stencil_a <- make_stencil_dag(5, 10)
  laplace_a <- make_laplace_dag(4, 3)
  rbind(
    run_benchmark(ge_a,      "Gaussian"),
    run_benchmark(stencil_a, "Stencil"),
    run_benchmark(laplace_a, "Laplace")
  )
}

# ── 6. Scalability / runtime ──────────────────────────────────────────────────

time_padis_once <- function(n, m) {
  deadlines <- runif(n, 1, 40)
  exec      <- matrix(0, m, n)
  for (i in 1:m) exec[i, ] <- runif(n, 1, 20)
  as.numeric(system.time({
    pp <- prep_dag(toadj(randDAG(n, n / 3)))
    simulate_padis(pp, exec, deadlines)
  })["elapsed"])
}

timing_sweep <- function(ns   = c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048),
                          m    = 10,
                          reps = 20) {
  df <- data.frame(
    n       = ns,
    m       = m,
    seconds = sapply(ns, function(n)
      round(mean(replicate(reps, time_padis_once(n, m))), 6))
  )
  print(df)
  invisible(df)
}

# ── 7. Paper sanity check ─────────────────────────────────────────────────────

sanity_check <- function() {
  a        <- matrix(0, 3, 3); a[1, 2] <- 1   # T1->T2, T3 independent
  pp       <- prep_dag(a)
  deadline <- c(1.5, 3, 10)
  exec     <- rbind(c(1, 3, 5), c(3, 1, 2.5))
  cat(sprintf("Paper example — PADIS: %d  (expected 3)\n",
              run_padis(pp, exec, deadline)))
}
