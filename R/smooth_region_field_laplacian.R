#' Smooth region values over an adjacency graph
#'
#' @param region_id_mat matrix of region identifiers.
#' @param region_values named numeric vector or data frame with `region_id` and `value`.
#' @param outside_ids unused compatibility placeholder for non-region values.
#' @param adjacency neighborhood definition: `"rook"` or `"queen"`.
#' @param lambda smoothing strength.
#' @param keep_na_outside logical; keep non-region pixels as `NA`.
#' @return A list with the interpolated pixel matrix and region-level predictions.
#' @export
smooth_region_field_laplacian <- function(
    region_id_mat,
    region_values,              # named numeric vector OR data.frame(region_id, value)
    outside_ids = c(NA_integer_, 0L),
    adjacency = c("rook", "queen"),
    lambda = 5,                 # higher = smoother
    keep_na_outside = TRUE
) {
  adjacency <- match.arg(adjacency)
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }

  # --- normalize inputs ---
  idmat <- region_id_mat
  idmat[idmat == 0L] <- NA_integer_

  # region values -> data.frame
  if (is.data.frame(region_values)) {
    stopifnot(all(c("region_id", "value") %in% names(region_values)))
    df <- region_values[, c("region_id", "value")]
    df$region_id <- as.integer(df$region_id)
    df$value <- as.numeric(df$value)
  } else {
    if (is.null(names(region_values))) stop("region_values must be a named numeric vector.")
    df <- data.frame(
      region_id = as.integer(names(region_values)),
      value = as.numeric(region_values)
    )
  }

  # regions present in matrix
  region_ids_in_mat <- sort(unique(as.integer(idmat[is.finite(idmat)])))
  df <- df[df$region_id %in% region_ids_in_mat, , drop = FALSE]
  if (nrow(df) < 2) stop("Need at least 2 regions with values.")

  # universe of regions to smooth over = all regions in matrix
  region_ids <- region_ids_in_mat
  n <- length(region_ids)
  id_to_idx <- setNames(seq_len(n), region_ids)

  # --- build adjacency edges from pixel lattice ---
  ny <- nrow(idmat); nx <- ncol(idmat)

  collect_pairs <- function(A, B) {
    E <- cbind(as.vector(A), as.vector(B))
    ok <- is.finite(E[,1]) & is.finite(E[,2]) & (E[,1] != E[,2])
    E <- E[ok, , drop = FALSE]
    a <- pmin(E[,1], E[,2])
    b <- pmax(E[,1], E[,2])
    unique(cbind(a, b))
  }

  edges <- collect_pairs(idmat[, -nx, drop=FALSE], idmat[, -1, drop=FALSE])
  edges <- unique(rbind(edges, collect_pairs(idmat[-ny, , drop=FALSE], idmat[-1, , drop=FALSE])))

  if (adjacency == "queen") {
    edges <- unique(rbind(
      edges,
      collect_pairs(idmat[-ny, -nx, drop=FALSE], idmat[-1, -1, drop=FALSE]),
      collect_pairs(idmat[-ny, -1,  drop=FALSE], idmat[-1, -nx, drop=FALSE])
    ))
  }

  # restrict edges to known region_ids
  keep <- edges[,1] %in% region_ids & edges[,2] %in% region_ids
  edges <- edges[keep, , drop = FALSE]
  if (nrow(edges) == 0) stop("No adjacency edges found (are regions disconnected?).")

  i <- id_to_idx[as.character(edges[,1])]
  j <- id_to_idx[as.character(edges[,2])]

  # --- build Laplacian L = D - W (unweighted) ---
  W <- Matrix::sparseMatrix(
    i = c(i, j),
    j = c(j, i),
    x = 1,
    dims = c(n, n)
  )
  deg <- Matrix::rowSums(W)
  D <- Matrix::Diagonal(x = deg)
  L <- D - W

  # --- observed y (can be missing for some regions) ---
  y <- rep(NA_real_, n)
  y[id_to_idx[as.character(df$region_id)]] <- df$value

  obs <- which(is.finite(y))
  if (length(obs) < 2) stop("Need at least 2 finite region values to smooth.")

  # Solve only for all nodes but anchor to observations via data term:
  # minimize ||M(z - y)||^2 + lambda z^T L z
  # => (M + lambda L) z = M y, where M has 1 on observed nodes, 0 otherwise.
  M <- Matrix::Diagonal(n = n, x = as.numeric(is.finite(y)))
  A <- M + lambda * L
  b <- as.numeric(M %*% y)

  zhat <- as.numeric(Matrix::solve(A, b))

  pred_by_region <- setNames(zhat, region_ids)

  # paint back to pixels
  out <- matrix(if (keep_na_outside) NA_real_ else 0, ny, nx)
  okpix <- is.finite(idmat)
  out[okpix] <- pred_by_region[as.character(idmat[okpix])]

  list(
    interpolated_matrix = out,
    predicted_region_values = pred_by_region,
    lambda = lambda,
    adjacency = adjacency
  )
}
