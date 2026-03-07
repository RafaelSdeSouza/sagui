#' Reindex positive region labels to a compact sequence
#'
#' @param labels Matrix of integer labels.
#' @return Matrix with positive labels remapped to `1:n`.
#' @export
reindex_regions <- function(labels) {
  stopifnot(is.matrix(labels))
  region_ids <- sort(unique(labels[is.finite(labels) & labels > 0]))
  if (!length(region_ids)) {
    return(labels)
  }

  remap <- seq_along(region_ids)
  names(remap) <- region_ids

  out <- labels
  keep <- is.finite(out) & out > 0
  out[keep] <- as.integer(remap[as.character(out[keep])])
  out
}

#' Filter regions by pixel count
#'
#' @param labels Matrix of integer labels.
#' @param min_size Minimum region size to keep.
#' @param max_size Maximum region size to keep.
#' @return Relabeled matrix containing only regions within the requested size range.
#' @export
filter_regions_by_size <- function(labels, min_size = 30L, max_size = Inf) {
  stopifnot(is.matrix(labels))
  if (!any(labels > 0, na.rm = TRUE)) {
    return(labels)
  }

  counts <- table(labels[labels > 0])
  keep <- as.integer(names(counts)[counts >= min_size & counts <= max_size])

  out <- labels
  out[!(out %in% keep)] <- 0L
  reindex_regions(out)
}

#' Fill zero-valued holes inside each labeled region
#'
#' @param labels Matrix of integer labels.
#' @return Matrix with enclosed background holes filled per region.
#' @export
fill_region_holes <- function(labels) {
  stopifnot(is.matrix(labels))

  height <- nrow(labels)
  width <- ncol(labels)
  out <- labels

  for (region_id in sort(unique(labels[labels > 0]))) {
    idx <- which(labels == region_id)
    if (!length(idx)) {
      next
    }

    yy <- ((idx - 1L) %% height) + 1L
    xx <- ((idx - 1L) %/% height) + 1L
    y1 <- min(yy)
    y2 <- max(yy)
    x1 <- min(xx)
    x2 <- max(xx)

    sub <- out[y1:y2, x1:x2, drop = FALSE]
    bg <- sub == 0L
    h <- nrow(sub)
    w <- ncol(sub)
    visited <- matrix(FALSE, h, w)
    queue <- integer(h * w)
    head <- 1L
    tail <- 0L

    push <- function(y, x) {
      if (y >= 1L && y <= h && x >= 1L && x <= w && bg[y, x] && !visited[y, x]) {
        visited[y, x] <<- TRUE
        tail <<- tail + 1L
        queue[tail] <<- (x - 1L) * h + y
      }
    }

    for (x in seq_len(w)) {
      push(1L, x)
      push(h, x)
    }
    for (y in seq_len(h)) {
      push(y, 1L)
      push(y, w)
    }

    while (head <= tail) {
      p <- queue[head]
      head <- head + 1L
      y <- ((p - 1L) %% h) + 1L
      x <- ((p - 1L) %/% h) + 1L

      push(y - 1L, x)
      push(y + 1L, x)
      push(y, x - 1L)
      push(y, x + 1L)
    }

    holes <- bg & !visited
    sub[holes] <- region_id
    out[y1:y2, x1:x2] <- sub
  }

  out
}
