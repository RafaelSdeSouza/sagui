#' Clean a logical segmentation support mask
#'
#' Optionally applies binary close/open morphology and removes small connected
#' components from a 2-D support mask.
#'
#' @param mask Logical or numeric matrix. Positive/`TRUE` finite pixels are kept.
#' @param min_mask_area Minimum connected-component area to retain.
#' @param close_size Size of the optional binary closing brush. Use `1` to skip.
#' @param open_size Size of the optional binary opening brush. Use `1` to skip.
#' @param keep_largest Logical; keep only the largest connected component after
#'   area filtering.
#' @param connectivity Pixel connectivity for component labeling. Use `8` for
#'   diagonal connectivity or `4` for edge-only connectivity.
#'
#' @return A logical matrix with the same dimensions as `mask`.
#' @export
clean_segmentation_mask <- function(mask,
                                    min_mask_area = 1L,
                                    close_size = 1L,
                                    open_size = 1L,
                                    keep_largest = FALSE,
                                    connectivity = c(8L, 4L)) {
  stopifnot(is.matrix(mask))
  connectivity <- as.integer(match.arg(as.character(connectivity), c("8", "4")))

  out <- is.finite(mask) & mask > 0
  before <- sum(out, na.rm = TRUE)

  close_size <- as.integer(close_size)
  open_size <- as.integer(open_size)
  if (!is.finite(close_size) || close_size < 1L) {
    close_size <- 1L
  }
  if (!is.finite(open_size) || open_size < 1L) {
    open_size <- 1L
  }

  if (close_size > 1L || open_size > 1L) {
    out <- .binary_close_open(out, close_size = close_size, open_size = open_size)
  }

  min_mask_area <- as.integer(min_mask_area)
  if (!is.finite(min_mask_area) || min_mask_area < 1L) {
    min_mask_area <- 1L
  }

  if (min_mask_area > 1L || isTRUE(keep_largest)) {
    labels <- .label_mask_components(out, connectivity = connectivity)
    ids <- sort(unique(labels[labels > 0L]))
    if (length(ids)) {
      sizes <- tabulate(labels[labels > 0L], nbins = max(ids))
      keep <- ids[sizes[ids] >= min_mask_area]
      if (isTRUE(keep_largest) && length(keep)) {
        keep <- keep[which.max(sizes[keep])]
      }
      out <- matrix(labels %in% keep, nrow = nrow(labels), ncol = ncol(labels))
    }
  }

  attr(out, "cleanup") <- list(
    before_pixels = before,
    after_pixels = sum(out, na.rm = TRUE),
    min_mask_area = min_mask_area,
    close_size = close_size,
    open_size = open_size,
    keep_largest = isTRUE(keep_largest),
    connectivity = connectivity
  )
  out
}

.label_mask_components <- function(mask, connectivity = 8L) {
  mask <- as.matrix(mask)
  mask <- is.finite(mask) & mask > 0
  nr <- nrow(mask)
  nc <- ncol(mask)
  labels <- matrix(0L, nrow = nr, ncol = nc)

  if (!any(mask, na.rm = TRUE)) {
    return(labels)
  }

  if (connectivity == 4L) {
    neigh <- rbind(c(-1L, 0L), c(1L, 0L), c(0L, -1L), c(0L, 1L))
  } else {
    neigh <- as.matrix(expand.grid(dr = -1L:1L, dc = -1L:1L))
    neigh <- neigh[rowSums(abs(neigh)) > 0L, , drop = FALSE]
  }

  component_id <- 0L
  queue_r <- integer(nr * nc)
  queue_c <- integer(nr * nc)

  for (r0 in seq_len(nr)) {
    for (c0 in seq_len(nc)) {
      if (!mask[r0, c0] || labels[r0, c0] != 0L) {
        next
      }

      component_id <- component_id + 1L
      head <- 1L
      tail <- 1L
      queue_r[tail] <- r0
      queue_c[tail] <- c0
      labels[r0, c0] <- component_id

      while (head <= tail) {
        r <- queue_r[head]
        c <- queue_c[head]
        head <- head + 1L

        for (i in seq_len(nrow(neigh))) {
          rr <- r + neigh[i, 1]
          cc <- c + neigh[i, 2]
          if (rr < 1L || rr > nr || cc < 1L || cc > nc) {
            next
          }
          if (!mask[rr, cc] || labels[rr, cc] != 0L) {
            next
          }
          tail <- tail + 1L
          queue_r[tail] <- rr
          queue_c[tail] <- cc
          labels[rr, cc] <- component_id
        }
      }
    }
  }

  labels
}

.apply_mask_cleanup <- function(mask_info,
                                clean_mask = FALSE,
                                min_mask_area = 1L,
                                close_size = 1L,
                                open_size = 1L,
                                keep_largest = FALSE) {
  if (!isTRUE(clean_mask)) {
    return(mask_info)
  }

  cleaned <- clean_segmentation_mask(
    mask_info$mask,
    min_mask_area = min_mask_area,
    close_size = close_size,
    open_size = open_size,
    keep_largest = keep_largest
  )
  mask_info$mask <- cleaned
  mask_info$mask_cleanup <- attr(cleaned, "cleanup")
  mask_info
}
