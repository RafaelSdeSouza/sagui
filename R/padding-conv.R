# nocov start
# Internal: 1D reflect padding for a numeric vector
.reflect_pad_vec <- function(v, pad) {
  n <- length(v)
  if (pad <= 0L || n == 0L) return(v)
  pad <- max(0L, min(pad, n - 1L))     # <- cap pad
  if (pad == 0L) return(v)
  left  <- rev(v[seq_len(pad)])
  right <- rev(v[seq.int(n - pad + 1L, n)])
  c(left, v, right)
}


# Internal: row-wise convolution with reflection
.conv_reflect_rows <- function(mat, k) {
  nrow_m <- nrow(mat); ncol_m <- ncol(mat)
  p <- floor(length(k) / 2)
  out <- matrix(0, nrow_m, ncol_m)
  for (r in seq_len(nrow_m)) {
    v  <- mat[r, ]
    pv <- .reflect_pad_vec(v, p)       # <- dotted
    for (x in seq_len(ncol_m)) out[r, x] <- sum(k * pv[x:(x + length(k) - 1L)])
  }
  out
}

.conv_reflect_cols <- function(mat, k) {
  nrow_m <- nrow(mat); ncol_m <- ncol(mat)
  p <- floor(length(k) / 2)
  out <- matrix(0, nrow_m, ncol_m)
  for (c in seq_len(ncol_m)) {
    v  <- mat[, c]
    pv <- .reflect_pad_vec(v, p)       # <- dotted
    for (y in seq_len(nrow_m)) out[y, c] <- sum(k * pv[y:(y + length(k) - 1L)])
  }
  out
}


# Internal: separable smoothing (X then Y)
.smooth_sep_raw <- function(img, k) {
  tmp <- .conv_reflect_rows(img, k)
  .conv_reflect_cols(tmp, k)
}

# Masked/normalized separable smoothing: (K*(img*mask))/(K*mask)
.smooth_sep_masked <- function(img, k, mask) {
  num <- .smooth_sep_raw(img * mask, k)
  den <- .smooth_sep_raw(mask,     k)
  out <- num / pmax(den, 1e-12)        # avoid divide-by-zero
  out[den < 1e-12] <- NA_real_ # undefined where kernel sees no valid pixels
  out
}
# nocov end

.smooth_sep <- function(img, k) {
  tmp <- .conv_reflect_rows(img, k)  # along x
  .conv_reflect_cols(tmp, k)         # along y
}
