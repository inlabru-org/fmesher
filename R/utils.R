#' @title Call stack utility functions
#' @description
#' Helper functions for displaying call stack information
#'
#' @param which The number of frames to go back from the caller
#' @param override character; Overrides the automated function name logic
#' @returns `fm_caller_name` returns a string with the the name of a calling
#'   function
#' @export
#' @name call-stack
#' @rdname call-stack
#' @keywords internal
#' @examples
#' fun <- function() {
#'   print(fm_caller_name())
#'   nm <- fm_caller_name()
#'   print(nm)
#' }
#' fun()
fm_caller_name <- function(which = 0L, override = NULL) {
  if (is.null(override)) {
    which <- -abs(which) - 1L
    if (abs(which) > sys.nframe()) {
      name <- ""
    } else {
      fun <- sys.call(which)
      if (is.null(fun)) {
        name <- ""
      } else {
        name <- as.character(fun)[[1]]
      }
    }
  } else {
    name <- override
  }
  name
}

#' @describeIn call-stack
#'
#' @param start The stack starting point
#' @param end The stack end point
#' @param with_numbers INclude call stack location numbers
#' @param \dots Currently unused
#' @returns `fm_call_stack` returns a character vector
#' @export
fm_call_stack <- function(start = 0L, end = 0L, with_numbers = TRUE, ...) {
  stack <- sys.calls()
  stack <- lapply(
    as.list(stack),
    function(x) as.character(deparse(x))
  )[
    start + seq_len(max(0, length(stack) - (abs(end) + 1L) - start))
  ]
  if (length(stack) > 0) {
    msg <-
      paste0(
        if (with_numbers) {
          paste0(seq_along(stack), ": ")
        } else {
          ""
        },
        lapply(
          stack,
          function(x) {
            paste0(
              vapply(
                x,
                function(x) {
                  if (nchar(x) > 80) {
                    paste0(
                      substr(x, 1, 74),
                      " [...]"
                    )
                  } else {
                    x
                  }
                },
                ""
              ),
              collapse = paste0("\n   ")
            )
          }
        )
      )
  } else {
    msg <- "Empty"
  }
  msg
}


#' @describeIn call-stack Inspired by `berryFunctions::tryStack`
#'
#' @param expr An `expression` to evaluate
#' @returns `fm_try_callstack` If successful, returns (invisibly) the value from
#'   the evaluated expression, otherwise an error object with call stack
#'   information attached to the error message.
#' @export
fm_try_callstack <- function(expr) {
  try_envir <- new.env()
  assign("error_stack", value = NULL, envir = try_envir)
  error_fun <- function(e) {
    # Get whole stack except the handlers
    stack <- fm_call_stack(start = 0, end = 2, with_numbers = FALSE)
    # Remove the fm_try_callstack tryCatch calls part(s),
    # There are 6 of them. First find the fm_try_callstack call (or multiple
    # calls for nested use, which should theoretically (almost) never happen,
    # since the inner call shouldn't fail!)
    self <- which(
      vapply(stack, function(x) {
        grepl("^fm_try_callstack\\(", x)
      }, TRUE) |
        vapply(stack, function(x) {
          grepl("^fmesher::fm_try_callstack\\(", x)
        }, TRUE) |
        vapply(stack, function(x) {
          grepl("^fmesher:::fm_try_callstack\\(", x)
        }, TRUE)
    )
    for (idx in rev(self)) {
      stack <- stack[-(idx + seq_len(6))]
      stack[idx] <- "fm_try_callstack(...)"
    }
    stack <- paste0(seq_len(length(stack)), ": ", stack, collapse = "\n")
    assign("error_stack", value = stack, envir = try_envir)
  }
  result <- try(
    withCallingHandlers(
      expr,
      error = error_fun
    ),
    silent = TRUE
  )
  if (inherits(result, "try-error")) {
    result[length(result) + 1] <- paste0(
      try_envir$error_stack,
      collapse = "\n"
    )
  }
  invisible(result)
}



fm_require_message <- function(pkg, msg = NULL, override = NULL) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(TRUE)
  }
  name <- fm_caller_name(1L, override = override)
  if (is.null(msg)) {
    msg <- paste0("Please install '", pkg, "'.")
  }
  message(
    paste0(
      "The function `",
      name,
      "()` requested the package '",
      pkg,
      "' but it is unavailable.",
      "\n",
      msg
    )
  )
  return(FALSE)
}
fm_require_stop <- function(pkg, msg = NULL, override = NULL) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(TRUE)
  }
  name <- fm_caller_name(1L, override = override)
  msg0 <- paste0("Please install '", pkg, "'.")
  if (is.null(msg)) {
    msg <- msg0
  } else {
    msg <- paste0(msg, "\n", msg0)
  }
  stop(
    paste0(
      "The function `",
      name,
      "()` requested the package '",
      pkg,
      "' but it is unavailable.",
      "\n",
      msg
    )
  )
}



#' @title Conversion between sparse matrix types
#' @rdname fmesher_sparse
#' @param x Object to be converted
#' @importFrom methods as
# Explicit import of something from Matrix to appease automated checks:
#' @importFrom Matrix as.matrix
#' @keywords internal
#' @returns `fm_as_dgCMatrix` returns a [Matrix::dgCMatrix-class] object.
#' @export
#' @examples
#' library(Matrix)
#' str(A <- fm_as_dgCMatrix(matrix(c(1, 2, 0, 0, 0, 3, 4, 0, 5), 3, 3)))
#' str(fm_as_dgTMatrix(A))
#' str(fm_as_unpackedMatrix(A))
#' str(fm_as_fmesher_sparse(A))
fm_as_dgCMatrix <- function(x) {
  UseMethod("fm_as_dgCMatrix")
}

#' @param unique logical; if `TRUE`, ensures that the sparse triplet
#' representation has a single entry for each non-zero matrix element.
#' @rdname fmesher_sparse
#' @returns `fm_as_dgTMatrix` returns a [Matrix::dgTMatrix-class] object.
#' @export
fm_as_dgTMatrix <- function(x, unique = TRUE, ...) {
  UseMethod("fm_as_dgTMatrix")
}

#' @rdname fmesher_sparse
#' @returns `fm_as_unpackedMatrix` returns an object of virtual class
#'  [Matrix::unpackedMatrix-class].
#' @export
fm_as_unpackedMatrix <- function(x) {
  UseMethod("fm_as_unpackedMatrix")
}

#' @rdname fmesher_sparse
#' @returns `fm_as_fmesher_sparse` returns an `fmesher_sparse` object.
#' @export
fm_as_fmesher_sparse <- function(x) {
  x <- fm_as_dgTMatrix(x, unique = TRUE)
  y <- structure(
    list(
      i = slot(x, name = "i"),
      j = slot(x, name = "j"),
      x = slot(x, name = "x"),
      dims = slot(x, name = "Dim")
    ),
    class = "fmesher_sparse"
  )
  y
}

#' @rdname fmesher_sparse
#' @export
fm_as_dgCMatrix.default <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    x
  } else {
    as(as(as(x, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  }
}

#' @rdname fmesher_sparse
#' @export
fm_as_dgCMatrix.fmesher_sparse <- function(x) {
  Matrix::sparseMatrix(
    i = x[["i"]] + 1L,
    j = x[["j"]] + 1L,
    x = x[["x"]],
    dims = x[["dims"]],
    repr = "C"
  )
}

#' @rdname fmesher_sparse
#' @export
fm_as_dgTMatrix.default <- function(x, unique = TRUE, ...) {
  if (unique) {
    as(fm_as_dgCMatrix(x), "TsparseMatrix")
  } else {
    if (inherits(x, "dgTMatrix")) {
      x
    } else {
      as(as(as(x, "dMatrix"), "generalMatrix"), "TsparseMatrix")
    }
  }
}

#' @rdname fmesher_sparse
#' @export
fm_as_unpackedMatrix.default <- function(x) {
  as(x, "unpackedMatrix")
}

#' @rdname fmesher_sparse
#' @export
fm_as_unpackedMatrix.fmesher_sparse <- function(x) {
  as(
    Matrix::sparseMatrix(
      i = x[["i"]] + 1L,
      j = x[["j"]] + 1L,
      x = x[["x"]],
      dims = x[["dims"]],
      repr = "C"
    ),
    "unpackedMatrix"
  )
}


#' @rdname fmesher_sparse
#' @export
fm_as_dgTMatrix.fmesher_sparse <- function(x, unique = TRUE, ...) {
  Matrix::sparseMatrix(
    i = x[["i"]] + 1L,
    j = x[["j"]] + 1L,
    x = x[["x"]],
    dims = x[["dims"]],
    repr = "T"
  )
}



#' Row-wise Kronecker products
#'
#' Takes two Matrices and computes the row-wise Kronecker product.  Optionally
#' applies row-wise weights and/or applies an additional 0/1 row-wise Kronecker
#' matrix product.
#'
#' @param M1 A matrix that can be transformed into a sparse Matrix.
#' @param M2 A matrix that can be transformed into a sparse Matrix.
#' @param repl An optional index vector.  For each entry, specifies which
#' replicate the row belongs to, in the sense used in
#' `INLA::inla.spde.make.A`
#' @param n.repl The maximum replicate index, in the sense used in
#' `INLA::inla.spde.make.A()`.
#' @param weights Optional scaling weights to be applied row-wise to the
#' resulting matrix.
#' @return A `Matrix::sparseMatrix` object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @export
#' @examples
#' fm_row_kron(rbind(c(1, 1, 0), c(0, 1, 1)), rbind(c(1, 2), c(3, 4)))
#'
fm_row_kron <- function(M1, M2, repl = NULL, n.repl = NULL, weights = NULL # ,
                        # method. = 1
) {
  if (!inherits(M1, "Matrix")) {
    M1 <- as(M1, "Matrix")
  }
  if (!inherits(M2, "Matrix")) {
    M2 <- as(M2, "Matrix")
  }
  n1 <- nrow(M1)
  n2 <- nrow(M2)
  if ((n1 == 1) && (n2 > 1)) {
    M1 <- Matrix::kronecker(rep(1, n2), M1)
    n <- n2
  } else if ((n1 > 1) && (n2 == 1)) {
    M2 <- Matrix::kronecker(rep(1, n1), M2)
    n <- n1
  } else if (n1 != n2) {
    stop(paste0("Size mismatch for row.kron, (n1, n2) = (", n1, ", ", n2, ")"))
  } else {
    n <- n1
  }
  if (is.null(repl)) {
    repl <- rep(1L, n)
  }
  if (is.null(n.repl)) {
    n.repl <- max(repl)
  }
  if (is.null(weights)) {
    weights <- rep(1, n)
  } else if (length(weights) == 1L) {
    weights <- rep(weights[1], n)
  }

  ## OK: Checked robustness for all-zero rows 2022-10-20, matrix 1.5-2
  ## TODO: Maybe move big sparseMatrix call outside the loop.
  ## TODO: Automatically choose M1 or M2 for looping.

  M1 <- fm_as_dgTMatrix(M1, unique = TRUE)
  M2 <- fm_as_dgTMatrix(M2, unique = TRUE)
  n1 <- (as.vector(Matrix::sparseMatrix(
    i = 1L + M1@i, j = rep(1L, length(M1@i)),
    x = 1L, dims = c(n, 1)
  )))
  n2 <- (as.vector(Matrix::sparseMatrix(
    i = 1L + M2@i, j = rep(1L, length(M2@i)),
    x = 1L, dims = c(n, 1)
  )))

  # if (identical(method., 1)) {
  #   M <- (Matrix::sparseMatrix(
  #     i = integer(0), j = integer(0), x = numeric(0),
  #     dims = c(n, ncol(M2) * ncol(M1) * n.repl)
  #   ))
  # } else {
  iii <- integer(0)
  jjj <- integer(0)
  xxx <- numeric(0)
  # }

  n1 <- n1[1L + M1@i]
  for (k in unique(n1)) {
    sub <- which(n1 == k)
    n.sub <- length(sub)

    i.sub <- 1L + M1@i[sub]
    j.sub <- 1L + M1@j[sub]
    o1 <- order(i.sub, j.sub)
    jj <- rep(seq_len(k), times = n.sub / k)
    i.sub <- i.sub[o1]
    j.sub <- (Matrix::sparseMatrix(
      i = i.sub,
      j = jj,
      x = j.sub[o1],
      dims = c(n, k)
    ))
    x.sub <- (Matrix::sparseMatrix(
      i = i.sub,
      j = jj,
      x = weights[i.sub] * M1@x[sub][o1],
      dims = c(n, k)
    ))
    sub2 <- which(is.element(1L + M2@i, i.sub))

    if (length(sub2) > 0) {
      i <- 1L + M2@i[sub2]
      ii <- rep(i, times = k)
      repl.i <- repl[ii]

      # if (identical(method., 1)) {
      #   M <- (M +
      #           Matrix::sparseMatrix(
      #             i = ii,
      #             j = (1L + rep(M2@j[sub2], times = k) +
      #                    ncol(M2) * (as.vector(j.sub[i, ]) - 1L) +
      #                    ncol(M2) * ncol(M1) * (repl.i - 1L)),
      #             x = (rep(M2@x[sub2], times = k) *
      #                    as.vector(x.sub[i, ])),
      #             dims = c(n, ncol(M2) * ncol(M1) * n.repl)
      #           ))
      # } else {
      iii <- c(iii, ii)
      jjj <- c(
        jjj,
        (1L + rep(M2@j[sub2], times = k) +
          ncol(M2) * (as.vector(j.sub[i, ]) - 1L) +
          ncol(M2) * ncol(M1) * (repl.i - 1L))
      )
      xxx <- c(
        xxx,
        (rep(M2@x[sub2], times = k) *
          as.vector(x.sub[i, ]))
      )
      # }
    }
  }

  #  if (!identical(method., 1)) {
  M <- Matrix::sparseMatrix(
    i = iii, j = jjj, x = xxx,
    dims = c(n, ncol(M2) * ncol(M1) * n.repl)
  )
  #  }

  return(M)
}


# @title Find S3 method supported classes
# @description Calls `utils::.S3Methods` and extracts the class information
# as a character vector
# @param f character; the name of an S3 generic
# @keyword internal
method_classes <- function(f) {
  gsub(
    pattern = paste0("^", f, "\\.([^*]*)\\*?"),
    replacement = "\\1",
    x = format(utils::.S3methods(f))
  )
}
