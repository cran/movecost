#' Defunct functions in movecost
#'
#' The user-facing functions of the movecost 2.x series were removed in
#' movecost 3.0.0 as part of the architectural redesign (a single cost surface
#' is now built once with \code{\link{mc_surface}} and reused by every analysis
#' function). The 2.x names below are retained only as defunct stubs: each one
#' raises an error naming its replacement, so that code written for the 2.x
#' series fails with an informative message rather than an opaque
#' "could not find function" error. See the 2.x-to-3.0 correspondence table in
#' the \code{NEWS} file for the full mapping.
#'
#' @param ... Arguments of the former 2.x functions; accepted but ignored.
#' @return These functions are defunct and return nothing; each one stops with
#'   an error that directs you to its 3.0 replacement.
#' @name movecost-defunct
#' @keywords internal
NULL

# Each stub calls .Defunct(), which raises an error of class "defunctError".
# Where a single 3.0 function replaces the old one, its name is passed as the
# first argument so that help redirection points at the right page; movecost()
# mapped to two steps, so it carries a bespoke message instead.

#' @rdname movecost-defunct
#' @export
movecost <- function(...) {
  .Defunct(msg = paste("movecost() was removed in movecost 3.0.0; build the cost",
                       "surface once with mc_surface(), then call mc_accum() or",
                       "mc_paths(). See the 2.x-to-3.0 table in NEWS.md."))
}

#' @rdname movecost-defunct
#' @export
movecorr <- function(...) {
  .Defunct("mc_corridor",
           msg = "movecorr() was removed in movecost 3.0.0; use mc_corridor(). See NEWS.md.")
}

#' @rdname movecost-defunct
#' @export
movealloc <- function(...) {
  .Defunct("mc_alloc",
           msg = "movealloc() was removed in movecost 3.0.0; use mc_alloc(). See NEWS.md.")
}

#' @rdname movecost-defunct
#' @export
movebound <- function(...) {
  .Defunct("mc_boundary",
           msg = "movebound() was removed in movecost 3.0.0; use mc_boundary(). See NEWS.md.")
}

#' @rdname movecost-defunct
#' @export
movenetw <- function(...) {
  .Defunct("mc_network",
           msg = "movenetw() was removed in movecost 3.0.0; use mc_network(). See NEWS.md.")
}

#' @rdname movecost-defunct
#' @export
moverank <- function(...) {
  .Defunct("mc_rank",
           msg = "moverank() was removed in movecost 3.0.0; use mc_rank(). See NEWS.md.")
}

#' @rdname movecost-defunct
#' @export
movecomp <- function(...) {
  .Defunct("mc_comp",
           msg = "movecomp() was removed in movecost 3.0.0; use mc_comp(). See NEWS.md.")
}
