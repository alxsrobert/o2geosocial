## This function was contributed by Rich Fitzjohn. It modifies default arguments
## using user-provided values. The argument 'strict' triggers and error
## behaviour: if strict==TRUE: all new values need to be part of the defaults.

modify_defaults <- function (defaults, x, strict = TRUE) {
  extra <- setdiff(names(x), names(defaults))
  if (strict && (length(extra) > 0L)) {
    stop("Additional invalid options: ", paste(extra, collapse = ", "))
  }
  utils::modifyList(defaults, x, keep.null = TRUE)
}
