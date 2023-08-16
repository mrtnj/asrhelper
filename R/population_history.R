#' Turn a data frame of population sizes generations at generations ago into
#' a vector of population sizes for forward simulation.
#'
#' @param history A data frame with columns Ne and generations_ago.
#'
#' @return A vector of population sizes.
#' @export
generations_ago_to_pop_size <- function(history) {
  last_ix <- nrow(history)

  result <- vector(mode = "list",
                   length = last_ix)
  k <- 1
  for (event_ix in last_ix:1) {
    if (event_ix > 1) {
      repeat_generations <-
        history$generations_ago[event_ix] -
        history$generations_ago[event_ix - 1]
    } else {
      repeat_generations <- history$generations_ago[event_ix]
    }

    result[[k]] <- rep(history$Ne[event_ix],
                       repeat_generations)
    k <- k + 1
  }

  unlist(result)
}
