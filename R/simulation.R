


#' Forward-simulate a population with a given sequence of population sizes.
#'
#' @param pop Starting population (AlphaSimR population object)
#' @param pop_sizes Vector of population sizes for each generation
#' @param simparam Simulation parameters (AlphaSimR SimParam object)
#'
#' @return A list of populations, representing each generation.
forward_simulate_population <- function(pop,
                                        pop_sizes,
                                        simparam) {

  n_gen <- length(pop_sizes)
  generations <- vector(mode = "list",
                        length = n_gen)
  generations[[1]] <- pop

  for (gen in 2:n_gen) {
      generations[[gen]] <- AlphaSimR::randCross(pop = generations[[gen - 1]],
                                                 nCrosses = pop_sizes[gen],
                                                 nProgeny = 1,
                                                 simParam = simparam)

      generations[[gen]] <- AlphaSimR::mutate(generations[[gen]],
                                              simParam = simparam)

  }

  generations

}

