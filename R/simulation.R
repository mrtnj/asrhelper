


#' Simulate a founder genome based on a genome table.
#'
#' @param genome_table Data frame of genome parameters, must contain columns
#' "length" with chromosome length in bp and "genetic lengths" with genetic
#' length in cM
#' @param n_ind Number of individual
#' @param seg_sites_per_bp Optional, number of segregating sites to be kept per
#' basepair. NULL for all
#' @param final_Ne Final effective population size. Passed to AlphaSimR.
#' @param historical_Ne Optional historical population sizes. Passed to AlphaSimR.
#' @param historical_Ne_time Optional timings for population size changes.
#' Passed to AlphaSimR.
#' @param split Optional split parameter. Passed to AlphaSimR.
#'
#' @return An AlphaSimR PopMap object to be used as founder of a simulation
#' @export
simulate_founder_genome <- function(genome_table,
                                    n_ind,
                                    seg_sites_per_bp,
                                    final_Ne,
                                    historical_Ne,
                                    historical_Ne_time,
                                    split) {

  if (!is.null(seg_sites_per_bp)) {
    seg_sites <- genome_table$length /1e6 * seg_sites_per_bp
  } else {
    seg_sites <- NULL
  }

  n_chr <- nrow(genome_table)

  pop_chr <- vector(mode = "list",
                    length = n_chr)

  for (chr_ix in 1:n_chr) {

    pop_chr[[chr_ix]] <-
      AlphaSimR::runMacs2(bp = genome_table$length[chr_ix],
                          genLen = genome_table$genetic_length[chr_ix] / 100,
                          nInd = n_ind,
                          nChr = 1,
                          segSites = seg_sites,
                          Ne = final_Ne,
                          histNe = historical_Ne,
                          histGen = historical_Ne_time,
                          split = split)

    ## Fix names
    sim_map <- pop_chr[[chr_ix]]@genMap
    old_names <- names(sim_map[[1]])
    names_split <- strsplit(old_names, split = "_")
    names_position <- unlist(lapply(names_split, "[", 2))
    names(pop_chr[[chr_ix]]@genMap) <- chr_ix
    names(pop_chr[[chr_ix]]@genMap[[1]]) <- paste(chr_ix, names_position, sep = "_")
  }

  pop <- Reduce(AlphaSimR::cChr, pop_chr)

  pop

}


#' Forward-simulate a population with a given sequence of population sizes.
#'
#' @param pop Starting population (AlphaSimR population object)
#' @param pop_sizes Vector of population sizes for each generation
#' @param simparam Simulation parameters (AlphaSimR SimParam object)
#'
#' @return A list of populations, representing each generation.
#' @export
forward_simulate_population <- function(pop,
                                        pop_sizes,
                                        simparam) {

  pop_sizes <- c(pop@nInd, pop_sizes)
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

  generations[-1]

}

