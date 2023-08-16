founders <- AlphaSimR::quickHaplo(nInd = 1000,
                                  nChr = 10,
                                  segSites = 1000)
simparam <- AlphaSimR::SimParam$new(founders)
pop <- AlphaSimR::newPop(founders, simparam)
pop_sizes <- rep(100, 10)
pop_forward <- forward_simulate_population(pop,
                                           pop_sizes = pop_sizes,
                                           simparam = simparam)


test_that("returns list of right length", {

  expect_vector(pop_forward,
                ptype = list(),
                size = 10)
})


test_that("list contains pops", {

  for (gen_ix in 1:length(pop_forward)) {
    expect_s4_class(pop_forward[[gen_ix]], "Pop")
  }

})


test_that("generation has right size", {

  for (gen_ix in 1:length(pop_forward)) {
    expect_equal(pop_forward[[gen_ix]]@nInd, pop_sizes[gen_ix])
  }

})






genome_table <- data.frame(chr = 1:3,
                           length = rep(1e6, 3),
                           genetic_length = rep(1, 3))

founders_genome <- simulate_founder_genome(genome_table,
                                           n_ind = 10,
                                           final_Ne = 1000,
                                           seg_sites_per_bp = 30,
                                           historical_Ne = NULL,
                                           historical_Ne_time = NULL,
                                           split = NULL)

test_that("right number of chromosomes and sites", {

  expect_equal(founders_genome@nChr, 3)
  expect_equal(founders_genome@nLoci, c(30, 30, 30))

})
