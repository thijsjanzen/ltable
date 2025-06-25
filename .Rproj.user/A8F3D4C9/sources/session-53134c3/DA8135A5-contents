#' Function to simulate a tree, conditional on observing all states.
#'
#' @param parameters vector with two parameters consisting of per species rates
#' of: speciating, extinction.
#' @param crown_age age of the crown
#' @param max_spec condition on tree not being larger than this
#' @param min_spec condition on tree not being smaller than this
#' @param max_species_extant should the maximum number of species only include
#' extant species (if TRUE), or also extinct species?
#' @param tree_size_hist returns a histogram of tree sizes, including failed
#' trees (e.g. trees that went extinct or were outside conditioning criteria)
#' @param non_extinction condition on non-extinction if TRUE
#' @param verbose provide verbose output
#' @param max_tries maximum number of tries before giving up
#' @param drop_extinct should extinct lineages be returned?
#' @param start_at_crown should the simulation start at the root, or the crown?
#' @param seed random seed
#' By default, the algorithm keeps simulating until it generates a tree where
#' both crown lineages survive to the present - this is to ensure that the tree
#' has a crown age that matches the used crown age. You can modify
#' 'non-extinction' to deviate from this behaviour.
#' @export
#'
#' @rawNamespace useDynLib(ltable, .registration = TRUE)
#' @rawNamespace import(Rcpp)
sim_bd <- function(parameters,
                crown_age,
                max_spec = 1e5,
                min_spec = 2,
                max_species_extant = TRUE,
                tree_size_hist = FALSE,
                non_extinction = TRUE,
                verbose = FALSE,
                max_tries = 1e6,
                drop_extinct = FALSE,
                start_at_crown = TRUE,
                seed = NULL) {

  if (is.null(seed)) seed <- -1

  res <- sim_bd_cpp(parameters,
                         crown_age,
                         max_spec,
                         max_species_extant,
                         min_spec,
                         non_extinction,
                         verbose,
                         max_tries,
                         seed,
                         tree_size_hist,
                         start_at_crown)

  if (length(res) < 1) { # this happens upon a throw
    return(list(phy = "ds",
                traits = 0))
  }

  Ltable        <- res$ltable

  out_hist <- 0
  if (tree_size_hist == TRUE) out_hist <- res$hist_tree_size

  if (start_at_crown == FALSE && sum(Ltable[, 4] == -1) == 1) {
    # fake phy
    phy <- ape::rphylo(n = 2, birth = 0.2, death = 0)
    phy$edge.length[-2]
    phy$tip.label <- phy$tip.label[-2]
    phy$edge <- phy$edge[-2, ]
    phy$edge <- matrix(data = phy$edge, nrow = 1) # important!
    phy$edge.length <- crown_age # this is now root age
    return(list(phy = phy,
                traits = res$traits[[1]],
                extinct = res$tracker[2],
                overshoot = res$tracker[3],
                conditioning = res$tracker[4],
                small = res$tracker[6],
                size_hist = out_hist))


  } else if (sum(Ltable[, 4] == -1) < 2) {
    warning("crown lineages died out")
    return(list(phy = "ds",
                traits = 0,
                extinct = res$tracker[2],
                overshoot = res$tracker[3],
                conditioning = res$tracker[4],
                small = res$tracker[6],
                size_hist = out_hist))
  }

  if (sum(res$tracker) >= max_tries) {
    warning("Couldn't simulate a tree in enough tries,
            try increasing max_tries")

    return(list(phy = "ds",
                traits = 0,
                extinct = res$tracker[2],
                overshoot = res$tracker[3],
                conditioning = res$tracker[4],
                small = res$tracker[6],
                size_hist = out_hist))
  }


  Ltable[, 1]   <- crown_age - Ltable[, 1] # simulation starts at 0,
  # not at crown age
  notmin1 <- which(Ltable[, 4] != -1)
  Ltable[notmin1, 4] <- crown_age - c(Ltable[notmin1, 4])
  Ltable[which(Ltable[, 4] == crown_age + 1), 4] <- -1

  phy <- treestats::l_to_phylo(Ltable, drop_extinct = drop_extinct)

  if (sum(Ltable[, 4] < 0)) {
    return(list(phy = phy,
                extinct = res$tracker[2],
                overshoot = res$tracker[3],
                conditioning = res$tracker[4],
                small = res$tracker[6],
                size_hist = out_hist,
                ltable = Ltable))
  } else {
    warning("simulation did not meet minimal requirements")
    return(list(phy = "ds",
                traits = 0,
                extinct = res$tracker[2],
                overshoot = res$tracker[3],
                conditioning = res$tracker[4],
                small = res$tracker[6],
                size_hist = out_hist))
  }
}
