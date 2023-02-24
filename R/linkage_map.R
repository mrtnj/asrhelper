##################################################################

## Using real genetic maps to adjust simulated maps



## Linear interpolation of genetic marker position

interpolate_cM <- function(marker_pos,
                           bp_before,
                           bp_after,
                           cM_before,
                           cM_after) {

  delta_bp <- bp_after - bp_before
  delta_cM <- cM_after - cM_before

  cM_per_bp <- delta_cM/delta_bp

  cM_before +  cM_per_bp * (marker_pos - bp_before)
}



## Get positions in cM based on a real map for simulated map, which is taken
## to show the relative position in basepairs

get_cM_positions_chr <- function(sim_map_chr,
                                 real_map_chr) {

  chr_length_bp <- max(real_map_chr$position_bp)
  chr_length_morgan <- max(real_map_chr$position_cM) / 100
  n_markers <- length(sim_map_chr)

  cM <- numeric(n_markers)

  for (marker_ix in 1:n_markers) {

    marker_pos <- chr_length_bp * sim_map_chr[marker_ix] / chr_length_morgan

    ix_before <- max(which(real_map_chr$position_bp < marker_pos))
    ix_after <- min(which(real_map_chr$position_bp > marker_pos))

    ## First marker on chromosome
    if (ix_after == 1) {
      cM_position <- interpolate_cM(marker_pos,
                                    bp_before = 0,
                                    bp_after = real_map_chr$position_bp[ix_after],
                                    cM_before = 0,
                                    cM_after = real_map_chr$position_cM[ix_after])
    } else {
      ## Most cases
      cM_position <- interpolate_cM(marker_pos,
                                    bp_before = real_map_chr$position_bp[ix_before],
                                    bp_after = real_map_chr$position_bp[ix_after],
                                    cM_before = real_map_chr$position_cM[ix_before],
                                    cM_after = real_map_chr$position_cM[ix_after])
    }

    cM[marker_ix] <- cM_position

  }

  cM
}


#' Adjust a simulated linkage map by interpolating positions on a real
#' linkage map
#'
#' @param sim_map List of maps for each chromosome in Morgans
#' @param real_map List of maps for each chromosome in Morgans
#'
#' @return A list containing the modified map
#' @export
make_adjusted_map <- function(sim_map,
                              real_map) {

  new_map <- vector(length = length(sim_map),
                    mode = "list")

  n_chr <- length(new_map)

  for (chr_ix in 1:n_chr) {
    new_map[[chr_ix]] <- get_cM_positions_chr(sim_map[[chr_ix]],
                                              real_map[[chr_ix]])
    ## Fix names
    old_names <- names(sim_map[[chr_ix]])
    names_split <- strsplit(old_names, split = "_")
    names_position <- unlist(lapply(names_split, "[", 2))
    names(new_map[[chr_ix]]) <- paste(chr_ix, names_position, sep = "_")
  }

  names(new_map) <- 1:n_chr

  new_map
}


