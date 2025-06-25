#' @keywords internal
get_y0 <- function(to_plot, t, parent) {

  index <- min(which(to_plot$ID == parent))
  entry <- to_plot[index, ]
  # now we have the line
  slope <- (entry$y1[[1]] - entry$y0[[1]]) / (entry$x1[[1]] - entry$x0[[1]])
  yval <- entry$y0 + slope * (t - entry$x0[[1]])
  return(yval)
}


#' @keywords internal
count_daughters <- function(ltab, id) {
  daughters <- ltab[which(ltab[, 2] == id), 3]
  num_daughters <- 0
  if (length(daughters) > 0) {
    for (id2 in daughters) {
      num_daughters <- num_daughters + 1 + count_daughters(ltab, id2)
    }
  }
  return(num_daughters)
}

update_pos <- function(ltable, ID, pos_table) {
  num_tips_connected <- 1 + count_daughters(ltable, ID)
  parent <- ltable[abs(ID), 2]
  # we now need to claim tip space
  available <- which(pos_table[, 1] == parent)

  parent_direction <- ltable[abs(parent), ncol(ltable)]

  to_claim <- available[1:num_tips_connected]
  if (parent_direction == -1) {
    to_claim <- tail(available, num_tips_connected)
  }
  pos_table[to_claim, 1] <- ID

  chosen_pos <- min(to_claim)
  if (parent_direction == -1) chosen_pos <- max(to_claim)

  y1 <- pos_table[chosen_pos, 2]
  pos_table[chosen_pos, 1] <- sign(ID) * 1e6 + sign(ID) * ID
  return(list("table" = pos_table,
              y1 = y1,
              parent_direction = parent_direction))
}

#' plot ltable as tree
#' @param ltable ltable
#' @export
plot_ltable <- function(ltable,
                        plot_hybrid = FALSE,
                        use_factor_ID = FALSE,
                        modify_theme = TRUE) {

  # make time increase again
  if (ltable[1, 1] == min(ltable[, 1])) {
    ltable[, 1] <- -1 * ltable[, 1]
    extinct_lineages <- which(ltable[, 4] != -1)
    ltable[extinct_lineages, 4] <- -1 * ltable[extinct_lineages, 4]
  }

  max_x <- 0
  min_x <- ltable[1, 1]
  ltable <- cbind(ltable,  sign(ltable[, 3]))

  num_tips <- length(ltable[, 1])

  tip_dist <- 10

  total_y <- num_tips * tip_dist

  num_tips_right <- length(which(ltable[, 3] > 0))
  num_tips_left <- length(which(ltable[, 3] < 0))

  max_y <- (num_tips_right - 1) * tip_dist
  min_y <- -1 * (tip_dist + (num_tips_left - 1) * tip_dist)

  root_pos <- (max_y + min_y) / 2

  right_pos <- cbind(rep(2, num_tips_right), seq(0, max_y, by = tip_dist))
  left_pos  <- cbind(rep(-1, num_tips_left), seq(min_y, -tip_dist, by = tip_dist))

  # create list of lines
  to_plot <- tibble::tibble(
    x0 = numeric(),
    x1 = numeric(),
    y0 = numeric(),
    y1 = numeric(),
    ID = numeric()
  )

  for (i in 1:2) {
    entry <- ltable[i, ]
    ID <- entry[3]
    parent <- entry[2]
    if (parent %in% c(0, -1)) { # crown
      x0 <- entry[1]
      x1 <- 0.0
      if (entry[4] != -1) x1 <- entry[4]

      y0 <- root_pos
      y1 <- tip_dist

      if (ID < 0) {
        potential_pos <- which(left_pos[, 1] == -1)
        chosen_pos <- min(potential_pos)
        y1 <- left_pos[chosen_pos, 2]
        left_pos[chosen_pos, 1] <- -1e6 - 1
      } else {
        # right tips
        potential_pos <- which(right_pos[, 1] == 2)
        chosen_pos <- max(potential_pos)
        y1 <- right_pos[chosen_pos, 2]
        right_pos[chosen_pos, 1] <- 1e6 + 2
      }

      to_add <- tibble::tibble(
        x0 = x0,
        x1 = x1,
        y0 = y0,
        y1 = y1,
        ID = ID
      )

      to_plot <- rbind(to_plot, to_add)
    }
  }

  if (nrow(ltable) > 2) {

    for (i in 3:nrow(ltable)) {
      entry <- ltable[i, ]
      ID <- entry[3]
      parent <- entry[2]
      x0 <- entry[1]
      x1 <- 0.0
      if (entry[4] != -1) x1 <- entry[4]

      y0 <- get_y0(to_plot, x0, parent)

      parent_direction <- 1
      if (ID < 0) {
        res <- update_pos(ltable, ID, left_pos)
        left_pos <- res$table
        y1 <- res$y1
        parent_direction <- res$parent_direction
      } else {
        res <- update_pos(ltable, ID, right_pos)
        right_pos <- res$table
        y1 <- res$y1
        parent_direction <- res$parent_direction
      }

      # this is a daughter branch, we need to figure out how many daughters there
      # are
      ltable[abs(ID), ncol(ltable)] <- -parent_direction

      to_add <- tibble::tibble(
        x0 = x0,
        x1 = x1,
        y0 = y0,
        y1 = y1,
        ID = ID
      )
      to_plot <- rbind(to_plot, to_add)

      if (plot_hybrid == TRUE) {
        parent2 <- entry[5]
        # check if there is another parent
        if (parent2 != parent) {
         y0 <- get_y0(to_plot, x0, parent2)
         to_add <- tibble::tibble(
           x0 = x0,
           x1 = x1,
           y0 = y0,
           y1 = y1,
           ID = ID
         )
         to_plot <- rbind(to_plot, to_add)
        }
      }

    }
  }

  p <- ggplot2::ggplot(to_plot)

  if (use_factor_ID) {
    p <- p + ggplot2::geom_segment(ggplot2::aes(x = .data[["x0"]],
                                                y = .data[["y0"]],
                                                xend = .data[["x1"]],
                                                yend = .data[["y1"]],
                                                col = as.factor(.data[["ID"]])))
  } else {
    p <- p + ggplot2::geom_segment(ggplot2::aes(x = .data[["x0"]],
                                                y = .data[["y0"]],
                                                xend = .data[["x1"]],
                                                yend = .data[["y1"]],
                                                col = .data[["ID"]]))
  }
  if (modify_theme) {
    p <- p + ggplot2::xlim(min_x, max_x) +
      ggplot2::coord_flip() +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.line.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()) +
      ggplot2::ylab("") +
      ggplot2::xlab("Time before present")
  }
  return(p)
}
