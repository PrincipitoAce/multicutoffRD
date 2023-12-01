#' Creates visualiation of new and old cutoff values based on output of optimalcutoffs()
#'
#' @param data A data frame with the following values: Original Cutoffs, New Cutoffs, & Difference
#' @param color A palette of 2 colors
#' @return Graph of new and old cutoff values
#'
#' @examples
#' cutoffsviz(NOT DONE)
#'
#' @export

cutoffsviz = function(data, color){
  # Tidying data for plotting
  toPlot <-
    data |>
    dplyr::select(-`Difference`) |>
    dplyr::mutate(Group = row_number()) |>
    tidyr::gather(`Cutoff Type`, Value, -Group)

  # Plotting
  toPlot |>
    ggplot2::ggplot(mapping = aes(x = Group, y = Value, color = `Cutoff Type`)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = c("New Cutoffs" = color[1], "Original Cutoffs" = color[2]))
}
