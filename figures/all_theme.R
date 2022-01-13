all_theme <- function() {
  theme_tufte(base_size = 16) %+replace%
    theme(
      plot.title = element_text(
        size = 12,
        family = "Arial",
        face = "bold"),
      axis.title = element_text(
        size = 8,
        family = "Arial"),
      axis.text = element_text(
        size = 7,
        family = "Arial"),
      strip.text.x = element_text(
        size = 7,
        family = "Arial"),
      legend.title = element_text(
        size = 8,
        family = "Arial"),
      legend.text = element_text(
        size = 7,
        family = "Arial"),
      axis.line = element_line(
        size = 0.3,
        linetype = "solid"),
      axis.ticks = element_line(
        size = 0.3,
        linetype = "solid"),
      complete = TRUE
    )
}
