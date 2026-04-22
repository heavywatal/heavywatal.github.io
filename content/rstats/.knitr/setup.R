options(
  mc.cores = parallel::detectCores(),
  wtl.printdf.summarize = FALSE,
  wtl.printdf.classes = FALSE,
  pillar.print_max = 8L,
  tibble.print_max = 8L,
  width = 9999L,
  cli.width = 60L,
  dplyr.summarise.inform = FALSE,
  readr.num_columns = 0L,
  readr.show_progress = FALSE,
  readr.show_col_types = FALSE
)
library(tibble) # nolint: unused_import_linter.
library(ggplot2)
theme_set(wtl::theme_wtl())
registerS3method("print", "tbl", wtl::printdf)
registerS3method("print", "tbl_df", wtl::printdf)
knitr::opts_chunk$set(comment = "")
knitr::opts_chunk$set(dev = "ragg_png")
knitr::opts_chunk$set(dpi = 108)
knitr::opts_chunk$set(fig.process = wtl::oxipng)
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE)
litedown::reactor(
  comment = "",
  message = NA,
  warning = NA,
  fig.path = "../figure/",
  dev = ragg::agg_png,
  dev.args = list(res = 108),
  cache = TRUE,
  attr.source = "r",
  print = NA
)
set.seed(24601)

src_alt_fig_chunk = function(label, ext = "png", number = 1L) {
  src = knitr::fig_chunk(label, ext, number)
  paste0('src="', src, '" alt="plot of chunk ', label, '"')
}
