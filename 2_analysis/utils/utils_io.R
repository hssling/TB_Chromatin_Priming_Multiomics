# utils_io.R
suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
})

read_config <- function(path = "config/config.yaml") {
  yaml::read_yaml(path)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

write_table <- function(df, path) {
  ensure_dir(dirname(path))
  data.table::fwrite(df, path)
}
