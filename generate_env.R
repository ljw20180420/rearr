library(rix)

rix(
  r_ver = "frozen-edge",
  r_pkgs = c("languageserver", "devtools", "tidyverse"),
  system_pkgs = NULL,
  git_pkgs = NULL,
  ide = "none",
  project_path = ".",
  overwrite = TRUE,
  print = TRUE
)
