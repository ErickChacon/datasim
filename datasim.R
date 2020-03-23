
roxygen2::roxygenize("..", roclets = c("collate", "namespace", "rd"))
devtools::check_built("..")

devtools::load_all("..")
