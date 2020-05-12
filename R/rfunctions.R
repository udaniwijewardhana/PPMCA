#' Runs the SpatialEpiApp Shiny web application.
#' @export
run_app <- function() {
  shiny::runApp(system.file('PPMCA', package='PPMCA'))
}
