#' Launch the SurprisalAnalysis Shiny App
#'
#' @param port port to run the app on (passed to shiny::runApp)
#' @param host host to listen on
#' @param launch.browser should launch a browser? set to TRUE by default
#' @param ... Further arguments passed along to shiny::runApp
#' @export
runSurprisalApp <- function(port      = getOption("shiny.port", 3838),
                            host      = getOption("shiny.host", "127.0.0.1"),
                            launch.browser = getOption("shiny.launch.browser", TRUE),
                            ...) {
  app_dir <- system.file("shiny", package = "SurprisalAnalysis")
  if (app_dir == "" || !dir.exists(app_dir)) {
    stop("Cannot find the Shiny app directory. Reinstall the package?", call. = FALSE)
  }
  shiny::runApp(appDir = app_dir,
                port      = port,
                host      = host,
                launch.browser = launch.browser,
                ...)
}
