#' Launch the SurprisalAnalysis Shiny App
#'
#' @param port port to run the app on (passed to shiny::runApp)
#' @param host host to listen on
#' @param launch.browser should launch a browser? set to TRUE by default
#' @param run boolean value, is set to TRUE by default. If set to FALSE it will not launch the graphical user interface
#' @param ... Further arguments passed along to shiny::runApp
#' @return no return value, running the function will launch an application with graphical user interface
#' @importFrom shiny fluidPage navbarPage tabPanel sidebarLayout sidebarPanel
#' @importFrom shiny fileInput selectInput checkboxInput sliderInput actionButton
#' @importFrom shiny strong em mainPanel div p tags br fluidRow column plotOutput renderPlot
#' @importFrom shiny downloadButton downloadHandler req observeEvent updateSliderInput updateTabsetPanel icon HTML
#' @importFrom shinythemes shinytheme
#' @importFrom shinyjs useShinyjs disable enable
#' @importFrom shinycssloaders withSpinner
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_x_continuous scale_y_continuous
#'   labs theme_minimal theme element_blank element_line element_text geom_bar coord_flip
#'   scale_fill_gradient xlab ylab
#' @importFrom latex2exp TeX
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#' @importFrom clusterProfiler enrichGO
#' @import patchwork
#' @import shinyWidgets
#' @import dplyr
#' @import tidyr
#' @import tidyverse
#' @import httpuv
#' @importFrom stats quantile
#' @importFrom utils head
#' @examples
#'
#' runSurprisalApp(port = httpuv::randomPort(), run = FALSE)
#'
#' @export runSurprisalApp
runSurprisalApp <- function(port      = getOption("shiny.port", 3838),
                            host      = getOption("shiny.host", "127.0.0.1"),
                            launch.browser = getOption("shiny.launch.browser", TRUE),
                            run = TRUE,
                            ...) {
  app_dir <- system.file("shiny", package = "SurprisalAnalysis")

  if (app_dir == "" || !dir.exists(app_dir)) {
    stop("Cannot find the Shiny app directory. Reinstall the package?", call. = FALSE)
  }

  app <- shiny::shinyAppDir(app_dir)
  if (isTRUE(run)) {
  app <- shiny::runApp(appDir = app_dir,
                port      = port,
                host      = host,
                launch.browser = launch.browser,
                ...)
  }
  invisible(app)
}
