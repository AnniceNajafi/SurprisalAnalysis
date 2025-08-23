options(shiny.maxRequestSize = 2000*1024^2)

library(shiny)
library(shinythemes)
library(shinyjs)
library(ggplot2)
library(latex2exp)
library(DT)
library(shinycssloaders)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(patchwork)

find_demo_file <- function(fname, pkg = "SurprisalAnalysis") {
  candidates <- c(
    file.path("data", fname),
    file.path("inst", "shiny", "data", fname),
    system.file("shiny", "data", fname, package = pkg)
  )
  hit <- candidates[file.exists(candidates)][1]
  if (length(hit) == 0 || is.na(hit)) return(NULL)
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

simulate_fileinput_display <- function(id, filename, duration_ms = 700) {
  js <- sprintf("
    (function(){
      var $ct   = $('#%s').closest('.shiny-input-container');
      var $text = $ct.find('input.form-control[type=\"text\"]');
      var $prog = $('#%s_progress');
      var $bar  = $prog.find('.progress-bar');
      if ($text.length) { $text.val('%s'); }
      if ($prog.length) {
        $bar.css('width','0%%').attr('aria-valuenow', 0);
        $prog.removeClass('hidden');
        var val = 0, step = 100 / Math.max(1, Math.floor(%d/30));
        var timer = setInterval(function(){
          val = Math.min(100, val + step);
          $bar.css('width', val + '%%').attr('aria-valuenow', val);
          if (val >= 100) { clearInterval(timer); setTimeout(function(){ $prog.addClass('hidden'); }, 250); }
        }, 30);
      }
    })();
  ", id, id, filename, as.integer(duration_ms))
  shinyjs::runjs(js)
}

detect_id_type <- function(id_vec, sample_n = 500, min_conf = 0.20) {
  x <- unique(toupper(as.character(id_vec)))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return("Gene Symbol")
  if (length(x) > sample_n) { set.seed(1L); x <- sample(x, sample_n) }

  pct <- function(m) mean(m, na.rm = TRUE)
  ens_like    <- pct(grepl("^ENSG\\d+", x))
  entrez_like <- pct(grepl("^\\d+$", x))
  sym_like    <- pct(grepl("^[A-Z0-9\\-\\.]+$", x) & !grepl("^ENSG\\d+|^\\d+$", x))

  if (ens_like    >= 0.80) return("Ensembl")
  if (entrez_like >= 0.80) return("Entrez Id")
  regex_hint <- if (sym_like >= 0.80) "Gene Symbol" else NA_character_

  keytypes_try <- c(ENSEMBL = "Ensembl", ENTREZID = "Entrez Id",
                    SYMBOL  = "Gene Symbol", PROBEID = "Probe Id")

  scores <- vapply(names(keytypes_try), function(kt) {
    m <- tryCatch(
      AnnotationDbi::mapIds(org.Hs.eg.db, keys = x, column = "ENTREZID",
                            keytype = kt, multiVals = "first"),
      error = function(e) rep(NA_character_, length(x))
    )
    mean(!is.na(m))
  }, numeric(1))

  best_kt <- names(scores)[which.max(scores)]
  best_ui <- keytypes_try[[best_kt]]
  if (!is.na(regex_hint) && scores[best_kt] < min_conf) return(regex_hint)
  if (scores[best_kt] >= min_conf) best_ui else "Gene Symbol"
}

# Core linear algebra for surprisal
turning_alphas <- function(log.mat) {
  mat <- as.matrix(log.mat)
  # don't touch zeros here; upstream transform ensures valid logs
  C <- crossprod(mat)
  ed <- eigen(C)
  P <- ed$vectors
  D <- diag(sqrt(ed$values))
  Y <- mat
  n <- ncol(Y)
  alph_lst <- lapply(seq_len(n), function(i) colSums(solve(D))[i] * Y %*% P[, i])
  alph_all <- do.call(cbind, alph_lst)
  hold_lst <- lapply(seq_len(n), function(i) P[, i] * D[i, i])
  holds <- do.call(cbind, hold_lst)
  list(holds, alph_all)
}

# ---- NEW: data-scale detection & transformation helpers ----------------------

# Heuristic: counts-like if nonnegative, large-ish dynamic range, mostly integers.
is_counts_like <- function(mat) {
  x <- as.numeric(mat)
  x <- x[is.finite(x)]
  if (!length(x)) return(FALSE)
  if (any(x < 0, na.rm = TRUE)) return(FALSE)
  int_frac <- mean(abs(x - round(x)) < 1e-8)
  dyn_range <- quantile(x, 0.99, na.rm = TRUE) - quantile(x, 0.01, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  (int_frac >= 0.70 && (max_val >= 50 || dyn_range >= 50))
}

# Apply selected (or auto) transform; returns log-scale numeric matrix
apply_transform <- function(mat, method = c("Auto", "log1p", "pseudocount"), pseudocount = 1e-6) {
  method <- match.arg(method)
  if (method == "Auto") {
    method <- if (is_counts_like(mat)) "log1p" else "pseudocount"
  }
  m <- as.matrix(mat)
  storage.mode(m) <- "double"
  if (method == "log1p") {
    # Works well for counts (including zeros)
    if (any(m < -1, na.rm = TRUE)) stop("log1p transform invalid for values < -1.")
    return(log1p(m))
  } else {
    # Pseudocount log; protect only zeros to preserve scale
    m[m <= 0] <- 0  # avoid taking log of negatives; you can adjust if needed
    return(log(m + pseudocount))
  }
}

# -----------------------------------------------------------------------------

ui <- fluidPage(
  theme = shinytheme("yeti"),
  useShinyjs(),

  tags$style(HTML("
    .navbar-default { background-color: #76ABAE !important; }
    .navbar-default .navbar-nav > li > a { color: white !important; }
    .irs-bar, .irs-bar-edge, .irs-single {
      background: #76ABAE !important;
      border-top: 1px solid #76ABAE !important;
      border-bottom: 1px solid #76ABAE !important;
    }
  ")),

  navbarPage(
    "Surprisal Analysis", id = "sur",

    tabPanel(
      "Instructions",
      sidebarLayout(
        sidebarPanel(
          fileInput("data_in", "Please upload a gene expression matrix csv file."),
          actionButton(
            "load_demo",
            strong(em("Load demo dataset")),
            style = "color:#000;background-color:#76ABAE;border-color:#76ABAE;width:100%;margin-bottom:8px;"
          ),
          helpText("Tip: If you don't have a file, click the button above to auto-load the demo data."),
          div(style="margin-top:6px; font-style:italic;",
              "Ready file: ", textOutput("ready_name", inline = TRUE)),
          tags$hr(),
          selectInput("data_type", "Select biomarker input type.",
                      choices = c("Ensembl","Entrez Id","Gene Symbol","Probe Id")),
          # ---- NEW: transformation control
          selectInput(
            "transform_method",
            "Transformation",
            choices = c("Auto", "log1p", "Log with pseudocount"),
            selected = "Auto"
          ),
          helpText(HTML(
            "Auto: picks <b>log1p</b> for counts-like data; otherwise <b>Log with pseudocount</b>."
          )),
          checkboxInput("lambda1_check", "Flip lambda 1 value"),
          sliderInput("percentile_GO", "GO transcript percentile",
                      min = 80, max = 95, value = 95, step = 1, ticks = FALSE),
          actionButton("submit", strong(em("Submit")),
                       style = "color:#000;background-color:#76ABAE;border-color:#76ABAE")
        ),
        mainPanel(
          div(
            style = "border:1px solid #78A2A5;background:#E4EEEF;padding:20px;text-align:justify;",
            p(strong("Welcome!")),
            p("I. Upload a csv file on the left, or click 'Load demo dataset' to auto-load ",
              strong("helper_T_cell_0_test.csv"), "."),
            p("II. The first column of the csv will be used as gene names and all other columns as samples, in the order provided."),
            p("III. An example of the input data is provided below and under the 'Sample Files' tab."),
            tags$img(src = 'sample_file.png', height = 200, width = 500),
            p("IV. Please ensure that you select the correct identifier type according to your transcript column."),
            p("V. The application performs GO enrichment analysis on the genes that contribute most strongly to each surprisal pattern. Gene “importance” is quantified by their weights in the corresponding transcript weight vector (G₁ for pattern 1, G₂ for pattern 2).

You can use the percentile slider on the left to select the cutoff. For example, if the slider is set to 95%, only genes whose weights fall within the top 5% (95th percentile and above) of the chosen pattern are included in the GO analysis.

If the “Flip λ₁” option is selected, the weights for pattern 1 are inverted before percentile selection, ensuring the “top” genes correspond to the same direction seen in the λ plots."),
            p(strong("Surprisal Analysis Applications:")),
            tags$a(href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0108549",
                   "Kravchenko-Balasha N et al. PLoS ONE. 2014;9(11):e108549."), br(), br(),
            tags$a(href="https://www.pnas.org/doi/full/10.1073/pnas.1414714111",
                   "Zadran S et al. PNAS. 2014;111(36):13235-40."), br(), br(),
            tags$a(href="https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007034",
                   "Su Y et al. PLoS Comput Biol. 2019;15(6):e1007034."), br(), br(),
            tags$a(href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5903653/",
                   "Bogaert KA et al. PLoS ONE. 2018;13(4):e0195142."), br(), br()
          )
        )
      )
    ),

    tabPanel(
      "Results",
      h4("Lambda values over time"),
      sliderInput("tick_size","Tick label size",min=2,max=30,value=12,step=1),
      fluidRow(
        column(4, sliderInput("y_range0","λ₀ Y-axis range",min=0,max=1,value=c(0,1),step=0.1)),
        column(4, sliderInput("y_range1","λ₁ Y-axis range",min=0,max=1,value=c(0,1),step=0.1)),
        column(4, sliderInput("y_range2","λ₂ Y-axis range",min=0,max=1,value=c(0,1),step=0.1))
      ),
      fluidRow(
        column(6, checkboxInput("sort1","Sort λ₁ points by value")),
        column(6, checkboxInput("sort2","Sort λ₂ points by value"))
      ),
      plotOutput("LinePlot") %>% withSpinner(color="#76ABAE"),
      tabsetPanel(
        id = "myGrid",
        tabPanel("Lambdas",
                 DT::dataTableOutput("alphas") %>% withSpinner(color="#76ABAE"),
                 downloadButton("downloadData","Download lambda values")),
        tabPanel("Transcript Weight",
                 DT::dataTableOutput("Gs") %>% withSpinner(color="#76ABAE"),
                 downloadButton("downloadData2","Download transcript weights"))
      ),
      h4("GO analysis on the first pattern"),
      plotOutput("GOPlot") %>% withSpinner(color="#76ABAE"),
      sliderInput("terms_1","Top terms:",min=10,max=30,value=15),
      h4("GO analysis on the second pattern"),
      plotOutput("GOPlot_2") %>% withSpinner(color="#76ABAE"),
      sliderInput("terms_2","Top terms:",min=10,max=30,value=15),
      h4("Stable state transcript stability"),
      plotOutput("stabPlot") %>% withSpinner(color="#76ABAE")
    ),

    tabPanel(
      "Sample Files",
      p("Please use the links below to download example files:", br(),
        tags$a(href="https://drive.google.com/file/d/1exoPw_Cnn_vNJACea68oSMTJ4Fg7DNN3/view?usp=sharing","Th0 In-vitro polarization"),
        br(),
        tags$a(href="https://drive.google.com/file/d/1y_5rGlPK52JAeiQnX6YqxNGDvpQqvAFq/view","Th1 In-vitro polarization")
      )
    ),

    tabPanel(
      "Contact Us",
      fluidRow(
        column(12, align="center", h1(icon("envelope", lib="font-awesome"))),
        column(12, align="center",
               p("Please contact annicenajafi27@gmail.com for any related issues."), br())
      ), br(), fluidRow(column(12, br()))
    )
  )
)

server <- function(input, output, session) {
  disable("submit")
  data_path <- reactiveVal(NULL)

  output$ready_name <- renderText({
    x <- data_path(); if (is.null(x)) "none" else basename(x)
  })

  auto_set_controls <- function(dt_head) {
    if (!is.null(dt_head) && nrow(dt_head) > 0) {
      detected <- detect_id_type(rownames(dt_head))
      updateSelectInput(session, "data_type", selected = detected)
      showNotification(paste("Detected identifier type:", detected), type = "message", duration = 4)

      # ---- NEW: auto-select transform
      tf <- if (is_counts_like(as.matrix(dt_head))) "log1p" else "Log with pseudocount"
      updateSelectInput(session, "transform_method", selected = "Auto")
      showNotification(paste("Auto will use:", tf), type = "message", duration = 4)
    }
  }

  observeEvent(input$data_in, {
    validate(need(!is.null(input$data_in$datapath), "No file found."))
    data_path(input$data_in$datapath)
    enable("submit")

    dt_head <- tryCatch(read.csv(input$data_in$datapath, row.names = 1,
                                 check.names = FALSE, nrows = 2000),
                        error = function(e) NULL)
    auto_set_controls(dt_head)
  }, ignoreInit = TRUE)

  observeEvent(input$load_demo, {
    demo_abs <- find_demo_file("helper_T_cell_0_test.csv")
    if (is.null(demo_abs)) {
      showModal(modalDialog(
        title = "Demo file not found",
        HTML("I looked for <code>helper_T_cell_0_test.csv</code> in:<br>
             • <code>data/</code> (next to app.R)<br>
             • <code>inst/shiny/data/</code> (package root)<br>
             • installed package via <code>system.file('shiny/data', ...)</code>"),
        easyClose = TRUE, footer = modalButton("OK")
      ))
      return()
    }
    simulate_fileinput_display("data_in", "helper_T_cell_0_test.csv", duration_ms = 700)
    data_path(demo_abs); enable("submit")

    dt_head <- tryCatch(read.csv(demo_abs, row.names = 1,
                                 check.names = FALSE, nrows = 2000),
                        error = function(e) NULL)
    auto_set_controls(dt_head)
    showModal(modalDialog("Demo dataset loaded and ready. Click Submit to analyze.",
                          easyClose = TRUE))
  }, ignoreInit = TRUE)

  read_current <- reactive({
    req(data_path())
    read.csv(data_path(), row.names = 1, check.names = FALSE)
  })

  dt_submitted <- eventReactive(input$submit, {
    read_current()
  })

  # ---- NEW: one source of truth for log-scale matrix -------------------------
  log_mat <- reactive({
    dt <- dt_submitted()
    req(nrow(dt) > 0, ncol(dt) > 0)
    method <- switch(input$transform_method,
                     "Auto" = "Auto",
                     "log1p" = "log1p",
                     "Log with pseudocount" = "pseudocount")
    apply_transform(as.matrix(dt), method = method, pseudocount = 1e-6)
  })
  # ----------------------------------------------------------------------------

  observeEvent(dt_submitted(), {
    h <- turning_alphas(log_mat())[[1]]
    for (i in 0:2) {
      rng   <- range(h[, i + 1])
      nice  <- pretty(rng, n = 3)
      diff  <- nice[length(nice)] - nice[1]
      extra <- diff * 0.1
      min_i <- nice[1] - extra
      max_i <- nice[length(nice)] + extra
      updateSliderInput(session, paste0("y_range", i),
                        min = min_i, max = max_i,
                        value = nice[c(1, length(nice))])
    }
    updateTabsetPanel(session, "sur", selected = "Results")
  }, ignoreInit = TRUE)

  output$LinePlot <- renderPlot({
    dt  <- dt_submitted()
    tp  <- seq_len(ncol(dt)); sn <- colnames(dt)
    h   <- turning_alphas(log_mat())[[1]]
    if (input$lambda1_check) h[,2] <- -h[,2]

    df1 <- data.frame(x = tp, y = h[,1])
    p1  <- ggplot(df1, aes(x, y)) +
      geom_point(colour = "#76ABAE", shape = 8, size = 1, stroke = 2) +
      geom_line(colour = "#76ABAE", size = 1) +
      scale_x_continuous(breaks = tp, labels = sn) +
      scale_y_continuous(limits = input$y_range0) +
      labs(x = "Sample", y = TeX("$\\lambda_0$"), tag = 'A') +
      theme_minimal() +
      theme(panel.grid = element_blank(), axis.line = element_line(colour = "black"),
            axis.text = element_text(size = input$tick_size),
            axis.title = element_text(size = input$tick_size + 2))

    df2 <- data.frame(sample = sn, y = h[,2])
    if (input$sort1) df2 <- df2[order(df2$y), ]
    df2$x <- seq_len(nrow(df2))
    p2    <- ggplot(df2, aes(x, y)) +
      geom_point(colour = "#5E898B", shape = 8, size = 1, stroke = 2) +
      geom_line(colour = "#5E898B", size = 1) +
      scale_x_continuous(breaks = df2$x, labels = df2$sample) +
      scale_y_continuous(limits = input$y_range1) +
      labs(x = "Sample", y = TeX("$\\lambda_1$"), tag = 'B') +
      theme_minimal() +
      theme(panel.grid = element_blank(), axis.line = element_line(colour = "black"),
            axis.text = element_text(size = input$tick_size),
            axis.title = element_text(size = input$tick_size + 2))

    df3 <- data.frame(sample = sn, y = h[,3])
    if (input$sort2) df3 <- df3[order(df3$y), ]
    df3$x <- seq_len(nrow(df3))
    p3    <- ggplot(df3, aes(x, y)) +
      geom_point(colour = "#426061", shape = 8, size = 1, stroke = 2) +
      geom_line(colour = "#426061", size = 1) +
      scale_x_continuous(breaks = df3$x, labels = df3$sample) +
      scale_y_continuous(limits = input$y_range2) +
      labs(x = "Sample", y = TeX("$\\lambda_2$"), tag = 'C') +
      theme_minimal() +
      theme(panel.grid = element_blank(), axis.line = element_line(colour = "black"),
            axis.text = element_text(size = input$tick_size),
            axis.title = element_text(size = input$tick_size + 2))

    p1 + p2 + p3
  })

  output$alphas <- DT::renderDataTable({
    dt  <- dt_submitted()
    h   <- turning_alphas(log_mat())[[1]]
    rownames(h) <- colnames(dt)
    data.frame(Sample = rownames(h), h, stringsAsFactors = FALSE)
  })

  output$Gs <- DT::renderDataTable({
    dt <- dt_submitted()
    a  <- turning_alphas(log_mat())[[2]]
    rownames(a) <- rownames(dt)
    data.frame(Gene = rownames(a), a, stringsAsFactors = FALSE)
  })

  output$GOPlot <- renderPlot({
    dt  <- dt_submitted()
    a   <- turning_alphas(log_mat())[[2]]
    if (input$lambda1_check) a[,2] <- -a[,2]
    kt  <- switch(input$data_type,
                  "Ensembl"     = "ENSEMBL",
                  "Entrez Id"   = "ENTREZID",
                  "Gene Symbol" = "SYMBOL",
                  "Probe Id"    = "PROBEID")
    ids <- toupper(rownames(a))[a[,2] > quantile(a[,2], input$percentile_GO/100)]
    entrez <- tryCatch(
      mapIds(org.Hs.eg.db, keys = ids, column = "ENTREZID", keytype = kt, multiVals = "first"),
      error = function(e) NULL
    )
    entrez <- entrez[!is.na(entrez)]
    if (length(entrez) > 0) {
      GOres <- enrichGO(gene = entrez, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "BP")
      if (nrow(GOres@result) > 0) {
        df_top <- head(GOres@result, input$terms_1)
        ggplot(df_top, aes(x = Description, y = Count, fill = p.adjust)) +
          geom_bar(stat = "identity") + coord_flip() +
          scale_fill_gradient(low = "#790915", high = "#062c5c") +
          theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.title.y = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 20),
                text = element_text(size = 18))
      }
    }
  })

  output$GOPlot_2 <- renderPlot({
    dt  <- dt_submitted()
    a   <- turning_alphas(log_mat())[[2]]
    if (input$lambda1_check) a[,2] <- -a[,2]
    kt  <- switch(input$data_type,
                  "Ensembl"     = "ENSEMBL",
                  "Entrez Id"   = "ENTREZID",
                  "Gene Symbol" = "SYMBOL",
                  "Probe Id"    = "PROBEID")
    ids <- toupper(rownames(a))[a[,3] > quantile(a[,3], input$percentile_GO/100)]
    entrez <- tryCatch(
      mapIds(org.Hs.eg.db, keys = ids, column = "ENTREZID", keytype = kt, multiVals = "first"),
      error = function(e) NULL
    )
    entrez <- entrez[!is.na(entrez)]
    if (length(entrez) > 0) {
      GOres <- enrichGO(gene = entrez, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "BP")
      if (nrow(GOres@result) > 0) {
        df_top <- head(GOres@result, input$terms_2)
        ggplot(df_top, aes(x = Description, y = Count, fill = p.adjust)) +
          geom_bar(stat = "identity") + coord_flip() +
          scale_fill_gradient(low = "#790915", high = "#062c5c") +
          theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.title.y = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 20),
                text = element_text(size = 18))
      }
    }
  })

  output$stabPlot <- renderPlot({
    res <- turning_alphas(log_mat())
    h   <- res[[1]]; a <- res[[2]]
    if (input$lambda1_check) h[,2] <- -h[,2]
    df  <- data.frame(x = a[,1]*h[,1], y = a[,2]*h[,2])

    ggplot(df, aes(x, y)) + geom_point(fill = "#76ABAE", color = "#76ABAE", shape = 8) +
      xlab(TeX("$G_{0i}\\lambda_0$")) +
      ylab(TeX("$G_{1i}\\lambda_1$")) +
      theme(panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5, size = 20),
            text = element_text(size = 18)) +
      labs(tag = "F", title = TeX("Stable state stability"))
  })

  output$downloadData <- downloadHandler(
    filename = function() "lambdas.csv",
    content  = function(file) {
      write.csv(turning_alphas(log_mat())[[1]], file)
    }
  )
  output$downloadData2 <- downloadHandler(
    filename = function() "weights.csv",
    content  = function(file) {
      write.csv(turning_alphas(log_mat())[[2]], file)
    }
  )
}

shinyApp(ui, server)
