options(shiny.maxRequestSize = 2000*1024^2)

turning_alphas <- function(log.mat) {
  mat <- as.matrix(log.mat)
  mat[mat == 0] <- 1e-6
  C <- crossprod(mat)
  ed <- eigen(C)
  P <- ed$vectors
  D <- diag(sqrt(ed$values))
  Y <- mat
  n <- ncol(Y)
  alph_lst <- lapply(seq_len(n), function(i) colSums(solve(D))[i] * Y %*% P[,i])
  alph_all <- do.call(cbind, alph_lst)
  hold_lst <- lapply(seq_len(n), function(i) P[,i] * D[i,i])
  holds <- do.call(cbind, hold_lst)
  list(holds, alph_all)
}

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

    tabPanel("Instructions",

             sidebarLayout(
               sidebarPanel(
                 fileInput("data_in","Please upload a gene expression matrix csv file."),
                 selectInput("data_type","Select biomarker input type.",
                             choices=c("Ensembl","Entrez Id","Gene Symbol","Probe Id")),
                 checkboxInput("lambda1_check","Flip lambda 1 value"),
                 sliderInput("percentile_GO","GO transcript percentile",
                             min=80,max=95,value=95,step=1,ticks=FALSE),
                 actionButton("submit",strong(em("Submit")),
                              style="color:#000;background-color:#76ABAE;border-color:#76ABAE")
               ),

               mainPanel(

                 div(style="border:1px solid #78A2A5;background:#E4EEEF;padding:20px;text-align:justify;",
                     p(strong("Welcome!")),
                     p("I. Please upload a csv file using the slot on the left."),
                     p("II. The first column of the csv will be used as gene names and all other columns as samples, in the order provided."),
                     p("III. An example of the input data is provided below and under the 'Sample Files' tab."),
                     tags$img(src='sample_file.png',height=200,width=500),
                     p("IV. Please ensure that you select the correct identifier type according to your transcript column."),
                     p(strong("Surprisal Analysis Applications:")),
                     tags$a(href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0108549",
                            "Kravchenko-Balasha N et al. PLoS ONE. 2014;9(11):e108549."),
                     br(),br(),
                     tags$a(href="https://www.pnas.org/doi/full/10.1073/pnas.1414714111",
                            "Zadran S et al. PNAS. 2014;111(36):13235-40."),
                     br(),br(),
                     tags$a(href="https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007034",
                            "Su Y et al. PLoS Comput Biol. 2019;15(6):e1007034."),
                     br(),br(),
                     tags$a(href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5903653/",
                            "Bogaert KA et al. PLoS ONE. 2018;13(4):e0195142."),
                     br(),br()
                 )
               )
             )
    ),

    tabPanel("Results",

             h4("Lambda values over time"),
             sliderInput("tick_size","Tick label size",min=6,max=30,value=12,step=1),
             fluidRow(
               column(4, sliderInput("y_range0","λ₀ Y-axis range",min=0,max=1,value=c(0,1),step=0.1)),
               column(4, sliderInput("y_range1","λ₁ Y-axis range",min=0,max=1,value=c(0,1),step=0.1)),
               column(4, sliderInput("y_range2","λ₂ Y-axis range",min=0,max=1,value=c(0,1),step=0.1))
             ),
             fluidRow(
               column(6, checkboxInput("sort1","Sort λ₁ points by value")),
               column(6, checkboxInput("sort2","Sort λ₂ points by value"))
             ),
             plotOutput("LinePlot")%>%withSpinner(color="#76ABAE"),
             tabsetPanel(id="myGrid",
                         tabPanel("Lambdas",
                                  DT::dataTableOutput("alphas")%>%withSpinner(color="#76ABAE"),
                                  downloadButton("downloadData","Download lambda values")
                         ),
                         tabPanel("Transcript Weight",
                                  DT::dataTableOutput("Gs")%>%withSpinner(color="#76ABAE"),
                                  downloadButton("downloadData2","Download transcript weights")
                         )
             ),
             h4("GO analysis on the first pattern"),
             plotOutput("GOPlot")%>%withSpinner(color="#76ABAE"),
             sliderInput("terms_1","Top terms:",min=10,max=30,value=15),
             h4("GO analysis on the second pattern"),
             plotOutput("GOPlot_2")%>%withSpinner(color="#76ABAE"),
             sliderInput("terms_2","Top terms:",min=10,max=30,value=15),
             h4("Stable state transcript stability"),
             plotOutput("stabPlot")%>%withSpinner(color="#76ABAE")

    ),
    tabPanel("Sample Files",

             p("Please use the links below to download example files:",br(),
               tags$a(href="https://drive.google.com/file/d/1exoPw_Cnn_vNJACea68oSMTJ4Fg7DNN3/view?usp=sharing","Th0 In-vitro polarization"),
               br(),
               tags$a(href="https://drive.google.com/file/d/1y_5rGlPK52JAeiQnX6YqxNGDvpQqvAFq/view","Th1 In-vitro polarization")
             )
    ),
    tabPanel("Contact Us",
             fluidRow(
               column(12,align="center",h1(icon("envelope",lib="font-awesome"))),
               column(12,align="center",
                      p("Please contact annicenajafi27@gmail.com for any related issues."),br()
               )
             ),br(),fluidRow(column(12,br()))

    )
  )
)

server <- function(input, output, session) {

  disable("send")

  observeEvent(input$data_in,{
    dt <- read.csv(input$data_in$datapath,row.names=1,check.names=FALSE)
    mat <- as.matrix(dt); mat[mat==0] <- 1e-6
    h   <- turning_alphas(log(mat))[[1]]
    for(i in 0:2) {
      rng    <- range(h[,i+1])
      nice   <- pretty(rng, n=3)
      diff   <- nice[length(nice)] - nice[1]
      extra  <- diff * 0.1
      min_i  <- nice[1] - extra
      max_i  <- nice[length(nice)] + extra
      updateSliderInput(session,paste0("y_range",i),
                        min   = min_i,
                        max   = max_i,
                        value = nice[c(1,length(nice))]
      )
    }
  })


  output$LinePlot <- renderPlot({
    req(input$data_in)
    dt  <- read.csv(input$data_in$datapath,row.names=1,check.names=FALSE)
    tp  <- seq_len(ncol(dt)); sn <- colnames(dt)
    mat <- as.matrix(dt); mat[mat==0] <- 1e-6; lm <- log(mat)
    h   <- turning_alphas(lm)[[1]]
    if(input$lambda1_check) h[,2] <- -h[,2]

    df1 <- data.frame(x=tp,y=h[,1])
    p1  <- ggplot(df1,aes(x,y))+
      geom_point(colour="#76ABAE", shape=8, size=1, stroke = 2)+
      geom_line(colour="#76ABAE", size=1)+
      scale_x_continuous(breaks=tp,labels=sn)+
      scale_y_continuous(limits=input$y_range0)+
      labs(x="Sample",y=TeX("$\\lambda_0$"), tag='A')+
      theme_minimal()+
      theme(panel.grid=element_blank(),axis.line=element_line(colour="black"),
            axis.text=element_text(size=input$tick_size),
            axis.title=element_text(size=input$tick_size+2))

    df2 <- data.frame(sample=sn,y=h[,2])
    if(input$sort1) df2 <- df2[order(df2$y),]
    df2$x <- seq_len(nrow(df2))
    p2    <- ggplot(df2,aes(x,y))+
      geom_point(colour="#5E898B", shape=8, size=1, stroke = 2)+
      geom_line(colour="#5E898B", size=1)+
      scale_x_continuous(breaks=df2$x,labels=df2$sample)+
      scale_y_continuous(limits=input$y_range1)+
      labs(x="Sample",y=TeX("$\\lambda_1$"), tag='B')+
      theme_minimal()+
      theme(panel.grid=element_blank(),axis.line=element_line(colour="black"),
            axis.text=element_text(size=input$tick_size),
            axis.title=element_text(size=input$tick_size+2))

    df3 <- data.frame(sample=sn,y=h[,3])

    if(input$sort2) df3 <- df3[order(df3$y),]

    df3$x <- seq_len(nrow(df3))

    p3    <- ggplot(df3,aes(x,y))+
      geom_point(colour="#426061", shape=8, size=1, stroke = 2)+
      geom_line(colour="#426061", size=1)+
      scale_x_continuous(breaks=df3$x,labels=df3$sample)+
      scale_y_continuous(limits=input$y_range2)+
      labs(x="Sample",y=TeX("$\\lambda_2$"), tag='C')+
      theme_minimal()+
      theme(panel.grid=element_blank(),axis.line=element_line(colour="black"),
            axis.text=element_text(size=input$tick_size),
            axis.title=element_text(size=input$tick_size+2))

    p1 + p2 + p3
  })

  output$alphas <- DT::renderDataTable({

    req(input$data_in)
    dt  <- read.csv(input$data_in$datapath,row.names=1,check.names=FALSE)
    mat <- as.matrix(dt); mat[mat==0] <- 1e-6; lm <- log(mat)
    h   <- turning_alphas(lm)[[1]]
    rownames(h) <- colnames(dt)
    data.frame(Sample=rownames(h),h,stringsAsFactors=FALSE)

  })

  output$Gs <- DT::renderDataTable({

    req(input$data_in)
    dt <- read.csv(input$data_in$datapath,row.names=1,check.names=FALSE)
    a  <- turning_alphas(log(as.matrix(dt)))[[2]]
    rownames(a) <- rownames(dt)
    data.frame(Gene=rownames(a),a,stringsAsFactors=FALSE)

  })

  output$GOPlot <- renderPlot({

    req(input$data_in)
    dt  <- read.csv(input$data_in$datapath,row.names=1,check.names=FALSE)
    a   <- turning_alphas(log(as.matrix(dt)))[[2]]
    if(input$lambda1_check) a[,2] <- -a[,2]
    kt  <- switch(input$data_type,
                  "Ensembl"     = "ENSEMBL",
                  "Entrez Id"   = "ENTREZID",
                  "Gene Symbol" = "SYMBOL",
                  "Probe Id"    = "PROBEID")
    ids <- toupper(rownames(a))[a[,2] > quantile(a[,2],input$percentile_GO/100)]
    entrez <- tryCatch(mapIds(org.Hs.eg.db,keys=ids,column="ENTREZID",keytype=kt,multiVals="first"),error=function(e)NULL)
    entrez <- entrez[!is.na(entrez)]

    if(length(entrez)>0) {

      GOres <- enrichGO(gene=entrez,OrgDb="org.Hs.eg.db",keyType="ENTREZID",ont="BP")

      if(nrow(GOres@result)>0) {

        df_top <- head(GOres@result,input$terms_1)

        ggplot(df_top,aes(x=Description,y=Count,fill=p.adjust))+
          geom_bar(stat="identity")+coord_flip()+
          scale_fill_gradient(low = "#790915", high = "#062c5c")+
          theme(
            # Remove panel border
            panel.border=element_blank(),
            #plot.border = element_blank(),
            # Remove panel grid lines
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # Add axis line
            axis.line = element_line(colour = "black"),
            #axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            #axis.text = element_blank(),
            #legend.position = "none",
            plot.title = element_text(hjust = 0.5, size=20),
            #axis.text = element_text(size = 15),

            text = element_text(size=18)
          )

      }

    }

  })

  output$GOPlot_2 <- renderPlot({

    req(input$data_in)
    dt  <- read.csv(input$data_in$datapath,row.names=1,check.names=FALSE)
    a   <- turning_alphas(log(as.matrix(dt)))[[2]]

    if(input$lambda1_check) a[,2] <- -a[,2]

    kt  <- switch(input$data_type,
                  "Ensembl"     = "ENSEMBL",
                  "Entrez Id"   = "ENTREZID",
                  "Gene Symbol" = "SYMBOL",
                  "Probe Id"    = "PROBEID")
    ids <- toupper(rownames(a))[a[,3] > quantile(a[,3],input$percentile_GO/100)]
    entrez <- tryCatch(mapIds(org.Hs.eg.db,keys=ids,column="ENTREZID",keytype=kt,multiVals="first"),error=function(e)NULL)
    entrez <- entrez[!is.na(entrez)]

    if(length(entrez)>0) {

      GOres <- enrichGO(gene=entrez,OrgDb="org.Hs.eg.db",keyType="ENTREZID",ont="BP")

      if(nrow(GOres@result)>0) {
        df_top <- head(GOres@result,input$terms_2)
        ggplot(df_top,aes(x=Description,y=Count,fill=p.adjust))+
          scale_fill_gradient(low = "#790915", high = "#062c5c")+
          geom_bar(stat="identity")+coord_flip()+
          theme(
            # Remove panel border
            panel.border=element_blank(),
            #plot.border = element_blank(),
            # Remove panel grid lines
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # Add axis line
            axis.line = element_line(colour = "black"),
            #axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            #axis.text = element_blank(),
            #legend.position = "none",
            plot.title = element_text(hjust = 0.5, size=20),
            #axis.text = element_text(size = 15),

            text = element_text(size=18)
          )

      }

    }

  })

  output$stabPlot <- renderPlot({
    req(input$data_in)

    dt  <- read.csv(input$data_in$datapath,row.names=1,check.names=FALSE)
    res <- turning_alphas(log(as.matrix(dt)))
    h   <- res[[1]]; a <- res[[2]]

    if(input$lambda1_check) h[,2] <- -h[,2]

    df  <- data.frame(x=a[,1]*h[,1],y=a[,2]*h[,2])

    ggplot(df,aes(x,y))+geom_point(fill="#76ABAE", color="#76ABAE", shape=8)+
      xlab(TeX("$G_{0i}\\lambda_0$"))+
      ylab(TeX("$G_{1i}\\lambda_1$"))+
      theme(
        # Remove panel border
        panel.border=element_blank(),
        #plot.border = element_blank(),
        # Remove panel grid lines
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Add axis line
        axis.line = element_line(colour = "black"),
        #axis.title.x = element_blank(),
        #axis.text = element_blank(),
        #legend.position = "none",
        plot.title = element_text(hjust = 0.5, size=20),
        #axis.text = element_text(size = 15),

        text = element_text(size=18)
      ) +labs(tag="F", title=TeX("Stable state stability"))

  })

  observeEvent(input$submit,{
    updateTabsetPanel(session,"sur",selected="Results")
  })

  output$downloadData <- downloadHandler(
    filename = function() "lambdas.csv",
    content  = function(file) {
      dt <- read.csv(input$data_in$datapath,row.names=1,check.names=FALSE)
      write.csv(turning_alphas(log(as.matrix(dt)))[[1]],file)
    }
  )

  output$downloadData2 <- downloadHandler(
    filename = function() "weights.csv",
    content  = function(file) {
      dt <- read.csv(input$data_in$datapath,row.names=1,check.names=FALSE)
      write.csv(turning_alphas(log(as.matrix(dt)))[[2]],file)
    }
  )
}

shinyApp(ui, server)
