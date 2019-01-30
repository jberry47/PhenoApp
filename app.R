library(reshape2)
library(grid)
library(stringr)
library(scales)
library(gridExtra)
library(shiny)
library(shinyBS)
library(shinydashboard)
library(DT)
library(shinyjs)
library(plyr)
library(ggplot2)
library(lme4)
library(RSQLite)
library(FactoMineR)
library(factoextra)
library(ggridges)
library(car)
library(survival)
library(vegan)
library(shinycssloaders)


ui <- dashboardPage(skin="black", title="Phenotyper Analysis Tool",
                    dashboardHeader(
                      title = tagList(
                        tags$span(
                          class = "logo-mini", "PAT"
                        ),
                        tags$span(
                          class = "logo-lg", "Phenotyper Analysis Tool"
                        )
                      ),
                      titleWidth = 450
                    ),
                    dashboardSidebar(
                      tags$script(HTML("$('body').addClass('sidebar-mini');")),
                      width = 150,
                      sidebarMenu(
                        menuItem("Overview", tabName = "overview"),
                        menuItem("Get Started",tabName = "get_started")
                      )
                    ),
                    dashboardBody(
                      useShinyjs(),
                      fluidRow(
                        tags$head(tags$style("#container * {display: inline;}")),
                        tags$style(HTML("
                                        .shiny-progress-notification .progress-text {
                                        font-size: 17pt;
                                        }
                                        .tabbable > .nav > li[class=active]    > a {background-color: #444444; color:white}
                                        .multicol{
                                        -webkit-column-count: 4; /* Chrome, Safari, Opera */
                                        -moz-column-count: 4; /* Firefox */
                                        column-count: 4;
                                        }
                                        .twocol{
                                        -webkit-column-count: 2; /* Chrome, Safari, Opera */
                                        -moz-column-count: 2; /* Firefox */
                                        column-count: 2;
                                        }
                                        .warning { 
                                        color: red;
                                        }"
                      )),
                      tabItems(
                        tabItem(tabName = "overview",
                                box(width=10,title = "Welcome",solidHeader = T,status = 'success',collapsible = TRUE,
                                    p("This analysis tool is designed to take input from either PlantCV or PhenotyperCV that was used to process images from 
                                      the Bellweather Phenotyping Facility at Donald Danforth Plant Science Center. Depending on what program was used to 
                                      analyze the images, there is a different import process. In both cases, after importing, four more boxes will appear: Outlier 
                                      Detection and Removal, Shapes Analysis, VIS Analysis, and NIR Analysis. These analyses are intended to be a first pass look 
                                      at data from this facility and does not perform any statistical inferences. If an effect is observed, it is
                                      required that proper statistical is testing is done outside of this framework.")
                                    ),
                                box(width=10,title = "Design File",solidHeader = T,status = 'success',collapsible = TRUE,
                                    p("In both cases, importing PhenotyperCV or PlantCV output, a common design file is needed. This file is the link that takes the barcode of a plant and assigns the
                                      the design parameters. So at most there should be 1141 rows in this file. See last paragraph in this box for exceptions. There are a couple of rules that 
                                      need to be followed for this tool to work properly."),
                                    tags$ul(
                                      tags$li("There must be a column called ",tags$b("Barcodes"),"with exactly this
                                              spelling."),
                                      tags$li("The other columns in this file shouldn't have any entries that have a ",code("."),"in them. For example sorghum has a genotype 
                                              commonly called B.Az9504 and if this is how it's labeled in this file, things will turn out funny. Please edit entries such 
                                              as this."),
                                      tags$li("While this tool is setup to handle more than two design columns, many of the plots are
                                              coded to only accept two columns (such as Genotype and Treatment). If your design is more complicated than that (like the second example design file)
                                              , text will appear indicating that it's averaging over specific design parameters.")
                                      ),
                                    p("The first row of the design file should be: Barcodes, Var1, Var2, ..., VarM. Where Var1, Var2, ..., VarM could be called anything. The first 10 rows of two example design files are
                                      shown here. "),
                                    column(6,
                                           tableOutput("design_ex1")),
                                    column(6,
                                           tableOutput("design_ex2")),
                                    p("If you only want to analyze a subset of your data, outside of this framework, select only those things you which to look 
                                      at in the design file and upload that. The merging process will remove plants from the dataset where there is no matching barcode 
                                      in the design file. ")
                                    ),
                                box(width=5,title = "Importing from PhenotyperCV",solidHeader = T,status = 'success',collapsible = TRUE,
                                    p("When processing images using this program, there are 5 files, 4 of which are required: design, snapshot, shapes, color, and nir-color (optional). 
                                      The design file is shown above and the snapshot file is the one that came with the image download (SnapshotInfo.csv). For the 
                                      shapes, color and nir-color files, these are created when running PhenotyperCV and no alterations are required to load the
                                      data into this framework. For example: "),
                                    code("find Images/ -name 'VIS_SV*' | xargs -P8 -I{} ./PhenotyperCV VIS_CH {} vis_background_image.png shapes.txt color.txt"),
                                    p("created shapes.txt and color.txt and"),
                                    code("find Images/ -name 'NIR_SV*' | xargs -P8 -I{} ./PhenotyperCV NIR {} nir_background_image.png nir_color.txt"),
                                    p("created nir_color.txt")
                                    ),
                                box(width=5,title = "Importing from PlantCV",solidHeader = T,status = 'success',collapsible = TRUE,
                                    p("When processing images using this program, there are 3 files that are required: design, snapshot, and sqlite3. The design file is shown 
                                      above,the snapshot file is the one that came with the image download (SnapshotInfo.csv), and the sqlite3 file is created when running PlantCV with the",code("-s"), "flag. An example bash script to run image
                                      analysis using this program and getting a sqlite3 file out is shown here:"),
                                    code("#!/bin/bash",br(),"
                                         /home/jberry/plantcv/plantcv-pipeline.py \\",br(),"
                                         -d /home/jberry/Phenotyper/Exp2 \\",br(),"
                                         -a phenofront \\",br(),"
                                         -p /home/jberry/Phenotyper/vis_sv_z1_L0.py \\",br(),"
                                         -i /home/jberry/Phenotyper/plantcv_SV_output \\",br(),"
                                         -T 40 \\",br(),"
                                         -c \\",br(),"
                                         -s PlantCV_SV_drought.sqlite3 \\",br(),"
                                         -f imgtype_camera_frame_zoom_lifter_gain_exposure_other_id\\",br(),"
                                         -M imgtype:VIS,camera:SV \\",br(),"
                                         -t png \\",br(),"
                                         -C NIR"),
                                    p("which creates the PlantCV_SV_drought.sqlite3 file."),
                                    p("PlantCV allows for co-processing the NIR images along side of the VIS images and that is done with the",code("-C"),"flag. This import
                                      method is designed to handle with and without the NIR data."),
                                    div(id="container",p("For more information about PlantCV, click "),tags$a(href="https://plantcv.readthedocs.io/en/latest/","here",target="_blank"))
                                  ),
                                  box(width=10,title = "Contributing and error reporting",solidHeader = T,status = 'success',collapsible = TRUE,collapsed=F,
                                    p("If you'd like to contribute to this app, you can! This application is on Github at",tags$a(href="https://github.com/jberry47/Shiny-PhenoAnalyzer","github.com/jberry47/Shiny-PhenoAnalyzer",target="_blank"),
                                      "and you just have to follow a few simple steps:"),
                                    tags$ol(
                                      tags$li("Fork the repository to your personal repository"),
                                      tags$li("Cloned the forked repository to your local machine"),
                                      tags$li("Create a local branch to do your work"),
                                      tags$li("Make the edits/additions/deletions you'd like to make"),
                                      tags$li("Push all the changes to your fork"),
                                      tags$li("Create a pull request from your fork to the master branch of this app")
                                    ),
                                    p("After the pull request is made, admins of the page will review it and accept it if there aren't conflicts."),
                                    p("If you're using the app and something breaks or it doesn't work as expected, then please make a 'New Issue' in the github
                                      repository and we'll look into it as soon as we can.")
                                  )
                        ),
                        tabItem(tabName = "get_started",
                                box(width=10,title = "Merging Files",solidHeader = T,status = 'success',collapsible = TRUE,
                                    h5("How many days before the first day of imaging were the plants planted?"),
                                    numericInput("dap_offset", "DAP Offset", value=2,width=80),
                                    br(),
                                    tabsetPanel(
                                      tabPanel(title="PhenotyperCV",
                                               fileInput("phenocv_design_file", "Choose design file",
                                                         multiple = F,
                                                         accept = c(".csv")),
                                               fileInput("phenocv_snapshot_file", "Choose snapshot file",
                                                         multiple = F,
                                                         accept = c(".csv")),
                                               fileInput("phenocv_shapes_file", "Choose shapes file",
                                                         multiple = F,
                                                         accept = c(".txt")),
                                               fileInput("phenocv_color_file", "Choose color file",
                                                         multiple = F,
                                                         accept = c(".txt")),
                                               radioButtons("pheno_nir_q", "Analyze NIR data",choices = c("Yes"="Yes","No"="No")),
                                               uiOutput("phenocv_nir_q_ui"),
                                               uiOutput("phenocv_go_ui"),
                                               br(),
                                               uiOutput("phenocv_download_merged_button")
                                      ),
                                      tabPanel(title = "PlantCV",
                                               fileInput("plantcv_design_file", "Choose design file",
                                                         multiple = F,
                                                         accept = c(".csv")),
                                               fileInput("plantcv_snapshot_file", "Choose snapshot file",
                                                         multiple = F,
                                                         accept = c(".csv")),
                                               fileInput("plantcv_sql_path", "Choose sqlite3 database",
                                                         multiple = F,
                                                         accept = c(".sqlite3")),
                                               uiOutput("plantcv_go_ui")
                                      )
                                    )
                                ),
                                uiOutput("summary_ui"),
                                uiOutput("outlier_removal"),
                                uiOutput("shapes_ui"),
                                uiOutput("vis_ui"),
                                uiOutput("nir_ui")
                                )
                        )
                      )
                      )
                    )

server <- function(input, output){
  report_error <- function(err){
    locs <- extractStackTrace(conditionStackTrace(err))$loc
    loc_lines <- max(as.numeric(str_sub(na.omit(unlist(lapply(strsplit(locs,"#"),function(i) i[2]))),end = -2)))
    showModal(modalDialog(
      p("Error in code block: ",tags$b(imp_error_step$data)),
      p("Error caught: ",code(err)),
      p("Line number ",as.numeric(loc_lines),": ",code(scan("app.R", '', skip = as.numeric(loc_lines)-1, nlines = 1, sep = '\n')))
    ))
  }
  
  output$phenocv_nir_q_ui <- renderUI({
    if(input$pheno_nir_q == "Yes"){
      fileInput("phenocv_nir_file", "Choose nir file",
                multiple = F,
                accept = c(".txt"))
    }
  })
  
  output$design_ex1 <- renderTable({
    head(read.csv("data/design_ex1.csv",stringsAsFactors = F),n=10)
  },bordered=TRUE,spacing = "xs",align = "l")
  output$design_ex2 <- renderTable({
    head(read.csv("data/design_ex2.csv",stringsAsFactors = F),n=10)
  },bordered=TRUE,spacing = "xs",align = "l")
  
  
  get_color <- function(file_name,snapshot1,design1,start,stop){
    color_data <- read.table(file_name,header = F,stringsAsFactors = F,sep = " ")
    color_data <- color_data[,-ncol(color_data)]
    color_data$id <- unlist(lapply(strsplit(color_data$V1,"/"),function(i) strsplit(i[str_detect(i,"snapshot")],"snapshot")[[1]][2]))
    color_data$imgname <- unlist(lapply(strsplit(color_data$V1,"/"),function(i) strsplit(i[str_detect(i,"png")],"[.]")[[1]][1]))    
    color_data <- join(color_data,snapshot1[,c("id","Barcodes","timestamp")],by="id")
    color_data <- join(color_data,design1,by="Barcodes")
    color_data$timestamp <- strptime(color_data$timestamp,format = "%Y-%m-%d %H:%M:%S")+as.numeric(input$dap_offset)
    beg <- min(color_data$timestamp)
    color_data$DAP <- floor(as.numeric((color_data$timestamp - beg)/60/60/24))+as.numeric(input$dap_offset)
    color_data[,start:stop] <- t(apply(color_data[,start:stop],1,function(i){i/(sum(i,na.rm = T)+1)}))*100
    color_data$hr <- as.POSIXlt(color_data$timestamp)$hour
    return(color_data)
  }
  
  options(shiny.maxRequestSize=2000*1024^2)
  
  #***********************************************************************************************
  # Merging Files Box
  #***********************************************************************************************
  output$phenocv_go_ui <- renderUI({
    b <- c(input$phenocv_design_file$name,input$phenocv_snapshot_file$name,input$phenocv_shapes_file$name,input$phenocv_color_file$name,input$phenocv_nir_file$name)
    #b <- c(input$phenocv_design_file$name,input$phenocv_snapshot_file$name,input$phenocv_shapes_file$name,input$phenocv_color_file$name)
    if((input$pheno_nir_q == "Yes" & length(b) == 5)|(input$pheno_nir_q == "No" & length(b) == 4)){
      actionButton("phenocv_merge","Merge Data")
    }
  })
  
  output$plantcv_go_ui <- renderUI({
    b <- c(input$plantcv_sql_path$name,input$plantcv_design_file$name,input$plantcv_snapshot_file$name)
    if(length(b) == 3){
      actionButton("plantcv_merge","Merge Data")
    }
  })
  
  output$phenocv_download_merged_button <- renderUI({
    if(!is.null(merged$data)){
      downloadButton("phenocv_merged_table","Download Merged Data (tsv)")
    }
  })
  
  output$phenocv_merged_table <- downloadHandler(
    filename = function() {"phenocv_merged_data.tsv"},
    content = function(file){
      write.table(merged$data,file,row.names = FALSE, quote = FALSE,sep = "\t")
    }
  )
  
  #***********************************************************************************************
  # Merging Files Action
  #***********************************************************************************************
  merged <- reactiveValues(data=NULL)
  design <- reactiveValues(data=NULL)
  shapes <- reactiveValues(data=NULL)
  vis <- reactiveValues(data=NULL)
  nir <- reactiveValues(data=NULL)
  empties1 <- reactiveValues(data=NULL)
  from <- reactiveValues(data=NULL)
  snapshot <- reactiveValues(data=NULL)
  imp_error_step <- reactiveValues(data=NULL)

  observeEvent(input$phenocv_merge,{
    merged$data <- NULL; design$data <- NULL; shapes$data <- NULL; vis$data <- NULL; nir$data <- NULL; empties1$data <- NULL; from$data <- NULL; snapshot$data <- NULL; nir_ready_checker$data <- FALSE; vis_ready_checker$data <- FALSE; nir_caps$data <- NULL; vis_caps$data <- NULL; outlier_check$data <- FALSE; cooksd$data <- NULL; outlier_fmla$data <- NULL; imp_error_step$data <- NULL
    from$data <- "phenocv"
    res <- try(withCallingHandlers(withLogErrors({
      if(input$pheno_nir_q == "Yes"){
        n <- 12
      }else{
        n <- 9
      }
      
      withProgress(message = '', value = 0, {
        imp_error_step$data <- "PhenoCV - Reading design file"
        incProgress(1/n, detail = "Reading design file...")
        assoc <- read.csv(input$phenocv_design_file$datapath,header=T,stringsAsFactors = F)
        assoc_empty <- assoc[rowSums(sapply(colnames(assoc),function(i) !(assoc[,i] %in% c("Blank","Empty","blank","empty"))))==ncol(assoc),]
        design$data <- assoc

        imp_error_step$data <- "PhenoCV - Reading snapshot file"      
        incProgress(1/n, detail = "Reading snapshot file...")
        img_to_barcode <- read.csv(input$phenocv_snapshot_file$datapath,header = T,stringsAsFactors = F)
        img_to_barcode$timestamp <- as.POSIXct(strptime(img_to_barcode$timestamp,format = "%Y-%m-%d %H:%M:%S"))
        colnames(img_to_barcode)[3] <- "Barcodes"
        snapshot1 <- join(assoc_empty,img_to_barcode, by = "Barcodes")
        snapshot$data <- snapshot1[-snapshot1$weight.before < 0,]
        img_to_barcode <- img_to_barcode[img_to_barcode$tiles != "",]
        img_to_barcode <- img_to_barcode[,c("id","Barcodes","timestamp")]

        imp_error_step$data <- "PhenoCV - Reading shapes file"
        incProgress(1/n, detail = "Reading shapes file...")
        sv_shapes <- read.table(input$phenocv_shapes_file$datapath,header = F,stringsAsFactors = F,sep = " ")
        sv_shapes <- sv_shapes[,which(!as.logical(apply(sv_shapes,2,FUN=function(i) all(is.na(i)))))]
        
        if(ncol(sv_shapes)==21){
          colnames(sv_shapes) <- c("meta","area","hull_area","solidity","perimeter","width","height","cmx","cmy","hull_verticies","ex","ey","emajor","eminor","angle","eccen","circ","round","ar","fd","oof")
          shapes$data <- sv_shapes
        }else{
          colnames(sv_shapes) <- c("meta","area","hull_area","solidity","perimeter","width","height","cmx","cmy","hull_verticies","ex","ey","emajor","eminor","angle","eccen","circ","round","ar","fd","oof", "det")
          shapes$data <- sv_shapes
        }
        
        sv_shapes$id <- unlist(lapply(strsplit(sv_shapes$meta,"/"),function(i) strsplit(i[str_detect(i,"snapshot")],"snapshot")[[1]][2]))
        sv_shapes$imgname <- unlist(lapply(strsplit(sv_shapes$meta,"/"),function(i) strsplit(i[str_detect(i,"png")],"[.]")[[1]][1]))
        removeNotification(id)
        
        imp_error_step$data <- "PhenoCV - Joining shapes and snapshot"
        incProgress(1/n, detail = "Joining shapes and snapshot...")
        sv_shapes <- join(sv_shapes,img_to_barcode[,c("id","Barcodes","timestamp")],by="id")
        removeNotification(id)
        
        imp_error_step$data <- "PhenoCV - Joining shapes and design"
        incProgress(1/n, detail = "Joining shapes and design...")
        sv_shapes <- join(sv_shapes,assoc_empty,by="Barcodes")
        sv_shapes <- sv_shapes[rowSums(sapply(colnames(assoc),function(i) !is.na(sv_shapes[,i])))==ncol(assoc),]
        removeNotification(id)
        
        imp_error_step$data <- "PhenoCV - Adding time columns"
        incProgress(1/n, detail = "Adding time columns...")
        sv_shapes$timestamp <- strptime(sv_shapes$timestamp,format = "%Y-%m-%d %H:%M:%S")
        beg <- min(sv_shapes$timestamp)
        sv_shapes$DAP <- floor(as.numeric((sv_shapes$timestamp - beg)/60/60/24))+as.numeric(input$dap_offset)
        sv_shapes$hour <- lubridate::hour(sv_shapes$timestamp)
        removeNotification(id)
        
        imp_error_step$data <- "PhenoCV - Removing empty pots"
        incProgress(1/n, detail = "Removing empty pots...")
        empties <-names(which(lapply(split(sv_shapes,sv_shapes$Barcodes),function(i) all(i$area < 10))==T))
        empties1$data <- data.frame("Barcodes" = empties, stringsAsFactors = F)
        sv_shapes <- sv_shapes[!(sv_shapes$Barcodes %in% empties),]
        sv_shapes[which(sv_shapes == Inf,arr.ind = T)] <- NaN    
        removeNotification(id)
        
        merged$data <- sv_shapes
        
        imp_error_step$data <- "PhenoCV - Reading VIS color data"
        incProgress(1/n, detail = "Reading VIS color data...")
        vis$data <- get_color(input$phenocv_color_file$datapath,img_to_barcode,assoc_empty,2,181)
        removeNotification(id)
        incProgress(1/n, detail = "Removing empty pots...")
        vis$data <- vis$data[!(vis$data$Barcodes %in% empties),]
        vis$data <- vis$data[rowSums(sapply(colnames(assoc),function(i) !is.na(vis$data[,i])))==ncol(assoc),]
        removeNotification(id)
        
        if(input$pheno_nir_q == "Yes"){
          imp_error_step$data <- "PhenoCV - Reading NIR color data"
          incProgress(1/n, detail = "Reading NIR color data...")
          nir$data <- get_color(input$phenocv_nir_file$datapath,img_to_barcode,assoc_empty,2,256)
          removeNotification(id)
          incProgress(1/n, detail = "Removing empty pots...")
          nir$data <- nir$data[!(nir$data$Barcodes %in% empties),]
          nir$data <- nir$data[rowSums(sapply(colnames(assoc),function(i) !is.na(nir$data[,i])))==ncol(assoc),]
          removeNotification(id)
          incProgress(1/n, detail = "Calculating NIR average...")
          nir$data$intensityAVG <- apply(nir$data[,2:256],1,function(i){sum((i/100)*(2:256),na.rm = T)})
          removeNotification(id)
        }
        id <- showNotification(h3("Done!"), duration = 1)
      })

    }),warning=function(war){},error=function(err){
      removeNotification(id)
      report_error(err)
    }))
  })
  
  #Data import
  observeEvent(input$plantcv_merge,{
    merged$data <- NULL; design$data <- NULL; shapes$data <- NULL; vis$data <- NULL; nir$data <- NULL; empties1$data <- NULL; from$data <- NULL; snapshot$data <- NULL; nir_ready_checker$data <- FALSE; vis_ready_checker$data <- FALSE; nir_caps$data <- NULL; vis_caps$data <- NULL; outlier_check$data <- FALSE; cooksd$data <- NULL; outlier_fmla$data <- NULL; imp_error_step$data <- NULL
    from$data <- "plantcv"
    res <- try(withCallingHandlers(withLogErrors({
    n <- 17
      withProgress(message = '', value = 0, {
        imp_error_step$data <- "PlantCV - Connecting to db"
        incProgress(1/n, detail = "Connecting to db...")
        db <- input$plantcv_sql_path$datapath
        drv <- dbDriver("SQLite")
        conn <- dbConnect(drv, dbname = db)
        removeNotification(id)
        
        imp_error_step$data <- "PlantCV - Reading design file"
        incProgress(1/n, detail = "Reading design file...")
        assoc <- read.csv(input$plantcv_design_file$datapath,header=T,stringsAsFactors = F)
        assoc_empty <- assoc[rowSums(sapply(colnames(assoc),function(i) !(assoc[,i] %in% c("Blank","Empty","blank","empty"))))==ncol(assoc),]
        design$data <- assoc
        removeNotification(id)
        
        imp_error_step$data <- "PlantCV - Reading snapshot file"
        incProgress(1/n, detail = "Reading snapshot file...")
        img_to_barcode <- read.csv(input$plantcv_snapshot_file$datapath,header = T,stringsAsFactors = F)
        img_to_barcode$timestamp <- as.POSIXct(strptime(img_to_barcode$timestamp,format = "%Y-%m-%d %H:%M:%S"))
        colnames(img_to_barcode)[3] <- "Barcodes"
        snapshot1 <- join(assoc_empty,img_to_barcode, by = "Barcodes")
        snapshot$data <- snapshot1[-snapshot1$weight.before < 0,]
        removeNotification(id)
        
        imp_error_step$data <- "PlantCV - Querying db for shapes data"
        incProgress(1/n, detail = "Querying db for shapes data...")
        meta <- colnames(dbGetQuery(conn = conn,'SELECT * FROM metadata'))
        shapes.df <- dbGetQuery(conn = conn,'SELECT * FROM metadata NATURAL JOIN features')
        colnames(shapes.df)[colnames(shapes.df) == "plantbarcode"] <- "Barcodes"
        shapes.df <- shapes.df[shapes.df$imgtype == "VIS",]
        shapes.df <- shapes.df[,which(!apply(shapes.df == 0, 2, all))]
        removeNotification(id)
        
        imp_error_step$data <- "PlantCV - Joining design file"
        incProgress(1/n, detail = "Joining design file...")
        colnames(shapes.df)[colnames(shapes.df) == "plantbarcode"] <- "Barcodes"
        sv_shapes <- join(shapes.df,assoc_empty,by="Barcodes")
        sv_shapes <- sv_shapes[rowSums(sapply(colnames(assoc),function(i) !is.na(sv_shapes[,i])))==ncol(assoc),]
        removeNotification(id)
        
        imp_error_step$data <- "PlantCV - Adding time columns"
        incProgress(1/n, detail = "Adding time columns...")
        sv_shapes$timestamp <- strptime(sv_shapes$timestamp,format = "%Y-%m-%d %H:%M:%S")
        beg <- min(sv_shapes$timestamp)
        sv_shapes$DAP <- floor(as.numeric((sv_shapes$timestamp - beg)/60/60/24))+as.numeric(input$dap_offset)
        sv_shapes$hour <- lubridate::hour(sv_shapes$timestamp)
        sv_shapes$timestamp <- as.character(sv_shapes$timestamp)
        removeNotification(id)
        
        imp_error_step$data <- "PlantCV - Removing empty pots"
        incProgress(1/n, detail = "Removing empty pots...")
        empties <- sv_shapes[sv_shapes$DAP == (max(sv_shapes$DAP)-1) & sv_shapes$area <10,"Barcodes"]
        sv_shapes <- sv_shapes[!(sv_shapes$Barcodes %in% empties),]
        colnames(sv_shapes) <- gsub("-","_",colnames(sv_shapes))
        sv_shapes[which(sv_shapes == Inf,arr.ind = T)] <- NaN
        removeNotification(id)
        
        #VIS data organization
        imp_error_step$data <- "PlantCV - Querying db for VIS data"
        incProgress(1/n, detail = "Querying db for VIS data...")
        vis.df <- dbGetQuery(conn = conn, 'SELECT * FROM metadata NATURAL JOIN signal WHERE channel_name = "hue"')
        if(nrow(vis.df)!= 0){
          colnames(vis.df)[colnames(vis.df) == "plantbarcode"] <- "Barcodes"
          vis.df <- cbind(data.frame("meta"=rep("meta",nrow(vis.df))),data.frame(do.call(rbind,lapply(strsplit(vis.df$values,", "),function(i)100*(as.numeric(i)/sum(as.numeric(i)))))),vis.df[,which(!(colnames(vis.df)%in%(c("bin_values","values"))))])
          removeNotification(id)
          
          incProgress(1/n, detail = "Joining with design...")
          vis.df <- join(vis.df,assoc_empty,by="Barcodes")
          vis.df <- vis.df[rowSums(sapply(colnames(assoc),function(i) !is.na(vis.df[,i])))==ncol(assoc),]
          removeNotification(id)
          
          incProgress(1/n, detail = "Adding time columns...")
          vis.df$timestamp <- strptime(vis.df$timestamp,format = "%Y-%m-%d %H:%M:%S")
          vis.df$DAP <- floor(as.numeric((vis.df$timestamp - beg)/60/60/24))+as.numeric(input$dap_offset)
          vis.df$hour <- lubridate::hour(vis.df$timestamp)
          vis.df$timestamp <- as.character(vis.df$timestamp)
          removeNotification(id)
          
          incProgress(1/n, detail = "Removing empty pots...")
          vis.df <- vis.df[!(vis.df$Barcodes %in% empties),]
          removeNotification(id)
          
          incProgress(1/n, detail = "Merging shapes data...")
          vis_shapes <- join(vis.df, sv_shapes, by = c("image_id","zoom","camera","cartag","exposure","frame","gain","id","imgtype","lifter","Line_name","measurementlabel","other","run_id","treatment","Treatment","Barcodes","timestamp","hour","DAP"))
          sv_shapes <- vis_shapes[,c(colnames(sv_shapes))]
          vis$data <- vis_shapes[,c(colnames(vis.df))]
          removeNotification(id)
        }else{
          removeNotification(id)
          id <- showNotification(h3("No VIS data detected"), duration = 1)
          vis$data <- NULL
        }
        
        #NIR data organization
        imp_error_step$data <- "PlantCV - Querying db for NIR data"
        incProgress(1/n, detail = "Querying db for NIR data...")
        nir.df <- dbGetQuery(conn = conn, 'SELECT * FROM metadata NATURAL JOIN signal WHERE channel_name = "nir"')
        if(nrow(nir.df)!= 0){
          colnames(nir.df)[colnames(nir.df) == "plantbarcode"] <- "Barcodes"
          nir.df <- cbind(data.frame("meta"=rep("meta",nrow(nir.df))),data.frame(do.call(rbind,lapply(strsplit(nir.df$values,", "),function(i)100*(as.numeric(i)/sum(as.numeric(i)))))),nir.df[,which(!(colnames(nir.df)%in%(c("bin_values","values"))))])
          removeNotification(id)
          
          incProgress(1/n, detail = "Joining with design...")
          nir.df <- join(nir.df,assoc_empty,by="Barcodes")
          nir.df <- nir.df[rowSums(sapply(colnames(assoc),function(i) !is.na(nir.df[,i])))==ncol(assoc),]
          removeNotification(id)
          
          incProgress(1/n, detail = "Adding time columns...")
          nir.df$timestamp <- strptime(nir.df$timestamp,format = "%Y-%m-%d %H:%M:%S")
          nir.df$DAP <- floor(as.numeric((nir.df$timestamp - beg)/60/60/24))+as.numeric(input$dap_offset)
          nir.df$hour <- lubridate::hour(nir.df$timestamp)
          removeNotification(id)
          
          incProgress(1/n, detail = "Removing empty pots...")
          nir.df <- nir.df[!(nir.df$Barcodes %in% empties),]
          nir.df$intensityAVG <- apply(nir.df[,2:256],1,function(i){sum((i/100)*(2:256),na.rm = T)})
          nir$data <- nir.df
          removeNotification(id)
        }else{
          removeNotification(id)
          id <- showNotification(h3("No NIR data detected"), duration = 1)
          nir$data <- NULL
        }
        
        incProgress(1/n, detail = "Finalizing...")
        imp_error_step$data <- "PlantCV - Finalizing"
        merged$data <- sv_shapes
        shapes$data <- merged$data[,colnames(merged$data)[colnames(merged$data) %in% c('image','image_id','area','hull_area','solidity','perimeter','width','height','longest_axis','center_of_mass_x','center_of_mass_y','hull_vertices','in_bounds','ellipse_center_x','ellipse_center_y','ellipse_major_axis','ellipse_minor_axis','ellipse_angle','ellipse_eccentricity','y_position','height_above_bound','height_below_bound','above_bound_area','percent_above_bound_area','below_bound_area','percent_below_bound_area')]]
        empties1$data <- data.frame("Barcodes" = empties, stringsAsFactors = F)
        
        id <- showNotification(h3("Done!"), duration = 1)
        dbDisconnect(conn)
      })

    }),warning=function(war){},error=function(err){
      removeNotification(id)
      dbDisconnect(conn)
      report_error(err)
    }))
  })
  
  #***********************************************************************************************
  # Outlier detection box
  #***********************************************************************************************
  output$outlier_removal <- renderUI({
    if(!is.null(merged$data)){
      box(width=10,title = "Outlier Detection and Removal",solidHeader = T,status = 'success',collapsible = TRUE,collapsed = TRUE,
          p("This step is not required"),
          if(!outlier_check$data){
            actionButton("detect_outliers","Detect Outliers")
            },
          actionButton("outlier_about",label = NULL,icon("question-circle"),style="background-color: white; border-color: white"),
          textOutput("outliers_model"),
          textOutput("num_outliers"),
          plotOutput("cooksd_plot"),
          uiOutput("remove_outliers_ui"),
          br(),
          uiOutput("download_cooks_ui")
      ) 
    }
  })
  
  cooksd <- reactiveValues(data=NULL)
  outlier_fmla <- reactiveValues(data=NULL)
  observeEvent(input$detect_outliers,{
    res <- try(withCallingHandlers(withLogErrors({
      imp_error_step$data <- "Detect outliers"
      disable("detect_outliers")
      id <- showNotification(h3("Calculating Cook's Distance..."), duration = NULL)
      des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
      outlier_fmla$data <- paste("as.numeric(area) ~",paste(c(des,"as.factor(DAP)"),collapse = ":"))
      fmla <- as.formula(outlier_fmla$data)
      cooksd$data <- cooks.distance(glm(data=merged$data,fmla))
      removeNotification(id)
    }),warning=function(war){},error=function(err){
      removeNotification(id)
      report_error(err)
    }))
  })
  
  output$num_outliers <- renderText({
    if(!is.null(cooksd$data)){
      n <- sum(cooksd$data >= 3*mean(cooksd$data),na.rm = T)
      paste("Outliers detected:",n,"(Approximately",round(100*(n/length(cooksd$data))),"% of dataset)")
    }else{""}
  })
  
  output$outliers_model <- renderText({
    if(!is.null(outlier_fmla$data)){
      paste0("Model: cooks.distance(glm(",outlier_fmla$data,"))")
    }else{""}
  })
  
  cooks_plot <- reactive({
    df <- data.frame("index"=1:length(cooksd$data),"cooksd"=cooksd$data)
    ggplot(df,aes(index,cooksd))+
      geom_point()+
      geom_hline(yintercept = 3*mean(cooksd$data),color="blue",linetype="dashed",size=2)+
      theme_light()+
      theme(axis.text = element_text(size = 14),
            axis.title= element_text(size = 18))+
      theme(strip.background=element_rect(fill="gray50"),
            strip.text.x=element_text(size=14,color="white"),
            strip.text.y=element_text(size=14,color="white"))
  })
  
  output$cooksd_plot <- renderPlot({
    if(!is.null(cooksd$data)){
      cooks_plot()
    }
  })
  
  output$cooksd_download <- downloadHandler(
    filename = function() {"pheno_cooksd_plot.png"},
    content=function(file){
      ggsave(file,cooks_plot(),device = "png",width = 10,height = 4,dpi = 300)
    })
  
  output$download_cooks_ui <- renderUI({
    if(!is.null(cooksd$data)){
      downloadButton("cooksd_download","Download Plot")
    }
  })
  
  output$remove_outliers_ui <- renderUI({
    if(!is.null(cooksd$data) & !outlier_check$data){
      actionButton("remove_outliers","Remove Outliers")
    }
  })
  
  outlier_check <-reactiveValues(data=FALSE)
  
  observeEvent(input$remove_outliers,{
    res <- try(withCallingHandlers(withLogErrors({
      imp_error_step$data <- "Removing outliers"
      id <- showNotification(h3("Removing from shapes, VIS, and NIR files..."), duration = NULL)
      merged$data <- merged$data[cooksd$data < 3*mean(cooksd$data),]
      if(from$data == "plantcv"){
        if(!is.null(vis$data)){
          vis$data <- vis$data[cooksd$data < 3*mean(cooksd$data),]
        }
        if(!is.null(nir$data)){
          nir$data <- nir$data[cooksd$data < 3*mean(cooksd$data),]
        }
      }else{
        vis$data <- vis$data[cooksd$data < 3*mean(cooksd$data),]
        outliers <- merged$data[cooksd$data >= 3*mean(cooksd$data),]
        outliers$camera_angle <- unlist(lapply(strsplit(outliers$meta,"_"),function(i) i[3]))
        outliers$unique_id <- paste(outliers$Barcodes,outliers$DAP,outliers$camera_angle,sep="_")
        if(input$pheno_nir_q == "Yes"){
          nir$data$camera_angle <- unlist(lapply(strsplit(as.character(nir$data$V1),"_"),function(i) i[3]))
          nir$data$unique_id <- paste(nir$data$Barcodes,nir$data$DAP,nir$data$camera_angle,sep="_")
          nir$data <- nir$data[!(nir$data$unique_id %in% outliers$unique_id),] 
        }
      }
      removeNotification(id)
      outlier_check$data <- TRUE
    }),warning=function(war){},error=function(err){
      removeNotification(id)
      report_error(err)
    }))
  })
  
  observeEvent(input$remove_outliers,{
    disable("remove_outliers")
    disable("detect_outliers")
  })
  
  #***********************************************************************************************
  # Shapes Box
  #***********************************************************************************************
  output$shapes_ui <- renderUI({
    des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
    s <- colnames(shapes$data)[!(colnames(shapes$data) %in% c("meta","image","image_id","in_bounds", "oof"))]
    if(!is.null(merged$data)){
      box(width=10,title = "Shapes Analysis",solidHeader = T,status = 'success',collapsible = TRUE,collapsed = TRUE,
          tabsetPanel(
            tabPanel(title="Shapes ANOVA",
                     selectInput("which_day","Which Day",sort(unique(merged$data$DAP)),max(unique(merged$data$DAP))),
                     div(id="container", actionButton("make_anova","Calculate ANOVA"),
                         actionButton("shapes_anova_about",label = NULL,icon("question-circle"),style="background-color: white; border-color: white")
                     ),
                     withSpinner(plotOutput("anova_plot"), type = 5),
                     uiOutput("download_shapes_anova_ui")
            ),
            tabPanel(title="Temporal ANOVA",
                    selectInput("anova_ts_shape","Which Shape",s,"area"),
                    div(id="container", actionButton("make_anova_ts","Calculate ANOVA"),
                        actionButton("anova_ts_about",label = NULL,icon("question-circle"),style="background-color: white; border-color: white")
                    ),
                    withSpinner(plotOutput("anova_ts_plot"), type = 5),
                    uiOutput("download_anova_ts_ui")
            ),
            tabPanel(title="Trends",
                     selectInput("dep_var","Y-axis",s,"area"),
                     selectInput("color_by","Color By",des,des[1]),
                     selectInput("facet_by","Facet By",des,des[2]),
                     textOutput("trends_collapsed_over"),
                     plotOutput("trends_plot"),
                     div(id="container", uiOutput("download_shapes_trends_ui"),
                        actionButton("trends_about",label = NULL,icon("question-circle"),style="background-color: white; border-color: white")
                     )
            ),
            tabPanel(title="Heatmap",
                     selectInput("h_color_by","Color By",s,"area"),
                     selectInput("h_group_by","Group By",des,des[1]),
                     selectInput("h_facet_by","Facet By",des,des[2]),
                     textOutput("h_collapsed_over"),
                     plotOutput("trends_heatmap"),
                     div(id="container", uiOutput("download_shapes_heatmap_ui"),
                        actionButton("heatmap_about",label = NULL,icon("question-circle"),style="background-color: white; border-color: white")
                     )
            ),
            tabPanel(title="Boxplots",
                     selectInput("box_dep_var","Y-axis",s,"area"),
                     selectInput("box_which_day","Which Day",sort(unique(merged$data$DAP)),max(unique(merged$data$DAP))),
                     selectInput("box_group_by","Group By",des,des[1]),
                     selectInput("box_facet_by","Facet By",des,des[2]),
                     textOutput("box_collapsed_over"),
                     plotOutput("boxplot_shapes"),
                     div(id="container", uiOutput("download_shapes_boxplot_ui"),
                        actionButton("boxplot_about",label = NULL,icon("question-circle"),style="background-color: white; border-color: white")
                     )
            )
          )
      ) 
    }
  })
  
  output$trends_collapsed_over <- renderText({
    des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
    left <- des[!(des %in% c(input$color_by,input$facet_by))]
    if(!length(left) == 0){
      paste0("Trend lines are collapsed over: ",paste(left,collapse=" ")) 
    }
  })
  
  output$box_collapsed_over <- renderText({
    des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
    left <- des[!(des %in% c(input$box_group_by,input$box_facet_by))]
    if(!length(left) == 0){
      paste0("Boxes are collapsed over: ",paste(left,collapse=" ")) 
    }
  })
  
  output$h_collapsed_over <- renderText({
    des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
    left <- des[!(des %in% c(input$h_group_by,input$h_facet_by))]
    if(!length(left) == 0){
      paste0("Heatmap values are collapsed over: ",paste(left,collapse=" ")) 
    }
  })
  

  #***********************************************************************************************
  # Shapes ANOVA
  #***********************************************************************************************
  anova_dat <- reactiveValues(data=NULL)
  observeEvent(input$make_anova,{
    res <- try(withCallingHandlers(withLogErrors({
      imp_error_step$data <- "Shapes Anova"
      id <- showNotification(h3("Calculating variances..."), duration = NULL)
      s <- colnames(shapes$data)[!(colnames(shapes$data) %in% c("meta","image","image_id","in_bounds", "oof"))]
      des <- sort(colnames(design$data)[!(colnames(design$data) %in% "Barcodes")])
      
      ext <- FALSE
      if(length(des)==2){
        ind_fmla <- paste0("(1|",des[1],")+(1|",des[2],")+(1|",des[1],":",des[2],")")
      }else{
        ind_fmla <- paste(paste0("(1|",des,")"),collapse = "+")
        ext <- TRUE
      }
      
      dat <- na.omit(merged$data[merged$data$DAP==as.numeric(input$which_day),])
      H2 <- c()
      for(e in s){
        fmla <- as.formula(paste0("as.numeric(",e,") ~ ",ind_fmla))
        model <- lmer(fmla,data = dat)
        re<- VarCorr(model)
        res<-attr(VarCorr(model), "sc")^2
        
        if(!ext){
          interaction.var <- as.numeric(attr(re[[which(str_detect(names(re),":"))]],"stddev"))^2
          des1.var <- as.numeric(attr(re[[des[1]]],"stddev"))^2
          des2.var <- as.numeric(attr(re[[des[2]]],"stddev"))^2
          
          tot.var<-sum(as.numeric(re),res)
          unexp <- 1-sum(as.numeric(re))/sum(as.numeric(re),res)
          
          h2 <- c((des1.var/tot.var),
            (des2.var/tot.var),
            (interaction.var/tot.var),
            unexp)
          H2 <- rbind(H2,h2) 
        }else{
          var <- lapply(des,function(i){as.numeric(attr(re[[i]],"stddev"))^2})
          
          tot.var <- sum(as.numeric(re),res)
          unexp <- 1-sum(as.numeric(re))/sum(as.numeric(re),res)
          
          h2 <- c(unlist(var)/tot.var,unexp)
          H2 <- rbind(H2,h2)
        }
      }
      H2 <- data.frame(H2,row.names = s)
      H2$Shape <- rownames(H2)
      rownames(H2) <- NULL
      if(!ext){
        colnames(H2) <- c(des[1],des[2],"Interaction","Unexplained","Shape")
      }else{
        colnames(H2) <- c(des,"Unexplained","Shape")
      }
      H2$Shape <-  ordered(H2$Shape,levels=H2$Shape[order(H2$Unexplained)])
      H2_melt <- melt(H2,id=c("Shape"))
      
      if(!ext){
        H2_melt$variable <- ordered(H2_melt$variable,levels=c("Unexplained",des[1],des[2],"Interaction"))
      }else{
        H2_melt$variable <- ordered(H2_melt$variable,levels=c("Unexplained",des))
      }
      anova_dat$data <- H2_melt
      removeNotification(id)
    }),warning=function(war){},error=function(err){
      removeNotification(id)
      report_error(err)
    }))
  })
  
  shapes_anova <- reactive({
    if(!is.null(anova_dat$data)){
      des <- sort(colnames(design$data)[!(colnames(design$data) %in% "Barcodes")])
      p <- ggplot(data=anova_dat$data,aes(Shape,value*100))+
        geom_bar(stat = "identity",aes(fill=variable))+
        ylab("Variance Explained (%)")+
        theme_bw()+
        theme(strip.background=element_rect(fill="gray50"),
              strip.text.x=element_text(size=14,color="white"),
              strip.text.y=element_text(size=14,color="white"))+
        theme(axis.text = element_text(size = 14),
              axis.title.y= element_text(size = 18),
              axis.title.x = element_blank())+
        theme(axis.ticks.length=unit(0.2,"cm"),
              plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
        theme(panel.border = element_rect(colour = "gray60", fill=NA, size=1,linetype = 1))+
        theme(legend.position = "top")+
        guides(fill = guide_legend(title = ""))+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      if(length(des) == 2){
        p <- p+scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))
      }else{
        hues <- seq(15, 375, length = length(des) + 1)
        cols <- hcl(h = hues, l = 65, c = 100)[1:length(des)]
        p <- p+scale_fill_manual(values = c("gray60",cols))
      }
      p
    }
  })
  
  output$anova_plot <- renderPlot({
    if(!is.null(anova_dat$data)){
      shapes_anova()
    }
  })
  
  output$shapes_anova_download <- downloadHandler(
    filename = function() {"shapes_anova.png"},
    content=function(file){
      ggsave(file,shapes_anova(),device = "png",width = 8,height = 4,dpi = 300)
    })
  
  output$download_shapes_anova_ui <- renderUI({
    if(!is.null(anova_dat$data)){
      downloadButton("shapes_anova_download","Download Plot")
    }
  })
  
  #***********************************************************************************************
  # Temporal ANOVA
  #***********************************************************************************************
  anova_ts_dat <- reactiveValues(data=NULL)
  observeEvent(input$make_anova_ts,{
    res <- try(withCallingHandlers(withLogErrors({
      imp_error_step$data <- "Temporal Anova"
      id <- showNotification(h3("Calculating variances..."), duration = NULL)
      which_shape <- input$anova_ts_shape
      des <- sort(colnames(design$data)[!(colnames(design$data) %in% "Barcodes")])
      
      ext <- FALSE
      if(length(des)==2){
        ind_fmla <- paste0("(1|",des[1],")+(1|",des[2],")+(1|",des[1],":",des[2],")")
      }else{
        ind_fmla <- paste(paste0("(1|",des,")"),collapse = "+")
        ext <- TRUE
      }
      
      H2 <- c()
      fmla <- as.formula(paste0("as.numeric(",which_shape,") ~ ",ind_fmla))
      
      ctab <- lapply(split(merged$data,merged$data$DAP),function(i) any(!apply(i[,des],2,function(j) length(unique(j))>1)))
      for(day in as.numeric(as.character(names(ctab)[!unlist(ctab)]))){
        dat <- na.omit(merged$data[merged$data$DAP == as.numeric(day),])
        model <- lmer(fmla,data = dat)
        re<- VarCorr(model)
        res<-attr(VarCorr(model), "sc")^2
        
        if(!ext){
          interaction.var <- as.numeric(attr(re[[which(str_detect(names(re),":"))]],"stddev"))^2
          des1.var <- as.numeric(attr(re[[des[1]]],"stddev"))^2
          des2.var <- as.numeric(attr(re[[des[2]]],"stddev"))^2
          
          tot.var<-sum(as.numeric(re),res)
          unexp <- 1-sum(as.numeric(re))/sum(as.numeric(re),res)
          
          h2 <- c((des1.var/tot.var),
            (des2.var/tot.var),
            (interaction.var/tot.var),
            unexp,
            day)
          H2 <- rbind(H2,h2) 
        }else{
          var <- lapply(des,function(i){as.numeric(attr(re[[i]],"stddev"))^2})
          
          tot.var <- sum(as.numeric(re),res)
          unexp <- 1-sum(as.numeric(re))/sum(as.numeric(re),res)
          
          h2 <- c(unlist(var)/tot.var,unexp,day)
          H2 <- rbind(H2,h2)
        }
      }
      H2 <- data.frame(H2)
      rownames(H2) <- NULL
      if(!ext){
        colnames(H2) <- c(des[1],des[2],"Interaction","Unexplained","Day")
      }else{
        colnames(H2) <- c(des,"Unexplained","Day")
      }
      
      H2_melt <- melt(H2,id=c("Day"))
      
      if(!ext){
        H2_melt$variable <- ordered(H2_melt$variable,levels=c("Unexplained",des[1],des[2],"Interaction"))
      }else{
        H2_melt$variable <- ordered(H2_melt$variable,levels=c("Unexplained",des))
      }
      H2_melt$shape <- which_shape
      anova_ts_dat$data <- H2_melt
      removeNotification(id)
    }),warning=function(war){},error=function(err){
      removeNotification(id)
      report_error(err)
    }))
  })
  
  anova_ts <- reactive({
    if(!is.null(anova_ts_dat$data)){
      des <- sort(colnames(design$data)[!(colnames(design$data) %in% "Barcodes")])
      which_shape <- anova_ts_dat$data$shape[1]
      p <- ggplot(data=anova_ts_dat$data,aes(as.numeric(Day),as.numeric(value)*100))+
        geom_line(aes(color=variable),size=1)+
        geom_point(aes(color=variable),size=3)+
        ggtitle(which_shape)+
        ylab("Variance Explained (%)")+
        xlab("DAP")+
        theme_bw()+
        theme(strip.background=element_rect(fill="gray50"),
          strip.text.x=element_text(size=14,color="white"),
          strip.text.y=element_text(size=14,color="white"))+
        theme(axis.text = element_text(size = 14),
          axis.title.y= element_text(size = 18))+
        theme(axis.ticks.length=unit(0.2,"cm"),
          plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
        theme(panel.border = element_rect(colour = "gray60", fill=NA, size=1,linetype = 1))+
        theme(legend.position = "right")+
        guides(color = guide_legend(title = ""))
      if(length(des) == 2){
        p <- p+scale_color_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))
      }else{
        hues <- seq(15, 375, length = length(des) + 1)
        cols <- hcl(h = hues, l = 65, c = 100)[1:length(des)]
        p <- p+scale_color_manual(values = c("gray60",cols))
      }
      p
    }
  })
  
  output$anova_ts_plot <- renderPlot({
    if(!is.null(anova_ts_dat$data)){
      anova_ts()
    }
  })
  
  output$anova_ts_download <- downloadHandler(
    filename = function() {"anova_ts.png"},
    content=function(file){
      ggsave(file,anova_ts(),device = "png",width = 8,height = 3.5,dpi = 300)
    })
  
  output$download_anova_ts_ui <- renderUI({
    if(!is.null(anova_ts_dat$data)){
      downloadButton("anova_ts_download","Download Plot")
    }
  })
  
  
  #***********************************************************************************************
  # Trends plots box
  #***********************************************************************************************
  shapes_trends <- reactive({
    ggplot(merged$data,aes_string("DAP",paste("as.numeric(",input$dep_var,")",collapse = "")))+
      facet_grid(~eval(parse(text=input$facet_by)))+
      geom_smooth(aes_string(color=input$color_by))+
      ylab(input$dep_var)+
      theme_light()+
      theme(axis.text = element_text(size = 14),
            axis.title= element_text(size = 18))+
      theme(strip.background=element_rect(fill="gray50"),
            strip.text.x=element_text(size=14,color="white"),
            strip.text.y=element_text(size=14,color="white"))
  })
  
  output$trends_plot <- renderPlot({
    shapes_trends()
  })
  
  output$shapes_trends_download <- downloadHandler(
    filename = function() {"shapes_trends.png"},
    content=function(file){
      ggsave(file,shapes_trends(),device = "png",width = 5.5,height = 4,dpi = 300)
    })
  
  output$download_shapes_trends_ui <- renderUI({
    if(!is.null(merged$data)){
      downloadButton("shapes_trends_download","Download Plot")
    }
  })
  
  
  #***********************************************************************************************
  # Growth heatmap box
  #***********************************************************************************************
  shapes_heatmap <- reactive({
    des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
    fmla <- as.formula(paste("as.numeric(",input$h_color_by,")","~",paste(c(des,"DAP"),collapse = "+")))
    df <- aggregate(data=merged$data,fmla,FUN = "mean")
    colnames(df)[ncol(df)] <- input$h_color_by
    ggplot(df,aes_string("DAP",input$h_group_by))+
      facet_grid(~eval(parse(text=input$h_facet_by)))+
      geom_tile(aes_string(fill=input$h_color_by))+
      ylab(input$h_group_by)+
      theme_light()+
      theme(axis.text = element_text(size = 14),
            axis.title= element_text(size = 18))+
      theme(strip.background=element_rect(fill="gray50"),
            strip.text.x=element_text(size=14,color="white"),
            strip.text.y=element_text(size=14,color="white"))
  })
  
  output$trends_heatmap <- renderPlot({
    shapes_heatmap()
  })
  
  output$shapes_heatmap_download <- downloadHandler(
    filename = function() {"shapes_heatmap.png"},
    content=function(file){
      ggsave(file,shapes_heatmap(),device = "png",width = 8,height = 4,dpi = 300)
    })
  
  output$download_shapes_heatmap_ui <- renderUI({
    if(!is.null(design$data)){
      downloadButton("shapes_heatmap_download","Download Plot")
    }
  })
  
  #***********************************************************************************************
  # Shapes Boxplots
  #***********************************************************************************************
  shapes_boxplot <- reactive({
    ggplot(merged$data[merged$data$DAP == input$box_which_day,],aes_string(input$box_group_by,paste("as.numeric(",input$box_dep_var,")",collapse = "")))+
      facet_grid(~eval(parse(text=input$box_facet_by)))+
      geom_violin(fill="gray50",alpha=.2)+
      geom_boxplot(width=.25)+
      ylab(input$box_dep_var)+
      theme_light()+
      theme(axis.text = element_text(size = 12),
            axis.title= element_text(size = 18))+
      theme(plot.title = element_text(hjust = 0.5),
            strip.background=element_rect(fill="gray50"),
            strip.text.x=element_text(size=14,color="white"),
            strip.text.y=element_text(size=14,color="white"))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$boxplot_shapes <- renderPlot({
    shapes_boxplot()
  })
  
  output$shapes_boxplot_download <- downloadHandler(
    filename = function() {"shapes_boxplot.png"},
    content=function(file){
      ggsave(file,shapes_boxplot(),device = "png",width = 8,height = 4,dpi = 300)
    })
  
  output$download_shapes_boxplot_ui <- renderUI({
    if(!is.null(merged$data)){
      downloadButton("shapes_boxplot_download","Download Plot")
    }
  })
  
  #***********************************************************************************************
  # Color Helpers
  #***********************************************************************************************
  hist_avg <- function(data,start,stop){
    des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
    sub <- data
    test <- data.frame(do.call("rbind",lapply(split(sub,sub[,des[1]]),function(t){
      data.frame(do.call("rbind",lapply(split(t,t[,des[2]]),function(g){
        data.frame(do.call("rbind",lapply(split(g,g$DAP),function(m){
          colMeans(m[,start:stop],na.rm = T)
        }
        )))
      }
      )))
    })))
    return(test)
  }
  
  hist_sd <- function(data,day,start,stop){
    sub <- data
    des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
    test <- data.frame(do.call("rbind",lapply(split(sub,sub[,des[1]]),function(t){
      data.frame(do.call("rbind",lapply(split(t,t[,des[2]]),function(g){
        data.frame(do.call("rbind",lapply(split(g,g$DAP),function(m){
          apply(m[,start:stop],2,function(i){sd(i,na.rm = T)})
        }
        )))
      }
      )))
    })))
    return(test)
  }

  
  #***********************************************************************************************
  # VIS box
  #***********************************************************************************************
  output$vis_ui <- renderUI({
    if(!is.null(vis$data)){
      des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
      box(width=10,title = "VIS Analysis",solidHeader = T,status = 'success',collapsible = TRUE,collapsed = TRUE,
          tabsetPanel(
            tabPanel("CAPS",
                     br(),
                     column(width=4,
                            br(),
                            selectInput("vis_caps_main", width = 300,
                                        label = "Main effect: ",
                                        choices = c("--",des),
                                        selected = "--"),
                            tags$b("Partialled out variables:    "),
                            uiOutput("vis_caps_partial"),
                            selectInput("vis_caps_dist",width=300,
                                        label="Distance Type: ",
                                        choices = c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis"),
                                        selected = "euclidean"),
                            selectInput("vis_caps_which_day","Which Day:",sort(unique(vis$data$DAP)),max(unique(vis$data$DAP,na.rm = T)),width = 300),
                            div(id="container", actionButton("make_vis_caps", "Go"),
                                 actionButton("vis_caps_about",label = NULL,icon("question-circle"),style="background-color: white; border-color: white")
                             ),
                            br(),
                            textOutput("vis_caps_warning"),
                            br(),
                            uiOutput("download_vis_caps")
                     ),
                     column(width=7,
                            plotOutput("vis_caps_out")
                     )
            ),
            tabPanel(title="Joyplot",
                     br(),
                     column(width = 4,
                        selectInput("vis_joyplot_which_day","Which Day",sort(unique(vis$data$DAP)),max(unique(vis$data$DAP,na.rm = T)))
                     ),
                     column(width = 4,
                         sliderInput("hue_range","HUE Degree Range", 0, 360, c(0,150), 1)   
                     ),
                     plotOutput("vis_joyplot"),
                     br(),
                     br(),
                     br(),
                     br(),
                     br(),
                     div(id="container", uiOutput("download_vis_joyplot_ui"),
                        actionButton("vis_joyplot_about",label = NULL,icon("question-circle"),style="background-color: white; border-color: white")
                     )
            )
          ), style='height: 625px'
      ) 
    }
  })
  
  vis_not_main <- reactiveValues(data=NULL)
  vis_ready_checker <- reactiveValues(data=FALSE)
  vis_caps <- reactiveValues(data=NULL)
  
  output$vis_caps_partial <- renderUI({
    des <- sort(colnames(design$data)[!(colnames(design$data) %in% "Barcodes")])
    if(length(des)==1){
      paste("--")
    }else if(length(des)>1){
      p(paste(vis_not_main$data,collapse = ", "))
    }
  })
  
  output$vis_caps_warning <- renderText({
    paste("This may take a few minutes to calculate.")
  })
  
  observeEvent(input$vis_caps_main,{
    vis_ready_checker$data <- FALSE
    des <- sort(colnames(design$data)[!(colnames(design$data) %in% "Barcodes")])
    main <- input$vis_caps_main
    pos <- des
    vis_not_main$data <- pos[!(pos %in% main)]
  })
  
  observeEvent({ 
    input$vis_caps_main
    input$vis_caps_dist
    input$vis_caps_which_day
  },{enable("make_vis_caps")})
  
  observeEvent(input$make_vis_caps,{
    res <- try(withCallingHandlers(withLogErrors({
      imp_error_step$data <- "VIS CAPS"
      disable("make_vis_caps")
      id <- showNotification(h3("Subsetting VIS data..."), duration = NULL)
      des <- sort(colnames(design$data)[!(colnames(design$data) %in% "Barcodes")])
      sub <- vis$data[vis$data$DAP == input$vis_caps_which_day,]
      sub <- sub[which(rowSums(sub[,2:181])!=0),]
      y <- sub[,2:181]
      y <- floor(y[,2:ncol(y)])
      y <- y[,which(colSums(y)!=0)]
      x <- setNames(data.frame(sub[,des]),des)
      if(length(des)==1){
        form <- as.formula(paste0("y ~ ",input$vis_caps_main))
      }else if(length(des)>1){
        form <- as.formula(paste0("y ~ ",input$vis_caps_main,"+Condition(",paste0(vis_not_main$data,collapse="*"),")",collapse = ""))
      }
      removeNotification(id)
      
      id <- showNotification(h3("Calculating CAPS..."), duration = NULL)
      cap <- capscale(form,x,dist=input$vis_caps_dist,sqrt.dist = T)
      removeNotification(id)
      
      id <- showNotification(h3("Done!"), duration = NULL)
      test <- data.frame(summary(cap)$sites)
      test <- cbind(test,x)
      vis_caps$data <- test
      vis_ready_checker$data <- TRUE
      removeNotification(id)
    }),warning=function(war){},error=function(err){
      removeNotification(id)
      report_error(err)
    }))

  })

  vis_caps_plot <- reactive({
    if(vis_ready_checker$data){
      ggplot(vis_caps$data,aes(eval(parse(text=colnames(vis_caps$data)[1])),eval(parse(text=colnames(vis_caps$data)[2]))))+
        geom_point(aes(color=factor(eval(parse(text=input$vis_caps_main)))))+
        stat_ellipse(aes(fill=factor(eval(parse(text=input$vis_caps_main))),color=factor(eval(parse(text=input$vis_caps_main)))),geom = "polygon",alpha=0.25)+
        xlab(colnames(vis_caps$data)[1])+
        ylab(colnames(vis_caps$data)[2])+
        geom_hline(yintercept=0,linetype="dashed",color="gray20",size=1)+
        geom_vline(xintercept=0,linetype="dashed",color="gray20",size=1)+
        theme_light()+
        theme(strip.background=element_rect(fill="gray50"),
              strip.text.x=element_text(size=14,color="white"),
              strip.text.y=element_text(size=14,color="white"))+
        theme(axis.title= element_text(size = 18))+
        theme(axis.text = element_text(size = 14))+
        theme(axis.ticks.length=unit(0.2,"cm"))+
        guides(color = guide_legend(title = input$vis_caps_main))+
        guides(fill = guide_legend(title = input$vis_caps_main))
    }
  })
  
  output$vis_caps_out <- renderPlot({
    if(input$vis_caps_main != "--"){
      vis_caps_plot()
    }else{
      ggplot()
    }
  })
  
 output$vis_caps_download <- downloadHandler(
   filename = function() {"vis_caps.png"},
   content=function(file){
     ggsave(file,vis_caps_plot(),device = "png",width = 5.2,height = 4,dpi = 300)
   })

 output$download_vis_caps <- renderUI({
   if(vis_ready_checker$data){
     downloadButton("vis_caps_download","Download Plot")
   }
 })
  
  vis_joyplot <- reactive({
    res <- try(withCallingHandlers(withLogErrors({
      imp_error_step$data <- "VIS Joyplot"
      sub <- vis$data[vis$data$DAP==input$vis_joyplot_which_day,]
      test_avg <- hist_avg(sub,start = 2,stop = 181)
      test_sd <- hist_sd(sub,start = 2,stop = 181)
      test_avg <- data.frame(melt(t(test_avg)))
      test_avg$sd <- data.frame(melt(t(test_sd)))[,3]
      test_avg$bin <- (2*(as.numeric(str_sub(test_avg$Var1,2,4))))
      test_avg$meta1 <- unlist(lapply(strsplit(as.character(test_avg$Var2),"[.]"),function(i)i[1]))
      test_avg$meta2 <- unlist(lapply(strsplit(as.character(test_avg$Var2),"[.]"),function(i)i[2]))
    }),warning=function(war){},error=function(err){
       removeNotification(id)
       report_error(err)
    }))
    if(class(res) != "try-error"){
      ggplot(data=test_avg,aes(x=bin,y=meta1, height=value))+
        facet_grid(~meta2)+
        geom_density_ridges_gradient(stat = "identity", aes(fill=bin),alpha=0.5, scale = 1)+
        scale_fill_gradientn(colors=hue_pal(l=65)(360))+
        scale_x_continuous(limits=c(input$hue_range[1],input$hue_range[2]),oob = rescale_none)+
        scale_y_discrete(expand = c(0.01, 0)) +
        ylab("")+
        xlab("Hue Channel")+
        #theme_ridges(grid=T,center_axis_labels = T)+
        theme_light()+
        theme(legend.position='none')+
        theme(axis.text = element_text(size = 12),
          axis.title= element_text(size = 18))+
        theme(plot.title = element_text(hjust = 0.5),
          strip.background=element_rect(fill="gray50"),
          strip.text.x=element_text(size=12,color="white"),
          strip.text.y=element_text(size=14,color="white"))+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  })
  
  output$vis_joyplot <- renderPlot({
    vis_joyplot()
  })
  
  output$vis_joyplot_download <- downloadHandler(
    filename = function() {"vis_joyplot.png"},
    content=function(file){
      ggsave(file,vis_joyplot(),device = "png",width = 8,height = 4,dpi = 300)
    })
  
  output$download_vis_joyplot_ui <- renderUI({
    if(!is.null(vis$data)){
      downloadButton("vis_joyplot_download","Download Plot")
    }
  })
  
  
  #***********************************************************************************************
  # NIR Analysis
  #***********************************************************************************************
  output$nir_ui <- renderUI({
    if(!is.null(nir$data)){
      des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
      box(width=10,title = "NIR Analysis",solidHeader = T,status = 'success',collapsible = TRUE,collapsed = TRUE,
          tabsetPanel(
            tabPanel("CAPS",
                     column(width=4,
                            br(),
                            selectInput("nir_caps_main", width = 300,
                                        label = "Main effect: ",
                                        choices = c("--",des),
                                        selected = "--"),
                            tags$b("Partialled out variables:    "),
                            uiOutput("nir_caps_partial"),
                            selectInput("nir_caps_dist",width=300,
                                        label="Distance Type: ",
                                        choices = c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis"),
                                        selected = "euclidean"),
                            selectInput("nir_caps_which_day","Which Day:",sort(unique(nir$data$DAP)),max(unique(nir$data$DAP,na.rm = T)),width = 300),
                            div(id="container", actionButton("make_nir_caps", "Go"),
                               actionButton("nir_caps_about",label = NULL,icon("question-circle"),style="background-color: white; border-color: white")
                            ),
                            textOutput("nir_caps_warning"),
                            br(),
                            uiOutput("download_nir_caps")
                     ),
                     column(width=7,
                            withSpinner(plotOutput("nir_caps_out"), type = 5)
                     )
            ),
            tabPanel(title="Heatmap",
                     selectInput("nir_day_start", "Day Start",sort(unique(nir$data$DAP)),min(unique(nir$data$DAP),na.rm = T)),
                     selectInput("nir_collapse_by", "Collapse By",des,des[1]),
                     plotOutput("nir_heatmap_nofacet"),
                     div(id="container", uiOutput("download_nir_heatmap_nofacet_ui"),
                         actionButton("nir_heatmap_nofacet_about",label = NULL,icon("question-circle"),style="background-color: white; border-color: white")
                     ),
                     br(),
                     plotOutput("nir_heatmap_withfacet"),
                     div(id="container", uiOutput("download_nir_heatmap_facet_ui"),
                        actionButton("nir_heatmap_withfacet_about",label = NULL,icon("question-circle"),style="background-color: white; border-color: white")
                     )
            )
          )
      )
    }
  })
  
  nir_not_main <- reactiveValues(data=NULL)
  nir_ready_checker <- reactiveValues(data=FALSE)
  nir_caps <- reactiveValues(data=NULL)
  
  
  output$nir_caps_partial <- renderUI({
    des <- sort(colnames(design$data)[!(colnames(design$data) %in% "Barcodes")])
    if(length(des)==1){
      paste("--")
    }else if(length(des)>1){
      p(paste(nir_not_main$data,collapse = ", "))
    }
  })
  
  output$nir_caps_warning <- renderText({
    paste("This may take a few minutes to calculate.")
  })
  
  observeEvent(input$nir_caps_main,{
    nir_ready_checker$data <- FALSE
    des <- sort(colnames(design$data)[!(colnames(design$data) %in% "Barcodes")])
    main <- input$nir_caps_main
    pos <- des
    nir_not_main$data <- pos[!(pos %in% main)]
  })
  
  observeEvent({ 
    input$nir_caps_main
    input$nir_caps_dist
    input$nir_caps_which_day
  },{enable("make_nir_caps")})
  
  observeEvent(input$make_nir_caps,{
    res <- try(withCallingHandlers(withLogErrors({
      imp_error_step$data <- "NIR CAPS"
      disable("make_nir_caps")
      id <- showNotification(h3("Subsetting NIR data..."), duration = NULL)
      des <- sort(colnames(design$data)[!(colnames(design$data) %in% "Barcodes")])
      sub <- nir$data[nir$data$DAP == input$nir_caps_which_day,]
      sub <- sub[which(rowSums(sub[,2:256])!=0),]
      y <- sub[,2:256]
      y <- floor(y[,2:ncol(y)])
      y <- y[,which(colSums(y)!=0)]
      x <- setNames(data.frame(sub[,des]),des)
      if(length(des)==1){
        form <- as.formula(paste0("y ~ ",input$nir_caps_main))
      }else if(length(des)>1){
        form <- as.formula(paste0("y ~ ",input$nir_caps_main,"+Condition(",paste0(nir_not_main$data,collapse="*"),")",collapse = ""))
      }
      removeNotification(id)
      
      id <- showNotification(h3("Calculating CAPS..."), duration = NULL)
      cap <- capscale(form,x,dist=input$nir_caps_dist,sqrt.dist = T)
      removeNotification(id)
      
      id <- showNotification(h3("Done!"), duration = NULL)
      test <- data.frame(summary(cap)$sites)
      test <- cbind(test,x)
      nir_caps$data <- test
      nir_ready_checker$data <- TRUE
      removeNotification(id)
    }),warning=function(war){},error=function(err){
      removeNotification(id)
      report_error(err)
    }))
  })
  
  nir_caps_plot <- reactive({
    if(nir_ready_checker$data){
      ggplot(nir_caps$data,aes(eval(parse(text=colnames(nir_caps$data)[1])),eval(parse(text=colnames(nir_caps$data)[2]))))+
        geom_point(aes(color=factor(eval(parse(text=input$nir_caps_main)))))+
        stat_ellipse(aes(fill=factor(eval(parse(text=input$nir_caps_main))),color=factor(eval(parse(text=input$nir_caps_main)))),geom = "polygon",alpha=0.25)+
        xlab(colnames(nir_caps$data)[1])+
        ylab(colnames(nir_caps$data)[2])+
        geom_hline(yintercept=0,linetype="dashed",color="gray20",size=1)+
        geom_vline(xintercept=0,linetype="dashed",color="gray20",size=1)+
        theme_light()+
        theme(strip.background=element_rect(fill="gray50"),
              strip.text.x=element_text(size=14,color="white"),
              strip.text.y=element_text(size=14,color="white"))+
        theme(axis.title= element_text(size = 18))+
        theme(axis.text = element_text(size = 14))+
        theme(axis.ticks.length=unit(0.2,"cm"))+
        guides(color = guide_legend(title = input$nir_caps_main))+
        guides(fill = guide_legend(title = input$nir_caps_main))
    }
  })
  
  output$nir_caps_out <- renderPlot({
    if(input$nir_caps_main != "--"){
      nir_caps_plot()
    }else{
      ggplot()
    }
  })
  
  output$nir_caps_download <- downloadHandler(
    filename = function() {"nir_caps.png"},
    content=function(file){
      ggsave(file,nir_caps_plot(),device = "png",width = 5.2,height = 4,dpi = 300)
    })
  
  output$download_nir_caps <- renderUI({
    if(nir_ready_checker$data){
      downloadButton("nir_caps_download","Download Plot")
    }
  })
  
  nir_heatmap_nofacet <- reactive({
    test <- aggregate(data=nir$data[nir$data$intensityAVG != 0 & nir$data$DAP >= as.numeric(input$nir_day_start),],as.formula(paste("intensityAVG~",input$nir_collapse_by,"+DAP",collapse="")),FUN = function(i)mean(i,na.rm=T))
    ggplot(test,aes_string("DAP",paste("as.factor(",input$nir_collapse_by,")",collapse = "")))+
      ylab(input$nir_collapse_by)+
      geom_tile(aes(fill=intensityAVG))+
      scale_fill_gradient2(midpoint = mean(test$intensityAVG),high ="gray10",low= "#56B1F7",mid = "#d7e4ef")+
      theme_light()+
      theme(axis.text = element_text(size = 12),
            axis.title= element_text(size = 18))+
      theme(plot.title = element_text(hjust = 0.5),
            strip.background=element_rect(fill="gray50"),
            strip.text.x=element_text(size=14,color="white"),
            strip.text.y=element_text(size=14,color="white"))
  })
  
  output$nir_heatmap_nofacet <- renderPlot({
    nir_heatmap_nofacet()
  })
  
  output$nir_heatmap_nofacet_download <- downloadHandler(
    filename = function() {"nir_heatmap_nofacet.png"},
    content=function(file){
      ggsave(file,nir_heatmap_nofacet(),device = "png",width = 8,height = 4,dpi = 300)
    })
  
  output$download_nir_heatmap_nofacet_ui <- renderUI({
    if(!is.null(nir$data)){
      downloadButton("nir_heatmap_nofacet_download","Download Plot")
    }
  })
  
  nir_heatmap_facet <- reactive({
    des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
    test <- aggregate(data=nir$data[nir$data$intensityAVG != 0 & nir$data$DAP >= as.numeric(input$nir_day_start),],as.formula(paste("intensityAVG~",des[1],"+",des[2],"+DAP")),FUN = function(i)mean(i,na.rm=T))
    ggplot(test,aes_string("DAP",des[1]))+
      facet_grid(~eval(parse(text=des[2])))+
      geom_tile(aes(fill=intensityAVG))+
      scale_fill_gradient2(high ="gray10",low= "#56B1F7",midpoint = mean(test$intensityAVG))+
      theme_light()+
      theme(axis.text = element_text(size = 12),
            axis.title= element_text(size = 18))+
      theme(plot.title = element_text(hjust = 0.5),
            strip.background=element_rect(fill="gray50"),
            strip.text.x=element_text(size=14,color="white"),
            strip.text.y=element_text(size=14,color="white"))
  })
  
  output$nir_heatmap_withfacet <- renderPlot({
    nir_heatmap_facet()
  })
  
  output$nir_heatmap_facet_download <- downloadHandler(
    filename = function() {"nir_heatmap_facet.png"},
    content=function(file){
      ggsave(file,nir_heatmap_facet(),device = "png",width = 8,height = 4,dpi = 300)
    })
  
  output$download_nir_heatmap_facet_ui <- renderUI({
    if(!is.null(nir$data)){
      downloadButton("nir_heatmap_facet_download","Download Plot")
    }
  })
  
  #***********************************************************************************************
  # "Summary" Box
  #***********************************************************************************************
  output$summary_ui <- renderUI({
    if(!is.null(snapshot$data)){
      des <- sort(colnames(design$data)[!(colnames(design$data) %in% "Barcodes")])
      box(width=10,title = "Summary",solidHeader = T,status = 'success',collapsible = TRUE,collapsed = TRUE,
          tabsetPanel(
            tabPanel(title = "Empties",
                    textOutput("num_empties"),
                    tableOutput("empties_table")
            ),
            tabPanel(title = "Image Quality",
                     br(),
                     textOutput("no_det"),       
                     column(4, 
                            numericInput("iqv_ylim_up", "Upper y-axis limit:", value = 10, width = 180)
                     ),
                     column(4, 
                            numericInput("iqv_ylim_low", "Lower y-axis limit:", value = -100, width= 180)
                     ),
                     br(),
                     plotOutput("iqv_plot", width = 400),
                     br(),
                     br(),
                     br(),
                     div(id="container", uiOutput("iqv_download_ui"),
                          actionButton("iqv_about",label = NULL,icon("question-circle"),style="background-color: white; border-color: white")
                      ),
                     br()
            ),
            tabPanel(title = "Water",
                   selectInput("water_facet_by", "Facet By:", des, des[1], width = 180),
                   selectInput("water_color_by", "Color By:", des, des[2], width = 180),
                   selectInput("water_var", "Water Measure:", c("weight.before","weight.after","water.amount"), "weight.before", width = 180),
                   plotOutput("water_plot"),
                   div(id="container", uiOutput("water_download_ui"),
                       actionButton("water_about",label = NULL,icon("question-circle"),style="background-color: white; border-color: white")
                   )
                   
            ),
            tabPanel(title = "OOF",
                   textOutput("oof_warn"),
                   withSpinner(plotOutput("oof_plot", height = 650), type = 5),
                   div(id="container", uiOutput("oof_download_ui"),
                       actionButton("oof_about",label = NULL,icon("question-circle"),style="background-color: white; border-color: white")
                   )
                   
            ),
            tabPanel(title = "Emergence Rate",
                    textOutput("er_warn"),
                    withSpinner(plotOutput("er_plot", height = 650), type = 5),
                    div(id="container", uiOutput("er_download_ui"),
                        actionButton("er_about",label = NULL,icon("question-circle"),style="background-color: white; border-color: white")
                    )
                    
            )
          )
      )
    }
  })
  
  output$num_empties <- renderText({
    if(!is.null(shapes$data)){
      paste0("Number of empty/blank pots: ", length(unique(empties1$data$Barcodes)))
    }
  })
  
  output$empties_table <- renderTable({
    des <- sort(colnames(design$data)[!(colnames(design$data) %in% "Barcodes")])
    fmla <- as.formula(paste("Barcodes~",paste(des,collapse = "+")))
    aggregate(data = design$data[which(design$data$Barcodes %in% empties1$data$Barcodes),], fmla, FUN = function(i) length(unique(i)))
  })
  
  iqv <- reactive({
    ggplot(merged$data, aes(x = hour, y = det))+
      geom_jitter()+
      scale_x_continuous(breaks = seq(0, 24, 2), limits = c(0, 24))+
      ylim(input$iqv_ylim_low, input$iqv_ylim_up)+
      ylab("Deviance (D)")+
      theme_light()+
      theme(axis.text = element_text(size = 12),
            axis.title= element_text(size = 18))+
      theme(plot.title = element_text(hjust = 0.5),
            strip.background=element_rect(fill="gray50"),
            strip.text.x=element_text(size=14,color="white"),
            strip.text.y=element_text(size=14,color="white"))
  })
  
  output$iqv_plot <- renderPlot({
    if(!is.null(merged$data$det)){
      iqv()
    }
  })
  
  output$iqv_download <- downloadHandler(
    filename = function() {"iqv_plot.png"},
    content = function(file){
      ggsave(file,iqv(),device = "png",width = 5,height = 4,dpi = 300)
    }
  )
  
  output$iqv_download_ui <- renderUI({
    if(!is.null(merged$data)){
      downloadButton("iqv_download","Download Plot")
    }
  })
  
  output$no_det <- renderText({
    if(is.null(merged$data$det)){
    paste("Image quality data was not found.")
    }
  })
  
  water <- reactive({
    ggplot(snapshot$data, aes_string(x = "timestamp", y = input$water_var))+
      geom_point(aes_string(color = paste0("as.factor(",input$water_color_by,")")))+
      facet_grid(~eval(parse(text=input$water_facet_by)))+
      theme_light()+
      theme(axis.text = element_text(size = 12),
            axis.title= element_text(size = 18))+
      theme(plot.title = element_text(hjust = 0.5),
            strip.background=element_rect(fill="gray50"),
            strip.text.x=element_text(size=14,color="white"),
            strip.text.y=element_text(size=14,color="white"))+
      guides(color = guide_legend(title = input$water_color_by))
  })
  
  output$water_plot <- renderPlot({
    water()
  })
  
  output$water_download <- downloadHandler(
    filename = function() {"water_plot.png"},
    content=function(file){
      ggsave(file,water(),device = "png",width = 9,height = 4,dpi = 300)
    })
  
  output$water_download_ui <- renderUI({
    if(!is.null(merged$data)){
      downloadButton("water_download","Download Plot")
    }
  })

  #*************************************************************************************************
  # OOF Survival
  #*************************************************************************************************
  output$oof_warn <- renderText({
    des <- sort(colnames(design$data)[!(colnames(design$data) %in% "Barcodes")])
    if(length(des)>3){
      "Design is too complex to show in one plot. Please subset your design file."
    }
  }) 
  
  oof_fig <- reactive({
    res <- try(withCallingHandlers(withLogErrors({
      imp_error_step$data <- "Survival Plot"
      des <- sort(colnames(design$data)[!(colnames(design$data) %in% "Barcodes")])
      if(from$data == "phenocv"){
        dat <- do.call("rbind",lapply(split(merged$data,merged$data$Barcodes),function(i) if(any(i$oof == 1)){
          sub <- i[i$oof == 1,]
          sub[order(sub$DAP),][1,]
        }else{
          i[i$DAP == max(i$DAP),][1,]
        }))
        dat$srv <- with(dat,Surv(time=DAP,event=oof))
      }else{
        dat <- do.call("rbind",lapply(split(merged$data,merged$data$Barcodes),function(i) if(any(i$in_bounds=="False")){
          sub <- i[i$in_bounds=="False",]
          sub[order(sub$DAP),][1,]
        }else{
          i[i$DAP == max(i$DAP),][1,]
        }))
        dat$srv <- with(dat,Surv(time=DAP,event=(in_bounds=="False")))
      }
      fmla <- as.formula(paste0("srv~",paste(des,collapse="+")))
      mod1 <- summary(survfit(fmla, data = dat, conf.type = "log-log"),time=min(merged$data$DAP):max(merged$data$DAP))
      mod_df <- data.frame("DAP"=mod1$time,"strata"=as.character(mod1$strata),"surv"=mod1$surv,"low"=mod1$lower,"high"=mod1$upper,stringsAsFactors = F)
      mod_df <- cbind(mod_df,setNames(data.frame(sapply(des,function(m){unlist(lapply(str_split(mod_df$strata,", "),function(i) trimws(str_split(i[str_detect(i,m)],"=")[[1]][2])))}),stringsAsFactors = F),des))
    }),warning=function(war){},error=function(err){
      removeNotification(id)
      report_error(err)
    }))
    if(class(res)!="try-error"){
      p <- ggplot(mod_df,aes(DAP,surv))
      if(length(des)>3){
        ggplot()
      }else if(length(des)==3){
        p <- p+facet_grid(as.formula(paste0(des[1],"~",des[2])))+
          geom_line(aes_string(color=des[3]))
      }else if(length(des)==2){
        p <- p+facet_grid(as.formula(paste0("~",des[1])))+
          geom_line(aes_string(color=des[2]))
      }else{
        p <- p+geom_line(aes_string(color=des[1]))
      }
      
      p <- p+ylab("Out Of Frame Risk")+
        scale_y_continuous(limits = c(0,1),breaks = seq(0,1,.2))+
        theme_light()+
        theme(axis.text = element_text(size = 14),
              axis.title= element_text(size = 18))+ 
        theme(strip.background=element_rect(fill="gray50"),
              strip.text.x=element_text(size=14,color="white"),
              strip.text.y=element_text(size=14,color="white"))
      p 
    }
  })
  
  output$oof_plot <- renderPlot({
    oof_fig()
  })
  
  output$oof_download <- downloadHandler(
    filename = function() {"oof_plot.png"},
    content=function(file){
      ggsave(file,oof_fig(),device = "png",width = 8,height = 4,dpi = 300)
    })
  
  output$oof_download_ui <- renderUI({
    if(!is.null(merged$data)){
      downloadButton("oof_download","Download Plot")
    }
  })
  
#*************************************************************************************************
# Emergence Rate (Germination/Survival)
#*************************************************************************************************

  output$er_warn <- renderText({
    des <- sort(colnames(design$data)[!(colnames(design$data) %in% "Barcodes")])
    if(length(des)>3){
      "Design is too complex to show in one plot. Please subset your design file."
    }
  })   
  
  er_fig <- reactive({
    res <- try(withCallingHandlers(withLogErrors({
      imp_error_step$data <- "Emergence Plot"
      des <- sort(colnames(design$data)[!(colnames(design$data) %in% "Barcodes")])
      #if(from$data == "phenocv"){
        dat <- do.call("rbind",lapply(split(merged$data,merged$data$Barcodes),function(i) if(any(i$area >= 10)){
          sub <- i[i$area >= 10,]
          sub[order(sub$DAP,decreasing=F),][1,]
        }else{
          sub <- i[i$DAP == max(i$DAP),][1,]
          sub[,"DAP"] <- sub[,"DAP"]+1
          sub[1,]
        }))
        dat$srv <- with(dat,Surv(time=DAP,event=(!DAP==(max(DAP)+1))))
        fmla <- as.formula(paste0("srv~",paste(des,collapse="+")))
        mod1 <- summary(survfit(fmla, data = dat, conf.type = "log-log"),time=min(merged$data$DAP):(max(merged$data$DAP)))
        mod_df <- data.frame("DAP"=mod1$time,"strata"=as.character(mod1$strata),"surv"=mod1$surv,"low"=mod1$lower,"high"=mod1$upper,stringsAsFactors = F)
        mod_df <- cbind(mod_df,setNames(data.frame(sapply(des,function(m){unlist(lapply(str_split(mod_df$strata,", "),function(i) trimws(str_split(i[str_detect(i,m)],"=")[[1]][2])))}),stringsAsFactors = F),des))
    }),warning=function(war){},error=function(err){
        removeNotification(id)
        report_error(err)
    }))
    if(class(res)!="try-error"){
      p <- ggplot(mod_df,aes(DAP,surv))
      if(length(des)>3){
        ggplot()
      }else if(length(des)==3){
        p <- p+facet_grid(as.formula(paste0(des[1],"~",des[2])))+
          geom_line(aes_string(color=des[3]))
      }else if(length(des)==2){
        p <- p+facet_grid(as.formula(paste0("~",des[1])))+
          geom_line(aes_string(color=des[2]))
      }else{
        p <- p+geom_line(aes_string(color=des[1]))
      }
      
       p <- p+ylab("Emergence Risk")+
         scale_y_continuous(limits = c(0,1),breaks = seq(0,1,.2))+
         theme_light()+
         theme(axis.text = element_text(size = 14),
               axis.title= element_text(size = 18))+ 
         theme(strip.background=element_rect(fill="gray50"),
               strip.text.x=element_text(size=14,color="white"),
               strip.text.y=element_text(size=14,color="white"))
       p
    }
    })
  
  output$er_plot <- renderPlot({
    er_fig()
  })
  
  output$er_download <- downloadHandler(
    filename = function() {"er_plot.png"},
    content=function(file){
      ggsave(file,er_fig(),device = "png",width = 8,height = 4,dpi = 300)
    })
  
  output$er_download_ui <- renderUI({
    if(!is.null(merged$data)){
      downloadButton("er_download","Download Plot")
    }
  })
  
  isolate({source("data/documentation.R",local=T)})
  
}

shinyApp(ui, server)