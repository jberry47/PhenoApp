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
                      fluidRow(
                        tags$head(tags$style("#container * {display: inline;}")),
                        tags$style(HTML("
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
                                      tags$li("While this tool is setup to handle more than two design columns, many of the analyses are
                                              coded to only accept two columns (such as Genotype and Treatment). If your design is more complicated than that, consider doing
                                              all analyses outside of this framework.")
                                      ),
                                    p("The first row of the design file should be: Barcodes, Var1, Var2. Where Var1 and Var2 could be called anything. The first 10 rows of two example design files are
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
                                    p("When processing images using this program, there are 2 files that are required: design, and sqlite3. The design file is shown 
                                      above and the sqlite3 file is created when running PlantCV with the",code("-s"), "flag. An example bash script to run image
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
                                    )
                                    ),
                        tabItem(tabName = "get_started",
                                box(width=10,title = "Merging Files",solidHeader = T,status = 'success',collapsible = TRUE,
                                    h5("How many days before the first day of imaging were the plants planted?"),
                                    textInput("dap_offset", "DAP Offset", value = 2,width=80),
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
                                               uiOutput("phenocv_go_ui")
                                      ),
                                      tabPanel(title = "PlantCV",
                                               fileInput("plantcv_design_file", "Choose design file",
                                                         multiple = F,
                                                         accept = c(".csv")),
                                               fileInput("plantcv_sql_path", "Choose sqlite3 database",
                                                         multiple = F,
                                                         accept = c(".sqlite3")),
                                               uiOutput("plantcv_go_ui")
                                      )
                                    )
                                ),
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
    b <- c(input$plantcv_sql_path$name,input$plantcv_design_file$name)
    if(length(b) == 2){
      actionButton("plantcv_merge","Merge Data")
    }
  })
  
  #***********************************************************************************************
  # Merging Files Action
  #***********************************************************************************************
  merged <- reactiveValues(data=NULL)
  design <- reactiveValues(data=NULL)
  shapes <- reactiveValues(data=NULL)
  vis <- reactiveValues(data=NULL)
  nir <- reactiveValues(data=NULL)
  empties <- reactiveValues(data=NULL)
  from <- reactiveValues(data=NULL)
  
  observeEvent(input$phenocv_merge,{
    merged$data <- NULL; design$data <- NULL; shapes$data <- NULL; vis$data <- NULL; nir$data <- NULL; empties$data <- NULL; from$data <- NULL;
    from$data <- "phenocv"
    id <- showNotification(h3("Reading snapshot file..."), duration = NULL)
    img_to_barcode <- read.csv(input$phenocv_snapshot_file$datapath,header = T,stringsAsFactors = F)
    img_to_barcode <- img_to_barcode[img_to_barcode$tiles != "",]
    colnames(img_to_barcode)[3] <- "Barcodes"
    img_to_barcode <- img_to_barcode[,c("id","Barcodes","timestamp")]
    removeNotification(id)
    
    id <- showNotification(h3("Reading shapes file..."), duration = NULL)
    sv_shapes <- read.table(input$phenocv_shapes_file$datapath,header = F,stringsAsFactors = F,sep = " ")
    sv_shapes <- sv_shapes[,which(!as.logical(apply(sv_shapes,2,FUN=function(i) all(is.na(i)))))]
    
    if(ncol(sv_shapes)==22){
      colnames(sv_shapes) <- c("meta","area","hull_area","solidity","perimeter","width","height","cmx","cmy","hull_verticies","ex","ey","emajor","eminor","angle","eccen","circ","round","ar","fd","oof","det")
      shapes$data <- sv_shapes
    }else if(ncol(sv_shapes)==21){
      colnames(sv_shapes) <- c("meta","area","hull_area","solidity","perimeter","width","height","cmx","cmy","hull_verticies","ex","ey","emajor","eminor","angle","eccen","circ","round","ar","oof","fd")
      shapes$data <- sv_shapes
    }else{
      colnames(sv_shapes) <- c("meta","area","hull_area","solidity","perimeter","width","height","cmx","cmy","hull_verticies","ex","ey","emajor","eminor","angle","eccen","circ","round","ar","oof")
      shapes$data <- sv_shapes
    }
    
    sv_shapes$id <- unlist(lapply(strsplit(sv_shapes$meta,"/"),function(i) strsplit(i[str_detect(i,"snapshot")],"snapshot")[[1]][2]))
    sv_shapes$imgname <- unlist(lapply(strsplit(sv_shapes$meta,"/"),function(i) strsplit(i[str_detect(i,"png")],"[.]")[[1]][1]))
    removeNotification(id)
    
    id <- showNotification(h3("Joining files..."), duration = NULL)
    sv_shapes <- join(sv_shapes,img_to_barcode[,c("id","Barcodes","timestamp")],by="id")
    removeNotification(id)
    
    id <- showNotification(h3("Reading design file..."), duration = NULL)
    assoc <- read.csv(input$phenocv_design_file$datapath,header=T,stringsAsFactors = F)
    design$data <- assoc
    removeNotification(id)
    
    id <- showNotification(h3("Joining files..."), duration = NULL)
    sv_shapes <- join(sv_shapes,assoc,by="Barcodes")
    sv_shapes <- sv_shapes[rowSums(sapply(colnames(assoc),function(i) !is.na(sv_shapes[,i])))==ncol(assoc),]
    removeNotification(id)
    
    id <- showNotification(h3("Adding time columns..."), duration = NULL)
    sv_shapes$timestamp <- strptime(sv_shapes$timestamp,format = "%Y-%m-%d %H:%M:%S")
    beg <- min(sv_shapes$timestamp)
    sv_shapes$DAP <- floor(as.numeric((sv_shapes$timestamp - beg)/60/60/24))+as.numeric(input$dap_offset)
    sv_shapes$hour <- lubridate::hour(sv_shapes$timestamp)
    removeNotification(id)
    
    id <- showNotification(h3("Removing empty pots..."), duration = NULL)
    sv_shapes <- sv_shapes[rowSums(sapply(colnames(assoc),function(i) !(sv_shapes[,i] %in% c("Blank","Empty","blank","empty"))))==ncol(assoc),]
    empties <- sv_shapes[sv_shapes$DAP == (max(sv_shapes$DAP)-1) & sv_shapes$area == 0,"Barcodes"]
    sv_shapes <- sv_shapes[!(sv_shapes$Barcodes %in% empties),]
    removeNotification(id)
    
    merged$data <- sv_shapes
    
    id <- showNotification(h3("Reading VIS color data..."), duration = NULL)
    vis$data <- get_color(input$phenocv_color_file$datapath,img_to_barcode,assoc,2,181)
    vis$data <- vis$data[rowSums(sapply(colnames(assoc),function(i) !(vis$data[,i] %in% c("Blank","Empty","blank","empty"))))==ncol(assoc),]
    vis$data <- vis$data[!(vis$data$Barcodes %in% empties),]
    vis$data <- vis$data[rowSums(sapply(colnames(assoc),function(i) !is.na(vis$data[,i])))==ncol(assoc),]
    removeNotification(id)
    
    if(input$pheno_nir_q == "Yes"){
      id <- showNotification(h3("Reading NIR color data..."), duration = NULL)
      nir$data <- get_color(input$phenocv_nir_file$datapath,img_to_barcode,assoc,2,256)
      nir$data <- nir$data[!(nir$data$Barcodes %in% empties),]
      nir$data <- nir$data[rowSums(sapply(colnames(assoc),function(i) !is.na(nir$data[,i])))==ncol(assoc),]
      nir$data$intensityAVG <- apply(nir$data[,2:256],1,function(i){sum((i/100)*(2:256),na.rm = T)})
      removeNotification(id)
    }
    id <- showNotification(h3("Done!"), duration = 1)
    
  })
  
  #Data import
  observeEvent(input$plantcv_merge,{
    merged$data <- NULL; design$data <- NULL; shapes$data <- NULL; vis$data <- NULL; nir$data <- NULL; empties$data <- NULL; from$data <- NULL;
    from$data <- "plantcv"
    id <- showNotification(h3("Connecting to db..."), duration = NULL)
    db <- input$plantcv_sql_path$datapath
    drv <- dbDriver("SQLite")
    conn <- dbConnect(drv, dbname = db)
    removeNotification(id)
    
    id <- showNotification(h3("Reading design file..."), duration = NULL)
    assoc <- read.csv(input$plantcv_design_file$datapath,header=T,stringsAsFactors = F)
    design$data <- assoc
    removeNotification(id)
    
    id <- showNotification(h3("Querying db for data..."), duration = NULL)
    meta <- colnames(dbGetQuery(conn = conn,'SELECT * FROM metadata'))
    shapes.df <- dbGetQuery(conn = conn,'SELECT * FROM metadata NATURAL JOIN features NATURAL JOIN signal WHERE signal.channel_name = "hue"')
#    shapes.df <- shapes.df[,as.numeric(which(colSums(shapes.df == "0") == 0))]
    removeNotification(id)
    
    id <- showNotification(h3("Extracting color data..."), duration = NULL)
    vis.df <- cbind(data.frame("meta"=rep("meta",nrow(shapes.df))),data.frame(do.call(rbind,lapply(strsplit(shapes.df$values,", "),function(i)100*(as.numeric(i)/sum(as.numeric(i)))))),shapes.df[,which((colnames(shapes.df)%in%meta))])
    shapes.df <- shapes.df[,(colnames(shapes.df)[!colnames(shapes.df) %in% c("bin-number","channel_name","values")])]
    colnames(shapes.df)[colnames(shapes.df) == "plantbarcode"] <- "Barcodes"
    colnames(vis.df)[colnames(vis.df) == "plantbarcode"] <- "Barcodes"
    shapes.df <- dbGetQuery(conn = conn,'SELECT * FROM metadata NATURAL JOIN features')
    shapes.df <- shapes.df[shapes.df$imgtype == "VIS",]
    shapes.df <- shapes.df[,which(!apply(shapes.df == 0, 2, all))]
    removeNotification(id)
    
    id <- showNotification(h3("Joining design file..."), duration = NULL)
    colnames(shapes.df)[colnames(shapes.df) == "plantbarcode"] <- "Barcodes"
    #shapes.df <- shapes.df[,(colnames(shapes.df)[!colnames(shapes.df) %in% c("bin-number","channel_name","values")])]
    sv_shapes <- join(shapes.df,assoc,by="Barcodes")
    sv_shapes <- sv_shapes[rowSums(sapply(colnames(assoc),function(i) !is.na(sv_shapes[,i])))==ncol(assoc),]
    removeNotification(id)
    
    id <- showNotification(h3("Adding time columns..."), duration = NULL)
    sv_shapes$timestamp <- strptime(sv_shapes$timestamp,format = "%Y-%m-%d %H:%M:%S")
    beg <- min(sv_shapes$timestamp)
    sv_shapes$DAP <- floor(as.numeric((sv_shapes$timestamp - beg)/60/60/24))+as.numeric(input$dap_offset)
    sv_shapes$hour <- lubridate::hour(sv_shapes$timestamp)
    sv_shapes$timestamp <- as.character(sv_shapes$timestamp)
    removeNotification(id)
    
    id <- showNotification(h3("Removing empty pots..."), duration = NULL)
    empties <- sv_shapes[sv_shapes$DAP == (max(sv_shapes$DAP)-1) & sv_shapes$area == 0,"Barcodes"]
    sv_shapes <- sv_shapes[rowSums(sapply(colnames(assoc),function(i) !(sv_shapes[,i] %in% c("Blank","Empty","blank","empty"))))==ncol(assoc),]
    sv_shapes <- sv_shapes[!(sv_shapes$Barcodes %in% empties),]
    colnames(sv_shapes) <- gsub("-","_",colnames(sv_shapes))
    removeNotification(id)
    
    #VIS data organization
    id <- showNotification(h3("Querying db for VIS data..."), duration = NULL)
    vis.df <- dbGetQuery(conn = conn, 'SELECT * FROM metadata NATURAL JOIN signal WHERE channel_name = "hue"')
    if(nrow(vis.df)!= 0){
#      vis.df <- vis.df[,which(!apply(vis.df == 0, 2, all))]
      colnames(vis.df)[colnames(vis.df) == "plantbarcode"] <- "Barcodes"
      vis.df <- cbind(data.frame("meta"=rep("meta",nrow(vis.df))),data.frame(do.call(rbind,lapply(strsplit(vis.df$values,", "),function(i)100*(as.numeric(i)/sum(as.numeric(i)))))),vis.df[,which(!(colnames(vis.df)%in%(c("bin_values","values"))))])
      removeNotification(id)
      
      id <- showNotification(h3("Joining with design..."), duration = NULL)
      vis.df <- join(vis.df,assoc,by="Barcodes")
      vis.df <- vis.df[rowSums(sapply(colnames(assoc),function(i) !is.na(vis.df[,i])))==ncol(assoc),]
      removeNotification(id)
      
      id <- showNotification(h3("Adding time columns..."), duration = NULL)
      vis.df$timestamp <- strptime(vis.df$timestamp,format = "%Y-%m-%d %H:%M:%S")
      vis.df$DAP <- floor(as.numeric((vis.df$timestamp - beg)/60/60/24))+as.numeric(input$dap_offset)
      vis.df$hour <- lubridate::hour(vis.df$timestamp)
      vis.df$timestamp <- as.character(vis.df$timestamp)
      removeNotification(id)
      
      id <- showNotification(h3("Removing empty pots..."), duration = NULL)
      vis.df <- vis.df[rowSums(sapply(colnames(assoc),function(i) !(vis.df[,i] %in% c("Blank","Empty","blank","empty"))))==ncol(assoc),]
      vis.df <- vis.df[!(vis.df$Barcodes %in% empties),]
      removeNotification(id)
      
      id <- showNotification(h3("Merging shapes data..."), duration = NULL)
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
    id <- showNotification(h3("Querying db for NIR data..."), duration = NULL)
    nir.df <- dbGetQuery(conn = conn, 'SELECT * FROM metadata NATURAL JOIN signal WHERE channel_name = "nir"')
    if(nrow(nir.df)!= 0){
#      nir.df <- nir.df[,as.numeric(which(colSums(nir.df == "0") == 0))]
      colnames(nir.df)[colnames(nir.df) == "plantbarcode"] <- "Barcodes"
      nir.df <- cbind(data.frame("meta"=rep("meta",nrow(nir.df))),data.frame(do.call(rbind,lapply(strsplit(nir.df$values,", "),function(i)100*(as.numeric(i)/sum(as.numeric(i)))))),nir.df[,which(!(colnames(nir.df)%in%(c("bin_values","values"))))])
      removeNotification(id)
      
      id <- showNotification(h3("Joining with design..."), duration = NULL)
      nir.df <- join(nir.df,assoc,by="Barcodes")
      nir.df <- nir.df[rowSums(sapply(colnames(assoc),function(i) !is.na(nir.df[,i])))==ncol(assoc),]
      removeNotification(id)
      
      id <- showNotification(h3("Adding time columns..."), duration = NULL)
      nir.df$timestamp <- strptime(nir.df$timestamp,format = "%Y-%m-%d %H:%M:%S")
      nir.df$DAP <- floor(as.numeric((nir.df$timestamp - beg)/60/60/24))+as.numeric(input$dap_offset)
      nir.df$hour <- lubridate::hour(nir.df$timestamp)
      removeNotification(id)
      
      id <- showNotification(h3("Removing empty pots..."), duration = NULL)
      nir.df <- nir.df[rowSums(sapply(colnames(assoc),function(i) !(nir.df[,i] %in% c("Blank","Empty","blank","empty"))))==ncol(assoc),]
      nir.df <- nir.df[!(nir.df$Barcodes %in% empties),]
      nir.df$intensityAVG <- apply(nir.df[,2:256],1,function(i){sum((i/100)*(2:256),na.rm = T)})
      nir$data <- nir.df
      removeNotification(id)
    }else{
      removeNotification(id)
      id <- showNotification(h3("No NIR data detected"), duration = 1)
      nir$data <- NULL
    }
    
    merged$data <- sv_shapes
    shapes$data <- merged$data[,colnames(merged$data)[colnames(merged$data) %in% c('image','image_id','area','hull_area','solidity','perimeter','width','height','longest_axis','center_of_mass_x','center_of_mass_y','hull_vertices','in_bounds','ellipse_center_x','ellipse_center_y','ellipse_major_axis','ellipse_minor_axis','ellipse_angle','ellipse_eccentricity','y_position','height_above_bound','height_below_bound','above_bound_area','percent_above_bound_area','below_bound_area','percent_below_bound_area')]]
    
    id <- showNotification(h3("Done!"), duration = 1)
    dbDisconnect(conn)
  })
  
  #***********************************************************************************************
  # Outlier detection box
  #***********************************************************************************************
  output$outlier_removal <- renderUI({
    if(!is.null(merged$data)){
      box(width=10,title = "Outlier Detection and Removal",solidHeader = T,status = 'success',collapsible = TRUE,collapsed = TRUE,
          p("This step is not required"),
          actionButton("detect_outliers","Detect Outliers"),
          textOutput("num_outliers"),
          plotOutput("cooksd_plot"),
          uiOutput("remove_outliers_ui"),
          br(),
          uiOutput("download_cooks_ui")
      ) 
    }
  })
  
  cooksd <- reactiveValues(data=NULL)
  observeEvent(input$detect_outliers,{
    id <- showNotification(h3("Calculating Cook's Distance..."), duration = NULL)
    des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
    fmla <- as.formula(paste("as.numeric(area) ~",paste(c(des,"as.factor(DAP)"),collapse = ":")))
    cooksd$data <- cooks.distance(glm(data=merged$data,fmla))
    removeNotification(id)
  })
  
  output$num_outliers <- renderText({
    if(!is.null(cooksd$data)){
      n <- sum(cooksd$data >= 3*mean(cooksd$data),na.rm = T)
      paste("Outliers detected:",n,"(Approximately",round(100*(n/length(cooksd$data))),"% of dataset)")
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
    if(!is.null(cooksd$data)){
      actionButton("remove_outliers","Remove Outliers")
    }
  })
  
  observeEvent(input$remove_outliers,{
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
  })
  
  #***********************************************************************************************
  # Shapes Box
  #***********************************************************************************************
  output$shapes_ui <- renderUI({
    des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
    s <- colnames(shapes$data)[!(colnames(shapes$data) %in% c("meta","image","image_id","in_bounds"))]
    if(!is.null(merged$data)){
      box(width=10,title = "Shapes Analysis",solidHeader = T,status = 'success',collapsible = TRUE,collapsed = TRUE,
          tabsetPanel(
            tabPanel(title="ANOVA",
                     selectInput("which_day","Which Day",sort(unique(merged$data$DAP)),max(unique(merged$data$DAP))),
                     actionButton("make_anova","Calculate ANOVA"),
                     plotOutput("anova_plot"),
                     uiOutput("download_shapes_anova_ui")
            ),
            tabPanel(title="Trends",
                     selectInput("dep_var","Y-axis",s,"area"),
                     selectInput("color_by","Color By",des,des[1]),
                     selectInput("facet_by","Facet By",des,des[2]),
                     plotOutput("trends_plot"),
                     uiOutput("download_shapes_trends_ui")
            ),
            tabPanel(title="Heatmap",
                     selectInput("h_color_by","Color By",s,"area"),
                     selectInput("h_group_by","Group By",des,des[1]),
                     selectInput("h_facet_by","Facet By",des,des[2]),
                     plotOutput("trends_heatmap"),
                     uiOutput("download_shapes_heatmap_ui")
            ),
            tabPanel(title="Boxplots",
                     selectInput("box_dep_var","Y-axis",s,"area"),
                     selectInput("box_which_day","Which Day",sort(unique(merged$data$DAP)),max(unique(merged$data$DAP))),
                     selectInput("box_group_by","Group By",des,des[1]),
                     selectInput("box_facet_by","Facet By",des,des[2]),
                     plotOutput("boxplot_shapes"),
                     uiOutput("download_shapes_boxplot_ui")
            )
          )
      ) 
    }
  })
  

  #***********************************************************************************************
  # ANOVA box
  #***********************************************************************************************
  anova_dat <- reactiveValues(data=NULL)
  observeEvent(input$make_anova,{
    id <- showNotification(h3("Calculating variances..."), duration = NULL)
    s <- colnames(shapes$data)[!(colnames(shapes$data) %in% c("meta","image","image_id","in_bounds"))]
    des <- sort(colnames(design$data)[!(colnames(design$data) %in% "Barcodes")])
    
    ext <- FALSE
    if(length(des)==2){
      ind_fmla <- paste0("(1|",des[1],")+(1|",des[2],")+(1|",des[1],":",des[2],")")
    }else{
      ind_fmla <- paste(paste0("(1|",des,")"),collapse = "+")
      ext <- TRUE
    }
    
    dat <- merged$data[merged$data$DAP==as.numeric(input$which_day),]
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
      geom_boxplot()+
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
    ggplot(merged$data[merged$data$DAP == input$box_which_day,],aes_string(input$box_group_by,paste("as.numeric(",input$box_dep_var,")",collapse = "")))+
      facet_grid(~eval(parse(text=input$box_facet_by)))+
      geom_boxplot()+
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
  
  makePCA <- function(data,day,start,stop,color_by){
    sub <- data[data$DAP==day,]
    channel.pca <- PCA(sub[,start:stop],graph = F)
    pca_df <- data.frame("Treatment"=as.character(sub[,color_by]),
                         "PC1"=channel.pca$ind$coord[,1],
                         "PC2"=channel.pca$ind$coord[,2])
    varexp <- signif(c(channel.pca$eig[1,2],channel.pca$eig[2,2]),4)
    
    p <- ggplot(data=pca_df, aes(PC1,PC2))+
      geom_point(data=aggregate(cbind(PC1,PC2)~Treatment,pca_df,mean),aes(color=Treatment),size=5)+
      stat_ellipse(aes(fill=Treatment,color=Treatment),geom = "polygon",alpha=0.25)+
      xlab(paste("PC1 (",varexp[1],"%)",sep = ""))+
      ylab(paste("PC2 (",varexp[2],"%)",sep = ""))+
      geom_vline(xintercept = 0,linetype="dashed")+
      geom_hline(yintercept = 0,linetype="dashed")+
      theme_minimal()+
      theme(legend.title=element_blank())+
      theme(axis.text = element_text(size = 18),
            axis.title= element_text(size = 24))+
      theme(panel.border = element_rect(colour = "gray60", fill=NA, size=1,linetype = 1))
    p
  }
  
  
  #***********************************************************************************************
  # VIS box
  #***********************************************************************************************
  output$vis_ui <- renderUI({
    if(!is.null(vis$data)){
      des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
      box(width=10,title = "VIS Analysis",solidHeader = T,status = 'success',collapsible = TRUE,collapsed = TRUE,
          tabsetPanel(
            tabPanel(title = "PCA",
                     selectInput("vis_which_day","Which Day",sort(unique(vis$data$DAP)),max(unique(vis$data$DAP),na.rm = T)),
                     selectInput("vis_color_by","Color By",des,des[1]),
                     plotOutput("vis_pca"),
                     uiOutput("download_vis_pca_ui")
            ),
            tabPanel(title="Joyplot",
                     selectInput("vis_joyplot_which_day","Which Day",sort(unique(vis$data$DAP)),max(unique(vis$data$DAP,na.rm = T))),
                     plotOutput("vis_joyplot"),
                     uiOutput("download_vis_joyplot_ui")
            )
          )
      ) 
    }
  })
  
  make_vis_pca <- reactive({
    makePCA(vis$data,input$vis_which_day,2,181,input$vis_color_by)
  })
  
  output$vis_pca <- renderPlot({
    makePCA(vis$data,input$vis_which_day,2,181,input$vis_color_by)
  })
  
  output$vis_pca_download <- downloadHandler(
    filename = function() {"vis_pca.png"},
    content=function(file){
      ggsave(file,make_vis_pca(),device = "png",width = 8,height = 4,dpi = 300)
    })
  
  output$download_vis_pca_ui <- renderUI({
    if(!is.null(vis$data)){
      downloadButton("vis_pca_download","Download Plot")
    }
  })
  
  vis_joyplot <- reactive({
    sub <- vis$data[vis$data$DAP==input$vis_joyplot_which_day,]
    test_avg <- hist_avg(sub,start = 2,stop = 181)
    test_sd <- hist_sd(sub,start = 2,stop = 181)
    test_avg <- data.frame(melt(t(test_avg)))
    test_avg$sd <- data.frame(melt(t(test_sd)))[,3]
    test_avg$bin <- (2*(as.numeric(str_sub(test_avg$Var1,2,4))))
    test_avg$meta1 <- unlist(lapply(strsplit(as.character(test_avg$Var2),"[.]"),function(i)i[1]))
    test_avg$meta2 <- unlist(lapply(strsplit(as.character(test_avg$Var2),"[.]"),function(i)i[2]))
    
    ggplot(data=test_avg,aes(x=bin,y=meta1, height=value))+
      facet_grid(~meta2)+
      geom_density_ridges(stat = "identity", aes(colour=meta2),alpha=0.5)+
      scale_x_continuous(breaks = c(0,90,180,270,360))+
      ylab("")+
      xlab("Hue Channel")+
      theme_ridges(grid=T,center_axis_labels = T)+
      theme(legend.position='none')+
      theme(axis.text = element_text(size = 12),
            axis.title= element_text(size = 18))+
      theme(plot.title = element_text(hjust = 0.5),
            strip.background=element_rect(fill="gray50"),
            strip.text.x=element_text(size=14,color="white"),
            strip.text.y=element_text(size=14,color="white"))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$vis_joyplot <- renderPlot({
    sub <- vis$data[vis$data$DAP==input$vis_joyplot_which_day,]
    test_avg <- hist_avg(sub,start = 2,stop = 181)
    test_sd <- hist_sd(sub,start = 2,stop = 181)
    test_avg <- data.frame(melt(t(test_avg)))
    test_avg$sd <- data.frame(melt(t(test_sd)))[,3]
    test_avg$bin <- (2*(as.numeric(str_sub(test_avg$Var1,2,4))))
    test_avg$meta1 <- unlist(lapply(strsplit(as.character(test_avg$Var2),"[.]"),function(i)i[1]))
    test_avg$meta2 <- unlist(lapply(strsplit(as.character(test_avg$Var2),"[.]"),function(i)i[2]))
    
    ggplot(data=test_avg,aes(x=bin,y=meta1, height=value))+
      facet_grid(~meta2)+
      geom_density_ridges(stat = "identity", aes(colour=meta2),alpha=0.5)+
      scale_x_continuous(breaks = c(0,90,180,270,360))+
      ylab("")+
      xlab("Hue Channel")+
      theme_ridges(grid=T,center_axis_labels = T)+
      theme(legend.position='none')+
      theme(axis.text = element_text(size = 12),
            axis.title= element_text(size = 18))+
      theme(plot.title = element_text(hjust = 0.5),
            strip.background=element_rect(fill="gray50"),
            strip.text.x=element_text(size=14,color="white"),
            strip.text.y=element_text(size=14,color="white"))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
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
            tabPanel(title = "PCA",
                     selectInput("nir_which_day","Which Day",sort(unique(nir$data$DAP)),max(unique(nir$data$DAP,na.rm = T))),
                     selectInput("nir_color_by","Color By",des,des[1]),
                     plotOutput("nir_pca"),
                     uiOutput("download_nir_PCA_ui")
            ),
            tabPanel(title="Heatmap",
                     selectInput("nir_day_start", "Day Start",sort(unique(nir$data$DAP)),min(unique(nir$data$DAP),na.rm = T)),
                     selectInput("nir_collapse_by", "Collapse By",des,des[1]),
                     plotOutput("nir_heatmap_nofacet"),
                     uiOutput("download_nir_heatmap_nofacet_ui"),
                     br(),
                     plotOutput("nir_heatmap_withfacet"),
                     uiOutput("download_nir_heatmap_facet_ui")
            )
          )
      )
    }
  })
  
  nir_make_pca <- reactive({
    makePCA(nir$data,input$nir_which_day,2,181,input$nir_color_by)
  })
  
  output$nir_pca <- renderPlot({
    makePCA(nir$data,input$nir_which_day,2,181,input$nir_color_by)
  })
  
  output$nir_pca_download <- downloadHandler(
    filename = function() {"nir_pca.png"},
    content=function(file){
      ggsave(file,nir_make_pca(),device = "png",width = 8,height = 4,dpi = 300)
    })
  
  output$download_nir_PCA_ui <- renderUI({
    if(!is.null(nir$data)){
      downloadButton("nir_pca_download","Download Plot")
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
  
}

shinyApp(ui, server)