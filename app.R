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


ui <- dashboardPage(skin="black", title="Phenotyping Analysis Tool",
                    dashboardHeader(
                      title = tagList(
                        tags$span(
                          class = "logo-mini", "PAT"
                        ),
                        tags$span(
                          class = "logo-lg", "Phenotyping Analysis Tool"
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
                                box(style = "overflow-y:scroll",width=10,title = "Welcome",solidHeader = T,status = 'success',collapsible = TRUE,
                                    p("test")
                                )
                        ),
                        tabItem(tabName = "get_started",
                                box(width=10,title = "Merging Files",solidHeader = T,status = 'success',collapsible = TRUE,
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
                                      fileInput("phenocv_nir_file", "Choose nir file",
                                        multiple = F,
                                        accept = c(".txt")),
                                      uiOutput("phenocv_go_ui")
                                    ),
                                    tabPanel(title = "PlantCV",
                                      fileInput("plantcv_sql_path", "Choose sqlite3 database",
                                        multiple = F,
                                        accept = c(".sqlite3")),
                                      fileInput("plantcv_design_file", "Choose design file",
                                        multiple = F,
                                        accept = c(".csv")),
                                      uiOutput("plantcv_go_ui")
                                    ),
                                    tabPanel(title = "Others?",
                                      p("DIRT and ImageJ are two example but I'm sure there are more")
                                    )
                                  )
                                ),
                                uiOutput("outlier_removal"),
                                uiOutput("shapes_anova"),
                                uiOutput("trends_plot_ui"),
                                uiOutput("heatmap_plot_ui"),
                                uiOutput("nir_heatmap_ui")
                        )
                      )
                    )
)

server <- function(input, output){
  
  get_color <- function(file_name,start,stop){
    color_data <- read.table(file_name,header = F,stringsAsFactors = F,sep = " ")[,-257]
    color_data$id <- as.character(sapply(color_data$V1,function(i) strsplit(strsplit(i,"/")[[1]][2],"snapshot")[[1]][2]))
    color_data$imgname <- as.character(sapply(color_data$V1,function(i) strsplit(strsplit(i,"/")[[1]][3],"[.]")[[1]][1]))
    color_data <- join(color_data,img_to_barcode[,c("id","Barcodes","timestamp")],by="id")
    color_data <- join(color_data,assoc,by="Barcodes")
    color_data$timestamp <- strptime(color_data$timestamp,format = "%Y-%m-%d %H:%M:%S")
    color_data$DAP <- floor(as.numeric((color_data$timestamp - beg)/60/60/24))+2
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
    if(length(b) == 5){
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
  observeEvent(input$phenocv_merge,{
    id <- showNotification(h3("Reading snapshot file..."), duration = NULL)
    img_to_barcode <- read.csv(input$phenocv_snapshot_file$datapath,header = T,stringsAsFactors = F)
    img_to_barcode <- img_to_barcode[img_to_barcode$tiles != "",]
    colnames(img_to_barcode)[3] <- "Barcodes"
    img_to_barcode <- img_to_barcode[,c("id","Barcodes","timestamp")]
    removeNotification(id)

    id <- showNotification(h3("Reading shapes file..."), duration = NULL)
    sv_shapes <- read.table(input$phenocv_shapes_file$datapath,header = F,stringsAsFactors = F,sep = " ")
    sv_shapes <- sv_shapes[,-ncol(sv_shapes)]
    colnames(sv_shapes) <- c("meta","area","hull_area","solidity","perimeter","width","height","cmx","cmy","hull_verticies","ex","ey","emajor","eminor","angle","eccen","circ","round","ar","fd","det")
    shapes$data <- sv_shapes
    sv_shapes$id <- substring(as.character(sapply(sv_shapes$meta,function(i) strsplit(i,"/")[[1]][2])),9)
    sv_shapes$imgname <- as.character(sapply(sv_shapes$meta,function(i) strsplit(strsplit(i,"/")[[1]][3],"[.]")[[1]][1]))
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
    sv_shapes$DAP <- floor(as.numeric((sv_shapes$timestamp - beg)/60/60/24))+2
    sv_shapes$hour <- lubridate::hour(sv_shapes$timestamp)
    removeNotification(id)
    
    id <- showNotification(h3("Removing empty pots..."), duration = NULL)
    sv_shapes <- sv_shapes[rowSums(sapply(colnames(assoc),function(i) !(sv_shapes[,i] %in% c("Blank","Empty","blank","empty"))))==ncol(assoc),]
    empties <- sv_shapes[sv_shapes$DAP == (max(sv_shapes$DAP)-1) & sv_shapes$area == 0,"Barcodes"]
    sv_shapes <- sv_shapes[!(sv_shapes$Barcodes %in% empties),]
    removeNotification(id)
    
    merged$data <- sv_shapes
    
    id <- showNotification(h3("Reading VIS color data..."), duration = NULL)
    vis$data <- get_color(input$phenocv_color_file$datapath,2,182)
    vis$data <- vis$data[!(vis$data$Barcodes %in% empties),]
    vis$data <- vis$data[rowSums(sapply(colnames(assoc),function(i) !is.na(vis$data[,i])))==ncol(assoc),]
    removeNotification(id)
    
    id <- showNotification(h3("Reading NIR color data..."), duration = NULL)
    nir$data <- get_color(input$phenocv_nir_file$datapath,2,256)
    nir$data <- nir$data[!(nir$data$Barcodes %in% empties),]
    nir$data <- nir$data[rowSums(sapply(colnames(assoc),function(i) !is.na(nir$data[,i])))==ncol(assoc),]
    nir$data$intensityAVG <- apply(nir$data[,3:255],1,function(i){sum((i/100)*(2:254),na.rm = T)})
    removeNotification(id)
    
    id <- showNotification(h3("Done!"), duration = 1)
    
  })
  
  observeEvent(input$plantcv_merge,{
    id <- showNotification(h3("Connecting to db..."), duration = NULL)
    db <- input$plantcv_sql_path$datapath
    drv <- dbDriver("SQLite")
    conn <- dbConnect(drv, dbname = db)
    removeNotification(id)
    
    id <- showNotification(h3("Querying db..."), duration = NULL)
    vis.df <- dbGetQuery(conn = conn, 'SELECT * FROM metadata NATURAL JOIN features WHERE imgtype = "VIS"')
    vis.df <- vis.df[,apply(vis.df[,seq(1, ncol(vis.df))], 2, function(x) unique(x)) != "0"]
    colnames(vis.df)[colnames(vis.df) == "plantbarcode"] <- "Barcodes"
    removeNotification(id)
    
    id <- showNotification(h3("Reading design file..."), duration = NULL)
    assoc <- read.csv(input$plantcv_design_file$datapath,header=T,stringsAsFactors = F)
    design$data <- assoc
    removeNotification(id)
    
    id <- showNotification(h3("Joining files..."), duration = NULL)
    sv_shapes <- join(vis.df,assoc,by="Barcodes")
    sv_shapes <- sv_shapes[rowSums(sapply(colnames(assoc),function(i) !is.na(sv_shapes[,i])))==ncol(assoc),]
    removeNotification(id)
    
    id <- showNotification(h3("Adding time columns..."), duration = NULL)
    sv_shapes$timestamp <- strptime(sv_shapes$timestamp,format = "%Y-%m-%d %H:%M:%S")
    beg <- min(sv_shapes$timestamp)
    sv_shapes$DAP <- floor(as.numeric((sv_shapes$timestamp - beg)/60/60/24))+2
    sv_shapes$hour <- lubridate::hour(sv_shapes$timestamp)
    removeNotification(id)
    
    id <- showNotification(h3("Removing empty pots..."), duration = NULL)
    sv_shapes <- sv_shapes[rowSums(sapply(colnames(assoc),function(i) !(sv_shapes[,i] %in% c("Blank","Empty","blank","empty"))))==ncol(assoc),]
    empties <- sv_shapes[sv_shapes$DAP == (max(sv_shapes$DAP)-1) & sv_shapes$area == 0,"Barcodes"]
    sv_shapes <- sv_shapes[!(sv_shapes$Barcodes %in% empties),]
    removeNotification(id)
    
    id <- showNotification(h3("Done!"), duration = 1)
    colnames(sv_shapes) <- gsub("-","_",colnames(sv_shapes))
    merged$data <- sv_shapes
    shapes$data <- merged$data[,c("image",'area','hull_area','solidity','perimeter','width','height','longest_axis','center_of_mass_x','center_of_mass_y','hull_vertices','in_bounds','ellipse_center_x','ellipse_center_y','ellipse_major_axis','ellipse_minor_axis','ellipse_angle','ellipse_eccentricity','y_position','height_above_bound','height_below_bound','above_bound_area','percent_above_bound_area','below_bound_area','percent_below_bound_area')]
    
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
        uiOutput("remove_outliers_ui")
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
  
  output$cooksd_plot <- renderPlot({
    if(!is.null(cooksd$data)){
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
    vis$data <- vis$data[cooksd$data < 3*mean(cooksd$data),]
    outliers <- merged$data[cooksd$data >= 3*mean(cooksd$data),]
    outliers$camera_angle <- unlist(lapply(strsplit(outliers$meta,"_"),function(i) i[3]))
    outliers$unique_id <- paste(outliers$Barcodes,outliers$DAP,outliers$camera_angle,sep="_")
    nir$camera_angle <- unlist(lapply(strsplit(as.character(nir$V1),"_"),function(i) i[3]))
    nir$data$unique_id <- paste(nir$data$Barcodes,nir$data$DAP,nir$data$camera_angle,sep="_")
    nir <- nir[!(nir$unique_id %in% outliers$unique_id),]
    removeNotification(id)
  })
  
  #***********************************************************************************************
  # ANOVA box
  #***********************************************************************************************
  output$shapes_anova <- renderUI({
    if(!is.null(merged$data)){
      box(width=10,title = "Shapes ANOVA",solidHeader = T,status = 'success',collapsible = TRUE,collapsed = TRUE,
        selectInput("which_day","Which Day",sort(unique(merged$data$DAP)),max(unique(merged$data$DAP))),
        actionButton("make_anova","Calculate ANOVA"),
        plotOutput("anova_plot")
      ) 
    }
  })
  
  anova_dat <- reactiveValues(data=NULL)
  observeEvent(input$make_anova,{
    id <- showNotification(h3("Calculating variances..."), duration = NULL)
    s <- colnames(shapes$data)[!(colnames(shapes$data) %in% c("meta","image","in_bounds"))]
    des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
    dat <- merged$data[merged$data$DAP==as.numeric(input$which_day),]
    H2 <- c()
    for(e in s){
      fmla <- as.formula(paste0("as.numeric(",e,")","~","(1|",des[1],")+(1|",des[2],")+(1|",des[1],":",des[2],")"))
      model <- lmer(fmla,data = dat)
      re<-as.numeric(VarCorr(model))
      res<-attr(VarCorr(model), "sc")^2
      interaction.var <- re[1]
      microbe.var<-re[2]
      drought.var<-re[3]
      tot.var<-sum(re,res)
      unexp <- 1-sum(re)/sum(re,res)
      
      h2 <- c((microbe.var/tot.var),
        (drought.var/tot.var),
        (interaction.var/tot.var),
        unexp)
      H2 <- rbind(H2,h2)
    }
    H2 <- data.frame(H2,row.names = s)
    H2$Shape <- rownames(H2)
    rownames(H2) <- NULL
    colnames(H2) <- c(des[1],des[2],"Interaction","Unexplained","Shape")
    H2$Shape <-  ordered(H2$Shape,levels=H2$Shape[order(H2$Unexplained)])
    H2_melt <- melt(H2,id=c("Shape"))
    H2_melt$variable <- ordered(H2_melt$variable,levels=c("Unexplained",des[2],des[1],"Interaction"))
    anova_dat$data <- H2_melt
    removeNotification(id)
  })
  
  output$anova_plot <- renderPlot({
    if(!is.null(anova_dat$data)){
      ggplot(data=anova_dat$data,aes(Shape,value*100))+
        geom_bar(stat = "identity",aes(fill=variable))+
        scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
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
    }
  })
  
  
  #***********************************************************************************************
  # Trends plots box
  #***********************************************************************************************
  output$trends_plot_ui <- renderUI({
    if(!is.null(merged$data)){
      s <- colnames(shapes$data)[!(colnames(shapes$data) %in% c("meta","image"))]
      des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
      box(width=10,title = " Growth Curves",solidHeader = T,status = 'success',collapsible = TRUE,collapsed = TRUE,
        selectInput("dep_var","Y-axis",s,"area"),
        selectInput("color_by","Color By",des,des[1]),
        selectInput("facet_by","Facet By",des,des[2]),
        plotOutput("trends_plot")
      ) 
    }
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
  
  
  #***********************************************************************************************
  # Growth heatmap box
  #***********************************************************************************************
  output$heatmap_plot_ui <- renderUI({
    if(!is.null(merged$data)){
      s <- colnames(shapes$data)[!(colnames(shapes$data) %in% c("meta","image"))]
      des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
      box(width=10,title = " Growth Heatmap",solidHeader = T,status = 'success',collapsible = TRUE,collapsed = TRUE,
        selectInput("h_color_by","Color By",s,"area"),
        selectInput("h_group_by","Group By",des,des[1]),
        selectInput("h_facet_by","Facet By",des,des[2]),
        plotOutput("trends_heatmap")
      ) 
    }
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
  
  #***********************************************************************************************
  # Color Helpers
  #***********************************************************************************************
  hist_avg <- function(data,start,stop){
    sub <- data
    test <- data.frame(do.call("rbind",lapply(split(sub,sub$Drought),function(t){
      data.frame(do.call("rbind",lapply(split(t,t$Microbes),function(g){
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
    test <- data.frame(do.call("rbind",lapply(split(sub,sub$Drought),function(t){
      data.frame(do.call("rbind",lapply(split(t,t$Microbes),function(g){
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
  # NIR Heatmap
  #***********************************************************************************************
  output$nir_heatmap_ui <- renderUI({
    if(!is.null(nir$data)){
      des <- colnames(design$data)[!(colnames(design$data) %in% "Barcodes")]
      box(width=10,title = "NIR Heatmap",solidHeader = T,status = 'success',collapsible = TRUE,collapsed = TRUE,
        plotOutput("nir_heatmap_nofacet"),
        plotOutput("nir_heatmap_withfacet")
      ) 
    }
  })
  
  output$nir_heatmap_nofacet <- renderPlot({
    test <- aggregate(data=nir$data[nir$data$intensityAVG != 0 & nir$data$DAP >8,],intensityAVG~Drought+Microbes+DAP,FUN = function(i)mean(i,na.rm=T))
    ggplot(test,aes(DAP,Drought))+
      #facet_grid(~Drought)+
      geom_tile(aes(fill=intensityAVG))+
      scale_fill_gradient2(limits=c(75,95),midpoint = 0.3+mean(nir$data$intensityAVG[nir$data$Drought == "AAA" & nir$data$intensityAVG != 0 & nir$data$DAP == 15]),high ="gray10",low= "#56B1F7",mid = "#d7e4ef")+
      theme_light()+
      theme(axis.text = element_text(size = 12),
        axis.title= element_text(size = 18))+
      theme(plot.title = element_text(hjust = 0.5),
        strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))
  })
  
  output$nir_heatmap_withfacet <- renderPlot({
    ggplot(nir$data[nir$data$intensityAVG != 0 & nir$data$DAP >8,],aes(DAP,Microbes))+
      facet_grid(~Drought)+
      geom_tile(aes(fill=intensityAVG))+
      scale_fill_gradient2(limits=c(75,95),midpoint = mean(nir$data$intensityAVG[nir$data$Drought == "AAA" & nir$data$intensityAVG != 0 & nir$data$DAP == 15]),high ="gray10",low= "#56B1F7",mid = "#d7e4ef")+
      theme_light()+
      theme(axis.text = element_text(size = 12),
        axis.title= element_text(size = 18))+
      theme(plot.title = element_text(hjust = 0.5),
        strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))
  })
  
}

shinyApp(ui, server)
