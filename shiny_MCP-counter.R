library(shiny)
library(MCPcounter)
library(mMCPcounter)
library(DT)
library(ComplexHeatmap)
library(dendextend)
library(circlize)
library(ggplot2)
library(reshape2)
library(ggsignif)
library(gridExtra)
library(dunn.test)

######################
# TO DO
#
# - Add warning for missing populations
# - Have all plots below heatmap
# - refine format testing function:
#     - test separator
#     - ...
# - Finish "what is MCP-counter" tab
# - Manage ENSEMBL IDs (after initial review of mMCP-counter paper)
#
######################





# Florent: local test data 
#gep <- read.table("~/Documents/Ligue/Shiny MCP-counter/20200421_humanTestData.csv",sep=";",header = TRUE, row.names = 1,stringsAsFactors = FALSE, check.names = FALSE)
#gep <- read.table("~/Documents/Ligue/Shiny MCP-counter/AD_norm.txt",sep=",",header = TRUE, row.names = 1,stringsAsFactors = FALSE, check.names = FALSE)
#estimates <- list(est = NULL, version = "h")
# Maxime : same
#gep <- read.csv("/Users/meylanmaxime/Documents/85 échantillons Nivoren data RNAseq (1).csv",sep=",",header = TRUE, row.names = 1,stringsAsFactors = FALSE, check.names = FALSE)
#estimates <- list(est = NULL, version = "h")

# set max upload size to 50Mb
options(shiny.maxRequestSize=50*1024^2)

# set max space for java to 2Gb (needs space to load excel spreadsheets)
options(java.parameters = "-Xmx2048m")

# user interface: specify all input and output that are present in the app
ui <- navbarPage(title = "webMCP",
                 
                 
  # Several tabs: 1 to run MCP-counter, the others for downstream analyses
  tabPanel("Step 1: run (m)MCP-counter",
           
           
           
           # Modify CSS style for position of Notification             
           tags$head(
             tags$style(
               HTML(".shiny-notification {
             position:fixed;
             top: calc(50%);
             left: calc(50%);
             }
             "
               )
             )
           ),
           
           
           
           
           
           img(src="three_stickers.png", width=100,height=100, align="right"),
           
           tags$h1("Welcome to webMCP"),
           
           tags$h4("Web app for MCP-counter and mMCP-counter"),
           
           sidebarLayout(
             
             # The sidebar contains all the input the use has to provide
             sidebarPanel(
               
               ## file input: user-supplied gene expression profiles
               fileInput(inputId = "geneExpressionProfiles",label = "Upload gene expression profiles file",accept = c("xlsx","txt","csv")),
               
               ## file type select box: specifies the file type (xlsx or text format, in this case also specifies the separator)
               selectInput(inputId = "fileType",label = "Select the file format",choices = c("xlsx","tab-separated text file (txt, tsv, csv)","comma-separated text file (txt, tsv, csv)","semi-colon-separated text file (txt, tsv, csv)")),
               
               ## radio buttons for the organism
               radioButtons(inputId = "organism",label = "Select organism", choices = c("Human (Homo sapiens)","Mouse (Mus musculus)")),
               
               ## run button
               actionButton(inputId = "runButton",label = "Run (m)MCP-counter"),
               
               ## text output for the format diagnostic results to be displayed
               textOutput("diagnosticResult")
               
             ),
             
             # the main panel contains the output table and download button
             mainPanel(
               
               ## conditional output for the download button
               uiOutput("downloadButton"),
               
               ## output: table of estimates
               dataTableOutput(outputId = "estimatesTable")
               
               
             )
             
           )
  ),
  
  
  tabPanel("Step 2: Downstream analyses - Clustering",
           
           sidebarLayout(
             
             sidebarPanel(
               
               # Slider for the number of clusters
               sliderInput("nrCluster",label = "Select number of clusters to display",min = 2, max = 10,value = 3),
               
               # Check boxes for choice of populations to include
               uiOutput("popCheckBoxes"),
               
               # All/non checkbox
               checkboxInput("allNonePop","All/None selected",value=TRUE),
               
               # Action button that runs the clustering analysis
               actionButton("runClustering","Update"),
               
               # Text output to print any necessary messages
               textOutput("clusteringMessage")
               
             ),
             
              mainPanel(
                
                plotOutput("heatmap"),
                
                ## conditional output for the download button
                uiOutput("downloadClustersButton"),
                
                
                
                plotOutput("MCPboxplots")
                
              )
             
           )
           ),
  
  
  tabPanel("What is MCP-counter?",
           
           tags$h2("What is MCP-counter"),
           tags$p("In this setion, we will rapidly describe what MCP-counter is. Do not hesitate to",tags$a(href="mailto:florent.petitprez@ligue-cancer.net","get in touch"),"if you have more questions."),
           tags$br(),
           
           tags$h3("Objective of the method"),
           tags$p("MCP-counter, and its murine counterpart method mMCP-counter, are deconvolution methods that use transcriptomic data to estimate the relative abundance of diverse immune and stromal population in heterogeneous bulk samples. The following table lists the populations that can be assessed using MCP-counter or mMCP-counter, depending on the organism."),
           tableOutput("populationsTable"),
           tags$br(),
           
           tags$h3("How does it work?"),
           tags$br(),
           tags$p("The global idea of MCP-counter is to seek genes that are expressed in one cell population (and all its sub-populations), and not expressed by all other cell types.",
                  "These genes are called",tags$i("transcriptomic markers."),"The plot below illustrates the expression pattern of one such transcriptomic marker for murine mast cells."),
           img(src="mMCPcounter_specificMarker.png", width=1000),
           tags$p(style="color:grey","Source: Petitprez et al., BioRXiv, 2020"),
           tags$p("The determination of a gene as a transcriptomic marker is based on three criteria: specific variablility (i.e. overexpression in the population of interest as compared to the expression in all other populations), aspecific variability (i.e. overexpression in the populations of interest as compared with the variability within all other populations), and sensitivity-specificity as measured with the area under the ROC curve."),
           tags$p("These selection criteria ensure sufficient specificity of the considered gene signatures. Therefore, the expression level of these genes is proportional the the abundance of the designated cell populations in the bulk samples."),
           
           tags$h3("How to interpret the scores?"),
           img(src="compCIBERSORT.png", width = 800),
           tags$p(style="color:grey","Source: Petitprez et al., Cancer Immunology Immunotherapy, 2017"),
           
           
           
           ),
  
  
  tabPanel("Citation",
           
           tags$p("If you use webMCP in a scientific publication, please cite the original article of the MCP-counter version you used (human or mouse), as well as:"),
           tags$p("CITATION INFO FOR THE SHINY APP"),
           tags$br(),
           tags$p("For the human MCP-counter: Becht, E., Giraldo, N.A., Lacroix, L. et al. ", tags$a(href="https://doi.org/10.1186/s13059-016-1070-5","Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression."),"Genome Biol 17, 218 (2016)."),
           tags$br(),
           tags$p("For the murine mMCP-counter: Petitprez, F., Lévy, S., Sun, C.-M., Meylan, M. et al. ", tags$a(href="https://doi.org/10.1101/2020.03.10.985176", "The murine Microenvironment Cell Population counter method to estimate abundance of tissue-infiltrating immune and stromal cell populations in murine samples using gene expression."), "bioRXiv (2020)"),
           tags$br(),
           tags$br(),
           tags$p(style="color:grey","webMCP was developped by",tags$a(href="https://florentpetitprez.netlify.app/","Florent Petitprez"),"and Maxime Meylan"),
           
           
           
           )
  
)





# server function: all computation using input
server <- function(input, output) {
  
  library(xlsx)
  
  
  
  ##############################
  # Tab 1: Running MCP-counter #
  ##############################
  
  ## static table for the "What is MCP-counter" tab
  output$populationsTable <- renderTable(data.frame('MCP-counter (human)' = c("T cells", "CD8+ T cells", "Cytotoxic lymphocytes", "NK cells", "B lineage", "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils", "Endothelial cells", "Fibroblasts","","","","","",""),
                                                    'mMCP-counter (mouse)' = c("T cells", "CD8+ T cells", "NK cells", "B-derived cells", "Memory B cells", "Monocytes/macrophages", "Monocytes", "Basophils", "Mast cells", "Eosinophils", "Granulocytes", "Neutrophils", "Vessels", "Endothelial cells", "Lymphatics","Fibroblasts"),check.names = FALSE))
  
  
  
  # initialize the estimates as a reactive value. This will be invalidated each time (m)MCP-counter in run.
  estimates <- reactiveValues(est = NULL, version = "h")
  
  # code to be run when the user clicks the "run (m)MCP-counter" button: read file path, determine format to read data, and run the appopriate version of MCP-counter depending on organism.
  observeEvent(input$runButton,{
    
    # get (tmp) path to user-provided file
    userFile <- (input$geneExpressionProfiles)$datapath
    
    
    # proceed only if file extension matches the declared file format
    if(formatFits(userFile,input$fileType)){
      
      
      # file manipulations needed if input provided as excel spreadsheet
      if(input$fileType=="xlsx"){
        gep <- read.xlsx2(userFile,sheetIndex = 1,check.names=FALSE,stringsAsFactors=FALSE,colClasses = c("character",rep("numeric",10000)))
        geneNames <- gep[,1]
        gep <- gep[,2:ncol(gep)]
        rownames(gep) <- geneNames
      }
      else{ # text format, only separator changes
        sepText <- substr(input$fileType,1,3)
        sep <- c("tab" = "\t", "com" = ",", "sem" = ";")[sepText]
        gep <- read.table(userFile,sep=sep,header = TRUE, row.names = 1,stringsAsFactors = FALSE, check.names = FALSE)
      }
      
      if(input$organism=="Mouse (Mus musculus)"){estimates$version <- "m"}
      
      diagnostic <- formatDiagnostic(gep,estimates$version)
      
      output$diagnosticResult <- renderText(diagnostic)
      
      # try to run MCP-counter only if file fits format requirements
      if(substr(diagnostic,1,50)=="File loaded succesfully and correctly interpreted."){
        
        # run the appropriate version of MCP-counter depending on the organism 
        if(estimates$version=="h"){
          MCPcountercall <- catchToList(MCPcounter.estimate(gep,featuresType = "HUGO_symbols"))
          estimates$est <- t(data.frame(MCPcountercall$value,check.names = FALSE))
          if(!is.null(MCPcountercall$warnings)){
            showNotification(MCPcountercall$warnings[1],type="error",duration = NULL)
          }
        }
        else{
          MCPcountercall <- catchToList(mMCPcounter.estimate(gep))
          estimates$est <- t(MCPcountercall$value)
          if(!is.null(MCPcountercall$warnings)){
            showNotification(MCPcountercall$warnings[1],type="error",duration = NULL)
          }
        }
        
      }
      
    }
    else{ # file extension does not match declared file format
      output$diagnosticResult <- renderText("File extension does not match the file format you declared. Please correct and try again.")
    }
    
   
  })
  
  
  output$downloadButton <- renderUI({
    if(is.null(estimates$est))return()
    else{downloadButton(outputId = "downloadEstimates",label = paste("Download ",c("h"="","m"="m")[estimates$version],"MCP-counter estimates",sep=""))}
  })
  
  
  # output the rounded estimates in a data table
  output$estimatesTable <- DT::renderDataTable(local.round(estimates$est, digits = 2),rownames = TRUE)
  
  # handling the download output button action
  output$downloadEstimates <- downloadHandler(filename = function(){"MCP-counterEstimates.tsv"},
                                              content = function(file){
                                                write.table(formatForExport(estimates$est,"Population"),file,row.names=FALSE,col.names=TRUE, sep = "\t")
                                                })
  
  
  
  ##############################
  # Tab 2: clustering analysis #
  ##############################
  
  clusterColorCode <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","#999999","#8dd3c7")
  names(clusterColorCode) <- paste("Cluster",1:10)
  
  #initialize variables
  clusters <- NULL
  populationsToUse <- reactiveValues(populations = NULL)
  
  
  output$popCheckBoxes <- renderUI(checkboxGroupInput("popChoice",label = "Check populations to be included in the analysis",
                                                      choices=colnames(estimates$est),selected = if(input$allNonePop)colnames(estimates$est)))
  
  
  # wait for the run button to be clicked
  observeEvent(input$runClustering,{
    
    
    
    if(is.null(estimates$est)){output$clusteringMessage <- "The MCP-counter estimates could not be found. Please ensure to first run the analysis of step 1."}
    else{
      
      populationsToUse$populations <- input$popChoice
      
      print(populationsToUse)
      
      # scale estimates for the heatmap
      est.norm <- t(apply(t(estimates$est[,populationsToUse$populations]),1,scale))
      colnames(est.norm) <- row.names(estimates$est[,populationsToUse$populations])
      
      # rows with sd=0 are set to NAs by the scale function. This chunk sets them to 0.
      est.norm[is.na(est.norm[,1]),] <- 0
      
      
      #dendrograms for the populations (mcp) and for the samples
      dend.mcp <- hclust(dist(est.norm,method = "euclidian"))
      dend.samples <- hclust(dist(t(est.norm),method = "euclidian"))
      dend.samples <- color_branches(dend.samples,k=input$nrCluster,col = clusterColorCode[1:input$nrCluster])
      
      # complexHeatmap definition of the heatmap
      hm <- Heatmap(as.matrix(est.norm), col = colorRamp2(c(-4,-2, 0, 2,4),c("#2166ac","#92c5de","#f7f7f7","#f4a582","#b2182b")),
                    cluster_rows = dend.mcp, cluster_columns = dend.samples, show_column_names = TRUE, show_row_names = TRUE, name = "Row Z-score")
      
      # plot the heatmap
      output$heatmap <- renderPlot(draw(hm))
      
      # get clusters
      clusters <- paste("Cluster",match(get_leaves_attr(dend.samples,"edgePar"),clusterColorCode))
      clusters_n <- paste(clusters,"\n(n=",table(clusters)[clusters],")",sep="")
      names(clusters_n) <- labels(dend.samples)
      clusters_n <- clusters_n[rownames(estimates$est)]
      
      # Get download button
      output$downloadClustersButton <- renderUI(downloadButton(outputId = "downloadClusters",label ="Download clusters information"))
      
      # handling the download output button action
      output$downloadClusters <- downloadHandler(filename = function(){"MCP-counterClusters.tsv"},
                                                  content = function(file){
                                                    write.table(data.frame(ID = colnames(est.norm),cluster = clusters),file,row.names=FALSE,col.names=TRUE, sep = "\t")
                                                  })
      
      
      
      # get boxplots 
      estimates_df <- data.frame(estimates$est[,populationsToUse$populations],check.names = F)
      estimates_df$clusters <- clusters_n
      melted_est <- melt_df(estimates_df,var_to_group = 'clusters')
      
      # MCPboxplots is a matrix that can be queried like the following
      # MCPboxplots["plot","T cells"]
      # MCPboxplots["plot",1:4]
      
      MCPboxplots <- sapply(colnames(estimates$est[,populationsToUse$populations]),function(x) plot_group_boxplot(data.m = melted_est,
                                                                                  variable=x,
                                                                                  violin = (min(table(clusters))>2), #get violin plot only if all clusters have at least 3 samples, otherwise it messes with the color code. If one cluster has 2 or less samples, simply plot boxplot
                                                                                  specify_col =  as.character(clusterColorCode[1:input$nrCluster]),
                                                                                  labs=c("","","MCP score") 
                                                                                ))
      
      output$MCPboxplots <- renderPlot({
        grid.arrange(grobs=MCPboxplots["plot",1:3],ncol=3)
      })
      
      # output$MCPboxplots <- renderPlot({
      #   par(mfrow = c(1,2))
      #   boxplot(estimates$est[,"T cells"]~clusters_n,xlab = "", ylab = "T cells", border = clusterColorCode,outline=FALSE, las=2)
      #   stripchart(estimates$est[,"T cells"]~clusters_n,vertical = TRUE, add=TRUE, method = "jitter",pch=16, col = clusterColorCode)
      #   boxplot(estimates$est[,"NK cells"]~clusters_n,xlab = "", ylab = "NK cells", border = clusterColorCode,outline=FALSE, las=2)
      #   stripchart(estimates$est[,"NK cells"]~clusters_n,vertical = TRUE, add=TRUE, method = "jitter",pch=16, col = clusterColorCode)
      # })
      
      
    }
    
    
    
    
  })
  
  
}
  


###################
# Local functions #
###################


# adapation of round: does nothing if applied to NULL. This is meant to avoid the error when no user data has been provided yet
local.round <- function(x,digits=2){
  if(!is.null(x))return(round(x,digits = digits))
}

# format a table for export. Adds the rownames as a first column with header specified as 2nd argument
formatForExport <- function(table,rowNameHeader){
  col1 <- rownames(table)
  names <- c(rowNameHeader,colnames(table))
  res <- cbind(col1,table)
  colnames(res) <- names
  return(res)
}


# formatFits test if there is an agreement between the file extension and the file format
formatFits <- function(filePath,fileFormat){
  
  extension = substr(filePath,nchar(filePath)-2, nchar(filePath))
  
  return(
    ((extension %in% c("txt","csv","tsv")) & (fileFormat !="xlsx")) |
      ((extension == "lsx") & (fileFormat =="xlsx"))
  )
  
}


# formatDiagnostic tests whether the table provided by the user fits the format required and returns a char that states if the format is correct or how to correct it.
formatDiagnostic <- function(table, version = c("h","m")[1]){
  
  # Test that table is not NULL
  if(is.null(table)){return("The file you provided could not be interpreted as a table. Please try again with an other file.")}
  
  # Test whether table is a data frame
  if(!is.data.frame(table)){return("The file you provided could not be transformed into a data set. Please pay attention to the formating guidelines and try again.")}
  
  # Test if all columns are numeric
  if(!all(apply(table,2,is.numeric))){return("Non-numeric columns have been detected. Please provide a numeric-only table, with the exception of sample and gene names.")}
  
  # Test if ENSEMBL gene IDs (to be supported in a future version)
  #if(length(grep("ENS",rownames(table))) > 50){return("ENSEMBL gene IDs detected. For now, this program only accepts HUGO Gene Symbols. Please provide the gene names as Gene Symbols. The support for ENSEMBL Gene IDs is planned in a future release.")}
  
  # Test if Human of Mouse based on the gene symbol case (Human: ALL CAPITAL, Mouse: First Letter Capital Only).
  firstLetterCapital = sum(toupper(substr(rownames(table),1,1))==substr(rownames(table),1,1))>0.8*nrow(table) # TRUE <=> at least 80% gene names start with capital letter
  if(!firstLetterCapital){return("Gene symbols have not been recognized. Please provide only human or murine data with genes identified through their Gene Symbols.")}
  allCapital = sum(toupper(rownames(table))==rownames(table))>0.8*nrow(table) # TRUE <=> at least 80% gene names are in full capital
  secondLetterCapital = sum(toupper(substr(rownames(table),2,2))==substr(rownames(table),2,2))>0.8*nrow(table) # TRUE <=> at least 80% gene names have a 2nd capital letter
  if(version=="h" & !(allCapital)){return("Gene symbols provided do not fit the human nomenclature. Please check the organism setting or provide human gene symbols.")}
  if(version=="m" & secondLetterCapital){return("Gene symbols provided do not fit the murine nomenclature. Please check the organism setting or provide murine gene symbols.")}
  
  return(paste("File loaded succesfully and correctly interpreted. ",c("h"="","m"="m")[version],"MCP-counter will be run shortly.",sep=""))
  
}

# boxplot function with appropriate tests and display of p-values
plot_group_boxplot <- function(data.m,
                               variable,
                               plot_only_signif=FALSE,
                               pval_thresh=0.05,
                               specify_col=NA,
                               palette="YlOrRd",
                               title="",
                               title_size=NA,
                               show_paired=F,
                               add_jitter=T,
                               color_points_id=NULL,
                               log_scale=F,
                               violin=F,
                               plot_outlier=T,
                               compare_groups=T,
                               add_hline=NA,
                               hide_test =F,
                               paired=F,
                               do_test=T,
                               ylim=NA,
                               labs=c("title","value","value")){
  #initialize variables
  return_plot <- FALSE
  dodge <- position_dodge(width=0.9)
  p <- NULL
  ylim1 <- NULL
  row_id <- which(variable == data.m[,"variable"])
  id_per_group <- lapply(unique(data.m$groups),function(x) which(data.m[row_id,"groups"] == x))
  names(id_per_group) <- unique(data.m$groups)
  data_to_sub <- data.m[row_id,]
  value_per_group <- lapply(id_per_group,function(x) data_to_sub[x,"value"])
  
  #Perfom global significance test
  if(do_test){class
    if(length(value_per_group) == 2 & !paired){
      test = wilcox.test(value_per_group[[1]],value_per_group[[2]])
      test_name <- "Mann-Whitney test"
    }
    if(length(value_per_group) > 2) {
      test <- kruskal.test(value_per_group)
      test_name <-"Kruskal-Wallis test"
    }
    if(paired){
      #add Ids for strips
      data.m$IDs <- c(id_per_group[[1]],id_per_group[[1]])
      test = wilcox.test(value_per_group[[1]],value_per_group[[2]],paired=T)
      test_name <- "Paired Mann-Whitney test"
    }
    pval <- ifelse(is.na(test$p.value),1,test$p.value)
  }else{
    pval <- 1
    test_name <- "No test"
  }
  
  if(plot_only_signif & pval <= pval_thresh) return_plot <-  TRUE
  if(plot_only_signif == FALSE ) return_plot <- TRUE
  if(return_plot){
    p <- ggplot(data = data.m[row_id,], aes(x=groups, y=value))
    p <- p + theme_bw()
    if(!is.na(title_size)){
      p <- p +   theme(plot.title = element_text(size = title_size))
    }
    if(all(!is.na(ylim))){
      p <- p + ylim(ylim)
    }
    if(any(is.na(specify_col))){
      p <- p + scale_fill_brewer(palette = palette)
    }else{
      p <- p + scale_fill_manual(values = specify_col)
    }
    if(violin){
      p <- p + geom_violin(aes(fill=groups), position = dodge,trim=T) + geom_boxplot(outlier.alpha = 0.003,width=0.1, fill="white")
    }else{
      p <- p + geom_boxplot(aes(fill=groups), position = dodge,outlier.shape=NA)
    }
    #display strips for paired data
    if(show_paired){
      p <- p + geom_line(aes(group=IDs),colour=alpha("grey50", alpha = 0.2),linetype="11")
      p <- p + geom_point(data=data.m[row_id,],aes(x =groups, y =value),alpha=0.03, colour ="black",show.legend = F)
    }
    if(add_jitter){
      p <- p + geom_jitter(position = position_jitter(width = .05, height = 0), alpha = 0.4)
    }
    if(!(is.na(add_hline))){
      p <- p + geom_hline(yintercept=add_hline, linetype="dashed", color = "red")
    }
    #display pvalue without test name
    if(hide_test){
      if(pval < 1e-16){
        p <- p + labs(subtitle=paste(" p < 10e-16",sep=" "),title=paste(variable,labs[1]),x=labs[2],y=labs[3]) +  guides(fill=FALSE)
        p <- p + theme(plot.subtitle=element_text(face="bold"))
      }else{
        p <- p + labs(subtitle=paste(" p = ",format(pval,digits = 3),sep=" "),title=paste(variable,labs[1]),x=labs[2],y=labs[3]) +  guides(fill=FALSE)
        if(pval < 0.05 ) p <- p + theme(plot.subtitle=element_text(face="bold"))
      }
    #display pvalue and test name
    }else{
      if(pval < 1e-16){
        p <- p + labs(subtitle=paste(test_name, " p < 10e-16",sep=" "),title=paste(variable,labs[1]),x=labs[2],y=labs[3]) +  guides(fill=FALSE)
        p <- p + theme(plot.subtitle=element_text(face="bold"))
      }else{
        p <- p + labs(subtitle=paste(test_name, " p = ",format(pval,digits = 3),sep=" "),title=paste(variable,labs[1]),x=labs[2],y=labs[3]) +  guides(fill=FALSE)
        if(pval < 0.05 ) p <- p + theme(plot.subtitle=element_text(face="bold"))
      }
    }
    p <- p + theme(plot.subtitle = element_text(size = 9))
    #color specified points
    if(!is.null(color_points_id)){
      p <- p + geom_point(data=data.m[num_ids,],aes(x =groups, y =value, colour ="red"))
    }
    if(log_scale){
      p <- p + scale_y_continuous(trans = "log10",breaks = base_breaks(), labels = prettyNum) + theme(panel.grid.minor = element_blank())
    }
    #if plot_outlier set to F, hide bottom 2% and top 5% values
    if(!plot_outlier){
      ylim1[1] <- quantile(data.m[row_id,"value"],0.02,na.rm=T)
      ylim1[2] <- quantile(data.m[row_id,"value"],0.95,na.rm=T)
      p <- p + coord_cartesian(ylim =ylim1)
    }
    #perform dunn.test to evaluate pairwize difference for group combination
    if(compare_groups & length(value_per_group) > 2 & sd(unlist(value_per_group))>0){
      #remove cat from dunn.test
      dunn_test <- dunn.test(x=data.m[row_id,"value"],g=data.m[row_id,'groups'],method = "bh")
      comparisons <- sapply(dunn_test$comparisons,function(x) strsplit(x,split =" - "))
      adjp_dunn <- format(dunn_test$`P.adjusted`,digits = 2,trim = T)
      to_bold <- ifelse(as.numeric(adjp_dunn) < 0.05, paste( "bold('" , adjp_dunn,"')"),adjp_dunn)
      p <- p + geom_signif(annotations = to_bold,parse = T, comparisons = comparisons, step_increase = 0.1)
      
      names(adjp_dunn) <- dunn_test$comparisons
      res <- list(plot=p,main_pval=pval,main_test=test_name,comparisons_pval=adjp_dunn)
      return(res)
    }
  }
  res <- list(plot=p,main_pval=pval,main_test=test_name,comparisons_pval=NA)
  return(res)
}

#Get nice axis ticks
#from Heather Turner @stackoverflow
base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = T, n = n)
  }
}
#Small edit to melt function
melt_df <- function(df,var_to_group){
  df.m <- melt(df,id = var_to_group,varnames = c("variable","value"))
  colnames(df.m)[1] <- "groups"
  df.m$value <- as.numeric(as.character(df.m$value))
  return(df.m)
}

# Function that allows to catch both value of expression and any warnings or errors
# from user "Aaron left Stack Overflow" @stackoverflow
catchToList <- function(expr) {
  val <- NULL
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, w$message)
    invokeRestart("muffleWarning")
  }
  myError <- NULL
  eHandler <- function(e) {
    myError <<- e$message
    NULL
  }
  val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
  list(value = val, warnings = myWarnings, error=myError)
} 




shinyApp(ui = ui, server = server)
