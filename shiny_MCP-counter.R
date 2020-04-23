library(shiny)
library(MCPcounter)
library(mMCPcounter)
library(DT)

# set max upload size to 50Mb
options(shiny.maxRequestSize=50*1024^2)

# set max space for java to 2Gb (needs space to load excel spreadsheets)
options(java.parameters = "-Xmx2048m")

# user interface: specify all input and output that are present in the app
ui <- navbarPage(title = "webMCP",
  
  # Several tabs: 1 to run MCP-counter, the others for downstream analyses
  tabPanel("Step 1: run (m)MCP-counter",
           
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
               #downloadButton(outputId = "downloadEstimates",label = "Download (m)MCP-counter estimates"),
               
               ## output: table of estimates
               dataTableOutput(outputId = "estimatesTable")
               
               
             )
             
           )
  ),
  
  
  tabPanel("Step 2: Downstream analyses",
           "Ongoing. Stay tuned!"),
  
  
  tabPanel("What is MCP-counter?",
           "Ongoing. Stay tuned!"),
  
  
  tabPanel("Citation",
           
           tags$p("If you use webMCP in a scientific publication, please cite the original article of the MCP-counter version you used (human or mouse), as well as:"),
           tags$p("CITATION INFO FOR THE SHINY APP"),
           tags$br(),
           tags$p("For the human MCP-counter: Becht, E., Giraldo, N.A., Lacroix, L. et al. ", tags$a(href="https://doi.org/10.1186/s13059-016-1070-5","Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression."),"Genome Biol 17, 218 (2016)."),
           tags$br(),
           tags$p("For the murine mMCP-counter: Petitprez, F., Lévy, S., Sun, C.-M., Meylan, M. et al. ", tags$a(href="https://doi.org/10.1101/2020.03.10.985176", "The murine Microenvironment Cell Population counter method to estimate abundance of tissue-infiltrating immune and stromal cell populations in murine samples using gene expression."), "bioRXiv (2020)")
           
           )
  
)


### Local functions

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
  if(length(grep("ENS",rownames(table))) > 50){return("ENSEMBL gene IDs detected. For now, this program only accepts HUGO Gene Symbols. Please provide the gene names as Gene Symbols. The support for ENSEMBL Gene IDs is planned in a future release.")}
  
  # Test if Human of Mouse based on the gene symbol case (Human: ALL CAPITAL, Mouse: First Letter Capital Only).
  firstLetterCapital = sum(toupper(substr(rownames(table),1,1))==substr(rownames(table),1,1))>0.8*nrow(table) # TRUE <=> at least 80% gene names start with capital letter
  if(!firstLetterCapital){return("Gene symbols have not been recognized. Please provide only human or murine data with genes identified through their Gene Symbols.")}
  allCapital = sum(toupper(rownames(table))==rownames(table))>0.8*nrow(table) # TRUE <=> at least 80% gene names are in full capital
  secondLetterCapital = sum(toupper(substr(rownames(table),2,2))==substr(rownames(table),2,2))>0.8*nrow(table) # TRUE <=> at least 80% gene names have a 2nd capital letter
  if(version=="h" & !(allCapital)){return("Gene symbols provided do not fit the human nomenclature. Please check the organism setting or provide human gene symbols.")}
  if(version=="m" & secondLetterCapital){return("Gene symbols provided do not fit the murine nomenclature. Please check the organism setting or provide murine gene symbols.")}
  
  
  return(paste("File loaded succesfully and correctly interpreted. ",c("h"="","m"="m")[version],"MCP-counter will be run shortly.",sep=""))
  
}








# server function: all computation using input
server <- function(input, output) {
  
  library(xlsx)
  
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
          estimates$est <- t(data.frame(MCPcounter.estimate(gep,featuresType = "HUGO_symbols"),check.names = FALSE))
        }
        else{
          estimates$est <- t(mMCPcounter.estimate(gep))
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
  
}
  


shinyApp(ui = ui, server = server)