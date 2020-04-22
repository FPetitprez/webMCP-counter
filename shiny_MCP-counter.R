library(shiny)
library(MCPcounter)
library(mMCPcounter)
library(DT)

# set max upload size to 50Mb
options(shiny.maxRequestSize=50*1024^2)

# set max space for java to 2Gb (needs space to load excel spreadsheets)
options(java.parameters = "-Xmx2048m")

# user interface: specify all input and output that are present in the app
ui <- fluidPage(
  
  "Welcome to the (m)MCP-counter server",
  
  ## file input: user-supplied gene expression profiles
  fileInput(inputId = "geneExpressionProfiles",label = "Upload gene expression profiles file",accept = c("xlsx","txt","csv")),
  
  ## file type select box: specifies the file type (xlsx or text format, in this case also specifies the separator)
  selectInput(inputId = "fileType",label = "Select the file format",choices = c("xlsx","tab-separated text file (txt, tsv, csv)","comma-separated text file (txt, tsv, csv)","semi-colon-separated text file (txt, tsv, csv)")),
  
  ## radio buttons for the organism
  radioButtons(inputId = "organism",label = "Select organism", choices = c("Human (Homo sapiens)","Mouse (Mus musculus)")),
  
  ## run button
  actionButton(inputId = "runButton",label = "Run (m)MCP-counter"),
  
  ## output: table of estimates
  dataTableOutput(outputId = "estimatesTable"),
  
  ## ask what format to use 
  
  ## download output button
  downloadButton(outputId = "downloadEstimates",label = "Download (m)MCP-counter estimates")
  
)

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

# server function: all computation using input
server <- function(input, output) {
  
  library(xlsx)
  
  # initialize the estimates as a reactive value. This will be invalidated each time (m)MCP-counter in run.
  estimates <- reactiveValues(est = NULL)
  
  # code to be run when the user clicks the "run (m)MCP-counter" button: read file path, determine format to read data, and run the appopriate version of MCP-counter depending on organism.
  observeEvent(input$runButton,{
    
    # get (tmp) path to user-provided file
    userFile <- (input$geneExpressionProfiles)$datapath
    
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
   
   # run the appropriate version of MCP-counter depending on the organism 
   if(input$organism=="Human (Homo sapiens)"){
      estimates$est <- t(data.frame(MCPcounter.estimate(gep,featuresType = "HUGO_symbols"),check.names = FALSE))
    }
    else{
      estimates$est <- t(mMCPcounter.estimate(gep))
    }
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