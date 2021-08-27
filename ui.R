library(shiny)
library(DT)


# set max upload size to 50Mb
options(shiny.maxRequestSize=50*1024^2)

# set max space for java to 2Gb (needs space to load excel spreadsheets)
options(java.parameters = "-Xmx2048m")

# user interface: specify all input and output that are present in the app
ui <- navbarPage(title = "webMCP-counter",
                 
                 
                 # Tracking code
                 tags$head(HTML(
                   "<script type='text/javascript'>
                     var _paq = window._paq = window._paq || [];
                     _paq.push(['trackPageView']);
                   _paq.push(['enableLinkTracking']);
                   (function() {
                     var u='http://134.157.229.105:3838/webMCP/';
                     _paq.push(['setTrackerUrl', u+'matomo.php']);
                     _paq.push(['setSiteId', '1']);
                     var d=document, g=d.createElement('script'), s=d.getElementsByTagName('script')[0];
                     g.type='text/javascript'; g.async=true; g.src=u+'matomo.js'; s.parentNode.insertBefore(g,s);
                   })();
                   </script>"
                 )),
                 
                 
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
                          
                          tags$h1("Welcome to webMCP-counter"),
                          
                          tags$h4("Web app for MCP-counter and mMCP-counter"),
                          
                          tags$p("To run MCP-counter or mMCP-counter using this web interface, you first need to prepare your data in a suitable format. 4 formats are accepted: Excel spreadsheet or text-based with tab, comma or semi-colon separator. Text-format are preferred as they are more memory-efficient. In all cases, the samples must be put in columns, and genes in rows. The first column must imperatively be composed of gene symbols or gene ENSEMBL IDs. The first line must be composed of the corresponding sample IDs, and the cell above gene symbols must be filled. If you have further questions or require assistance, please open a new issue on the",tags$a(href="https://github.com/FPetitprez/webMCP-counter/issues/new/choose","webMCP-counter Github page.",target="_blank")),
                          
                          
                          sidebarLayout(
                            
                            # The sidebar contains all the input the use has to provide
                            sidebarPanel(
                              
                              ## file input: user-supplied gene expression profiles
                              fileInput(inputId = "geneExpressionProfiles",label = "Upload gene expression profiles file",accept = c("xlsx","txt","csv")),
                              tags$p("or, alternatively,"),
                              checkboxInput("TCGA_SARC","try with the TCGA SARC example data:"),
                              tags$br(),
                              
                              ## file type select box: specifies the file type (xlsx or text format, in this case also specifies the separator)
                              selectInput(inputId = "fileType",label = "Select the file format",choices = c("xlsx","tab-separated text file (txt, tsv, csv)","comma-separated text file (txt, tsv, csv)","semi-colon-separated text file (txt, tsv, csv)")),
                              
                              ## radio buttons for the gene IDs
                              radioButtons(inputId = "geneIDs",label = "Select gene ID format", choices = c("Gene symbol","ENSEMBL")),
                              
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
                          tags$p("In this section, we will rapidly describe what MCP-counter is. If you have further questions or require assistance, please open a new issue on the",tags$a(href="https://github.com/FPetitprez/webMCP-counter/issues/new/choose","webMCP-counter Github page.",target="_blank")),
                          tags$br(),
                          
                          tags$h3("Objective of the method"),
                          tags$p("MCP-counter, and its murine counterpart method mMCP-counter, are deconvolution methods that use transcriptomic data to estimate the relative abundance of diverse immune and stromal population in heterogeneous bulk samples. The following table lists the populations that can be assessed using MCP-counter or mMCP-counter, depending on the organism."),
                          tableOutput("populationsTable"),
                          tags$br(),
                          
                          tags$h3("How does it work?"),
                          tags$p("The global idea of MCP-counter is to seek genes that are expressed in one cell population (and all its sub-populations), and not expressed by all other cell types.",
                                 "These genes are called",tags$i("transcriptomic markers."),"The plot below illustrates the expression pattern of one such transcriptomic marker for murine mast cells."),
                          img(src="mMCPcounter_specificMarker.png", width=1000),
                          tags$p(style="color:grey","Source: Petitprez et al., Genome Medicine, 2020"),
                          tags$p("The determination of a gene as a transcriptomic marker is based on three criteria: specific variablility (i.e. overexpression in the population of interest as compared to the expression in all other populations), aspecific variability (i.e. overexpression in the populations of interest as compared with the variability within all other populations), and sensitivity-specificity as measured with the area under the ROC curve. These selection criteria ensure sufficient specificity of the considered gene signatures. Therefore, the expression level of these genes is proportional the the abundance of the designated cell populations in the bulk samples."),
                          tags$br(),
                          
                          tags$h3("How to interpret the scores?"),
                          tags$p("The scores returned by MCP-counter or mMCP-counter are expressed in arbitrary units. These scores are proportional to the amount of the estimated cell populations in the sample. Each population having a different arbitrary unit. Therefore, it cannot be used to compare the abundance of different populations within one sample. However, these scores allow the comparison of the abundance of one cell population between samples in a cohort. This is a fundamental difference with other deconvolution methods such as CIBERSORT which estimate the relative composition within the overall immune infiltrate, and therefore allow to compare between populations within a sample, but not between samples. For more details about these difference and a benchmark of different deconvolution methods to estimate immune and stromal sample compositon, you can refer to",tags$a(href="https://doi.org/10.1093/bioinformatics/btz363", "Sturm et al., Bioinformatics, 2019")),
                          tags$p("The plot below illustrates this fundamental difference of approach. The left anel is a schematic representation of three possible cell mixtures, while the middle and right panels represent, respectively, the estimates that would be suggested by CIBERSORT and MCP-counter. We notice that the estimates of CIBERSORT for the first two mixes are similar, as they are expressed as percentages of cells among the screened populations only, regardless of the total infiltration in the sample. Conversely, MCP-counter scores are proportional to the amount of each cell population in the total sample, which allows inter-sample comparison for each population. However, these scores are expressed in a different arbitrary unit for each population, which prevents intra-sample comparison between populations. CIBERSORT allows this type of comparison"),
                          img(src="compCIBERSORT.png", width = 800),
                          tags$p(style="color:grey","Source: Petitprez et al., Cancer Immunology Immunotherapy, 2017")
                          
                          
                          
                 ),
                 
                 
                 tabPanel("Help for webMCP-counter",
                          
                          tags$p("This section describes the main steps to run webMCP-counter. If you need more assistance or wish to suggest an improvement to webMCP-counter, please ask open a new issue on the",tags$a(href="https://github.com/FPetitprez/webMCP-counter/issues/new/choose","webMCP-counter Github page.",target="_blank")),
                          
                          tags$h2("Prepare and load gene expression data"),
                          tags$p("To run MCP-counter or mMCP-counter using this web interface, you first need to prepare your data in a suitable format. 4 formats are accepted: Excel spreadsheet or text-based with tab, comma or semi-colon separator. Text-format are preferred as they are more memory-efficient. In all cases, the following formatting reules must be followed:"),
                          tags$p("- Samples must be put in columns, and genes in rows."),
                          tags$p("- The first column must imperatively be composed of gene symbols or gene ENSEMBL IDs."),
                          tags$p("- The first line must be composed of the corresponding sample IDs, and the cell above gene symbols must be filled."),
                          tags$p("Once your data is ready in a correct format, you can upload it using the upload button on the step 1 tab. Then you need to specify the format you have chosen (Excel spreadsheet or text file, in this case you need to specify the separator), the organism of origin (human or mouse) and the format for gene IDs (Gene symbol or ENSEMBL ID). Then simply click the `run (m)MCP-counter` button."),
                          tags$p("In case the file format is not correct, a text will appear below the button telling you what went wrong. If everything is correct, MCP-counter will be run on your data."),
                          tags$p("Alternatively, if you simply wish to see what webMCP-counter can do, you can select the TCGA SARC (Soft-tissue Sarcoma from The Cancer Genome Atlas project) dataset that is provided as an example. Simply tick the box and the run button. All parameters will be automatically set and MCP-counter will be run on the TCGA SARC dataset."),
                          tags$br(),
                          
                          tags$h2("Obtain the (m)MCP-counter scores"),
                          tags$p("Once MCP-counter is run, a table will appear on the right side of the step 1 tab. You can browse through it manually or use the search bar to find something of interest to you. A download button will also be available to download the full table. The download will be a .txt file, with columns separated by tabulations. This file format can be loaded and read with Microsoft Excel or LibreOffice Calc for instance."),
                          tags$br(),
                          
                          tags$h2("Downstream analysis"),
                          tags$p("After you have run (m)MCP-counter on your data, you can use the step 2 tab to conduct basic downstream analysis. In this tab, you will be able to plot a heatmap representing the scores on your data for the populations that are of interest to you. You will also be able to cluster your samples in as many groups as you wish (between 2 and 10), and visualize the differences between cell types abundances between the clusters."),
                          tags$p("To run this downstream analysis, specify the number of clusters you want, and select the populations you wish to be used for the clustering. Then click the `update` button. A heatmap appears on the right, representing normalized MCP-counter scores of the requested populations on your samples. The dendrogram (classification tree) above the heatmap is colored according to the clusters. Below the heatmap, violin plots or boxplots (if one cluster has 2 or less samples) are displayed for all requested populations, showing the differences between the clusters. Kruskal-Wallis tests, and pairwise comparisons using Dunn test are also performed to estimate significance of the inter-cluster differences. A download button also offers you to obtain the cluster to which all your samples belong."),
                          tags$p("If you wish to modify the number of clusters and/or the included populations, simply adapt the settings on the left panel and click on the update button."),
                          tags$br(),
                          
                          
                          tags$h2("Need more help?"),
                          tags$p("If you need more assistance or wish to suggest an improvement to webMCP-counter, please ask open a new issue on the",tags$a(href="https://github.com/FPetitprez/webMCP-counter/issues/new/choose","webMCP-counter Github page.",target="_blank")),
                          tags$br(),
                          
                          
                          img(src="diagram.png", width=700),
                          tags$p(style="color:grey","Diagram showing the 2 steps of webMCP-counter.")
                          
                          
                 ),
                 
                 
                 tabPanel("Citation",
                          
                          tags$p("If you use webMCP-counter in a scientific publication, please cite the original article of the MCP-counter version you used (human or mouse), as well as:"),
                          tags$p("Meylan, M., Becht, E., sautès-Fridman, C., de Reyniès, A., Fridman, W.H. and Petitprez F. ",tags$a(href="https://biorxiv.org/cgi/content/short/2020.12.03.400754v1","webMCP-counter: a web interface for transcriptomics-based quantification of immune and stromal cells in heterogeneous human or murine samples"), "bioRXiv (2020)."),
                          tags$br(),
                          tags$p("For the human MCP-counter: Becht, E., Giraldo, N.A., Lacroix, L. et al. ", tags$a(href="https://doi.org/10.1186/s13059-016-1070-5","Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression."),"Genome Biol 17, 218 (2016)."),
                          tags$br(),
                          tags$p("For the murine mMCP-counter: Petitprez, F., Lévy, S., Sun, C.-M., Meylan, M. et al. ", tags$a(href="https://doi.org/10.1186/s13073-020-00783-w", "The murine Microenvironment Cell Population counter method to estimate abundance of tissue-infiltrating immune and stromal cell populations in murine samples using gene expression."), "Genome Med 12, 86 (2020)."),
                          
                          
                 ),
                 
                 
                 tabPanel("Version",
                          
                          tags$p("This is webMCP-counter version 1.1 (July 2021)."),
                          tags$br(),
                          tags$p(paste0("This version runs the version ",packageVersion("MCPcounter")," of the MCPcounter package and version ", packageVersion("mMCPcounter"), " of the mMCPcounter package."))
                          
                          
                          
                 )
                 
)
