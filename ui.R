library(shiny)
library(shinycssloaders)
library(shinyWidgets)
library(colourpicker)
library(shinyalert)
library(shinydisconnect)
library(shinydashboard)
library(DT)
library(fresh)
library(shinyBS)
library(plotly)
source('R/figures.R')
source('R/pathway_module.R')


#--------------------------------------------------------------------------
# Author : Hadrien Sofr
# PEARL : Pathway Enrichment Analysis Reporting Launcher
# Description : This R Shiny application is designed to perform 
#               pathway enrichment analysis, volcano plots, and other
#               bioinformatics analyses. It provides interactive 
#               visualization and reporting functionalities.
# Date Created : 10 October 2024
# Last Modified : 22 November 2024
#--------------------------------------------------------------------------


my_theme = create_theme(adminlte_color(light_blue = "#D8BFD8", purple = "#D8BFD8"))

  dashboardPage(
    skin = 'purple',
    #========
    #headrer
    #========
    
    dashboardHeader(
      title = "PEARL",
      tags$li(
        class = "dropdown",
        style = "display: flex; align-items: center; height: 100%;",
        actionButton(
          inputId = "help",
          label = "Help",
          icon = icon(name = "question", lib = "font-awesome"),
          style = "margin: 0 px;"
        ),
        
        tags$div(
          selectInput(
            inputId = "language",
            label = NULL,
            choices = list(
              "English" = "en",
              "French" = "fr",
              "Spanish" = "es"
            ),
            selected = "en",
            width = "100px"
          ),
          style = "margin-left: 0px;margin-top: 15px" #regarder ici comment bien définir les onglets.
        )
      )
    ),
    
    #===========
    # End Header
    #===========
    
    
    #===========
    # Sidebar
    #===========
    
    dashboardSidebar(
      sidebarMenu(
        id = "sidebar",
        selected = "Home",
        
        menuItem("", tabName = "Home2"),
        
        menuItem(
          "Home",
          tabName = "Home",
          icon = icon(name = "home", lib = "font-awesome")
        ),
        
        sidebarSearchForm(
          label = "Enter a gene",
          textId = "searchGene",
          buttonId = "searchButton",
          icon = icon(name = 'dna', lib = 'font-awesome')
        ),
        
        menuItem(
          "Whole data inspection",
          icon = icon("dna", lib = 'font-awesome'),
          tabName = "whole_data_inspection",
          
          menuSubItem(
            "Volcano Plot",
            tabName = "volcano_plot",
            icon = icon(name = 'chart-area', lib = 'font-awesome')
          ),
          
          menuSubItem(
            "MA plot",
            tabName = "MA_plot",
            icon = icon(name = 'chart-bar', lib = 'font-awesome')
          )
        ),
        
        menuItem(
          "Gene Ontology enrichement",
          icon = icon("circle-nodes", lib = 'font-awesome'),
          tabName = "gene_ontology"
        ),
        
        menuItem(
          "Pathway enrichement",
          icon = icon("network-wired", lib = 'font-awesome'),
          tabName = "pathway_enrichment"
        ),
        
        menuItem(
          "Data Base",
          icon = icon("gg"),
          menuSubItem("KEGG", tabName = "kegg"),
          menuSubItem("Reactome", tabName = "Reactome")
        ),
        
        fileInput(
          inputId = 'file',
          label = 'select a CSV file',
          placeholder = "No file selected",
          accept = c(".csv")
        ),

        selectInput(
          'organisme',
          label = 'select organism name',
          choices = list('homo sapiens', 'dog', 'cat', 'palkia')
        ),
        
        div(
          style = "position: absolute; top: 58%; left: 50%; transform: translate(-50%, -50%); width: 100%;",
          actionButton("toggle_sidebar", "Ask an Agent!", class = "sidebar-toggle")
        )
      )
    ),
  
  #============
  # End Sidebar
  #============
  
  #============
  # Body
  #============
  
  dashboardBody(
    use_theme(my_theme),
    includeCSS("www/styles.css"),

    
    tags$head(
      tags$style(HTML("
        .right-sidebar {
          position: fixed;
          right: 0;
          top: 0;
          width: 300px;
          height: 100%;
          background-color: #f8f9fa;
          padding: 15px;
          box-shadow: -2px 0 5px rgba(0,0,0,0.1);
          z-index: 1000;
          overflow-y: auto;
          transition: transform 0.3s ease;
        }
        .right-sidebar.collapsed {
          transform: translateX(100%);
        }
        .sidebar-toggle {
          position: fixed;
          right: 10px;
          top: 10px;
          z-index: 1001;
        }
      "))
    ),
    #===============================
    # Define Tab Items for each page
    #===============================
    
    tabItems(
      #############
      # Home
      #############
      # Specific home page for the final part.
      
      
      tabItem(tabName = "Home2",
              
              h2("Welcome to PEARL"),
              h4("Pathway Enrichment Analysis Reporting Launcher"),
              
              ### Intro Panel
              box(title = "What is PEARL?", 
                  status = "primary", 
                  solidHeader = TRUE, 
                  width = 12,
                  p("PEARL is an interactive application designed to process and analyze RNA-Seq data."),
                  p("It integrates differential expression results with functional enrichment analysis, enabling intuitive exploration and publication-ready visualizations."),
                  p("Navigate through the menu on the left to access the available analysis modules.")
              ),
              
              ### Step-by-step workflow
              h3("How to use PEARL"),
              
              fluidRow(
                box(title = "Upload your data", 
                    status = "primary", 
                    solidHeader = TRUE, 
                    width = 4,
                    p("Upload a CSV file with your differential expression results."),
                    tags$ul(
                      tags$li("Go to: 'Whole Data Inspection'"),
                      tags$li("Generate Volcano or MA plots")
                    )
                ),
                
                box(title = "Visualize & Explore", 
                    status = "primary", 
                    solidHeader = TRUE, 
                    width = 4,
                    p("Interact with Volcano and MA plots:"),
                    tags$ul(
                      tags$li("Adjust thresholds and visualization options"),
                      tags$li("Zoom, filter and download plots and data")
                    )
                ),
                
                box(title = "Perform Enrichment Analysis", 
                    status = "primary", 
                    solidHeader = TRUE, 
                    width = 4,
                    p("Run enrichment analysis on your dataset:"),
                    tags$ul(
                      tags$li("Gene Ontology (GO)"),
                      tags$li("Pathway analysis (KEGG, Reactome)"),
                      tags$li("Compare ORA and GSEA methods")
                    )
                )
              ),
              
              ### Additional Info
              h3("Learn more about PEARL"),
              
              bsCollapse(id = "home_info_panel", open = NULL,
                         
                         bsCollapsePanel("How does PEARL work?", 
                                         HTML("
                    <p>PEARL integrates your differential expression data with functional enrichment analysis tools, including ORA and GSEA.</p>
                    <p>Users can:</p>
                    <ul>
                      <li>Customize thresholds and parameters</li>
                      <li>Choose between GO terms, pathways, and protein domain analyses</li>
                      <li>Interactively explore results with detailed visualizations</li>
                      <li>Export figures and tables for publications or reports</li>
                    </ul>
                    ")
                         ),
                         
                         bsCollapsePanel("What's new in this version?", 
                                         HTML("
                    <ul>
                      <li> Enhanced user interface with dynamic explanations</li>
                      <li> Full integration of ORA and GSEA analysis pipelines</li>
                      <li> Support for Gene Ontology, KEGG pathways, Reactome pathways, and protein domains</li>
                      <li> Improved figure export options (high-res PNG, PDF)</li>
                      <li> Better reproducibility of the entire analysis workflow</li>
                    </ul>
                    ")
                         )
              )
      ),
      
      
      
      
      #############
      # Fin Home
      #############
      
    
      ##############
      # Volcano plot 
      ##############
      
      tabItem(tabName = "volcano_plot", h2("Volcano Plot"),
        fluidRow(
        box(
          title = "Volcano Plot",
          solidHeader = TRUE,
          status = 'primary',
          color = "olive",
          plotlyOutput("volcano_plot"),
          
          # Download button for volcano plot
          downloadButton(
            outputId = 'download_volcano_plot',
            label = 'Download Volcano Plot',
            icon = shiny::icon('download')
          ),
          
          downloadButton(
            outputId = 'download_data_volcano',
            label = 'Download Data',
            icon = shiny::icon('download')
          ),
          
          checkboxInput(
            inputId = "show_DE_genes_volcano",
            label = "Show only differentially expressed genes",
            value = FALSE
          )
        ),
        
        # Volcano plot Box settings
        box(
          title = "Settings",
          solidHeader = TRUE,
          status = "primary",
          
          sliderInput(
            inputId = "slider_log2FC_volcano_plot",
            label = "Select log2 Fold Change (log2FC) threshold:",
            min = 0,
            max = 5,
            step = 0.5,
            value = 2,
            round = FALSE,
            animate = TRUE
          ),
          
          # p-adj threshold = the most common used
          sliderTextInput(
            inputId = "slider_padj_volcano_plot",
            label = "Select adjusted p-value (p-adj) threshold:",
            choices = c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1),
            selected = 0.05
          ),
          
          # Add color picker widgets for volcano plot colors
          colourInput(
            "color_up",
            "Color for 'UP' Differential Expression",
            value = "#D55E00"
          ),
          colourInput(
            "color_down",
            "Color for 'DOWN' Differential Expression",
            value = "#56B4E9"
          ),
          colourInput(
            "color_neutral",
            "Color for 'NO' Differential Expression",
            value = "#d3d3d3"
          ),
          
          # Add text inputs for title and subtitle
          textInput("title_text", "Title of the Volcano Plot", value = "Volcano plot"),
          textInput("subtitle_text", "Subtitle of the Volcano Plot", value = "Differential Gene Expression Analysis")
        )
      ), fluidRow(
        DTOutput(
          outputId = 'DTtable_volcano',
          width = "100%",
          height = "auto",
          fill = TRUE
        )
      )),
      
      #########
      # MA plot 
      #########
      
      tabItem(tabName = "MA_plot", h2("MA plot"), 
              fluidRow(
                box(
                  title = "MA plot",
                  solidHeader = TRUE,
                  status = 'primary',
                  color = "olive",
                  plotlyOutput("MA_plot"),
                  
                  # Download button for volcano plot
                  downloadButton(
                    outputId = 'download_MA_Plot',
                    label = 'Download MA Plot',
                    icon = shiny::icon('download')
                  ),
                  
                  downloadButton(
                    outputId = 'download_data_MA_plot',
                    label = 'Download Data',
                    icon = shiny::icon('download')
                  ),
                  
                  checkboxInput(
                    inputId = "show_DE_genes_MA",
                    label = "Show only differentially expressed genes",
                    value = FALSE
                  )
                ),
                
                # Volcano plot Box settings
                box(
                  title = "Settings",
                  solidHeader = TRUE,
                  status = "primary",
                  
                  sliderInput(
                    inputId = "slider_log2FC_MA_plot",
                    label = "Select log2 Fold Change for differential expression status",
                    min = 0,
                    max = 5,
                    step = 0.5,
                    value = 0.5,
                    round = FALSE,
                    animate = TRUE
                  ),
                  
                  # p-adj threshold = the most common used
                  sliderTextInput(
                    inputId = "slider_padj_MA_plot",
                    label = "Select adjusted p-value (p-adj) for differential expression statut",
                    choices = c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1),
                    selected = 0.05
                  ),
                  
                  # Add color picker widgets for volcano plot colors
                  colourInput(
                    "color_up_MA",
                    "Color for 'UP' Differential Expression",
                    value = "#D55E00"
                  ),
                  colourInput(
                    "color_down_MA",
                    "Color for 'DOWN' Differential Expression",
                    value = "#56B4E9"
                  ),
                  colourInput(
                    "color_neutral_MA",
                    "Color for 'NO' Differential Expression",
                    value = "#d3d3d3"
                  ),
                  
                  # Add text inputs for title and subtitle
                  textInput("title_text_MA_plot", "Title of the MA plot", value = "MA plot"),
                  textInput("subtitle_text_MA_plot", "Subtitle of the MA Plot", value = "Comparison of Gene Expression Levels")
                )
              ), fluidRow(
                DTOutput(
                  outputId = 'DTtable_MA_plot',
                  width = "100%",
                  height = "auto",
                  fill = TRUE
                )
              )
      ),
      #############
      # Fin MA plot 
      #############
      
      tabItem(
        tabName = "gene_ontology",
        h2("Gene Ontology Enrichment Analysis"),
        tags$head(
          tags$script(HTML("
      Shiny.addCustomMessageHandler('triggerDownload', function(message) {
        if (message === 'downloadAll') {
          document.querySelector('#download_all_plots').click();
        }
      });
    "))
        ),
        fluidRow(
          box(
            title = "Gene Ontology Settings",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            background = "light-blue",  
            
            fluidRow(
              column(
                width = 4,
                numericInput("pvalue_go", "P-value cutoff", value = 0.05),
                selectInput(
                  "Ontology",
                  label = h4("Ontology:"),
                  choices = c("BP", "CC", "MF", "All"),
                  selected = "CC"
                ),
                sliderInput(
                  inputId = "tresholdLog2FoldChange",
                  label = "Log2FoldChange (log2FC) cutoff",
                  min = 0,
                  max = 5,
                  step = 0.1,
                  value = 0
                ),
                sliderInput("pmc_years", "PMC Plot Period:", min = 2002, max = 2025, value = c(2000, 2022), step = 1),
              ),
              column(
                width = 4,
                radioButtons(
                  inputId = "method_go",
                  label = "Analysis method",
                  choices = list(
                    "Overrepresentation analysis (ORA)" = 1,
                    "Gene Set Enrichment Analysis (GSEA)" = 2
                  ),
                  selected = 1
                ),
                radioButtons(
                  "type_go",
                  label = h4("DEG type:"),
                  choices = list(
                    "Over expressed DEG only" = "over",
                    "Under expressed DEG only" = "under",
                    "Both" = "both"
                  ),
                  selected = "over"
                ),
                selectInput(
                  "espece_id",
                  label = h4("Organism:"),
                  choices = list(
                    "Human" = 1, "Mouse" = 2, "Rat" = 3, "Yeast" = 4, "Fly" = 5,
                    "Arabidopsis" = 6, "Zebrafish" = 7, "Bovine" = 8, "Worm" = 9,
                    "Chicken" = 10, "Canine" = 11, "Pig" = 12, "Rhesus" = 13,
                    "E coli strain K12" = 14, "Xenopus" = 15, "Chimp" = 16,
                    "Anopheles" = 17, "Malaria" = 18, "E coli strain Sakai" = 19
                  ),
                  selected = NULL
                ),
                conditionalPanel(
                  condition = "input.method_go == 2",
                  selectInput("showCategory_gseaplot", "Select Gene Set for GSEA Plot", choices = NULL)
                ),
              ),
              column(
                width = 4,
                numericInput("showCategory_dotplot", "Dotplot categories", value = 10),
              #  conditionalPanel(
              #    condition = "input.method_go == 2", #GSEA
              #  numericInput("showCategory_ridgeplot", "Ridgeplot categories", value = 10)
              #  ),
                #numericInput("showCategory_emapplot", "Emapplot categories", value = 30, min = 1),
                numericInput("showCategory_netplot", "Netplot categories", value = 5, min = 1),
                numericInput("showCategory_goNetplot", "GO Netplot categories", value = 10, min = 1),
                #conditionalPanel(
                #  condition = "input.method_go == 1", #ORA
                #numericInput("showCategory_ridge_go", "Show top N GO terms:", value = 10, min = 2, max = 50)
                #),
                
                numericInput("showCategory_barplot", label = "barplot categories", value = 10, min = 1, max = 50),
                numericInput("showCategory_goplot", "goplot categories", value = 5)
                
              )
            ),
            
            actionButton("Run_Annotation_go", "Run GO Enrichment", class = "btn-primary", style = "margin-top: 10px;"),
            actionButton("triggerDownload", "Download All Plots", class = "btn-primary", style = "margin-top: 10px;"),
            
            conditionalPanel(
              condition = "input.method_go == 2",
              tabsetPanel(
                tabPanel("GSEA Table", DTOutput("Table_go_GSEA")),
                tabPanel("GSEA Dotplot",
                         plotOutput("dotplot_gsea_go"),
                         downloadButton("download_dotplot_gsea_go", "Download Dotplot"),
                         
                         actionButton("show_read_plot_dotplot_gsea_go", "How to read this plot?"),
                         
                         bsModal("modal_read_plot_dotplot_gsea_go", "How to read Dotplot", "show_read_plot_dotplot_gsea_go", size = "large",
                                 htmlOutput("read_plot_text_dotplot_gsea_go")
                         )
                )
                ,
                tabPanel("GSEA Ridgeplot",
                         plotOutput("ridgeplot_go"),
                         downloadButton("download_ridgeplot_go", "Download Ridgeplot"),
                         
                         actionButton("show_read_plot_ridgeplot_go", "How to read this plot?"),
                         
                         bsModal("modal_read_plot_ridgeplot_go", "How to read Ridgeplot", "show_read_plot_ridgeplot_go", size = "large",
                                 htmlOutput("read_plot_text_ridgeplot_go")
                         )
                ),
                tabPanel("GSEA plot GO",
                         plotOutput("gsea_plot_go"),
                         downloadButton("download_gsea_plot_go", "Download GSEA Plot"),
                         
                         actionButton("show_read_plot_gsea_plot_go", "How to read this plot?"),
                         
                         bsModal("modal_read_plot_gsea_plot_go", "How to read GSEA Plot", "show_read_plot_gsea_plot_go", size = "large",
                                 htmlOutput("read_plot_text_gsea_plot_go")
                         )
                )
                ,
                tabPanel("Emap plot",
                         plotOutput("goEmapPlot"),
                         downloadButton("download_goEmapPlot", "Download Emap Plot"),
                         
                         actionButton("show_read_plot_goEmapPlot", "How to read this plot?"),
                         
                         bsModal("modal_read_plot_goEmapPlot", "How to read Emap Plot", "show_read_plot_goEmapPlot", size = "large",
                                 htmlOutput("read_plot_text_goEmapPlot")
                         )
                )
                ,
                tabPanel("GSEA Netplot",
                         plotOutput("goGseaNetplot"),
                         downloadButton("download_goGseaNetplot", "Download Netplot"),
                         
                         actionButton("show_read_plot_goGseaNetplot", "How to read this plot?"),
                         
                         bsModal("modal_read_plot_goGseaNetplot", "How to read Netplot", "show_read_plot_goGseaNetplot", size = "large",
                                 htmlOutput("read_plot_text_goGseaNetplot")
                         )
                ),
                tabPanel("PubTrend",
                         plotOutput("pmcplot_go_gsea"),
                         downloadButton("download_pmcplot_go_gsea", "Download PubTrend Plot"),
                         
                         actionButton("show_read_plot_pmcplot_go_gsea", "How to read this plot?"),
                         
                         bsModal("modal_read_plot_pmcplot_go_gsea", "How to read PubTrend Plot", "show_read_plot_pmcplot_go_gsea", size = "large",
                                 htmlOutput("read_plot_text_pmcplot_go_gsea")
                         )
                )
              )
            ),
            
            conditionalPanel(
              condition = "input.method_go == 1",
              tabsetPanel(
                tabPanel("ORA Table", DTOutput("Table_go_ORA")),
                tabPanel("GO Up setplot",
                         plotOutput("goUpsetplot"),
                         downloadButton("download_goUpsetplot", "Download Upset Plot"),
                         
                         actionButton("show_read_plot_goUpsetplot", "How to read this plot?"),
                         
                         bsModal("modal_read_plot_goUpsetplot", "How to read UpSet Plot", "show_read_plot_goUpsetplot", size = "large",
                                 htmlOutput("read_plot_text_goUpsetplot")
                         )
                ),
                tabPanel("Barplot",
                         plotOutput("barplot_ora_go"),
                         downloadButton("download_barplot_ora_go", "Download Barplot"),
                         
                         actionButton("show_read_plot_barplot_ora_go", "How to read this plot?"),
                         
                         bsModal("modal_read_plot_barplot_ora_go", "How to read Barplot", "show_read_plot_barplot_ora_go", size = "large",
                                 htmlOutput("read_plot_text_barplot_ora_go")
                         )
                ),
                tabPanel("Dotplot",
                         plotOutput("dotplot_sea_go"),
                         downloadButton("download_dotplot_sea_go", "Download Dotplot"),
                         
                         actionButton("show_read_plot_dotplot_sea_go", "How to read this plot?"),
                         
                         bsModal("modal_read_plot_dotplot_sea_go", "How to read Dotplot", "show_read_plot_dotplot_sea_go", size = "large",
                                 htmlOutput("read_plot_text_dotplot_sea_go")
                         )
                ),
                tabPanel("GO Plot",
                         plotOutput("goplot_sea"),
                         downloadButton("download_goplot_sea", "Download GO Plot"),
                         
                         actionButton("show_read_plot_goplot_sea", "How to read this plot?"),
                         
                         bsModal("modal_read_plot_goplot_sea", "How to read GO Plot", "show_read_plot_goplot_sea", size = "large",
                                 htmlOutput("read_plot_text_goplot_sea")
                         )
                ),
                tabPanel("GO Netplot",
                         plotOutput("goNetplot"),
                         downloadButton("download_goNetplot", "Download Netplot"),
                         
                         actionButton("show_read_plot_goNetplot", "How to read this plot?"),
                         
                         bsModal("modal_read_plot_goNetplot", "How to read GO Netplot", "show_read_plot_goNetplot", size = "large",
                                 htmlOutput("read_plot_text_goNetplot")
                         )
                ),
                tabPanel("ORA Emap Plot",
                         plotOutput("ORAEmapPlot"),
                         downloadButton("download_ORAEmapPlot", "Download Emap Plot"),
                         
                         actionButton("show_read_plot_ORAEmapPlot", "How to read this plot?"),
                         
                         bsModal("modal_read_plot_ORAEmapPlot", "How to read Emap Plot", "show_read_plot_ORAEmapPlot", size = "large",
                                 htmlOutput("read_plot_text_ORAEmapPlot")
                         )
                ),
                tabPanel("PubTrend",
                         plotOutput("pmcplot_go_ora"),
                         downloadButton("download_pmcplot_go_ora", "Download PubTrend Plot"),
                         
                         actionButton("show_read_plot_pmcplot_go_ora", "How to read this plot?"),
                         
                         bsModal("modal_read_plot_pmcplot_go_ora", "How to read PubTrend Plot", "show_read_plot_pmcplot_go_ora", size = "large",
                                 htmlOutput("read_plot_text_pmcplot_go_ora")
                         )
                )
              )
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Documentation and Help",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            
            actionButton("show_choose_method_go", "How to choose method"),
            actionButton("show_about_ora_go", "About ORA"),
            actionButton("show_about_gsea_go", "About GSEA"),
            actionButton("show_advanced_methodology_go", "Advanced Methodology")
          )
        ),
        
        # Modals
        bsModal("modal_choose_method_go", "How to choose method", "show_choose_method_go", size = "large",
                htmlOutput("read_choose_method_go")
        ),
        
        bsModal("modal_about_ora_go", "About ORA", "show_about_ora_go", size = "large",
                htmlOutput("read_about_ora_go")
        ),
        
        bsModal("modal_about_gsea_go", "About GSEA", "show_about_gsea_go", size = "large",
                htmlOutput("read_about_gsea_go")
        ),
        
        bsModal("modal_advanced_methodology_go", "Advanced Methodology", "show_advanced_methodology_go", size = "large",
                htmlOutput("read_advanced_methodology_go")
        ),
        
        
        
        
      
      ),
      
      
      tabItem(
        tabName = "pathway_enrichment",
        h2("Pathway Enrichment Analysis"),
        
        fluidRow(
          box(
            title = "Pathway Settings",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            
            selectInput("pathway_analysis_method", "Analysis Method", choices = c("ORA", "GSEA")),
            selectInput("deg_direction", "Gene Regulation", choices = c("Over", "Under", "Both")),
            selectInput("pathway_analysis_database", "Database", choices = c("KEGG", "Reactome")),
            numericInput("slider_pval_pathway_analysis", "p.adjust Threshold", value = 0.05, min = 0, max = 1, step = 0.01),
            actionButton("run_pathway_analysis", "Run Analysis")
          )
        ),
        
        fluidRow(
          box(
            title = "Documentation and Help",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            
            actionButton("show_choose_method", "How to choose method"),
            actionButton("show_about_ora", "About ORA"),
            actionButton("show_about_gsea", "About GSEA"),
            actionButton("show_advanced_methodology", "Advanced Methodology")
          )
        ),
        
        bsModal("modal_choose_method", "How to choose method", "show_choose_method", size = "large",
                htmlOutput("read_choose_method")
        ),
        
        bsModal("modal_about_ora", "About ORA", "show_about_ora", size = "large",
                htmlOutput("read_about_ora")
        ),
        
        bsModal("modal_about_gsea", "About GSEA", "show_about_gsea", size = "large",
                htmlOutput("read_about_gsea")
        ),
        
        bsModal("modal_advanced_methodology", "Advanced Methodology", "show_advanced_methodology", size = "large",
                htmlOutput("read_advanced_methodology")
        ),
        
        
        pathwayUI("pathway")
      ),
      
      tabItem(tabName = "kegg", h2("KEGG Database"), p("Content for KEGG.")),
      
      tabItem(tabName = "Reactome", h2("Reactome Database"), p("Content for Reactome."))
    ),


  # Right sidebar  
  ),
  
  tags$script(HTML("
      $(document).on('click', '#toggle_sidebar', function() {
        $('#right_sidebar').toggleClass('collapsed');
      });
    "))
  )

# End dashboardPage
  #============
  # End Body
  #============