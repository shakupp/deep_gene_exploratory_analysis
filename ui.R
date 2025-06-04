library(shiny)
library(shinycssloaders)
library(shinyWidgets)
library(colourpicker)
library(shinyalert)
library(shinydisconnect)
library(shinydashboard)
library(DT)
library(fresh)
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
              p("PEARL (Pathway Enrichment Analysis Reporting Launcher) is an interactive application designed to process and analyze RNA-Seq data. The application allows users to perform various analyses, visualize results, and explore gene expression data in a user-friendly interface."),
              p("Plots and different analysis options can be selected via the menu bar on the left. Choose from a range of tools, including Volcano plots, MA plots, and Gene Ontology enrichment, to tailor your analysis based on your data needs."),
              p("Navigate to 'Whole Genome Analyses' to select either a volcano plot or an MA plot."),
              p("Upload a CSV file by selecting 'select a csv file', then go to 'whole data inspection', Volcano plot, or MA plot."),
              p("You can then adjust the settings of your plot as desired. Once you're satisfied, you can download the plots via the 'Download ... Plot' button."),
              p("You can also download the data associated with your plot. The table will take into account any zoom on the plot."),
              p("A button allows you to select only the differentially expressed genes for your table."),
              p("Finally, switching from the Volcano plot to the MA plot will not erase your work! ;)") 
      ),
      tabItem(tabName = "Home",
              h2("Welcome to PEARL"),
              p("PEARL (Pathway Enrichment Analysis Reporting Launcher) is an interactive application designed to process and analyze RNA-Seq data. The application allows users to perform various analyses, visualize results, and explore gene expression data in a user-friendly interface."),
              p("Plots and different analysis options can be selected via the menu bar on the left. Choose from a range of tools, including Volcano plots, MA plots, and Gene Ontology enrichment, to tailor your analysis based on your data needs."),
              p("Navigate to 'Whole Genome Analyses' to select either a volcano plot or an MA plot."),
              p("Upload a CSV file by selecting 'select a csv file', then go to 'whole data inspection', Volcano plot, or MA plot."),
              p("You can then adjust the settings of your plot as desired. Once you're satisfied, you can download the plots via the 'Download ... Plot' button."),
              p("You can also download the data associated with your plot. The table will take into account any zoom on the plot."),
              p("A button allows you to select only the differentially expressed genes for your table."),
              p("Finally, switching from the Volcano plot to the MA plot will not erase your work! ;)") 
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
                         downloadButton("download_dotplot_gsea_go", "Download Dotplot")),
                tabPanel("GSEA Ridgeplot",
                         plotOutput("ridgeplot_go"),
                         downloadButton("download_ridgeplot_go", "Download Ridgeplot")),
                tabPanel("GSEA plot GO",
                         plotOutput("gsea_plot_go"),
                         downloadButton("download_gsea_plot_go", "Download GSEA Plot")),
                tabPanel("Emap plot",
                         plotOutput("goEmapPlot"),
                         downloadButton("download_goEmapPlot", "Download Emap Plot")),
                tabPanel("GSEA Netplot",
                         plotOutput("goGseaNetplot"),
                         downloadButton("download_goGseaNetplot", "Download Netplot")),
                tabPanel("PubTrend",
                         plotOutput("pmcplot_go_gsea"),
                         downloadButton("download_pmcplot_go_gsea", "Download PubTrend Plot"))
              )
            ),
            
            conditionalPanel(
              condition = "input.method_go == 1",
              tabsetPanel(
                tabPanel("ORA Table", DTOutput("Table_go_ORA")),
                tabPanel("GO Up setplot",
                         plotOutput("goUpsetplot"),
                         downloadButton("download_goUpsetplot", "Download Upset Plot")),
                tabPanel("Barplot",
                         plotOutput("barplot_ora_go"),
                         downloadButton("download_barplot_ora_go", "Download Barplot")),
                tabPanel("Dotplot",
                         plotOutput("dotplot_sea_go"),
                         downloadButton("download_dotplot_sea_go", "Download Dotplot")),
                tabPanel("GO Plot",
                         plotOutput("goplot_sea"),
                         downloadButton("download_goplot_sea", "Download GO Plot")),
                tabPanel("GO Netplot",
                         plotOutput("goNetplot"),
                         downloadButton("download_goNetplot", "Download Netplot")),
                tabPanel("ORA Emap Plot",
                         plotOutput("ORAEmapPlot"),
                         downloadButton("download_ORAEmapPlot", "Download Emap Plot")),
                tabPanel("PubTrend",
                         plotOutput("pmcplot_go_ora"),
                         downloadButton("download_pmcplot_go_ora", "Download PubTrend Plot"))
              )
            )
          )
        )
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
        
        pathwayUI("pathway")
      ),
      
      tabItem(tabName = "kegg", h2("KEGG Database"), p("Content for KEGG.")),
      
      tabItem(tabName = "Reactome", h2("Reactome Database"), p("Content for Reactome."))
    ),


  # Right sidebar  
  div(
    id = "right_sidebar",
    class = "right-sidebar",
    sidebarPanel(
      width = 12,
      # Add the icon/image at the top of the sidebar
      div(style = "text-align: center; margin-bottom: 15px;",
          img(src = "photo_iconAI.png", height = "290px"),
          h3("Your personalized AI expert is here for you!", style = "color: #6c757d;")
      ),
      textInput("user_question", "Ask about DEGs/biology:", 
                placeholder = "e.g., What's a gene?"),
      actionButton("ask_groq", "Ask all about it!"),
      tags$hr(),
      htmlOutput("groq_response")
    )
  ),
  
  tags$script(HTML("
      $(document).on('click', '#toggle_sidebar', function() {
        $('#right_sidebar').toggleClass('collapsed');
      });
    "))
  )
  ) # End dashboardPage
  #============
  # End Body
  #============