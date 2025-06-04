# pathway_module.R
plot_kegg_map <- function(df, selected_pathway) {
  if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    stop("org.Mm.eg.db is required but not installed")
  }
  
  # Create relative file if it needed
  tmp_dir <- file.path("www", "tmp")
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)
  
  converted_ids <- bitr(df$ID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  df <- merge(df, converted_ids, by.x = "ID", by.y = "ENSEMBL")
  gene_data <- df$log2FC
  names(gene_data) <- df$ENTREZID
  
  out_suffix <- paste0("pathview_", selected_pathway, "_", as.integer(Sys.time()))
  result <- try(
    pathview(
      gene.data = gene_data,
      pathway.id = selected_pathway,
      species = "mmu",
      out.suffix = out_suffix,
      kegg.native = TRUE,
      same.layer = FALSE
    ),
    silent = TRUE
  )

  # Find the good file
  files <- list.files(getwd(), pattern = paste0("^", selected_pathway, "\\.", out_suffix, ".*\\.png$"), full.names = TRUE)
  
  # generate unique name in www/tmp
  tmpfile_name <- paste0("kegg_", selected_pathway, "_", as.integer(Sys.time()), ".png")
  tmpfile_path <- file.path(tmp_dir, tmpfile_name)
  
  if (length(files) > 0 && file.exists(files[1])) {
    file.copy(files[1], tmpfile_path, overwrite = TRUE)
    message("File copy to :", tmpfile_path)
  } else {
    message("File not found for :", selected_pathway)
    png(tmpfile_path, width = 800, height = 600)
    plot.new()
    text(0.5, 0.5, "Pathview image not found.", cex = 1.5)
    dev.off()
  }
  
  list(
    src = file.path("tmp", tmpfile_name),  # relatif path www/tmp/
    contentType = 'image/png',
    width = 800,
    height = 600,
    alt = "KEGG pathway map"
  )
}


pathwayUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("pathway_selector")),
    textOutput(ns("gene_count")),
    tabsetPanel(
      tabPanel("Results", DTOutput(ns("table"))),
      tabPanel("Bar Plot", plotlyOutput(ns("barplot"))),
      tabPanel("Dot Plot", plotlyOutput(ns("dotplot"))),
      tabPanel("Ridge Plot", plotOutput(ns("ridgeplot"))),
      tabPanel("GSEA Plot", plotOutput(ns("gseaplot"))),
      tabPanel("Enrichment Map", plotOutput(ns("emapplot"))),
      tabPanel("Cnet Plot", plotOutput(ns("cnetplot"))),
      tabPanel("Upset Plot", plotOutput(ns("upsetplot"))),
      tabPanel("KEGG Mapping", uiOutput(ns("keggmap"))),
      tabPanel(
        "Tree Plot",
        plotOutput("pathway_gsea_treeplot")
      )
      ,
    ))
  }

pathwayServer <- function(id, input_data, method, direction, database, pval_thresh, run_trigger) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    filtered_data <- eventReactive(run_trigger(), {
      req(input_data())
      filter_input_data(input_data(), method(), direction())
    })
    
    pathway_result <- reactive({
      req(filtered_data())
      run_pathway_analysis(filtered_data(), method(), database())
    })
    
    output$gene_count <- renderText({
      df <- input_data()
      converted_ids <- bitr(df$ID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
      paste("Mapped genes:", nrow(converted_ids))
    })
    
    output$pathway_selector <- renderUI({
  res <- pathway_result()
  if (!is.null(res)) {
    paths <- res@result
    paths <- paths[paths$p.adjust < pval_thresh(), ]
    
    message("DEBUG - pathway_selector found n=", nrow(paths))
    message("DEBUG - pathway IDs: ", paste0(paths$ID, collapse = ", "))
    
    if (nrow(paths) == 0) return(NULL)
    choices <- setNames(paths$ID, paste0(paths$ID, " - ", paths$Description))
    selectInput(ns("selected_pathway"), "Select a pathway to map (KEGG only)", choices = choices)
  } else NULL
})
    
    output$table <- renderDT({
      req(pathway_result())
      res <- as.data.frame(pathway_result()@result)
      res <- res[res$p.adjust < pval_thresh(), ]
      datatable(res)
    })
    
    output$barplot <- renderPlotly({
      req(pathway_result())
      ggplotly(plot_barplot(pathway_result(), pval_thresh()))
    })
    
    output$dotplot <- renderPlotly({
      req(pathway_result())
      ggplotly(plot_dotplot(pathway_result(), pval_thresh()))
    })
    
    output$ridgeplot <- renderPlot({
      plot_ridgeplot(pathway_result(), pval_thresh())
    })
    
    output$gseaplot <- renderPlot({
      plot_gseaplot(pathway_result(), pval_thresh())
    })
    
    output$emapplot <- renderPlot({
      plot_emapplot(pathway_result(), pval_thresh())
    })
    
    output$cnetplot <- renderPlot({
      plot_cnetplot(pathway_result(), pval_thresh())
    })
    
    output$upsetplot <- renderPlot({
      plot_upsetplot(pathway_result(), pval_thresh())
    })
    
    output$keggmap <- renderUI({
      req(database() == "KEGG")
      req(filtered_data())
      req(!is.null(input$selected_pathway), input$selected_pathway != "")
      
      img_data <- plot_kegg_map(filtered_data(), input$selected_pathway)
      tags$img(src = img_data$src, width = "100%", alt = "KEGG pathway map")
    })
    

    
 

 
    output$pathway_gsea_treeplot <- renderPlot({
      message("Rendering Pathway Tree Plot")
      plot <- pathwayGSEATreePlot()
      if (is.null(plot)) {
        ggplot() + theme_void() + annotate("text", x = 1, y = 1, label = "An error occurred while generating the Tree Plot")
      } else {
        plot
      }
    })
    
    
    
    return(list(
      selected_pathway = reactive({ input$selected_pathway }),
      result = pathway_result
    ))
  })
}
