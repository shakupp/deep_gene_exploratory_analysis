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
      
      # Bar Plot
      tabPanel("Bar Plot", 
               plotlyOutput(ns("barplot")),
               actionButton(ns("show_read_plot_bar"), "How to read this plot?"),
               bsModal(ns("modal_read_plot_bar"), "How to read Bar Plot", ns("show_read_plot_bar"), size = "large",
                       htmlOutput(ns("read_plot_text_bar"))
               )
      ),
      
      # Dot Plot
      tabPanel("Dot Plot", 
               plotlyOutput(ns("dotplot")),
               actionButton(ns("show_read_plot_dot"), "How to read this plot?"),
               bsModal(ns("modal_read_plot_dot"), "How to read Dot Plot", ns("show_read_plot_dot"), size = "large",
                       htmlOutput(ns("read_plot_text_dot"))
               )
      ),
      
      # Ridge Plot
      tabPanel("Ridge Plot", 
               plotOutput(ns("ridgeplot")),
               actionButton(ns("show_read_plot_ridge"), "How to read this plot?"),
               bsModal(ns("modal_read_plot_ridge"), "How to read Ridge Plot", ns("show_read_plot_ridge"), size = "large",
                       htmlOutput(ns("read_plot_text_ridge"))
               )
      ),
      
      # GSEA Plot
      tabPanel("GSEA Plot", 
               plotOutput(ns("gseaplot")),
               actionButton(ns("show_read_plot_gsea"), "How to read this plot?"),
               bsModal(ns("modal_read_plot_gsea"), "How to read GSEA Plot", ns("show_read_plot_gsea"), size = "large",
                       htmlOutput(ns("read_plot_text_gsea"))
               )
      ),
      
      # Enrichment Map
      tabPanel("Enrichment Map", 
               plotOutput(ns("emapplot")),
               actionButton(ns("show_read_plot_emap"), "How to read this plot?"),
               bsModal(ns("modal_read_plot_emap"), "How to read Enrichment Map", ns("show_read_plot_emap"), size = "large",
                       htmlOutput(ns("read_plot_text_emap"))
               )
      ),
      
      # Cnet Plot
      tabPanel("Cnet Plot", 
               plotOutput(ns("cnetplot")),
               actionButton(ns("show_read_plot_cnet"), "How to read this plot?"),
               bsModal(ns("modal_read_plot_cnet"), "How to read Cnet Plot", ns("show_read_plot_cnet"), size = "large",
                       htmlOutput(ns("read_plot_text_cnet"))
               )
      ),
      
      # Upset Plot
      tabPanel("Upset Plot", 
               plotOutput(ns("upsetplot")),
               actionButton(ns("show_read_plot_upset"), "How to read this plot?"),
               bsModal(ns("modal_read_plot_upset"), "How to read Upset Plot", ns("show_read_plot_upset"), size = "large",
                       htmlOutput(ns("read_plot_text_upset"))
               )
      ),
      
      # KEGG Mapping
      tabPanel("KEGG Mapping", 
               uiOutput(ns("keggmap")),
               actionButton(ns("show_read_plot_kegg"), "How to read this plot?"),
               bsModal(ns("modal_read_plot_kegg"), "How to read KEGG Mapping", ns("show_read_plot_kegg"), size = "large",
                       htmlOutput(ns("read_plot_text_kegg"))
               )
      ),
      
      # Tree Plot
      tabPanel("Tree Plot", 
               plotOutput("pathway_gsea_treeplot"),
               actionButton(ns("show_read_plot_tree"), "How to read this plot?"),
               bsModal(ns("modal_read_plot_tree"), "How to read Tree Plot", ns("show_read_plot_tree"), size = "large",
                       htmlOutput(ns("read_plot_text_tree"))
               )
      )
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
    
    # Bar Plot
    output$read_plot_text_bar <- renderUI({
        HTML("
    <p><b>Bar plot:</b></p>
    <ul>
      <li><b>Y-axis:</b> Pathways (or terms)</li>
      <li><b>X-axis:</b> -log<sub>10</sub>(adjusted p-value)</li>
    </ul>
    <p><b>How to interpret:</b></p>
    <ul>
      <li>Higher bars indicate more statistically significant pathways.</li>
      <li>Use this plot to quickly compare significance between pathways.</li>
    </ul>
  ")
      })
    
    # Dot Plot
    output$read_plot_text_dot <- renderUI({
      HTML("
    <p><b>Dot plot:</b></p>
    <ul>
      <li><b>Y-axis:</b> Pathways (or terms)</li>
      <li><b>X-axis:</b> Gene ratio (significant genes / total pathway genes)</li>
      <li><b>Dot size:</b> Number of significant genes in the pathway</li>
      <li><b>Dot color:</b> Adjusted p-value</li>
    </ul>
    <p><b>How to interpret:</b></p>
    <ul>
      <li>Pathways with higher gene ratio and lower p-value are more relevant.</li>
      <li>Dot size gives context on gene counts.</li>
    </ul>
  ")
    })
    
    # Ridge Plot
    output$read_plot_text_ridge <- renderUI({
      HTML("
    <p><b>Ridge plot:</b></p>
    <ul>
      <li>Shows the distribution of enrichment scores across pathways.</li>
      <li>Each ridge corresponds to a pathway.</li>
      <li>Primarily used with GSEA results.</li>
    </ul>
    <p><b>How to interpret:</b></p>
    <ul>
      <li>Peak positions indicate where the enrichment signal is strongest.</li>
      <li>Compare the shape and height of ridges between pathways.</li>
    </ul>
  ")
    })
    
    # GSEA Plot
    output$read_plot_text_gsea <- renderUI({
      HTML("
    <p><b>GSEA Enrichment Plot:</b></p>
    <ul>
      <li><b>X-axis:</b> Genes ranked by expression metric</li>
      <li><b>Green curve:</b> Running enrichment score (ES)</li>
      <li><b>Vertical bars:</b> Positions of genes from the pathway in the ranked list</li>
      <li><b>Bottom track:</b> Ranked metric (e.g. log fold change)</li>
    </ul>
    <p><b>How to interpret:</b></p>
    <ul>
      <li>The peak of the green curve indicates maximum enrichment (ES).</li>
      <li>Genes concentrated on the left → enrichment in upregulated genes.</li>
      <li>Genes on the right → enrichment in downregulated genes.</li>
    </ul>
  ")
    })
    
    # Enrichment Map
    output$read_plot_text_emap <- renderUI({
      HTML("
    <p><b>Enrichment Map (Network plot):</b></p>
    <ul>
      <li><b>Nodes:</b> Enriched pathways or terms</li>
      <li><b>Node size:</b> Gene count or other attribute (depending on settings)</li>
      <li><b>Node color:</b> Adjusted p-value</li>
      <li><b>Edges:</b> Similarity between pathways (usually Jaccard index or overlap coefficient)</li>
    </ul>
    <p><b>How to interpret:</b></p>
    <ul>
      <li>Clusters of connected nodes represent functionally related pathways.</li>
      <li>Thicker edges = higher similarity between pathways.</li>
    </ul>
  ")
    })
    
    # Cnet Plot
    output$read_plot_text_cnet <- renderUI({
      HTML("
    <p><b>Cnet plot:</b></p>
    <ul>
      <li>Visualizes relationships between genes and enriched pathways.</li>
      <li><b>Nodes:</b> Genes (circles) and pathways (squares or other shapes)</li>
      <li><b>Edges:</b> Link genes to pathways where they are involved.</li>
    </ul>
    <p><b>How to interpret:</b></p>
    <ul>
      <li>Shows which genes contribute to multiple pathways.</li>
      <li>Helps identify shared mechanisms between pathways.</li>
    </ul>
  ")
    })
    
    # Upset Plot
    output$read_plot_text_upset <- renderUI({
      HTML("
    <p><b>Upset plot:</b></p>
    <ul>
      <li><b>Bottom:</b> Pathways or terms</li>
      <li><b>Connected dots:</b> Pathways sharing common significant genes</li>
      <li><b>Top barplot:</b> Size of each intersection (number of shared genes)</li>
    </ul>
    <p><b>How to interpret:</b></p>
    <ul>
      <li>Identifies pathways that share overlapping sets of genes.</li>
      <li>Useful to explore redundancy and relationships between pathways.</li>
    </ul>
  ")
    })
    
    # KEGG Mapping
    output$read_plot_text_kegg <- renderUI({
      HTML("
    <p><b>KEGG Pathway Map (Pathview):</b></p>
    <ul>
      <li>KEGG pathway diagram showing molecular interactions and processes.</li>
      <li><b>Node color:</b> Expression values or other mapped metrics (e.g. log fold change).</li>
      <li>Nodes usually represent genes, gene products, or metabolites.</li>
      <li>Edges represent known biological interactions (activation, inhibition, etc.).</li>
    </ul>
    <p><b>How to interpret:</b></p>
    <ul>
      <li>Color coding helps identify upregulated/downregulated genes within the pathway.</li>
      <li>Focus on sub-pathways with coordinated expression changes.</li>
    </ul>
  ")
    })
    
    # Tree Plot
    output$read_plot_text_tree <- renderUI({
      HTML("
    <p><b>Tree plot:</b></p>
    <ul>
      <li>Visualizes hierarchical relationships or similarity between pathways.</li>
      <li>Each branch represents a pathway cluster.</li>
      <li>Leaves (tips) correspond to individual pathways.</li>
    </ul>
    <p><b>How to interpret:</b></p>
    <ul>
      <li>Closely related pathways cluster together in the tree.</li>
      <li>Helps to understand functional grouping of enriched pathways.</li>
    </ul>
  ")
    })
    
    
    
    
    return(list(
      selected_pathway = reactive({ input$selected_pathway }),
      result = pathway_result
    ))
  })
}
