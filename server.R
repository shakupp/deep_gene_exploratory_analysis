library(shiny)
library(dplyr)
library(DT)
library(httr)
library(jsonlite)
library(plotly)
library(clusterProfiler)
library(org.Mm.eg.db)  # pour souris — a adapte
library(enrichplot)
library(pathview)
library(ReactomePA)
library(archive) 
source('R/figures.R')
source('R/verification.R')
source('R/formate.R')
source("R/pathway.R")
source("R/pathway_plot.R")
source("R/pathway_module.R")

# Define server
server <- function(input, output, session) {
  if (dir.exists("www/tmp")) {
    unlink("www/tmp/*", recursive = TRUE)
  }
  
  
  
  
  # download and verification
  data_table <- reactive({
    req(input$file)
    ext = tools::file_ext(input$file$name)
    
    # fill is a CSV?
    if (!is_csv_format(ext)) {
      return(NULL)
    }
    
    # If fill have the correct format, read it (to read csv files with sep = "," or ";")
    df_try1 <- try(read.csv(input$file$datapath, header = TRUE), silent = TRUE)
    if (inherits(df_try1, "try-error") || !is_colnames_true(colnames(df_try1))) {
      df_try2 <- try(read.csv(input$file$datapath, header = TRUE, sep = ";"), silent = TRUE)
      if (!inherits(df_try2, "try-error") && is_colnames_true(colnames(df_try2))) {
        df <- df_try2
      } else {
        return(NULL)
      }
    } else {
      df <- df_try1
    }
    
    return(df)
  })
  espece_id <- reactive({
    espece_id <- input$espece_id
    req(espece_id)
  })
  tresholdLog2FoldChange <- reactive({
    tresholdLog2FoldChange <- input$tresholdLog2FoldChange
    req(tresholdLog2FoldChange)
  })
  pvalue <- reactive({
    pvalue <- input$pvalue
    req(pvalue)
  })
  organisms_table <- reactive({
    read.csv("www/orga_translate_table.csv", sep = ";") # a changer le nom du file apres 
  })
  ##############
  # Volcano plot
  ##############
  data_volcano <- reactive({
    df = data_table()
    
    # add diferential expression columns
    df = set_differential_expression(
      data_ = df,
      log2FC_ = input$slider_log2FC_volcano_plot,
      padj_ = input$slider_padj_volcano_plot
    )
    return(df)
  })
  
  filtered_data_volcano <- reactive({
    df = data_volcano()
    # Check if the button has been clicked
    if (input$show_DE_genes_volcano) {
      df = df %>% filter(expressionStatus %in% c("UP", "DOWN"))
    }
    
    event_data_volcano = event_data("plotly_relayout", source = "volcano")
    
    if (is.null(event_data_volcano)){
      return(df) # Return all data or only filtered data
    }
    
    # Axes limits
    x_min = event_data_volcano[["xaxis.range[0]"]]
    x_max = event_data_volcano[["xaxis.range[1]"]]
    y_min = event_data_volcano[["yaxis.range[0]"]]
    y_max = event_data_volcano[["yaxis.range[1]"]]
    
    # filtered data with the plots limits
    if (!is.null(x_min) && !is.null(x_max) && !is.null(y_min) && !is.null(y_max)) {
      df = df %>% filter(log2FC >= x_min & log2FC <= x_max, -log10(padj) >= y_min & -log10(padj) <= y_max)
    }
    
    return(df)
  })
  
  # return datatable for UI
  output$DTtable_volcano  <- renderDT({
    filtered_data_volcano()
  })
  
  # making the volcanoPlot
  volcano_plot_figure = reactive({
    req(data_volcano())
    
    #Get colors & titles from the UI
    color_up = input$color_up
    color_down = input$color_down
    color_neutral = input$color_neutral
    title_text = input$title_text
    subtitle_text = input$subtitle_text
    
    #Generate the volcano plot
    volcano_plot(
      data_ = data_volcano(),
      title_ = title_text,
      subtitle_ = subtitle_text,
      log2FC_ = input$slider_log2FC_volcano_plot,
      padj_ = input$slider_padj_volcano_plot,
      color_ = c(color_down, color_neutral,color_up )
    )
  })
  
  # return volcano plot for UI
  output$volcano_plot <- renderPlotly({
    ggplotly(volcano_plot_figure(), tooltip = "label",source = "volcano")
  })
  
  # Download volcano plot
  output$download_volcano_plot <- downloadHandler(
    filename = function() {
      paste(
        "volcano_plot_log2FC_",
        input$slider_log2FC_volcano_plot,
        "_padj_",
        input$slider_padj_volcano_plot,
        "_",
        Sys.Date(),
        ".png",
        sep = ""
      )
    },
    
    content = function(file) {
      ggsave(
        file,
        plot = volcano_plot_figure(),
        device = "png",
        width = 8,
        height = 6
      )
    }
  )
  
  # Download data
  output$download_data_volcano <- downloadHandler(
    filename = function() {
      paste(
        "data_log2FC_",
        input$slider_log2FC_volcano_plot,
        "_padj_",
        input$slider_padj_volcano_plot,
        "_",
        if (input$show_DE_genes_volcano)
          "DE_genes"
        else
          "all_genes",
        "_",
        Sys.Date(),
        ".csv",
        sep = ""
      )
    },
    
    content = function(file) {
      # Save CSV
      write.csv(filtered_data_volcano(), file, row.names = FALSE)
    }
  )
  
  ###########
  # MA plot
  ###########
  
  data_MA <- reactive({
    df = data_table()
    
    # add differential expression columns
    df = set_differential_expression(
      data_ = df,
      log2FC_ = input$slider_log2FC_MA_plot,
      padj_ = input$slider_padj_MA_plot
    )
    return(df)
  })
  
  filtered_data_MA <- reactive({
    df = data_MA()
    # Check if the button has been clicked
    if (input$show_DE_genes_MA) {
      df <- df %>% filter(expressionStatus %in% c("UP", "DOWN"))
    }
    
    event_data_ma = event_data("plotly_relayout", source = 'ma')
    
    # Return all data if there is no zoom
    if (is.null(event_data_ma)) {
      return(df)
    }
    
    # ma plot limits
    x_min = event_data_ma[["xaxis.range[0]"]]
    x_max = event_data_ma[["xaxis.range[1]"]]
    y_min = event_data_ma[["yaxis.range[0]"]]
    y_max = event_data_ma[["yaxis.range[1]"]]
    
    if (!is.null(x_min) && !is.null(x_max) && !is.null(y_min) && !is.null(y_max)){
      df = df %>% filter(log2(baseMean) >= x_min & log2(baseMean) <= x_max,log2FC >= y_min & log2FC <= y_max)
    }
    
    return(df)
    
  })
  
  # return (MA) data table for UI
  output$DTtable_MA_plot  <- renderDT({
    filtered_data_MA()
  })
  
  # Making MA plot
  
  MA_plot_figure = reactive({
    req(data_MA())
    
    #Get colors & titles from the UI
    color_up_MA = input$color_up_MA
    color_down_MA = input$color_down_MA
    color_neutral_MA = input$color_neutral_MA
    title_text_MA_plot = input$title_text_MA_plot
    subtitle_text_MA_plot = input$subtitle_text_MA_plot
    
    #Generate the volcano plot
    MA_plot(
      data_ = data_MA(),
      title_ = title_text_MA_plot,
      subtitle_ = subtitle_text_MA_plot,
      color_ = c(color_down_MA,color_neutral_MA,color_up_MA)
    )
  })

  # return MA plot for UI
  output$MA_plot <- renderPlotly({
    ggplotly(MA_plot_figure(), tooltip = "label", source = 'ma')
  })

  # Download MA plot
  output$download_MA_Plot <- downloadHandler(
    filename = function() {
      paste(
        "MA_log2FC_",
        input$slider_log2FC_MA_plot,
        "_padj_",
        input$slider_padj_MA_plot,
        "_",
        Sys.Date(),
        ".png",
        sep = ""
      )
    },
    
    content = function(file) {
      ggsave(
        file,
        plot = MA_plot_figure(),
        device = "png",
        width = 8,
        height = 6
      )
    }
  )

  # Download data
  output$download_data_MA_plot <- downloadHandler(
    filename = function() {
      paste(
        "data_log2FC_",
        input$slider_log2FC_MA_plot,
        "_padj_",
        input$slider_padj_MA_plot,
        "_",
        if (input$show_DE_genes_MA)
          "DE_genes"
        else
          "all_genes",
        "_",
        Sys.Date(),
        ".csv",
        sep = ""
      )
    },
    
    content = function(file) {
      # Save CSV
      write.csv(filtered_data_MA(), file, row.names = FALSE)
    }
  )
  ### chargement of the species csv to fin correspondance with  ###
  #species_reference <- read.csv("species_reference.csv", sep = ";",header = TRUE)
  # function to get the appropriate species ID based on the user input
  #get_species_id <- function(species_name, database) {}
  
  ####################
  # Pathway analysis #
  ####################
    mod_pathway <- pathwayServer(
  id = "pathway",
  input_data = data_table,
  method = reactive(input$pathway_analysis_method),
  direction = reactive(input$deg_direction),
  database = reactive(input$pathway_analysis_database),
  pval_thresh = reactive(input$slider_pval_pathway_analysis),
  run_trigger = reactive(input$run_pathway_analysis)
)
  
  
  
  
  
  
  
  ####################
  # AI Agent #
  ####################
  
  ask_groq <- function(question) {
    # Replace with your actual Groq API key
    api_key <- "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"
    #Sys.setenv(GROQ_API_KEY = "your_api_key")
    #api_key <- Sys.getenv("GROQ_API_KEY") # peut etre apres faire comme ca 
    # System message to define AI behavior
    system_message <- list(
      role = "system",
      content = "You are an AI expert for analysis of DEGs and pathways and biology. Answer client questions concisely with relevant information."
    )
    
    # Prepare the request payload
    payload <- list(
      model = "meta-llama/llama-4-scout-17b-16e-instruct",
      messages = list(
        system_message,
        list(
          role = "user",
          content = question
        )
      ),
      max_tokens = 1000,
      temperature = 0.1
    )
    
    # Make the API request
    response <- httr::POST(
      url = "https://api.groq.com/openai/v1/chat/completions",
      httr::add_headers(
        "Authorization" = paste("Bearer", api_key),
        "Content-Type" = "application/json"
      ),
      body = jsonlite::toJSON(payload, auto_unbox = TRUE),
      encode = "json"
    )
    
    # Process the response
    if (httr::status_code(response) == 200) {
      content <- httr::content(response, "parsed")
      return(content$choices[[1]]$message$content)
    } else {
      return(paste("Error:", httr::content(response, "text")))
    }
  }
  
  # Event handler for the ask button
  observeEvent(input$ask_groq, {
    req(input$user_question)
    
    # Show loading message
    output$groq_response <- renderUI({
      HTML("<div style='color:blue;'>Thinking...</div>")
    })
    
    # Get response (wrap in tryCatch for error handling)
    tryCatch({
      answer <- ask_groq(input$user_question)
      answer <- gsub("*", "", answer, fixed = TRUE)
        output$groq_response <- renderUI({
          HTML(paste0("
    <div style='
      font-family: \"Segoe UI\", Roboto, sans-serif;
      background: linear-gradient(to right, #f9f9f9, #ffffff);
      padding: 20px;
      border-radius: 15px;
      box-shadow: 0 4px 12px rgba(0,0,0,0.1);
      margin-top: 20px;
      color: #333;
      font-size: 13px;
      line-height: 1.6;
    '>
      <div style='font-weight:600; color:#2c3e50; font-size:18px; margin-bottom:10px;'>
        🤖 AI Expert:
      </div>
      ", answer, "
    </div>
  "))
      })
    }, error = function(e) {
      output$groq_response <- renderUI({
        HTML(paste0("
    <div style='
      font-family: \"Segoe UI\", Roboto, sans-serif;
      background: linear-gradient(to right, #f9f9f9, #ffffff);
      padding: 20px;
      border-radius: 15px;
      box-shadow: 0 4px 12px rgba(0,0,0,0.1);
      margin-top: 20px;
      color: #333;
      font-size: 16px;
      line-height: 1.6;
    '>
      <div style='font-weight:600; color:#2c3e50; font-size:18px; margin-bottom:10px;'>
        🤖 AI Expert:
      </div>
      ", answer, "
    </div>
  "))
      })
    })
  })
  

  
  
  ####################
  # GSEA pour GO #
  ####################
  

  
  
  goGse_annotation <-  eventReactive(input$Run_Annotation_go, {
    req(input$method_go == 2)
    espece_id <- espece_id()
    organisms_table <- organisms_table()
    go <- data_table()
    original_gene_list <- go$log2FC
    names(original_gene_list) <- go$ID
    gene_list<-na.omit(original_gene_list)
    gene_list <- sort(gene_list, decreasing=TRUE)
    orgdb_pkg_name <- organisms_table[espece_id, "OrgDb"]
    if (!require(orgdb_pkg_name, character.only = TRUE)) {
      stop(paste("Could not load package:", orgdb_pkg_name))
    }
    gse <- gseGO(geneList=gene_list, 
                 ont = input$Ontology, 
                 keyType = 'ENSEMBL', 
                 pvalueCutoff = input$pvalue_go, 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 verbose = TRUE, 
                 #OrgDb = organisms_table[espece_id,2], not working 
                 #OrgDb <- get(orgdb_pkg_name), #try later 
                 OrgDb = org.Mm.eg.db,
                 pAdjustMethod = "BH")
    
  })
  
  output$Table_go_GSEA <- DT::renderDataTable({
    gse <- goGse_annotation()
    req(gse)
    
    data <- as.data.frame(gse) %>%
      mutate(GO_term = paste0("<a href='https://www.ebi.ac.uk/QuickGO/term/", ID,"' target='_blank'>", ID,"</a>"))
    col_a_afficher <- c('GO_term','setSize', 'enrichmentScore', 'NES', 'pvalue', 'p.adjust', 'rank')
    
    datatable(data[col_a_afficher], escape = FALSE)
  })
  
  output$dotplot_gsea_go <- renderPlot({
    gse <- goGse_annotation()
    req(gse)
    dotplot(gse, showCategory = input$showCategory_dotplot, title = "GSEA Dotplot", split = ".sign") + facet_grid(.~.sign)
  })
  
  output$ridgeplot_go <- renderPlot({
    gse <- goGse_annotation()
    req(gse)
    ridgeplot(gse, showCategory = 10)
  })
  
  
  observe({
    gse <- goGse_annotation()
    req(gse)
    ids <- gse@result$ID
    updateSelectInput(session, "showCategory_gseaplot", choices = ids, selected = ids[1])
  })
  
  output$gsea_plot_go <- renderPlot({
    gse <- goGse_annotation()
    req(gse)
    gseaplot2(gse, geneSetID = input$showCategory_gseaplot)
  })
  
  
  output$goEmapPlot <- renderPlot({
    gse <- goGse_annotation()
    req(gse)
    
    # Use pairwise_termsim on the gseaResult object itself
    # (only if it's of the right class, e.g. enrichResult or gseaResult)
    go_enrich_sim <- pairwise_termsim(gse)
    
    emapplot(go_enrich_sim, layout = "kk", showCategory = 15)
  })
  
  
  output$goGseaNetplot <- renderPlot({
    gse <- goGse_annotation()
    req(gse)
    
    n_cat <- input$showCategory_netplot
    req(n_cat)
    
    cnetplot(gse, showCategory = n_cat, foldChange = gse@geneList)
  })
  
  
  output$pmcplot_go_gsea <- renderPlot({
    res <-  goGse_annotation()
    req(!is.null(res))
    
    if (nrow(as.data.frame(res)) == 0) return(NULL)
    
    # Get top 5 significant terms
    top_terms <- head(res$Description, 5)
    
    # Add period argument (e.g., 2000:2024)
    tryCatch({
      pmcplot(top_terms, period = 2000:2024) + 
        ggtitle("PubMed Publication Trends") +
        theme(legend.position = "right")
    }, error = function(e) {
      showNotification(paste("PMC plot error:", e$message), type = "error")
      plot_empty_message("PMC plot failed")
    })
  })
  
  
  
  
  ####################################################
  ### ORA pour GO goGse_enrichement
  #####################################################
  
  
  
  goGse_enrichement <- eventReactive(input$Run_Annotation_go, {
    req(input$method_go == 1)
    
    # Validate required inputs
    req(espece_id())
    req(organisms_table())
    req(data_table())
    req(input$pvalue_go)
    req(input$type_go)
    
    tryCatch({
      go <- data_table()
      
      # Prepare gene list
      original_gene_list <- go$log2FC
      names(original_gene_list) <- go$ID
      gene_list <- na.omit(original_gene_list) %>% sort(decreasing = TRUE)
      
      # Get significant genes
      sig_genes_df <- subset(go, padj < input$pvalue_go)
      genes <- sig_genes_df$log2FC
      names(genes) <- sig_genes_df$ID
      genes <- na.omit(genes)
      
      # Filter by expression direction
      threshold <- tresholdLog2FoldChange()
      genes <- switch(input$type_go,
                      "over" = names(genes)[genes > threshold],
                      "under" = names(genes)[genes < -threshold],
                      "both" = names(genes)[abs(genes) > threshold])
      
      # Validate we have genes to analyze
      req(length(genes) > 0)
      go_enrich <- enrichGO(gene = genes,
                            universe = names(gene_list),
                            OrgDb = org.Mm.eg.db, #organisms_table[espece_id(), 2], 
                            keyType = 'ENSEMBL',
                            readable = T,
                            ont = input$Ontology,
                            pvalueCutoff = input$pvalue_go, 
                            #qvalueCutoff = 0.2, 
                            
      )  
      if (is.null(go_enrich) || nrow(as.data.frame(go_enrich)) == 0) {
        showNotification("No enrichment results found.", type = "warning")
        return(NULL)
      }
      
      # Check internal slot for geneSets: 
      print(str(go_enrich@geneSets))
      print(head(go_enrich@geneSets))
      
      go_enrich
    }, error = function(e) {
      showNotification(paste("GO enrichment failed:", e$message), type = "error")
      NULL
    })
  })
  
  output$Table_go_ORA <- DT::renderDataTable({
    res <- goGse_enrichement()
    
    # Check if the result is not NULL and has a $result slot
    if (!is.null(res) && !is.null(res@result)) {
      df <- as.data.frame(res@result)
      
      if (nrow(df) == 0) {
        return(DT::datatable(data.frame(Message = "No enrichment results found.")))
      }
      
      # Add clickable links to GO terms
      df$ID <- paste0(
        '<a href="https://www.ebi.ac.uk/QuickGO/term/',
        df$ID,
        '" target="_blank">',
        df$ID,
        '</a>'
      )
      
      # Optionally update input choices if needed
      updateSelectInput(session, "paths", choices = df$Description)
      
      # Show relevant columns
      cols_to_show <- c('ID', 'Description', 'GeneRatio', 'BgRatio', 'p.adjust', 'Count')
      df_subset <- df[, intersect(cols_to_show, colnames(df))]
      
      return(DT::datatable(df_subset, escape = FALSE, options = list(pageLength = 10)))
    } else {
      # If enrichment results are NULL
      return(DT::datatable(data.frame(Message = "GO enrichment not available.")))
    }
  })
  
  output$goUpsetplot <- renderPlot({
    # Get enrichment results
    res <- goGse_enrichement()
    
    # Validate results exist and are significant
    if (is.null(res)) {
      showNotification("GO enrichment not performed yet", type = "warning")
      return(plot_empty_message("Run GO enrichment first"))
    }
    
    if (!inherits(res, "enrichResult")) {
      showNotification("Invalid enrichment results format", type = "error")
      return(plot_empty_message("Enrichment failed"))
    }
    
    res_df <- as.data.frame(res)
    if (nrow(res_df) == 0) {
      showNotification("No significant GO terms found", type = "warning")
      return(plot_empty_message("No significant results\nAdjust parameters"))
    }
    
    # Create the plot with error handling
    tryCatch({
      enrichplot::upsetplot(
        res,
        showCategory = min(15, nrow(res_df)), # Show up to 15 or available terms
        n = 10, # Number of top terms to show
        matrix.color = "steelblue",
        main.bar.color = "darkred",
        sets.bar.color = "darkblue"
      )
    }, error = function(e) {
      showNotification(paste("Upset plot failed:", e$message), type = "error")
      plot_empty_message("Plot generation error")
    })
  })
  
  # Helper function for empty plots
  plot_empty_message <- function(message) {
    ggplot() + 
      annotate("text", x = 1, y = 1, label = message, size = 6) + 
      theme_void() +
      theme(plot.background = element_rect(fill = "white"))
  }
  
  
  output$barplot_ora_go <- renderPlot({
    gse <- goGse_enrichement()
    
    barplot(gse,
            drop = TRUE,
            showCategory = input$showCategory_barplot,
            title = "Barplot",
            font.size = 8)
  })
  
  
  observeEvent(input$Run_Annotation_go, {
    res <- goGse_enrichement()
    if (!is.null(res)) {
      print(head(as.data.frame(res)))
    }
  })
  
  
  output$dotplot_sea_go <-renderPlot({
    gse<-goGse_enrichement()
    dotplot(gse, showCategory = input$showCategory_dotplot)+ ggtitle("Dotpot for SEA")
  })
  
  
  output$goplot_sea <- renderPlot({
    gse <- goGse_enrichement()
    
    # Check if result is valid
    if (is.null(gse) || nrow(as.data.frame(gse)) == 0) {
      showNotification("No enrichment terms to plot in goplot.", type = "warning")
      return(NULL)
    }
    
    # Only use split if '.sign' exists in gse@result
    df <- as.data.frame(gse)
    if (!".sign" %in% colnames(df)) {
      goplot(gse, 
             showCategory = input$showCategory_goplot,
             font.size = 8)
    } else {
      goplot(gse, 
             showCategory = input$showCategory_goplot,
             font.size = 8,
             split = ".sign")
    }
  })
  
  
  output$goNetplot <- renderPlot({
    gse <- goGse_enrichement()
    
    # Defensive checks
    if (is.null(gse) || nrow(as.data.frame(gse)) == 0) {
      showNotification("No enrichment results to plot in Netplot.", type = "warning")
      return(NULL)
    }
    
    cnetplot(gse,
             showCategory = input$showCategory_goNetplot,
             foldChange = NULL,  # Optional: add foldChange vector if you have it
             circular = FALSE,
             colorEdge = TRUE)
  })
  
  
  
  
  
  output$ORAEmapPlot <- renderPlot({
    gse <- goGse_enrichement()
    req(gse)
    
    go_enrich_sim <- pairwise_termsim(gse)
    
    emapplot(go_enrich_sim, layout = "kk", showCategory = 15)
  })
  
  
  output$pmcplot_go_ora <- renderPlot({
    res <- goGse_enrichement()
    req(!is.null(res))
    
    if (nrow(as.data.frame(res)) == 0) return(NULL)
    
    # Get top 5 significant terms
    top_terms <- head(res$Description, 5)
    
    # Add period argument (e.g., 2000:2024)
    tryCatch({
      pmcplot(top_terms, period = 2000:2024) + 
        ggtitle("PubMed Publication Trends") +
        theme(legend.position = "right")
    }, error = function(e) {
      showNotification(paste("PMC plot error:", e$message), type = "error")
      plot_empty_message("PMC plot failed")
    })
  })
  
  
  #################download 
  
  
  
  
  
  generate_dotplot_gsea_go <- function() {
    # Replace this with your actual plot generation code
    plot <- ggplot(mtcars, aes(x = mpg, y = wt)) + geom_point()
    return(plot)
  }
  
  generate_ridgeplot_go <- function() {
    # Replace this with your actual plot generation code
    plot <- ggplot(mtcars, aes(x = cyl, y = hp)) + geom_boxplot()
    return(plot)
  }
  
  generate_gsea_plot_go <- function() {
    # Replace this with your actual plot generation code
    plot <- ggplot(mtcars, aes(x = disp, y = qsec)) + geom_line()
    return(plot)
  }
  
  generate_goEmapPlot <- function() {
    # Replace this with your actual plot generation code
    plot <- ggplot(mtcars, aes(x = hp, y = drat)) + geom_point()
    return(plot)
  }
  
  generate_goGseaNetplot <- function() {
    # Replace this with your actual plot generation code
    plot <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_smooth()
    return(plot)
  }
  
  generate_pmcplot_go_gsea <- function() {
    # Replace this with your actual plot generation code
    plot <- ggplot(mtcars, aes(x = qsec, y = disp)) + geom_point()
    return(plot)
  }
  
  generate_goUpsetplot <- function() {
    # Replace this with your actual plot generation code
    plot <- ggplot(mtcars, aes(x = mpg, y = wt)) + geom_point()
    return(plot)
  }
  
  generate_barplot_ora_go <- function() {
    # Replace this with your actual plot generation code
    plot <- ggplot(mtcars, aes(x = factor(cyl), y = hp)) + geom_bar(stat = "identity")
    return(plot)
  }
  
  generate_dotplot_sea_go <- function() {
    # Replace this with your actual plot generation code
    plot <- ggplot(mtcars, aes(x = mpg, y = wt)) + geom_point()
    return(plot)
  }
  
  generate_goplot_sea <- function() {
    # Replace this with your actual plot generation code
    plot <- ggplot(mtcars, aes(x = disp, y = qsec)) + geom_line()
    return(plot)
  }
  
  generate_goNetplot <- function() {
    # Replace this with your actual plot generation code
    plot <- ggplot(mtcars, aes(x = hp, y = drat)) + geom_point()
    return(plot)
  }
  
  generate_ORAEmapPlot <- function() {
    # Replace this with your actual plot generation code
    plot <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_smooth()
    return(plot)
  }
  
  generate_pmcplot_go_ora <- function() {
    # Replace this with your actual plot generation code
    plot <- ggplot(mtcars, aes(x = qsec, y = disp)) + geom_point()
    return(plot)
  }
  
  # Download handlers
  output$download_dotplot_gsea_go <- downloadHandler(
    filename = function() {
      "dotplot_gsea_go.png"
    },
    content = function(file) {
      plot <- generate_dotplot_gsea_go()
      png(file)
      print(plot)
      dev.off()
    }
  )
  
  output$download_ridgeplot_go <- downloadHandler(
    filename = function() {
      "ridgeplot_go.png"
    },
    content = function(file) {
      plot <- generate_ridgeplot_go()
      png(file)
      print(plot)
      dev.off()
    }
  )
  
  output$download_gsea_plot_go <- downloadHandler(
    filename = function() {
      "gsea_plot_go.png"
    },
    content = function(file) {
      plot <- generate_gsea_plot_go()
      png(file)
      print(plot)
      dev.off()
    }
  )
  
  output$download_goEmapPlot <- downloadHandler(
    filename = function() {
      "goEmapPlot.png"
    },
    content = function(file) {
      plot <- generate_goEmapPlot()
      png(file)
      print(plot)
      dev.off()
    }
  )
  
  output$download_goGseaNetplot <- downloadHandler(
    filename = function() {
      "goGseaNetplot.png"
    },
    content = function(file) {
      plot <- generate_goGseaNetplot()
      png(file)
      print(plot)
      dev.off()
    }
  )
  
  output$download_pmcplot_go_gsea <- downloadHandler(
    filename = function() {
      "pmcplot_go_gsea.png"
    },
    content = function(file) {
      plot <- generate_pmcplot_go_gsea()
      png(file)
      print(plot)
      dev.off()
    }
  )
  
  output$download_goUpsetplot <- downloadHandler(
    filename = function() {
      "goUpsetplot.png"
    },
    content = function(file) {
      plot <- generate_goUpsetplot()
      png(file)
      print(plot)
      dev.off()
    }
  )
  
  output$download_barplot_ora_go <- downloadHandler(
    filename = function() {
      "barplot_ora_go.png"
    },
    content = function(file) {
      plot <- generate_barplot_ora_go()
      png(file)
      print(plot)
      dev.off()
    }
  )
  
  output$download_dotplot_sea_go <- downloadHandler(
    filename = function() {
      "dotplot_sea_go.png"
    },
    content = function(file) {
      plot <- generate_dotplot_sea_go()
      png(file)
      print(plot)
      dev.off()
    }
  )
  
  output$download_goplot_sea <- downloadHandler(
    filename = function() {
      "goplot_sea.png"
    },
    content = function(file) {
      plot <- generate_goplot_sea()
      png(file)
      print(plot)
      dev.off()
    }
  )
  
  output$download_goNetplot <- downloadHandler(
    filename = function() {
      "goNetplot.png"
    },
    content = function(file) {
      plot <- generate_goNetplot()
      png(file)
      print(plot)
      dev.off()
    }
  )
  
  output$download_ORAEmapPlot <- downloadHandler(
    filename = function() {
      "ORAEmapPlot.png"
    },
    content = function(file) {
      plot <- generate_ORAEmapPlot()
      png(file)
      print(plot)
      dev.off()
    }
  )
  
  output$download_pmcplot_go_ora <- downloadHandler(
    filename = function() {
      "pmcplot_go_ora.png"
    },
    content = function(file) {
      plot <- generate_pmcplot_go_ora()
      png(file)
      print(plot)
      dev.off()
    }
  )
  # Download all plots as a zip file
  output$download_all_plots <- downloadHandler(
    filename = function() {
      "Gene_Ontology_Enrichment_Analysis.zip"
    },
    content = function(file) {
      temp_dir <- "data/figures_output"
      # If the directory doesn't exist in the working directory, create it
      if (!dir.exists(temp_dir)) {
        dir.create(temp_dir, recursive = TRUE)
      }
      
      # Create a subdirectory for the plots
      plots_dir <- file.path(temp_dir, "Gene_Ontology_Enrichment_Analysis")
      if (!dir.exists(plots_dir)) {
        dir.create(plots_dir, recursive = TRUE)
      }
      
      plots <- list(
        "dotplot_gsea_go" = generate_dotplot_gsea_go(),
        "ridgeplot_go" = generate_ridgeplot_go(),
        "gsea_plot_go" = generate_gsea_plot_go(),
        "goEmapPlot" = generate_goEmapPlot(),
        "goGseaNetplot" = generate_goGseaNetplot(),
        "pmcplot_go_gsea" = generate_pmcplot_go_gsea(),
        "goUpsetplot" = generate_goUpsetplot(),
        "barplot_ora_go" = generate_barplot_ora_go(),
        "dotplot_sea_go" = generate_dotplot_sea_go(),
        "goplot_sea" = generate_goplot_sea(),
        "goNetplot" = generate_goNetplot(),
        "ORAEmapPlot" = generate_ORAEmapPlot(),
        "pmcplot_go_ora" = generate_pmcplot_go_ora()
      )
      
      # Save each plot as a PNG file in the specified directory
      for (name in names(plots)) {
        file_path <- file.path(plots_dir, paste0(name, ".png"))
        png(file_path)
        print(plots[[name]])
        dev.off()
      }
      
      # Zip the directory containing the plots
      zip(zipfile = file, files = file.path(plots_dir, list.files(plots_dir, full.names = TRUE)))
    }
  )
  
  # Observe event to trigger download
  observeEvent(input$trigger_download, {
    session$sendCustomMessage(type = "triggerDownload", message = "downloadAll")
  })
  
}
  

