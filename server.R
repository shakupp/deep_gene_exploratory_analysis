library(shiny)
library(dplyr)
library(DT)
library(httr)
library(jsonlite)
library(plotly)
source('figures.R')
source('verification.R')
source('formate.R')

# Define server
server <- function(input, output) {
  # download and verification
  data_table <- reactive({
    req(input$file)
    ext = tools::file_ext(input$file$name)
    
    # fill is a CSV?
    if (!is_csv_format(ext)) {
      return(NULL)
    }
    
    # If fill have the correct format, read it
    df = read.csv(input$file$datapath, header = TRUE)
    
    # Verify column names
    if (!is_colnames_true(colnames(df))) {
      return(NULL)  # Exit early if columns are incorrect
    }
    
    # Verify if there are NA
    if (!is_NA(df)) {
      return(NULL)  # Exit early if there are NAs
    }
    
    # Here, add colnames format and value format to be sure value and colnames feet with other part of the script
    df = formating_dataframe(df)
    
    return(df)
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
  
  # return (MA) datatable for UI
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
  
  #######################
  # Gene Research Button
  #######################
  
}

