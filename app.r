library(shiny)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(pheatmap)
library(shinyjs)
library(ggplot2)
library(ggrepel)
library(gridExtra)

load("sorted_exprs_matrix.Rdata")
load("targets.Rdata")
load("summary_table.Rdata")
load("sig_genesets_up.Rdata")
load("sig_genesets_down.Rdata")
load("gsea_for_volc_r.Rdata")

cssdir <- glue::glue("{here::here()}/custom.css")

ui <- navbarPage(
  "Hypoxia Microarray Analysis",
  tags$head(
    includeCSS(cssdir),
  ),
  tabPanel(
    "Heatmap",
    tags$div(
      class = "tabpane-child",
      tags$div(
        class = "heatmap-row",
        tags$div(
          class = "heatmap-col",
          shinyWidgets::sliderTextInput("heat_thresh",
            "How many genes to display? \
          (Top differentially regulated genes ordered by ascending P value)",
            choices = 10 * 2^(0:10),
            selected = 40, grid = TRUE
          )
        ),
        tags$div(
          class = "heatmap-col",
          checkboxInput("srownames", "Show row names", FALSE),
          checkboxInput("logtransform", "Log transform values", FALSE),
          selectInput("select", "Select annotation",
            choices = c("none", colnames(targets)),
            selected = "Condition",
            multiple = TRUE,
            selectize = TRUE
          )
        ),
        tags$div(
          class = "heatmap-col",
          radioButtons("norm", "Scale by", choices = c("row", "column", "none"))
        )
      ),
      tags$hr(),
      tags$div(
        class = "heatmap-row",
        plotOutput("distPlot")
      )
    )
  ),
  tabPanel(
    "Volcano",
    tags$div(
      class = "tabpane-child",
      tags$div(
        class = "volcano-row",
        shinyWidgets::sliderTextInput("volc_thresh",
          "How many genes to display? \
        (Top differentially regulated genes ordered by ascending P value)",
          choices = seq(0, 200, by = 10),
          selected = 10, grid = TRUE
        )
      ),
      tags$div(
        class = "volcano-row",
        plotOutput("volcano")
      )
    )
  ),
  tabPanel(
    "Functional Enrichment Plots",
    tags$div(
      class = "tabpane-child",
      tags$div(
        id = "sigsets_div",
        #h2("Top Gene Sets"),
        tags$div(
          id = "upregulated_div",
          plotOutput("up_sigsets")
        ),
        tags$hr(),
        tags$div(
          id = "downregulated_div",
          plotOutput("down_sigsets")
        ),
        tags$hr()
      ),
      tags$div(
        id = "sigvolcano_div",
        h2("Exemplar Volcano Plots"),
        tags$div(
          class = "volcano-row",
          selectInput("selvol", "Select MSigDB Hallmark Gene Set",
            choices = c(
              "Hypoxia",
              "TNF_via_NFkB",
              "Oxidative_Phosphorylation",
              "Fatty_acid_metabolism"
            ),
            selected = "Hypoxia", multiple = FALSE,
            selectize = TRUE
          )
        ),
        tags$div(
          class = "volcano-row",
          plotOutput("gsevolcano")
        )
      )
    )
  ),
  tabPanel(
    "Table",
    class = "tablepane",
    mainPanel(dataTableOutput("table"))
  )
)

server <- function(input, output, session) {
  output$distPlot <- renderPlot({
    if (input$logtransform) {
      sorted_exprs_matrix <- log2(sorted_exprs_matrix + 1)
    }

    if (is.null(input$select)) {
      mysel <- NULL
    } else if (input$select[1] == "none") {
      mysel <- NULL
    } else if (length(input$select) == 1) {
      # if the data frame has one column it converts to a factor
      # force the type to be a data frame and restore row and column names
      mysel <- as.data.frame(targets[, input$select[1]])
      rownames(mysel) <- rownames(targets)
      colnames(mysel) <- input$select[1]
    } else {
      mysel <- targets[, input$select]
    }

    pheatmap(sorted_exprs_matrix[1:input$heat_thresh, ],
      fontsize_row = 10,
      fontsize_col = 10,
      show_rownames = input$srownames,
      scale = input$norm,
      annotation_col = mysel
    )
  })

  output$volcano <- renderPlot({
    ggplot(summary_table, aes(x = logFC, y = -log10(P.Value))) +
      geom_point(aes(colour = signif5pc, alpha = signif5pc)) +
      scale_alpha_manual(values = c(0.1, 0.5), guide = FALSE) +
      geom_text_repel(aes(
        x = logFC, y = -log10(P.Value),
        label = ifelse(
          summary_table$PROBEID == summary_table$PROBEID[1:input$volc_thresh],
          summary_table$SYMBOL[order(summary_table$P.Value)], ""
        )
      ), max.overlaps = Inf) +
      labs(title = "Hypoxia") +
      xlab(expression(paste(Log[2], "(Hypoxia/Normoxia)"))) +
      ylab(expression(paste(-Log[10], "(P-Value)"))) +
      ylim(0, 20) +
      theme(
        legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))
      ) +
      geom_vline(xintercept = c(-1, 1), linetype = "longdash", col = "black") +
      geom_hline(
        yintercept = -log10(0.0014),
        linetype = "longdash",
        col = "black"
      ) +
      scale_colour_manual(values = c("black", "red")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  })

  output$up_sigsets <- renderPlot({
    ggplot(sig_genesets_up, aes(
      x = rownames(sig_genesets_up),
      y = -log10(PValue)
    )) +
      geom_bar(stat = "identity") +
      labs(title = "\nSignificantly Up in Hypoxia") +
      xlab("MSigDB \nHallmark Gene sets") +
      ylab(expression(paste(-Log[10], "(P-Value)"))) +
      scale_x_discrete(
        labels = sub("HALLMARK_", "", rownames(sig_genesets_up))[
          order(sig_genesets_up$PValue,
            decreasing = TRUE
          )
        ],
        limits = rownames(sig_genesets_up)[
          order(sig_genesets_up$PValue,
            decreasing = TRUE
          )
        ]
      ) +
      coord_flip() +
      theme_bw()
  })

  output$down_sigsets <- renderPlot({
    ggplot(sig_genesets_down, aes(
      x = rownames(sig_genesets_down),
      y = -log10(PValue)
    )) +
      geom_bar(stat = "identity") +
      labs(title = "\nSignificantly Down in Hypoxia") +
      ylab(expression(paste(-Log[10], "(P-Value)"))) +
      scale_x_discrete(
        labels = sub("HALLMARK_", "", rownames(sig_genesets_down))[
          order(sig_genesets_down$PValue,
            decreasing = TRUE
          )
        ],
        limits = rownames(sig_genesets_down)[
          order(sig_genesets_down$PValue,
            decreasing = TRUE
          )
        ]
      ) +
      coord_flip() +
      theme_bw() +
      theme(
        plot.margin = unit(c(0, 0, 0, 2), "cm"),
        axis.title.y = element_blank()
      )
  })

  output$gsevolcano <- renderPlot({
    if (input$selvol == "Hypoxia") {
      k <- "Gene Set: Hypoxia"
      p <- hypoxin
    } else if (input$selvol == "TNF_via_NFkB") {
      k <- "Gene Set: TNF Via NFkB"
      p <- tnfin
    } else if (input$selvol == "Oxidative_Phosphorylation") {
      k <- "Gene Set: Oxidative Phosphorylation"
      p <- oxphin
    } else if (input$selvol == "Fatty_acid_metabolism") {
      k <- "Gene Set: Fatty Acid Metabolism"
      p <- fatin
    }

    ggplot(summary_table, aes(x = logFC, y = -log10(P.Value))) +
      geom_point(aes(colour = p, alpha = p)) +
      scale_alpha_manual(values = c(0.01, 1), guide = FALSE) +
      labs(title = k, colour = "Key:") +
      xlab(expression(paste(Log[2], "(Hypoxia/Normoxia)"))) +
      ylab(expression(paste(-Log[10], "(P-Value)"))) +
      ylim(0, 1 + max(-log10(summary_table$P.Value[p]))) +
      xlim(
        -0.5 - max(summary_table$logFC[p]),
        0.5 + max(summary_table$logFC[p])
      ) +
      theme(
        legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))
      ) +
      geom_vline(
        xintercept = c(-1, 1), linetype = "longdash",
        col = "black"
      ) +
      geom_hline(
        yintercept = -log10(0.0014), linetype = "longdash",
        col = "black"
      ) +
      scale_colour_manual(
        labels = c("Other Genes", input$selvol),
        values = c("black", "red")
      ) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  })

  output$table <- renderDataTable(summary_table)

  observeEvent(input$refresh, {
    session$invalidate
  })
}

shinyApp(ui = ui, server = server)
