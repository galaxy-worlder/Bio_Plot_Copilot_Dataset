# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(ggseqlogo)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# 读取数据
df_motif <- read.csv("./supplementary/supplemantarydata5.2.csv")
df_heatmap <- read.csv("./supplementary/supplemantarydata5.1.csv")

# Shiny UI
ui <- fluidPage(
    titlePanel("HLA-I motif and SHAP&LIME heatmap"),
    sidebarLayout(
        sidebarPanel(
            # Tabset for Motif and Heatmap
            tabsetPanel(
                id = "tabs", #添加一个 input$tab 变量来跟踪当前选中的选项卡
                tabPanel("Motif",
                         selectInput("hla", "Choose HLA Allele", choices = sort(unique(df_motif$hlatype))),
                         selectInput("len", "Choose Peptide Length", choices = sort(unique(df_motif$antigen_peptide_length)))
                ),
                tabPanel("Heatmap",
                         selectInput("mode", "Choose Mode", choices = sort(unique(df_heatmap$mode))),
                         selectInput("allele", "Choose HLA Allele", choices = sort(unique(df_heatmap$allele))),
                         # 默认的peptide选择为空，需要动态更新
                         selectInput("peptide", "Choose Peptide", choices = NULL)
                )
            )
        ),
        mainPanel(
            plotOutput("plotOutput", width = "550px", height = "300px")  # Both Motif and Heatmap will be rendered here
        )
    )
)

# Shiny Server
server <- function(input, output, session) { 
    # 动态生成Motif图
    output$plotOutput <- renderPlot({
    # 判断当前选中的是Motif界面还是Heatmap界面
    # 如果Motif标签被选中，input$hla 和 input$len会有值
    if (input$tabs == "Motif") {
        req(input$hla, input$len)  # 确保Motif界面的必要输入值已经选择
        # 绘制Motif图
        df_sub <- df_motif %>% filter(hlatype == input$hla & antigen_peptide_length == input$len)
        n_peptides <- nrow(df_sub)
        ggplot() +
            geom_logo(df_sub$antigen_peptide, method = "bits", seq_type = "aa") +
            ggtitle(paste0(input$hla, " ", input$len, "-mer (n=", n_peptides, ")")) +
            labs(y = "bits") +
            theme_minimal(base_size = 15) +
            theme(
                panel.grid = element_blank(),
                panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                panel.background = element_rect(fill = "white", color = "black"),
                plot.title = element_text(hjust = 0.5, size = 15),
                axis.text.x = element_text(angle = 0, hjust = 1, size = 15),
                axis.text.y = element_text(size = 15),
                legend.position = "none"
            )
    } else if (input$tabs == "Heatmap") {
        req(input$mode, input$allele, input$peptide)  # 确保Heatmap界面的必要输入值已经选择
        # 绘制Heatmap图
        data_sub <- df_heatmap %>% filter(mode == input$mode, allele == input$allele, peptide == input$peptide)
        order <- c("NetMHCpan4.1_EL", "NetMHCpan4.1_BA", "MHCflurry2.0_PS", "MHCflurry2.0_BA", "CapsNet-MHC_AN", "BigMHC", "STMHCpan", "TransPHLA")
        data_sub <- data_sub %>% mutate(tool = factor(tool, levels = order)) %>% arrange(tool)
        aa <- strsplit(unique(data_sub$peptide), split = "")[[1]]
        title <- paste0(unique(data_sub$mode), " for ", unique(data_sub$allele), "+", unique(data_sub$peptide))
        input_matrix <- as.matrix(data_sub[, grep("^pos_[1-9]$", colnames(data_sub))])
        # force all values to be numeric
        input_matrix <- matrix(as.numeric(unlist(input_matrix)),
                       nrow = nrow(input_matrix),
                       ncol = ncol(input_matrix),
                       dimnames = dimnames(input_matrix))
        # 准备行注释的数据和颜色
        anno_col <- data.frame(bind_result=factor(data_sub$bind_result_tool))
        anno_colors <- c("1" = "red", "0" = "blue")
        ha = rowAnnotation(bind = anno_simple(anno_col$bind_result, col = anno_colors))
        Heatmap(
            input_matrix,
            col = colorRamp2(c(-0.35, 0, 0.20), c("#A0C8E3", "#FDFDFD", "#FF8A8A")),
            border = TRUE,
            border_gp = gpar(col = "black"),
            cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", input_matrix[i, j]), x, y, gp = gpar(fontsize = 10))
            },
            column_title = title,
            column_title_gp = gpar(fontsize = 15),
            column_labels = aa,
            column_names_rot = 0,
            column_names_gp = gpar(fontsize = 15),
            row_labels = data_sub$tool,
            row_names_gp = gpar(fontsize = 15),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_heatmap_legend = TRUE,
            left_annotation = ha
            )
        }
    })
    # 动态更新peptide选择框
    observe({
        # 过滤出当前选中allele对应的peptides
        if (input$allele != "") {
            peptides_for_allele <- unique(df_heatmap %>% filter(allele == input$allele) %>% pull(peptide))
            # 更新peptide选择框的选项
            updateSelectInput(session, "peptide", choices = peptides_for_allele)
        }
    })
}

# 运行应用
shinyApp(ui = ui, server = server)
