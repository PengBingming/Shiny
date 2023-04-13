
rm(list = ls())

# 一、加载包 -------------------------------------------------------------------
library(shiny)
library(shinydashboard)

library(ggplot2)

library(forcats)
library(ggstance)
library(patchwork)
library(pathview)
library(ReactomePA)
library(enrichplot)

library(clusterProfiler)
library(org.Hs.eg.db) # 人基因数据库
library(org.Mm.eg.db)


# # ② 设置shiny界面参数 ui -----------------------------------------------------------
# # ②.1页面标题 ------------------------------------------------------------------
header <- dashboardHeader(title = "KEGG/GSEA")

geneName <- intersect(columns(org.Mm.eg.db),columns(org.Hs.eg.db) )
# # ②.2 侧边界面 ------------------------------------------------------------------
sidebar <- dashboardSidebar(
  fileInput(inputId = "file_deseq", "基因选择：.RDS 文件",
            multiple = TRUE),
  tags$hr(),
  radioButtons("species", "物种：",
               choices = c('人' = "hsa",
                           '鼠' = "mmu"),
               selected = "hsa"),
  tags$hr(),
  selectInput("from_geneName","基因名：form", geneName, selected = "ENTREZID"),
  selectInput("to_geneName","基因名：to", geneName , selected = "ENTREZID"),
  numericInput(inputId = "pvalue",
               label = "pvalue：",
               value = 1),
  numericInput(inputId = "padj",
               label = "padj：",
               value = 1),
  numericInput(inputId = "logFC",
               label = "logFC：",
               value = 1.5)
)


# # ②.3 网页呈现的内容  --------------------------------------------------------------
body <- dashboardBody( 
  
  # navbarPage 1 DEseq2
  tabsetPanel(
      tabPanel('KEGG enrich',
               fluidRow( 
                 column(width = 6, 
                        box(title = 'KEGG',
                            width = NULL, 
                            height = 10,
                            status = "primary",
                            solidHeader = T),
                        tags$br(),
                        box(
                          width = 4, 
                          height = 5, 
                          solidHeader = T,
                          selectInput(
                            "plot_kegg_id1",
                            "图形选择：",
                            c("p_dot",
                              "p_bar",
                              "p_emap") ) ),
                        box(
                          width = 4, 
                          height = 5, 
                          solidHeader = T,
                          numericInput(inputId = "kegg_num1", 
                                       label = '通路数',
                                       value = 10 ) )
                        ),
                 column(width = 6, 
                        box(
                          title = "Gene",
                          width = NULL, 
                          height = 10, 
                          status = "primary",
                          solidHeader = T) ,
                        tags$br(),
                        box(
                            width = 4, 
                            height = 5, 
                            solidHeader = T,
                            selectInput(
                              "plot_kegg_id2",
                              "图形选择：",
                              c("p1", 
                                "p2",
                                "p3",
                                "p4") ) ),
                        box(
                          width = 4, 
                          height = 5, 
                          solidHeader = T,
                          numericInput(inputId = "kegg_num2", 
                                       label = '通路',
                                       value = 3 ) ),
                        box(
                          width = 4, 
                          height = 5, 
                          solidHeader = T,
                          radioButtons("kegg_path_select", "展示选择：",
                                       choices = c('通路数目' = "num",
                                                   '特定通路' = "uni"),
                                       selected = "num")
                          )
                        
                        )
                 ),
               tags$br(),
               tags$br(),
               tags$br(),
               tags$br(),
               splitLayout(cellWidths = c("50%","50%"),
                           plotOutput("plot_kegg1") ,
                           plotOutput("plot_kegg2") ) 
               ), # tabPanel 1.1
     
       # 1.2 分析结果
      tabPanel('GSEA enrich',
               fluidRow( 
                 column(width = 6, 
                        box(title = 'GSEA',
                            width = NULL, 
                            height = 10,
                            status = "primary",
                            solidHeader = T),
                        tags$br(),
                        box(
                          width = 4, 
                          height = 5, 
                          solidHeader = T,
                          selectInput(
                            "plot_gsea_id1",
                            "图形选择：",
                            c(
                              "p_bar",
                              "p_emap") ) ),
                        box(
                          width = 4, 
                          height = 5, 
                          solidHeader = T,
                          numericInput(inputId = "gsea_num1", 
                                       label = '通路数',
                                       value = 10 ) )
                 ),
                 column(width = 6, 
                        box(
                          title = "Gene",
                          width = NULL, 
                          height = 10, 
                          status = "primary",
                          solidHeader = T) ,
                        tags$br(),
                        box(
                          width = 4, 
                          height = 5, 
                          solidHeader = T,
                          selectInput(
                            "plot_gsea_id2",
                            "图形选择：",
                            c("p_gsea",
                              "p1", 
                              "p2",
                              "p3",
                              "p4") ) ),
                        box(
                          width = 4, 
                          height = 5, 
                          solidHeader = T,
                          numericInput(inputId = "gsea_num2", 
                                       label = '通路',
                                       value = 3 ) ),
                        box(
                          width = 4, 
                          height = 5, 
                          solidHeader = T,
                          radioButtons("gsea_path_select", "展示选择：",
                                       choices = c('通路数目' = "num",
                                                   '特定通路' = "uni"),
                                       selected = "num")
                        )
                        
                 )
               ),
               tags$br(),
               tags$br(),
               tags$br(),
               tags$br(),
               splitLayout(cellWidths = c("50%","50%"),
                           plotOutput("plot_gsea1") ,
                           plotOutput("plot_gsea2") ) 
      ) # tabPanel 1.2
      
  ) # navbarPage
  
) # dashboardBody

# # ②.4 ui -----------------------------------------------------------------
ui <- dashboardPage(header, sidebar, body)

server <-  function(input, output, session){
  
  plot_kegg <- reactive({
    
    file_deseq <- input$file_deseq
    
    if(is.null(file_deseq)){return(NULL)}
     
    gene_list <- readRDS(file_deseq$datapath)
    
    gene_list$gene_up    -> gene_up
    gene_list$gene       -> gene 
    gene_list$gene_diff  -> gene_diff
    gene_list$gene_all   -> gene_all
    gene_list$gene_down  -> gene_down
    gene_list$geneList   -> geneList
    

        # 2.1 差异基因富集
    kk <- enrichKEGG(gene          = gene_diff,
                     organism      = input$species,
                     minGSSize     = 20,
                     maxGSSize     = 500,
                     pAdjustMethod = 'BH',
                     universe      = names(geneList),
                     pvalueCutoff  = input$pvalue,
                     qvalueCutoff  = input$padj )
    
    if(input$species=='hsa' ){
      
      kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      
    } else (
      
      kk <- setReadable(kk, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
      
    )
    

    # 2.2 气泡图
    p_dot <- dotplot(kk, showCategory = input$kegg_num1 )
    
    # 2.3 柱形图
    p_bar <- barplot(kk, showCategory = input$kegg_num1 )
    
    # 2.4 通路中的基因
    
    # y1 <- c('PI3K-Akt signaling pathway',
    #         'Fc gamma R-mediated phagocytosis',
    #         'Cell adhesion molecules',
    #         'Focal adhesion')
    
    # kk@result[["Description"]]
    
    
    y1 = input$kegg_num2

    if ( input$kegg_path_select =="num" ) {
      y1 = input$kegg_num2
    } else (
      y1 = kk@result[["Description"]][input$kegg_num2]
      )
    
    p1 = cnetplot(kk, foldChange = geneList,showCategory = y1)
    
    p2 = cnetplot(kk, foldChange = geneList,showCategory = y1, circular = T)
    p3 = cnetplot(kk, foldChange = geneList,showCategory = y1, colorEdge = T, node_label="gene")
    p4 = cnetplot(kk, foldChange = geneList,showCategory = y1, layout = "gem")
    # (p1 + p2) / (p3 + p4)
    
    # 2.5 通路间联系
    kk <- enrichplot::pairwise_termsim(kk)
    p_emap <- emapplot(kk, showCategory = 20, layout="kk", cex_category=1.5,min_edge = 0.8) 
  
    
    # 保存数据待用
    plot_kegg <- list()
    
    plot_kegg$p_dot  <- p_dot
    plot_kegg$p_bar  <- p_bar
    plot_kegg$p1  <- p1
    plot_kegg$p2  <- p2
    plot_kegg$p3  <- p3
    plot_kegg$p4  <- p4
    plot_kegg$p_emap  <- p_emap
    
    return(plot_kegg)
    
  } )
  
  
  plot_gsea <- reactive({
    
    file_deseq <- input$file_deseq
    
    if(is.null(file_deseq)){return(NULL)}
    
    gene_list <- readRDS(file_deseq$datapath)
    
    gene_list$gene_up    -> gene_up
    gene_list$gene       -> gene 
    gene_list$gene_diff  -> gene_diff
    gene_list$gene_all   -> gene_all
    gene_list$gene_down  -> gene_down
    gene_list$geneList   -> geneList
    
    
    # 2.1 差异基因富集
    kk_gse <- gseKEGG(geneList      = geneList,
                      organism      = input$species,
                      minGSSize     = 20,
                      maxGSSize     = 500,
                      pAdjustMethod = 'BH',
                      pvalueCutoff  = input$pvalue,
                      verbose       = FALSE)
    
    if(input$species =='hsa'){
     
       kk_gse <- setReadable(kk_gse, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
       
    } else (
      
      kk_gse <- setReadable(kk_gse, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
      
      )
   

    # 3.2 柱状图
   p_bar <- ggplot(kk_gse, aes(NES, fct_reorder(Description, NES), fill=qvalue), showCategory= input$gsea_num1) + 
      geom_bar(stat='identity') + 
      scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
      theme_minimal() + ylab(NULL)
    
    # 3.4 通路与相应基因关联图

    if ( input$gsea_path_select =="num" ) {
      
      y2 = input$gsea_num2
      y3 = 1:input$gsea_num2
  
      
    } else if( input$gsea_path_select =="uni"){
      
      y2 = kk_gse@result[["Description"]][input$gsea_num2]
      y3 = input$gsea_num2
      
    }
     
   # 3.3 GSEA: gseaplot 特定通路情况
    p_gsea <- gseaplot2(kk_gse, geneSetID = y3, pvalue_table = T)

    p1 = cnetplot(kk_gse, foldChange=geneList,showCategory = y2)
    p2 = cnetplot(kk_gse, foldChange=geneList,showCategory = y2, circular = T)
    p3 = cnetplot(kk_gse, foldChange=geneList,showCategory = y2, colorEdge = T, node_label="gene")
    p4 = cnetplot(kk_gse, foldChange=geneList,showCategory = y2, layout = "gem")
    
    # 3.5 通路与通路间的关联
    kk_gse <- enrichplot::pairwise_termsim(kk_gse)
    p_emap <- emapplot(kk_gse, showCategory = 15, layout="kk", min_edge = 0.1, color = "NES") 
    
   # 保存数据待用
    plot_gsea <- list()
    
    plot_gsea$p_gsea  <- p_gsea
    plot_gsea$p_bar  <- p_bar
    plot_gsea$p1  <- p1
    plot_gsea$p2  <- p2
    plot_gsea$p3  <- p3
    plot_gsea$p4  <- p4
    plot_gsea$p_emap <- p_emap
    
    
    return( plot_gsea )
    
  } )
  
  # # 三、参考数据 ----------------------------------------------------------------
  
  # 3.1 展示KEGG
  output$plot_kegg1 <- renderPlot({
    
    if (is.null(plot_kegg() ) ){ return() }
    
    plot_kegg <- plot_kegg()
    
    plot_kegg$p_dot  -> p_dot
    plot_kegg$p_bar  -> p_bar
    plot_kegg$p_emap -> p_emap
    
    if(input$plot_kegg_id1 == "p_dot") {
      return( p_dot )
    }
    else if(input$plot_kegg_id1 == "p_bar"){
      return( p_bar )
    }
    else { p_emap }
  })
    
  output$plot_kegg2 <-  renderPlot({
    
    if (is.null(plot_kegg() ) ){ return() }
    
    plot_kegg <- plot_kegg()
    
    plot_kegg$p1  -> p1
    plot_kegg$p2  -> p2
    plot_kegg$p3  -> p3
    plot_kegg$p4  -> p4
  
    if(input$plot_kegg_id2 == "p1"){
      return( p1 )
    }
    else if(input$plot_kegg_id2 == "p2"){
      return( p2 )
    }
    else if(input$plot_kegg_id2 == "p3"){
      return( p3 )
    }
    else if(input$plot_kegg_id2 == "p4"){
      return( p4 )
    }
  })
  
  
  # 3.1 展示GSEA
  output$plot_gsea1 <- renderPlot({
    
    if (is.null(plot_gsea() ) ){ return() }
    
    plot_gsea <- plot_gsea()
    
    plot_gsea$p_bar  -> p_bar
    plot_gsea$p_emap -> p_emap
    
  if(input$plot_gsea_id1 == "p_bar"){
      return( p_bar )
    }
    else { p_emap }
  })
  
  output$plot_gsea2 <-  renderPlot({
    
    if (is.null(plot_gsea() ) ){ return() }
    
    plot_gsea <- plot_gsea()
    
    plot_gsea$p_gsea  -> p_gsea
    plot_gsea$p1  -> p1
    plot_gsea$p2  -> p2
    plot_gsea$p3  -> p3
    plot_gsea$p4  -> p4
    
    
    if(input$plot_gsea_id2 == "p_gsea"){
      return( p_gsea )
    }
    else if(input$plot_gsea_id2 == "p1"){
      return( p1 )
    }
    else if(input$plot_gsea_id2 == "p2"){
      return( p2 )
    }
    else if(input$plot_gsea_id2 == "p3"){
      return( p3 )
    }
    else if(input$plot_gsea_id2 == "p4"){
      return( p4 )
    }
  })
  
  
} 

shinyApp(ui, server)


