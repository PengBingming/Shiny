
rm(list = ls())

if (!require('shinydashboard', quietly = TRUE)) BiocManager::install('shinydashboard')
if (!require('shiny', quietly = TRUE)) BiocManager::install('shiny')

library(clusterProfiler) # 转换 ID
if (!require('limma', quietly = TRUE)) BiocManager::install('limma')
if (!require("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!require("org.Mm.eg.db", quietly = TRUE)) BiocManager::install("org.Mm.eg.db")


# # ② 设置shiny界面参数 ui -----------------------------------------------------------
# # ②.1页面标题 ------------------------------------------------------------------
header <- dashboardHeader(title = "RNAseq 基因表达")

geneName <- intersect(columns(org.Mm.eg.db),columns(org.Hs.eg.db) )
# # ②.2 侧边界面 ------------------------------------------------------------------
sidebar <- dashboardSidebar(
  fileInput(inputId = "file_deseq", "差异分析的 .RDS 文件",
            multiple = TRUE),
  tags$hr(),
  radioButtons("type", "数据来源：",
               choices = c('DEseq2' = "deseq",
                           'Limma' = "limma") ),
  radioButtons("species", "物种：",
               choices = c('人' = "hsa",
                           '鼠' = "mmu"),
               selected = "hsa"),
  tags$hr(),
  selectInput("from_geneName","基因名：form", geneName, selected = "ENTREZID"),
  selectInput("to_geneName","基因名：to", geneName , selected = "ENTREZID"),
  numericInput(inputId = "pvalue",
               label = "pvalue：",
               value = 0.05),
  numericInput(inputId = "padj",
               label = "padj：",
               value = 0.05),
  numericInput(inputId = "logFC",
               label = "logFC：",
               value = 1.5)
  
  
)


# # ②.3 网页呈现的内容  --------------------------------------------------------------
body <- dashboardBody( 
  
  # navbarPage 1 DEseq2
  tabsetPanel(
    tabPanel('标准化的数据',
             tags$br(),
             fluidRow(
               column(width = 12, 
                      downloadButton("downloadSampleData", "参考数据"),
                      box(
                        title = "矩阵：",
                        width = NULL, 
                        height = 10, 
                        status = "primary",
                        solidHeader = T) 
               )
             ),
             tags$br(),
             splitLayout(cellWidths = c("90%"),
                         dataTableOutput("sample1_exp_deseq"),' ' ) 
    ), # tabPanel 1.1
    
    # 1.2 分析结果
    tabPanel('基因 ID 转换',
             tags$br(),
             fluidRow( column(width = 12, 
                              downloadButton("downloadData", "下载基因"),
                              box(title = "基因 ID 转换：",
                                  width = NULL, 
                                  height = 10, 
                                  status = "primary",
                                  solidHeader = T ) )
             ),
             tags$br(),
             splitLayout(cellWidths = c("90%"),
                         dataTableOutput("DEG_degseq") ) ) # tabPanel 1.2
    ) # navbarPage
  
) # dashboardBody

# # ②.4 ui -----------------------------------------------------------------
ui <- dashboardPage(header, sidebar, body)

server <-  function(input, output, session){
  
  plots_DEseq <- reactive({
    
    file_deseq <- input$file_deseq
    
    if(is.null(file_deseq)){return(NULL)}
    
    result <- readRDS(file_deseq$datapath)
    
    if(input$type == 'limma'){
      
      dat <- data.frame(result$DEG1)
     
      dat$pvalue <- dat$P.Value
      dat$padj <- dat$adj.P.Val
      
    } else if(input$type =='deseq' ){
      dat <- data.frame(result$DEG1)
      dat$logFC <- dat$log2FoldChange

    }
    
    if(input$from_geneName==input$to_geneName){
      
      dat1 <- dat
      dat1 <- cbind("ENTREZID"=rownames(dat1), dat1 )
      
      
    } else if(!input$from_geneName==input$to_geneName){
      
      dat$ID = rownames(dat)
      if(input$species == 'hsa' ){
        
        dfName <- bitr(unique(dat$ID), fromType <- input$from_geneName,
                       toType <- input$to_geneName,
                       OrgDb <- org.Hs.eg.db)
      }
      
      if(input$species == 'mmu' ){
        
        dfName <- bitr(unique(dat$ID), fromType <- input$from_geneName,
                       toType <- input$to_geneName,
                       OrgDb <- org.Mm.eg.db)
      }
      
      dat1 <- merge(dat,dfName,by.y=input$from_geneName, by.x='ID') 
      dat <- dat[,setdiff(colnames(dat), 'ID')] 
      
    } 
    
    # 保存数据待用
    df_exp <- list()
    
    df_exp$group <- result$group
    df_exp$DEG  <- dat
    df_exp$DEG_id <- dat1
    
    return(df_exp)
    
  } )
  
  # # 三、参考数据 ----------------------------------------------------------------
  
  # 3.1 展示参考数据
  output$sample1_exp_deseq <- renderDataTable({
    
    df_exp <- plots_DEseq()
    dat <- df_exp$DEG
    dat <- cbind('ID'=rownames(dat),dat)
    
    if (!is.null(plots_DEseq() ) ) { return(dat ) }
    sample1 <- readRDS("/srv/shiny-server/RNAseq/file/result.RDS")
    
    sample1 <- sample1$vsdmat
    
    sample1 <- cbind("ID"=rownames(sample1), sample1)
    
    return(sample1)
    
  }) 
  
  select_gene <- reactive({
    
    if (is.null(plots_DEseq() ) ) { return(NULL ) }
    
    df_exp <- plots_DEseq()
    DEG <- df_exp$DEG_id
    
    DEG$change <- as.factor(
      ifelse(
        DEG$pvalue < input$pvalue & abs(DEG$logFC) > input$logFC,
        ifelse(
          DEG$logFC> input$logFC,'UP','DOWN'),
        'NOT'))
    table(DEG$change) # 查看基因上、下调情况
    
    return(DEG)
    
  })
  
  
  # 4.2 转换 ID 后的结果
  output$DEG_degseq <- renderDataTable({
    
    if (is.null(select_gene() ) ) { return(NULL ) }
    
    DEG <- select_gene()
    
    return( DEG )
    
  } ) 
  

  # 提取基因数据
  # 1.4 确定基因为上调、下调或无明显改变

  
  output$downloadData <- downloadHandler(
    
    filename = function() {
      paste("gene_list.RDS")
    },
    
    content = function(file) {
      
      DEG <- select_gene()
      
      # 3.1 选出上调基因的'ENTREZID' ID
      gene_up <- DEG[DEG$change == 'UP','ENTREZID']
      
      # 3.2 下调基因ID
      gene_down <- DEG[DEG$change == 'DOWN','ENTREZID']
      
      # 3.3 差异基因ID
      gene_diff <- c(gene_up,gene_down)
      
      # 3.4 所有基因ID
      gene_all <- as.character(DEG[ ,'ENTREZID'] )
      
      # 3.4 geneList: LogFC
      geneList <- DEG$logFC # 把 DEG 数据logFC列值赋值给数据geneList
      names(geneList) <- DEG$ENTREZID # 把ID赋值给geneList数据的名字
      geneList <- sort(geneList, decreasing = T) # 把数据进行排序
      
      # 3.5 选取 logFC >= 2 的基因
      gene <- names(geneList)[abs(geneList) > input$logFC ]
      
      # 保存数据待用
      gene_list <- list()
      
        gene_list$gene_up <- gene_up
        gene_list$gene <- gene
        gene_list$gene_all <- gene_all
        gene_list$gene_diff <- gene_diff
        gene_list$gene_down <- gene_down
        gene_list$geneList <- geneList

        saveRDS(gene_list, file = file)
        
    } )
  
  
  # 3.2.1 下载参考数据  DEseq2
  output$downloadSampleData <- downloadHandler(
    filename = function() {
      paste('DEseq2.RDS')
    },
    content = function(file) {
      
      sample1 <- readRDS("/srv/shiny-server/RNAseq/file/result.RDS")
      saveRDS(sample1, file = file )
      
    }
  )
  
  
  
  
  
  
} 

shinyApp(ui, server)


