

# # ① 加载包 -------------------------------------------------------------------
if (!require('shinydashboard', quietly = TRUE)) BiocManager::install('shinydashboard')
if (!require('shiny', quietly = TRUE)) BiocManager::install('shiny')
if (!require('limma', quietly = TRUE)) BiocManager::install('limma')
if (!require("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!require("readxl", quietly = TRUE)) BiocManager::install("readxl")
if (!require("xlsx", quietly = TRUE)) BiocManager::install("xlsx")

library(clusterProfiler) # 转换 ID
if (!require("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!require("org.Mm.eg.db", quietly = TRUE)) BiocManager::install("org.Mm.eg.db")

if (!require("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
library("FactoMineR")#画主成分分析图需要加载这两个包
library("factoextra") 

load("./df_sample.Rdata")

# 基因表达箱图函数
gene_exp = function(data, name, group){         #自定义一个函数bp，函数为{}里的内容
  library(ggpubr)
  df = data.frame(gene = data[name,], stage = group)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter", ylab = name) +
  #  Add p-value
  # stat_compare_means()
  stat_compare_means(# method = "anova"
                     )

  return(p)
}


# # ② 设置shiny界面参数 ui -----------------------------------------------------------
# # ②.1页面标题 ------------------------------------------------------------------
header <- dashboardHeader(title = "RNAseq : DEseq2")

geneName <- intersect(columns(org.Mm.eg.db),columns(org.Hs.eg.db) )
# # ②.2 侧边界面 ------------------------------------------------------------------
sidebar <- dashboardSidebar(
  fileInput(inputId = "file", "输入Excel 文件：整数型矩阵",
            multiple = TRUE,),
  h5('支持：.xlsx .xls 文件'),
  actionButton("submit1", "展示示例/开始运算"),
  tags$hr(),
  numericInput(inputId = "rowSum_filter",
               label = h6("过滤低表达基因：基因在各样本的计数和不低于"),
               value = 1) ,
  tags$hr(),
  selectInput("untrt", "对照组：", c("") ),
  selectInput("trt",   "实验组：", c("") ),
  tags$hr(),
  radioButtons("species", "物种：",
               choices = c('人' = "hsa",
                           '鼠' = "mmu"),
               selected = "hsa"),
  tags$hr()
)


# # ②.3 网页呈现的内容  --------------------------------------------------------------
body <- dashboardBody( 
  
  # tabsetPanel 1 DEseq2
  tabsetPanel(
    # 1.1 未校正 参考数据 / 输入数据
    tabPanel('数据格式',
             fluidRow( 
               column(width = 2, 
                      box(title = '下载参考数据：',
                          width = NULL, 
                          height = 10,
                          status = "primary",
                          solidHeader = T)
               ),
               column(width = 7, 
                      box(
                        title = "矩阵：",
                        width = NULL, 
                        height = 10, 
                        status = "primary",
                        solidHeader = T),
                      
                      br(),
                      h5('Excel中，表1 为矩阵，基因列名需为：ID 。
                         样本名 需与 表2 group列对应。
                         condition 列为样本状态')
               ),
               column(width = 3 , 
                      box(
                        title = "分组：",
                        width = NULL, 
                        height = 10, 
                        status = "primary",
                        solidHeader = T),
                      br(),
                      h5('表2 列名分别为：group 与 condition')
                      )
               
             ),
             tags$br(),
             splitLayout(cellWidths = c("18%", "60%","2","20%"),
                         box(width = 8 ,
                             downloadButton("downloadSampleData", "参考数据"),
                             tags$br(),
                             tags$br(),
                             numericInput(inputId = "sample_exp_row",
                                          label = "下载矩阵行数：",
                                          value = 2000),
                             solidHeader = T),
                         dataTableOutput("sample1_exp"), 
                         '',
                         dataTableOutput("sample2_group") ) 
    ), # tabPanel 1.1
    
    # 1.2 分析结果
    tabPanel('分析结果',
             downloadButton("downloadData", "分析结果: .RDS文件"),
             downloadButton("downloadDEG", "分析结果: .csv文件"),
             tags$br(),
             tags$br(),
             fluidRow(
               column(width = 12 , 
                              box( 
                                title = "差异结果：DEG",
                                width = NULL, 
                                height = 10, 
                                status = "primary",
                                solidHeader = T ) )
             ),
             tags$br(),
             splitLayout(cellWidths = c("90%"),
                         dataTableOutput("DEG") ) ), # tabPanel 1.2
    # 1.3 标准化数据
    tabPanel('标准化数据',
             fluidRow(
               downloadButton("downloadnorm", "标准化数据"),
               tags$br(),
               tags$br(),
               column(width = 2,
                      box(title = "标准化方式：",
                          width = NULL, 
                          height = 5, 
                          solidHeader = T ,
                          status = "primary",
                          radioButtons("plots_id2", "标准化方式：",
                                       choices = c('norm' = "pca",
                                                   'vsd' = "vsd",
                                                   'rld' = "rld"),
                                       selected = "vsd") ) ),
               column( width = 10,
                       box(
                         title = "矩阵：",
                         width = NULL, 
                         height = 10, 
                         status = "primary",
                         solidHeader = T) 
               )
             ),
             tags$br(),
             splitLayout(cellWidths = c('15%','85%'),'',
                         dataTableOutput("sample1_exp_deseq") ) 
    ), # tabPanel 1.1,
    # 1.3 图形
    tabPanel('基础图形',
             fluidRow( column(width = 12, 
                              box( 
                                title = "各类图形：PCA、热图、火山图",
                                width = NULL, 
                                height = 10, 
                                status = "primary",
                                solidHeader = T ) )
             ),
             tags$br(),
             fluidRow( 
               column(width = 2,
                      actionButton("submit4", "开始画图") ),
               column(width = 2,
                      downloadButton("downloadPlots_png", "下载图形：.png 文件") ),
               column(width = 2,
                      downloadButton("downloadPlots", "下载图形：.Rdata 文件") )
               ),
             tags$br(),
             fluidRow( 
               column(width = 12,
                      box( width = 2, 
                           height = 5, 
                           solidHeader = T ,
                           radioButtons("plots_id1", "图形选择：",
                                        choices = c('PCA：pca' = "pca",
                                                    'PCA：vsd' = "vsd",
                                                    'PCA：rld' = "rld",
                                                    '热图：按Pvalue排序'     = "heatmap",
                                                    '热图：按logFC排序'     = "heatmap1",
                                                    '火山图'   = "volcano"),
                                        selected = "vsd") ),
                      box( 
                        width = 2, 
                        height = 5, 
                        solidHeader = T ,
                        numericInput(inputId = "heatmap_num",
                                     label = "热图基因数：",
                                     value = 50) ),
                      box(
                        width = 2, 
                        height = 5,
                        solidHeader = T ,
                        numericInput(inputId = "pvalue",
                                     label = "pvalue：热图、火山图",
                                     value = 0.05)),
                      box(
                        width = 2, 
                        height = 5, 
                        solidHeader = T ,
                        numericInput(inputId = "padj",
                                     label = "padj：热图、火山图",
                                     value = 0.1) ),
                      box(
                        width =2, 
                        height = 5,
                        solidHeader = T ,
                        numericInput(inputId = "logFC_cutoff",
                                     label = "logFC：热图、火山图",
                                     value = 1.5) )
                      
               ) 
             ),
             tags$br(),
             tags$br(),
             tags$br(),
             splitLayout(cellWidths = c("25%","50%","25%"),'' ,
                         plotOutput("plots"),'' ) ), # tabPanel 1.3

  # navbarPage 1 DEseq2

    # 1.3 图形
    tabPanel('基因表达',
             fluidRow( column(width = 5, 
                              box(
                                title = "标准化矩阵：基因 ID 转换",
                                width = NULL, 
                                height = 10, 
                                status = "primary",
                                solidHeader = T ) ),
                       column(width = 7, 
                              box(
                                title = "标准化矩阵：基因表达箱图：",
                                width = NULL, 
                                height = 10, 
                                status = "primary",
                                solidHeader = T ) )
             ),
             tags$br(),
             fluidRow( 
               column(width = 1,
                      br(),
                      actionButton("submit2", "开始转换")
                     ),
               column(width = 2,
                      selectInput("from_geneName","基因名：form", geneName, selected = "ENTREZID" )
                     ),
               column(width = 2,
                      selectInput("to_geneName","基因名：to", geneName , selected = "ENTREZID" )
                     ),
               column(width = 1,
                      br(),
                      actionButton("submit3", "开始画图")
                      ),
               column(width = 2,
                      selectInput("gene","基因",c("") )
                      )
             ),
             h5('from 为 ID 对应基因名类型，转换基因名相同时则不进行转换；转换后的基因名在最后一列'),
             splitLayout(cellWidths = c("42%","58%"),
                         dataTableOutput("DEG_degseq"),
                         plotOutput("plots_DEseq")
                           ) ) # tabPanel 1.3
  
  ) # tabsetPanel
) # dashboardBody

# # ②.4 ui -----------------------------------------------------------------
ui <- dashboardPage(header, sidebar, body)


# # ③ 设置 server -------------------------------------------------------------
server <- function(input, output, session) {
  
  # # 一、DEseq2 分析函数 ----------------------------------------------------------------
  
  
  
  # 1.1 DEseq2 差异分析 
  myDEseq1 <- reactive({
    
    file1 <- input$file
    if(is.null(file1)){
      df1 <- df_sample1
      df1 <- head(df1,2000)
      
      df2 <- df_sample2
    } 
    else{
      df1 <- data.frame( read_excel(file1$datapath,1) )
      df2 <- data.frame( read_excel(file1$datapath,2) )
      }

    df2 <- df2[order(df2$condition),]
    df1 <- df1[,c("ID",df2$group)]
    
    RNAseq <- list()
    RNAseq$exp <- df1
    RNAseq$group <- df2
    
    exp <- RNAseq$exp # 原始表达矩阵数据
    
    # colnames(exp) <- toupper(colnames(exp)) # 列名大写
    # library(limma)
    exp <- avereps(exp[ , setdiff(colnames(exp),'ID')], ID = exp$ID ) # 去重复，ID赋值到行名
    
    # 使用 DEseq 处理的数据必须是整数，设置数据为整数
    forceMatrixToInteger <- function(data){
      apply (data, c (1, 2), function (x) {
        (as.integer(x))
      })
    }
    
    exp <- forceMatrixToInteger(data = exp)
    
    # 1.2 分组信息
    
    group <- RNAseq$group # 分组信息
    
    colnames(group) <- tolower(colnames(group))
    
    colData <- data.frame(group$condition)
    rownames(colData) <- as.character(group$group)  # 设置行名为样本名信息
    colnames(colData) <- 'condition' # 列名
    
    # 二、构建 DEseq对象 -------------------------------------------------------------
    
    # 2.1 对表达矩阵数据进行筛选
    index1 <- c(rowSums(exp) > input$rowSum_filter ) # 判断
    exp2 <- exp[index1,] # 筛选
    
    # 2.2 构建分组数据框 colData
    colData$condition <- factor(colData$condition)
    
    # 三、DEseq2处理 --------------------------------------------------------------
    
    # 3.1 构建DEseq2对象
    # library("DESeq2")
    dds <- DESeqDataSetFromMatrix(countData = exp2,
                                  colData = colData,
                                  design = ~ condition)
    
    # 3.2 DEseq2
    # 差异分析。
    dds <- DESeq(dds)

    
    # 标准化
    norm_count <- counts(dds, normalized=T)
    norm_count_mad <- apply(norm_count, 1, mad)
    norm_count <- norm_count[order(norm_count_mad, decreasing = T),]
    # write.csv( norm_count,'airway_dds_norm_count.csv')
    
    # 4.2 方差稳定变换，The variance stabilizing transformation(vst)
    vsd<-vst(dds)
    
    # plotPCA(vsd, intgroup=c('condition'))
    vsdmat <- assay(vsd) # 提取转化后的矩阵
    vsdmat <- vsdmat[order(norm_count_mad,decreasing = T),]
    # write.csv( vsdmat,'airway_dds_norm_count_vsd.csv')
    
    # 4.3 正则化对数变换，The regularized-logarithm transformation(rlog)
    rld <- rlog(dds,blind = F)
    # plotPCA(rld,intgroup=c('condition'))
    
    rlogmat <- assay(rld)
    rlogmat <- rlogmat[order(norm_count_mad,decreasing = T),]
    # write.csv( rlogmat,'airway_dds_norm_count_rlog.csv')
    
    # 五、提取DEseq2分析的结果 ---------------------------------------------------------
    
    result <-list()
    
    result$norm_count <- norm_count
    
    result$exp     <- exp
    result$exp2     <- exp2
    result$group   <- group
    
    result$dds     <- dds
    
    result$vsd     <- vsd
    result$rld     <- rld
    result$vsdmat  <- vsdmat
    result$rlogmat <- rlogmat

    return(result)
    
    
  })
  
# 选择实验组与对照组
observe({
  
     file1 <- input$file
    if(is.null(file1)){
      # 显示参考数据，分组
      df2_deseq <- df_sample2
    } 
    else{ 
      # 显示输入数据
      df2_deseq <- data.frame(read_excel(file1$datapath,2 ) ) 
        }
    
    df2_deseq <- df2_deseq[order(df2_deseq$condition),]
    condition <- df2_deseq$condition
    # 实验组与对照组情况
    updateSelectInput(session, "untrt",label = '对照组', choices = unique(condition) , selected = unique(condition)[1] )
    updateSelectInput(session, "trt",  label = '实验组', choices = unique(condition) , selected = unique(condition)[2] )
  } )
  
 # 选择实验组与对照组后提取差异结果
observeEvent(input$submit1, {
      
  myDEseq <- reactive({
    
    if(is.null(myDEseq1() )){return(NULL)}
    
    result <- myDEseq1()
    
    dds <- result$dds
 
    res <- results(dds, contrast = c("condition", input$trt, input$untrt ) ) 
   
    DEG1 <- as.data.frame(res) # 数据框形式

    # 5.2 排序、筛选
    # order() 给出从小到大排序后的位置(默认升序)
    resOrdered <- res[order(res$pvalue),] 
    
    DEG2 <- as.data.frame(resOrdered)
    # write.csv(x = resSig, file = "results.csv")  # 输出2
    
    result$res     <- res
    
    result$DEG1    <- DEG1
    result$DEG2    <- DEG2
    
    return(result)
  })

  
 # 1.2、DEseq2 PCA、热图与火山图
observeEvent(input$submit4, {
 
  plots <- reactive({
    
    if(is.null(myDEseq() ) ){return(NULL)} 
    
    result <- myDEseq()
    df1    <- result$norm_count
    group  <- result$group
    group  <- group$condition
    
    df_pca <- t(df1 ) # 画PCA图时要求是行名是样本名，列名是探针名，因此此时需要转换 t()
    df_pca <- as.data.frame(df_pca ) # 将 matrix转换为data.frame
    df_pca <- cbind(df_pca, group ) # cbind横向追加，即将分组信息追加到最后一列
    
    dat.pca <- PCA(df_pca[,-ncol(df_pca)], graph = FALSE)
    p_pca <- fviz_pca_ind(dat.pca,
                          geom.ind = "point", # show points only (nbut not "text")
                          col.ind = df_pca$group, # color by groups
                          
                          addEllipses = TRUE, # Concentration ellipses
                          legend.title = "Groups")
    
    # 2、vsd
    p_vsd <- plotPCA(result$vsd, intgroup=c('condition'))
    
    # 3、rld
    p_rld <- plotPCA(result$rld,intgroup=c('condition'))
    
    # 4、热图
    library(pheatmap)
    DEG2 <- result$DEG2
    exp2 <- result$vsdmat
    group <- result$group

    
    # 选取存在差异的基因
    DEG2$change <- as.factor(
      ifelse(
        DEG2$pvalue < input$pvalue & DEG2$padj < input$padj & abs(DEG2$log2FoldChange) >= input$logFC_cutoff,
        ifelse(
          DEG2$log2FoldChange >= input$logFC_cutoff,'UP','DOWN'),
        'NOT'))
    
    DEG2 <- DEG2[which(!DEG2$change=='NOT'),]


    # 设置热图分组
    t <- table(group$condition)
    
    rep <- vector()
    for (i in 1:length(t) ) {
      rep <- c( rep, 1:t[i] )
    }
    
    annotation_col = data.frame(
      group = factor(  group$condition ) , 
      rep = rep )
    rownames(annotation_col)<-colnames(exp2)
    # 选择基因

    # 按 Pvalue 排序
    choose_gene <- head(rownames(DEG2), input$heatmap_num)
    choose_matrix1 <- exp2[choose_gene, ]
    
    p_heatmap <- pheatmap(choose_matrix1,
                          annotation_col = annotation_col,
                          scale = 'row',
                          cluster_cols = F,
                          show_colnames =T,
                          show_rownames = T) 
    
    # 按 logFC 排序
    DEG2$FC <- abs(DEG2$log2FoldChange) 
    DEG2 <-  DEG2[order( DEG2$FC,decreasing = T),]
    
    choose_gene <- head(rownames(DEG2), input$heatmap_num)
    choose_matrix1 <- exp2[choose_gene, ]
    
    p_heatmap1 <- pheatmap(choose_matrix1,
                          annotation_col = annotation_col,
                          scale = 'row',
                          cluster_cols = F,
                          show_colnames =T,
                          show_rownames = T) 
    
    # 5、火山图
    DEG1 <- result$DEG1
    DEG1 <- na.omit(DEG1) # 去除缺失值
    
    # 设置阈值
    
    # 确定基因为上调、下调或无明显改变
    DEG1$change <- as.factor(
      ifelse(
        DEG1$pvalue < input$pvalue & DEG1$padj < input$padj & abs(DEG1$log2FoldChange) > input$logFC_cutoff,
        ifelse(
          DEG1$log2FoldChange > input$logFC_cutoff,'UP','DOWN'),
        'NOT'))
    
    table(DEG1$change) # 查看基因上、下调情况
    
    # 设置火山图的标题
    this_tile=paste('Cutoff for logFC is ',round(input$logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG1[DEG1$change=='UP',]),
                    '\nThe number of down gene is ',nrow(DEG1[DEG1$change=='DOWN',]))
    
    # 画火山图
    p_volcano <- ggplot(data=DEG1,
                        aes(x=log2FoldChange,y=-log10(pvalue),   #这里将pvalue取负对数
                            color= change)) +
      geom_point(alpha=0.4,size=1.75) +     #绘制点图
      theme_set(theme_set(theme_bw(base_size=20))) +
      xlab("log2 fold change")+ylab("-log10 pvalue") +    #轴标签
      ggtitle(this_tile)+theme(plot.title=element_text(size=15,hjust=0.5)) +
      scale_color_manual(values=c('blue','black','red'))   #设定颜色
    
    
    plots <- list()
    plots$p_pca     <- p_pca
    plots$p_vsd     <- p_vsd
    plots$p_rld     <- p_rld
    plots$p_heatmap <- p_heatmap
    plots$p_heatmap1 <- p_heatmap1
    plots$p_volcano <- p_volcano
    
    return(plots)
    
  } )
  
  # 5.1 DEseq2 
  output$plots <-  renderPlot({
    if (is.null(plots() ) ) { return() }
    
    plots <- plots()
    
    plots <- plots()
    
    if(input$plots_id1 == "pca") {
      p <- plots$p_pca
    }
    else if(input$plots_id1 == "vsd"){
      p <- plots$p_vsd
    } 
    else if(input$plots_id1 == "rld"){
      p <-  plots$p_rld
    } 
    else if(input$plots_id1 == "heatmap"){
      p <- plots$p_heatmap
    }
    else if(input$plots_id1 == "heatmap1"){
      p <- plots$p_heatmap1
    }
    else if(input$plots_id1 == "volcano"){  
      p <- plots$p_volcano
    }
    
    return(p)
    
  })
  
  # 下载图像 .Rdata
  output$downloadPlots <- downloadHandler(
    
    filename = function() {
      paste("plot.Rdata")
    },
    
    content = function(file) {
      
      if (is.null(plots() ) ) { return() }
      
      plots <- plots()
      
      if(input$plots_id1 == "pca") {
        p <- plots$p_pca
      }
      else if(input$plots_id1 == "vsd"){
        p <- plots$p_vsd
      } 
      else if(input$plots_id1 == "rld"){
        p <-  plots$p_rld
      } 
      else if(input$plots_id1 == "heatmap"){
        p <- plots$p_heatmap
      }
      else if(input$plots_id1 == "heatmap1"){
        p <- plots$p_heatmap1
      }
      else if(input$plots_id1 == "volcano"){ 
        p <- plots$p_volcano
        }
      save(p, file = file)
      
    } )
  
  # 下载图像 .png
  output$downloadPlots_png <- downloadHandler(
    
    filename = function() {
      paste("plot.png")
    },
    
    content = function(file) {
      
     if (is.null(plots() ) ) { return() }
      
      png(file,width=10,height=8,unit="in",res=150)
      
      plots <- plots()
      
      if(input$plots_id1 == "pca") {
        p <- plots$p_pca
      }
      else if(input$plots_id1 == "vsd"){
        p <- plots$p_vsd
      } 
      else if(input$plots_id1 == "rld"){
        p <-  plots$p_rld
      } 
      else if(input$plots_id1 == "heatmap"){
        p <- plots$p_heatmap
      }
      else if(input$plots_id1 == "heatmap1"){
        p <- plots$p_heatmap1
      }
      else if(input$plots_id1 == "volcano"){ 
        p <- plots$p_volcano
      }
      
      print(p)
      
      dev.off()
      
    } )
  
  
})
  
  plots_DEseq <- reactive({
    
    if(is.null(myDEseq1() )){return(NULL)}
    
    result <- myDEseq1() 
    
      # 选择三种标准化数据
      if(input$plots_id2 == "pca") {
        dat <- data.frame(result$norm_count)
      }
      else if(input$plots_id2 == "vsd"){
        dat <- data.frame(result$vsdmat)
      } 
      else if(input$plots_id2 == "rld"){
        dat <- data.frame(result$rlogmat)
      } 
    
    
    if(input$from_geneName==input$to_geneName){
      dat1 <- dat
      dat1 <- cbind('ID'=rownames(dat1),dat1)
      
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
 # }  )
    
    # 保存数据待用
    df_exp <- list()
    
    df_exp$group <- result$group
    df_exp$dat  <- dat
    df_exp$dat1 <- dat1
    
    return(df_exp)
    
  } )

  
  # # 四、分析结果 ----------------------------------------------------------------
  
  # 4.1 展示分析结果：分组
  # 4.2 展示分析情况：差异结果
  
  output$DEG <- renderDataTable({
    
    if (is.null(myDEseq() ) ) { return() }
    
    result <- myDEseq()
    
    dds   <- result$dds
    
    res <- results(dds, contrast = c("condition",input$untrt, input$trt ) )
    
    DEG <- as.data.frame(res) 
    DEG <- cbind('ID'=rownames(DEG),DEG)
    
    return(DEG)
  })
  
  
  
  # 4.3 下载计算结果：DEseq2
  output$downloadData <- downloadHandler(
    
    filename = function() {
      paste("DEseq2_result.RDS")
    },
    
    content = function(file) {
      
      result1 <- myDEseq()
      
      result <- list()
      
      result$group <- result1$group 
      result$DEG1 <- result1$DEG1
      
      result$norm_count <- result1$norm_count
      result$vsdmat <- result1$vsdmat
      result$rlogmat <- result1$rlogmat
      
      saveRDS(result, file = file)
      
    } )
  
  # 下载差异结果
  output$downloadDEG <- downloadHandler(
    
    filename = function() {
      paste("DEG.csv")
    },
    
    content = function(file) {
      
      result1 <- myDEseq()
      DEG1 <- result1$DEG1
      DEG1 <- cbind('ID'=rownames(DEG1),DEG1)
      write.csv(DEG1, file,row.names = F, fileEncoding = "GB18030") 
    } )
  
  # 下载标准化结果
  output$downloadnorm <- downloadHandler(
    
    filename = function() {
      paste("norm.csv")
    },
    
    content = function(file) {
      
      if(is.null(myDEseq1() )){return(NULL)}
      
      result <- myDEseq1() 
      
      # 选择三种标准化数据
      if(input$plots_id2 == "pca") {
        dat <- data.frame(result$norm_count)
      }
      else if(input$plots_id2 == "vsd"){
        dat <- data.frame(result$vsdmat)
      } 
      else if(input$plots_id2 == "rld"){
        dat <- data.frame(result$rlogmat)
      } 
      
      dat <- cbind('ID'=rownames(dat), dat)
      
      write.csv(dat, file,row.names = F , fileEncoding = "GB18030") 
      
    } )
  


  # # 五、网页展示图形 --------------------------------------------------------------
  
  

  
  # # 三、参考数据 ----------------------------------------------------------------
  
  # 3.1 展示参考数据
  output$sample1_exp_deseq <- renderDataTable({
    
    df_exp <- plots_DEseq()
    dat <- df_exp$dat
    dat <- cbind('ID'=rownames(dat),dat)
    
    if (!is.null(plots_DEseq() ) ) { 
      return(dat ) 
    } else if(is.null( plots_DEseq() ) ){
      
      result <- readRDS("/srv/shiny-server/RNAseq/file/result.RDS")
      
      # 选择三种标准化数据
      if(input$plots_id2 == "pca") {
        dat <- data.frame(result$norm_count)
      }
      else if(input$plots_id2 == "vsd"){
        dat <- data.frame(result$vsdmat)
      } 
      else if(input$plots_id2 == "rld"){
        dat <- data.frame(result$rlogmat)
      } 
      
      return(dat)
}
    
  }) 
  

  # 4.2 转换 ID 后的结果
  observeEvent(input$submit2, {

    output$DEG_degseq <- renderDataTable({
      
      if (is.null(plots_DEseq() ) ) { return(NULL ) }
      
      df_exp <- plots_DEseq()
      dat1 <- df_exp$dat1
      
      return( dat1  )
      
    } ) 
    
    
  })
  
  
  # 基因表达箱图
  observe({
    if(is.null(plots_DEseq() ) ) {return(NULL) }  
    
    df_exp <-  plots_DEseq()
    dat1   <- df_exp$dat1

    if(input$from_geneName==input$to_geneName){
      
      dat1 <- dat1[, setdiff(colnames(dat1), "ID" ) ]
      dat1 <- as.matrix(dat1)
      
    } else(
      dat1   <- avereps(dat1[ , setdiff(colnames(dat1),c(input$to_geneName, "ID") ) ], ID = dat1[,input$to_geneName] ) 
    )
    
    updateSelectInput(session, "gene",label = '基因',choices = rownames(dat1)  )
  })
  
  observeEvent(input$submit3, {
    if(is.null(plots_DEseq() ) ) {return(NULL) }  
    
    df_exp <-  plots_DEseq()
    dat1   <- df_exp$dat1
    group  <- df_exp$group
    group  <- group$condition
    
    if(input$from_geneName==input$to_geneName){
      
      dat1 <- dat1[, setdiff(colnames(dat1), "ID" ) ]
      dat1 <- as.matrix(dat1)
      
    } else(
      dat1   <- avereps(dat1[ , setdiff(colnames(dat1),c(input$to_geneName, "ID") ) ], ID = dat1[,input$to_geneName] ) 
    )
    
    output$plots_DEseq <- renderPlot({
      
      p_exp <- gene_exp(data = dat1, name = input$gene , group = group) 
      return(p_exp)
      
    })
    
  })
  
  })


# # 三、参考数据 ----------------------------------------------------------------

# 3.1 展示参考数据

# 3.1.1、展示 DEseq2 参考数据

# 3.1.1.1 DEseq2 矩阵 
output$sample1_exp <- renderDataTable({
  file1 <- input$file
  if(is.null(file1)){
     df1 <- df_sample1
     return(head(df1, 2000) )
    
    }
  return(data.frame(read_excel(file1$datapath,1 ) ))

}) 

# 3.1.1.2 DEseq2 分组 
output$sample2_group <- renderDataTable({
  
  file1 <- input$file
  if(is.null(file1)){
    return( df_sample2 )
  }
  return( data.frame( read_excel(file1$datapath,2 ) ) )
  
}) 
  
  # 3.2、下载参考数据
  
  # 3.2.1 下载参考数据  DEseq2
  output$downloadSampleData <- downloadHandler(
    filename = function() {
      paste('DEseq2参考数据.xlsx')
    },
    content = function(file) {
      
      sample1 <- head( df_sample1, input$sample_exp_row)
      sample2 <- df_sample2
      
      write.xlsx(sample1, file, sheetName = '表达矩阵示例', row.names = F, showNA = F)
      write.xlsx(sample2, file, sheetName = '分组信息', append=TRUE, row.names = F, showNA = F)
    }
  )
  
    
}

# # ④ shinyApp -----------------------------------------------------------------
shinyApp(ui, server)


