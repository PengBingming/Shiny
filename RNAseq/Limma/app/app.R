
# # ① 加载包 -------------------------------------------------------------------
if (!require('shinydashboard', quietly = TRUE)) BiocManager::install('shinydashboard')
if (!require('shiny', quietly = TRUE)) BiocManager::install('shiny')

library(clusterProfiler) # 转换 ID
if (!require('limma', quietly = TRUE)) BiocManager::install('limma')
if (!require("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!require("org.Mm.eg.db", quietly = TRUE)) BiocManager::install("org.Mm.eg.db")

if (!require("readxl", quietly = TRUE)) BiocManager::install("readxl")
if (!require("xlsx", quietly = TRUE)) BiocManager::install("xlsx")

if (!require("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
library("FactoMineR")#画主成分分析图需要加载这两个包
library("factoextra") 

# 基因表达箱图函数
gene_exp = function(data, name, group){         #自定义一个函数bp，函数为{}里的内容
  library(ggpubr)
  df = data.frame(gene = data[name,], stage = group)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter", ylab = name) +
    #  Add p-value
    # stat_compare_means()
    stat_compare_means(method = "anova")
  
  return(p)
}

# # ② 设置shiny界面参数 ui -----------------------------------------------------------
# # ②.1页面标题 ------------------------------------------------------------------
header <- dashboardHeader(title = "RNAseq 数据分析")

geneName <- intersect(columns(org.Mm.eg.db),columns(org.Hs.eg.db) )
# # ②.2 侧边界面 ------------------------------------------------------------------
sidebar <- dashboardSidebar(
  fileInput(inputId = "file_limma", "标准化矩阵（Limma）：",
            multiple = TRUE,),
  tags$hr(),
  numericInput(inputId = "rowSum_filter",
               label = "过滤低表达参数：rowSums(exp) >",
               value = 1) ,
  tags$hr(),
  selectInput("untrt_limma", "对照组：", c(" ") ),
  selectInput("trt_limma",   "实验组：", c(" ") ),
  tags$hr(),
  radioButtons("species", "物种：",
               choices = c('人' = "hsa",
                           '鼠' = "mmu"),
               selected = "hsa"),
  tags$hr(),
  selectInput("from_geneName","基因名：form", geneName, selected = "ENTREZID" ),
  selectInput("to_geneName","基因名：to", geneName , selected = "ENTREZID" ),
  tags$hr()
)


# # ②.3 网页呈现的内容  --------------------------------------------------------------
body <- dashboardBody( 
  
  # tabsetPanel 1 limma
  tabsetPanel(
    tabPanel('数据格式',
             fluidRow( 
               column(width = 2, 
                      box(title = '下载参考数据',
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
                        solidHeader = T) 
                      ),
               column(width = 3 , 
                      box(
                        title = "分组：",
                        width = NULL, 
                        height = 10, 
                        status = "primary",
                        solidHeader = T) )
               ),
             tags$br(),
             splitLayout(cellWidths = c("18%", "60%","2","20%"),
                         box(width = 8 ,
                             downloadButton("downloadSampleData_limma", "参考数据"),
                             tags$br(),
                             tags$br(),
                             numericInput(inputId = "sample_exp_row_limma",
                                          label = "下载矩阵行数：",
                                          value = 5000),
                             solidHeader = T),
                         dataTableOutput("sample1_exp_limma"), 
                         '',
                         dataTableOutput("sample2_group_limma") ) 
             ), # tabPanel 2.3
    # 2.2 分析结果
    tabPanel('分析结果',
             downloadButton("downloadData_limma", "Limma 分析结果"),
             tags$br(),
             tags$br(),
             fluidRow( 
               column(width = 12, 
                      box( 
                        title = "差异结果：",
                        width = NULL, 
                        height = 10, 
                        status = "primary",
                        solidHeader = T ) )
               ),
             tags$br(),
             splitLayout(cellWidths = c("100%"),
                         dataTableOutput("DEG_limma") ) ) # tabPanel 2.2
    ,
    # 2.3 Limma 各类图形
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
               column(width = 12,
                      box( width = 2, 
                           height = 5, 
                           solidHeader = T ,
                           radioButtons("plots_id_limma", "图形选择：",
                                        choices = c('PCA'    ='pca',
                                                    '热图：Pvalue'   = "heatmap",
                                                    '热图：logFC'   = "heatmap1",
                                                    '火山图' = "volcano"),
                                        selected = "heatmap") ),
                      box( 
                        width = 2, 
                        height = 5, 
                        solidHeader = T ,
                        numericInput(inputId = "heatmap_num_limma",
                                     label = "热图基因数：",
                                     value = 50) ),
                      box(
                        width = 2, 
                        height = 5,
                        solidHeader = T ,
                        numericInput(inputId = "pvalue_limma",
                                     label = "火山图 pvalue：",
                                     value = 0.05)),
                      box(
                        width = 2, 
                        height = 5, 
                        solidHeader = T ,
                        numericInput(inputId = "padj_limma",
                                     label = "火山图 padj：",
                                     value = 0.1) ),
                      box(
                        width =2, 
                        height = 5,
                        solidHeader = T ,
                        numericInput(inputId = "logFC_cutoff_limma",
                                     label = "火山图 logFC：",
                                     value = 1.5) )
                      ) 
               ),
             tags$br(),
             tags$br(),
             tags$br(),
             splitLayout(cellWidths = c("25%","50%","25%"),
                         '',
                         plotOutput("plots_limma"),'' ) )  , # tabPanel 2.3

      tabPanel('标准化的数据',
               fluidRow(
                 column(width = 12, 
                        box(
                          title = "矩阵：",
                          width = NULL, 
                          height = 10, 
                          status = "primary",
                          solidHeader = T) 
                 )
               ),
               tags$br(),
               splitLayout(cellWidths = c("100%"),
                           dataTableOutput("sample1_exp_deseq")) 
      ), # tabPanel 1.1
      
      # 1.2 分析结果
      tabPanel('基因 ID 转换',
               fluidRow( column(width = 12, 
                                box(title = "基因 ID 转换：",
                                    width = NULL, 
                                    height = 10, 
                                    status = "primary",
                                    solidHeader = T ) )
               ),
               tags$br(),
               splitLayout(cellWidths = c("100%"),
                           dataTableOutput("DEG_degseq") ) ) # tabPanel 1.2
      ,
      # 1.3 图形
      tabPanel('基因表达情况',
               fluidRow( column(width = 12, 
                                box( 
                                  title = "基因表达情况：",
                                  width = NULL, 
                                  height = 10, 
                                  status = "primary",
                                  solidHeader = T ) )
               ),
               tags$br(),
               fluidRow( 
                 column(width = 12,
                        box( width = 2, 
                             height = 5, 
                             solidHeader = T,
                             selectInput("gene","基因",c(" ") ) )
                 ) 
               ),
               tags$br(),
               tags$br(),
               tags$br(),
               splitLayout(cellWidths = c("25%","75%"),
                           '' , plotOutput("plots_DEseq") ) ) # tabPanel 1.3
                  
) # tablesetpanel
  
) # dashboardBody

# # ②.4 ui -----------------------------------------------------------------
ui <- dashboardPage(header, sidebar, body)


# # ③ 设置 server -------------------------------------------------------------
server <- function(input, output, session) {
  
  # # 二、Limma 分析函数 ----------------------------------------------------------
  
  # 1.2 Limma 差异分析 
  
  # 4.2 展示分析情况：差异结果
  observe({
    file2 <- input$file_limma
    if(is.null(file2)){return(NULL)}
    
    df2 <- data.frame(read_excel(file2$datapath,2 ) )
    
    df2 <- df2[order(df2$condition),]
    group <- df2$condition
    
    updateSelectInput(session, "untrt_limma",label = '对照组', choices = levels(factor(group)), selected = levels(factor(group))[1] )
    updateSelectInput(session, "trt_limma",  label = '实验组', choices = levels(factor(group)), selected = levels(factor(group))[2] )
 
  myLimma <- reactive( {   
    file2 <- input$file_limma
    if(is.null(file2)){return(NULL)}
    
    df1 <- data.frame(read_excel(file2$datapath,1 ) )
    df2 <- data.frame(read_excel(file2$datapath,2 ) )
    
    df2 <- df2[order(df2$condition),]
    df1 <- df1[,c("ID",df2$group)]
    
    lis <- strsplit(df1$ID,'[.]')
    for (i in 1:length(lis) ) {
      
      df1$symbol[i] <- lis[[i]][1]
      
    } ; rm(i)
    
    df1 <- avereps(df1[,c(-1,-ncol(df1))],ID=df1$symbol)
    
    df1 <- df1[,colSums(df1) > 0 ]
    df1 <- df1[rowSums(df1) > input$rowSum_filter , ]
    df1 <- as.matrix(df1)
    
    group <- df2$condition
    # 差异分析走标准的limma流程 ---------------------------------------------------------
    
    #创建一个分组的矩阵
    design = model.matrix(~0+factor(group))#创建一个分组的矩阵
    colnames(design) = levels(factor(group))
    rownames(design) = colnames(df1)
    
    # dematrix <- as.data.frame(design)
    
    # 创建差异比较矩阵 ----------------------------------------------------------------
    
    # 这个矩阵声明，我们要 实验组 和 对照组 进行差异分析比较
    # contrast.matrix <- makeContrasts('HDM-Saline',levels = design)
    
    levels <- levels( factor(group) )
    
    contrast.matrix <- matrix( unique( design[, input$trt_limma ] - design[, input$untrt_limma ] ), 
                               dimnames = list(Levels = levels, 
                                               Contrasts = paste0(input$trt_limma , ' - ', input$untrt_limma) ) )
    # 第一步lmFit，#lmFit为每个基因给定一系列的阵列来拟合线性模型
    fit<-lmFit(df1,design)
    
    # 第二步eBayes，#eBayes给出了一个微阵列线性模型拟合，通过经验贝叶斯调整标准误差到一个共同的值来计算修正后的t统计量、修正后的f统计量和微分表达式的对数概率。
    fit1<- contrasts.fit(fit, contrast.matrix)
    fit1<- eBayes(fit1)
    
    # 第三步topTable,#topTable从线性模型拟合中提取出排名靠前的基因表。
    options(digits = 4) #设置全局的数字有效位数为4
    
    # topTable(fit1,coef=2,adjust='BH') 
    tempOutput <- topTable(fit1, coef=1, n=Inf) 
    tempOutput <- na.omit(tempOutput)#移除NA值
    
    # 上调基因和下调基因 ---------------------------------------------------------------
    
    
    
    # 
    # tempOutput$gene=ifelse(tempOutput$P.Value > 0.01,'stable', #if 判断：如果这一基因的P.Value>0.01，则为stable基因
    #                        ifelse( tempOutput$logFC > 1.2,'up', #接上句else 否则：接下来开始判断那些P.Value<0.01的基因，再if 判断：如果logFC >1.5,则为up（上调）基因
    #                                ifelse( tempOutput$logFC < -1.2,'down','stable') )#接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
    # )

    
    tem <- tempOutput[order(tempOutput$P.Value),]
    
    dat <- data.frame(df1)

    if(input$from_geneName==input$to_geneName){
      
      dat1 <- dat

      updateSelectInput(session, "gene",label = '基因',choices = rownames(dat1)  )
       
      dat1 <- as.matrix(dat1)

    } 
    else if(!input$from_geneName==input$to_geneName){
      
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
      dat <- as.matrix(dat1)  

      dat1   <- avereps(dat1[ , setdiff(colnames(dat1),c(input$to_geneName, "ID") ) ], ID = dat1[,input$to_geneName] )
      updateSelectInput(session, "gene",label = '基因',choices = rownames(dat1)  )
      
   
    } 
    
    # 输出结果
    result <- list()
    result$DEG1 <- tempOutput
    result$DEG2 <- tem
    result$group <- group
    result$exp <- df1
    
    result$df2 <- df2
    result$dat  <- dat
    result$dat1 <- dat1
    
    return(result)
    
  } )

  plots_limma <- reactive({

  result <- myLimma()
  if(is.null(myLimma() ) ){return(NULL)}

  df1 <- result$exp
  group <- result$group

  df_pca <- t(df1 ) # 画PCA图时要求是行名是样本名，列名是探针名，因此此时需要转换 t()
  df_pca <- as.data.frame(df_pca ) # 将 matrix转换为data.frame
  df_pca <- cbind(df_pca,group ) # cbind横向追加，即将分组信息追加到最后一列

  dat.pca <- PCA(df_pca[,-ncol(df_pca)], graph = FALSE)
  p_pca <- fviz_pca_ind(dat.pca,
                        geom.ind = "point", # show points only (nbut not "text")
                        col.ind = df_pca$group, # color by groups

                        addEllipses = TRUE, # Concentration ellipses
                        legend.title = "Groups"
  )

  # 3、热图
  library(pheatmap)
  DEG2 <- result$DEG2
  exp2 <- result$exp
  group <- result$df2

  # 选取存在差异的基因
  DEG2$change <- as.factor(
    ifelse(
      DEG2$P.Value < input$pvalue_limma & DEG2$adj.P.Val < input$padj_limma & abs(DEG2$logFC) >= input$logFC_cutoff_limma,
      ifelse(
        DEG2$logFC >= input$logFC_cutoff_limma,'UP','DOWN'),
      'NOT'))

  DEG2 <- DEG2[which(!DEG2$change=='NOT'),]

  # 设置热图分组
  t <- table(group$condition)

  rep <- vector()
  for (i in 1:length(t) ) {
    rep <- c( rep, 1:t[i] )
  }

  annotation_col = data.frame(
    group = factor( group$condition ) ,
    rep = rep )
  rownames(annotation_col)<-colnames(exp2)

  # 选择基因

  choose_gene <- head(rownames(DEG2), input$heatmap_num_limma)
  choose_matrix1 <- exp2[choose_gene, ]

  p_heatmap <- pheatmap(choose_matrix1,
                        annotation_col = annotation_col,
                        scale = 'row',
                        cluster_cols = F,
                        show_colnames =T,
                        show_rownames = T)

  # 按 logFC 排序
  DEG2$FC <- abs(DEG2$logFC) 
  DEG2 <-  DEG2[order( DEG2$FC,decreasing = T),]
  
  choose_gene <- head(rownames(DEG2), input$heatmap_num_limma)
  choose_matrix1 <- exp2[choose_gene, ]
  
  p_heatmap1 <- pheatmap(choose_matrix1,
                        annotation_col = annotation_col,
                        scale = 'row',
                        cluster_cols = F,
                        show_colnames =T,
                        show_rownames = T)
  
  # 4、火山图
  DEG1 <- result$DEG1
  DEG1 <- na.omit(DEG1) # 去除缺失值

  # 设置阈值

  # 确定基因为上调、下调或无明显改变
  DEG1$change <- as.factor(
    ifelse(
      DEG1$P.Value < input$pvalue_limma & DEG1$adj.P.Val < input$padj_limma & abs(DEG1$logFC) > input$logFC_cutoff_limma,
      ifelse(
        DEG1$logFC > input$logFC_cutoff_limma,'UP','DOWN'),
      'NOT'))

  table(DEG1$change) # 查看基因上、下调情况

  # 设置火山图的标题
  this_tile=paste('Cutoff for logFC is ',round(input$logFC_cutoff_limma,3),
                  '\nThe number of up gene is ',nrow(DEG1[DEG1$change=='UP',]),
                  '\nThe number of down gene is ',nrow(DEG1[DEG1$change=='DOWN',]))

  # 画火山图
  p_volcano <- ggplot(data=DEG1,
                      aes(x = logFC, y = -log10(P.Value ),   #这里将pvalue取负对数
                          color= change)) +
    geom_point(alpha=0.4,size=1.75) +     #绘制点图
    theme_set(theme_set(theme_bw(base_size=20))) +
    xlab("log2 fold change")+ylab("-log10 pvalue") +    #轴标签
    ggtitle(this_tile)+theme(plot.title=element_text(size=15,hjust=0.5)) +
    scale_color_manual(values=c('blue','black','red'))   #设定颜色


  plots <- list()
  plots$p_pca     <- p_pca
  plots$p_heatmap <- p_heatmap
  plots$p_heatmap1 <- p_heatmap1
  plots$p_volcano <- p_volcano

  return(plots)

} )


  output$DEG_limma <- renderDataTable({

    result <- myLimma()
    DEG <- result$DEG1
    DEG <- cbind(ID=rownames(DEG),DEG)

    return( DEG  )

  })


  output$plots_DEseq <- renderPlot({
      
      if(is.null(myLimma() ) ) {return(NULL) }
      
      result <-  myLimma()
      dat1   <- result$dat1
      group  <- result$group
      p_exp <- gene_exp(data = dat1, name = input$gene , group = group)
      
     return(p_exp)
  })


  output$sample1_exp_deseq <- renderDataTable({

    results <- myLimma()
    dat <- results$dat
    dat <- cbind('ID'=rownames(dat),dat)

    if (!is.null(myLimma() ) ) { return(dat ) }
    sample1 <- readRDS("/srv/shiny-server/RNAseq/file/result.RDS")

    sample1 <- sample1$vsdmat

    sample1 <- cbind("ID"=rownames(sample1), sample1)

    return(sample1)

  })


  # 4.2 转换 ID 后的结果
  output$DEG_degseq <- renderDataTable({

    if (is.null(myLimma() ) ) { return(NULL ) }

    result <- myLimma()
    dat <- result$dat
    dat <- cbind('ID'=rownames(dat),dat)
    
    return( dat )

  } )


  # 4.4 下载计算结果：Limma
  output$downloadData_limma <- downloadHandler(

    filename = function() {
      paste("Limma 处理结果.RDS")
    },

    content = function(file) {

      result1 <- myLimma()

      result <- list()

      result$DEG1 <- result1$DEG1
      result$DEG2 <- result1$DEG2
      result$exp <- result1$exp
      result$df2 <- result1$df2
      result$group <- result1$group

      saveRDS(result, file = file)

    } )

  # 5.2 Limma
  output$plots_limma <-  renderPlot({
    if (is.null(plots_limma() ) ){ return() }

    plots <- plots_limma()

    # 1、pca
    p_pca <- plots$p_pca

    # 2、heatmap
    p_heatmap <- plots$p_heatmap
    p_heatmap1 <- plots$p_heatmap1
    
    # 3、volcano
    p_volcano <- plots$p_volcano

    if(input$plots_id_limma == "pca") {
      return(p_pca)
    }
    else if(input$plots_id_limma == "heatmap"){
      return(p_heatmap)
    }
    else if(input$plots_id_limma == "heatmap1"){
      return(p_heatmap1)
    }
    else {p_volcano }
  
    })


  # 2.2 Limma PCA、热图与火山图


  # # 三、参考数据 ----------------------------------------------------------------

  # 3.1.2、展示 Limma 参考数据
  # 3.1.2.1 Limma 矩阵
  output$sample1_exp_limma <- renderDataTable({

    myLimma <- myLimma()
    df1 <- myLimma$exp
    sample1 <- cbind(rownames(df1), df1)

    if(!is.null(myLimma() ) ){ 
      return(sample1 ) 
      }
    else if(is.null(myLimma() ) ){
            library(xlsx)
            sample1 <- as.data.frame(read_excel('/srv/shiny-server/RNAseq/file/GSE112491_limma_log2normalized_geneCounts.xlsx',1) )
            return(sample1)
    }
    
  })

  # 3.1.2.2 Limma 分组
  output$sample2_group_limma <- renderDataTable({

    myLimma <- myLimma()
    sample2 <- myLimma$df2

    if(!is.null(sample2 ) ){ 
       return(sample2 ) 
      }
    else if(is.null(sample2 ) ){
          library(xlsx)
          sample2 <- as.data.frame( read_excel('/srv/shiny-server/RNAseq/file/GSE112491_limma_log2normalized_geneCounts.xlsx',2) )
          return(sample2)
    }
  })

  # 3.2.2 下载参考数据 limma
  output$downloadSampleData_limma <- downloadHandler(
    filename = function() {
      paste('Limma参考数据.xlsx')
    },
    content = function(file) {

      sample1 <- head( as.data.frame( read_excel('/srv/shiny-server/RNAseq/file/GSE112491_limma_log2normalized_geneCounts.xlsx',1) ), input$sample_exp_row_limma)
      sample2 <- as.data.frame( read_excel('/srv/shiny-server/RNAseq/file/GSE112491_limma_log2normalized_geneCounts.xlsx',2) )

      write.xlsx(sample1, file, sheetName = '表达矩阵示例', row.names = F, showNA = F)
      write.xlsx(sample2, file, sheetName = '分组信息', append=TRUE, row.names = F, showNA = F)
    }
  )
  })

 
  # 3.2.2 下载参考数据 limma
  output$downloadSampleData_limma <- downloadHandler(
    filename = function() {
      paste('Limma参考数据.xlsx')
    },
    content = function(file) {
      
      sample1 <- head( as.data.frame( read_excel('/srv/shiny-server/RNAseq/file/GSE112491_limma_log2normalized_geneCounts.xlsx',1) ), input$sample_exp_row_limma)
      sample2 <- as.data.frame( read_excel('/srv/shiny-server/RNAseq/file/GSE112491_limma_log2normalized_geneCounts.xlsx',2) )
      
      write.xlsx(sample1, file, sheetName = '表达矩阵示例', row.names = F, showNA = F)
      write.xlsx(sample2, file, sheetName = '分组信息', append=TRUE, row.names = F, showNA = F)
    }
  )
  
}

# # ④ shinyApp -----------------------------------------------------------------
shinyApp(ui, server)


