
library(shiny)
library(shinydashboard)

library(ggplot2)
library(ggpubr)

library(showtext) # 解决画图乱码
showtext_auto()

library(rstatix) # 添加 P 值

library(readxl)

# Create the UI
header <- dashboardHeader(title = h5("单因素方差分析/秩和检验")) # 页面标题
sidebar <- dashboardSidebar(
  h5("支持格式：.csv .xlsx .xls "),
  fileInput("file", "输入文件："),
  
  actionButton("submit", "Analyze Data"),
  
  hr(),
  selectInput(inputId = 'num',label = '样本组数：',
              c("两组" = 'two',
                "多组" = 'multi'),
              selected = 'multi'
  ),
  selectInput(inputId = 'method', label = '检验方法：',
              c("t检验/方差分析" = 't_anova',
                "秩和检验" = 'wilcox'),
              selected = 't_anova'
  ),

  selectInput(inputId = 'pair',label = '是否配对：', c("配对" = 'T', "独立" = 'F'),selected = 'F' ),
  
  selectInput(inputId = 'alt',label = '单边/双边检验：', c("two.sided" = 'two.sided', "less" = 'less',"greater"="greater"),
              selected = 'two.sided' ),
  
  selectInput(inputId = 'padj',label = '校正方法：', 
              c("Bonferroni" = 'bonferroni', 
                "Holm" = 'holm',
                'hochberg'='hochberg',
                'BH'='BH',
                'BY'='BY',
                'fdr'='fdr',
                'hommel'='hommel',
                "不校正"="none"),
              selected = 'bonferroni' ),
  h6("（多组存在差异时，两两比较的校正方法）"),
  hr()
  
)

body <- dashboardBody( 
  
  tabsetPanel(
    tabPanel(h4("数据"),
             fluidRow( 
               column(width =12,  
                      box( 
                        title = "格式 : value 列为对应数值，group 列对应分组",
                        width = NULL, 
                        height = 10, 
                        status = "primary",
                        solidHeader = T ),
                      box( 
                        h6("两组时只能有两个分组（ group 列中只有两种情况），
                            进行 t 检验或者秩和检验；
                            多组时（ group 列中有多种情况），
                           进行方差分析或者秩和检验"),
                        width = NULL, 
                        height = 0, 
                        solidHeader = T ) )
              
             ),
             tags$br(),
             fluidRow( 
               column(width = 2,
                      tags$br(),
                      downloadButton("downloadSampleData", "下载参考数据") ) 
             ),
             tags$br(),
             splitLayout(cellWidths = c("30","30%"),'',
                         dataTableOutput("sample_table")) ),
    tabPanel(h4("统计表"),
             fluidRow( 
               column(width = 12, 
                      box( 
                        title = "结果",
                        width = NULL, 
                        height = 10, 
                        status = "primary",
                        solidHeader = T ) )
             ),
             tags$br(),
             fluidRow( 
               column(width = 2,
                      tags$br(),
                      downloadButton("downloadcontent1", "下载检验结果") ),
               column(width = 4),
               column(width = 3,
                      tags$br(),
                      downloadButton("downloadcontent2", "下载两两比较结果") )
             ),
             tags$br(),
             splitLayout(cellWidths = c("49%","2%","49%"),
                         dataTableOutput("test_table"),"",
                         dataTableOutput("stat.test")) ),
    
    tabPanel(h4("统计图"),
             fluidRow(column(width =12, 
                             box( 
                               title = "图形",
                               width = NULL, 
                               height = 10, 
                               status = "primary",
                               solidHeader = T ) )
             ),
             tags$br(),
             fluidRow( 
               column(width = 1,
                      tags$br(),
                      downloadButton("downloadplot_pdf", ".pdf") 
                      ),
               column(width = 1,
                      tags$br(),
                      downloadButton("downloadplot_png", ".png") 
                      ),
               column(width = 1,
                      numericInput(inputId = 'w',label = '下载图形长',value = 10) ),
               column(width = 1,
                      numericInput(inputId = 'h',label = '下载图形高',value = 8) ),
               column(width = 1,
               ),
               column(width = 2,
                      textInput(inputId = 'group', label = 'x轴标签：', value ='group') ),
               column(width = 2,
                      textInput(inputId = 'value', label = 'y轴标签：', value ='value' ) ),
               column(width = 2,
                      textInput(inputId = 'title', label = '图形标题：',value ="title") )
             ),
             fluidRow( 
               column(width = 1,
                      numericInput(inputId = 'point',label = '散点：', value = 0.5,  min = 0,max = 2)
                      ),
               column(width = 1,
                      numericInput(inputId = 'errorbar', label = '误差条：', value = 0.5,   min = 0,  max = 2)
                      ),
               column(width = 1,
                      numericInput(inputId = 'box', label = '箱图线条：', value = 0.5,min = 0,max = 2 )
                      ),
               column(width = 2,
                      numericInput(inputId = 'violin', label = '小提琴图线条：', value = 0.5, min = 0, max = 2 )
               ),
               column(width = 2,
                      selectInput(inputId = 'star', label = '显示横杠',c('显示'='T','隐藏'="F"),selected = "T" )
               ),
               column(width = 2,
                      numericInput(inputId = 'star_size', label = '星号',max = 20,min = 0,value = 5 )
               )
               ),
             
             tags$br(),
             tags$br(),
             splitLayout(cellWidths = c("90%"),
                         plotOutput("plot")),
    )
  )
  
)

ui <- dashboardPage(header, sidebar, body)
# Create the server
server <- function(input, output,session) {
  
  # Load the data
  # 读取数据
  data <- reactive({
    file1 <- input$file
    if ( is.null(file1) ){
      
      if(input$num=='two'){
        data <- iris[,c(1,5)]
        colnames(data) <- c('value', 'group')
        data <- data[which(!data$group=='setosa'),]
        
      }
      else if(input$num=='multi'){
        data <- iris[,c(1,5)]
        colnames(data) <- c('value', 'group')
      }

      
    } 
    else{
      
      d <- tail( unlist(strsplit(file1$datapath,'[.]') ), 1)
      
      if( d=='csv' ){
        
        data <- data.frame( read.csv(file1$datapath,header=T, stringsAsFactors = FALSE, fileEncoding = 'GB18030') )
      } else{
        data <- data.frame( read_excel(file1$datapath,1) )
      } 
      
    } # else
    
    return( data )
  })
  
  
  observeEvent(input$submit, {
    
    result  <- reactive({
      # Perform ANOVA
      # 方差分析
      data = data() 
      data$group <- factor(data$group)
      
      if(input$num=='two'){
        
        if(input$method=='t_anova'){
          
        t.test <- t.test( value ~ group,data = data , piar = c(input$pair=='T') ,alternative = input$alt)
        test_table <- t(as.data.frame( unlist(t.test)[1:3]))
        stat.test <- test_table
        
        # 画图
        plot <- ggplot(data=data, aes(x= group, y= value, group = group, color= group ) ) +
          geom_boxplot(alpha = 1, size = input$box ) +
          stat_boxplot(geom = 'errorbar' , size = input$errorbar, width = .3) +
          geom_violin(aes( (fill= group),col=group), alpha = .1, size=input$violin ) +
          geom_jitter(size = input$point, alpha= .5 , position =  position_jitter(width = .1) ) +
          ggtitle(input$title) +
          labs(x= input$group ,y= input$value)+
          theme_classic(base_size =  18) +
          stat_compare_means(method = 't.test' ) + # 添加 P 值。 library(ggpubr)
          theme(legend.position = "none", plot.title = element_text(size = 13, hjust = 0.5) ) + 
          scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) #调整y轴
      }
        
        else if(input$method=='wilcox'){
          wilcox.test <- wilcox.test( value ~ group,data = data, piar= c(input$pair=='T'),alternative = input$alt )
          test_table <- t(as.data.frame( unlist(wilcox.test)[1:2]))
          stat.test <- test_table
          # 画图
          plot <- ggplot(data=data, aes(x= group, y= value, group = group, color= group ) ) +
            geom_boxplot(alpha = 1, size = input$box ) +
            stat_boxplot(geom = 'errorbar' , size = input$errorbar, width = .3) +
            geom_violin(aes( (fill= group),col=group), alpha = .1, size=input$violin ) +
            geom_jitter(size = input$point, alpha= .5 , position =  position_jitter(width = .1) ) +
            ggtitle(input$title) +
            labs(x= input$group ,y= input$value)+
            theme_classic(base_size =  18) +
            stat_compare_means( ) + # 添加 P 值。 library(ggpubr)
            theme(legend.position = "none", plot.title = element_text(size = 13, hjust = 0.5) ) + 
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) #调整y轴
        }
      }
      else if(input$num=='multi'){
        
        if(input$method=='t_anova'){
          
          aov1 <- aov( value ~ group,data = data )
          
          # 统计结果
          test_table <- summary(aov1)[[1]]
          test_table<- cbind('Item' =rownames(test_table), test_table)  
          
          stat.test <- pairwise_t_test(data,value~group, p.adjust.method =  input$padj, paired = (input$pair=="T")) %>% add_y_position()
          logi <- (min(stat.test$p.adj) < 0.05)
          
          # 画图
          plot <- ggplot(data=data, aes(x= group, y= value, group = group, color= group ) ) +
            geom_boxplot(alpha = 1, size = input$box ) +
            stat_boxplot(geom = 'errorbar' , size = input$errorbar, width = .3) +
            geom_violin(aes( (fill= group),col=group), alpha = .1, size=input$violin ) +
            geom_jitter(size = input$point, alpha= .5 , position =  position_jitter(width = .1) ) +
            ggtitle(input$title) +
            labs(x= input$group ,y= input$value)+
            theme_classic(base_size =  18) +
            stat_compare_means(method = 'anova' ) + # 添加 P 值。 library(ggpubr)
            theme(legend.position = "none", plot.title = element_text(size = 13, hjust = 0.5) ) + 
            stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0.01,
                               hide.ns =  logi ,remove.bracket = (input$star=='F') , size = input$star_size) + #隐藏无意义的p值
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) #调整y轴
        }
        
        else if(input$method=='wilcox'){
          
          test_table <- kruskal.test(data= data, value~group) %>% unlist() %>% as.data.frame()
          test_table<- cbind('Item' =rownames(test_table), test_table)  
          
          stat.test <- data %>% pairwise_wilcox_test(value~group, p.adjust.method = input$padj, paired = (input$pair=='T')) %>% add_y_position()
          logi <- (min(stat.test$p.adj) < 0.05)
          
          # 画图
          plot <- ggplot(data=data, aes(x= group, y= value, group = group, color= group ) ) +
            geom_boxplot(alpha = 1, size = input$box ) +
            stat_boxplot(geom = 'errorbar' , size = input$errorbar, width = .3) +
            geom_violin(aes( (fill= group),col=group), alpha = .1, size=input$violin ) +
            geom_jitter(size = input$point, alpha= .5 , position =  position_jitter(width = .1) ) +
            ggtitle(input$title) +
            labs(x= input$group ,y= input$value)+
            theme_classic(base_size =  18) +
            stat_compare_means() + # 添加 P 值。 library(ggpubr)
            theme(legend.position = "none", plot.title = element_text(size = 13, hjust = 0.5) ) + 
            stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0.01,
                               hide.ns =  logi , remove.bracket = (input$star=='F'),size = input$star_size ) + #隐藏无意义的p值
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) #调整y轴
        }
      
      }
      

      result <- list()
      result$plot <- plot

      result$stat.test <- stat.test
      result$test_table <-test_table
      
      return(result )
    
    } )
    
 
  output$stat.test <- renderDataTable({
      
      result <- result()
      stat.test <- result$stat.test 
      
      return( stat.test )
      
    } ) # 统计结果
    
    output$test_table <- renderDataTable({
      
      result <- result()
      test_table <- result$test_table
      
      return( test_table )
      
    } ) # 统计结果
    
    
    
    output$plot <- renderPlot( {
      
      result <- result()
      plot <- result$plot
      
      return( plot )
    })
    
    
    # 3.1、下载差异结果1
    output$downloadcontent1 <- downloadHandler(
      
      filename = function() {
        paste('analysis.csv')
      } ,
      
      content = function(file) {
        
        result <- result()
        test_table <-  result$test_table
        
        # 输出
        write.csv(test_table , file, row.names = F, fileEncoding = 'GB18030')
        
      }  ) 
    
    # 3.2、下载差异结果2
    output$downloadcontent2 <- downloadHandler(
      
      filename = function() {
        paste('contrast.csv')
      } ,
      
      content = function(file) {
        
        result <- result()
        stat.test <- as.matrix( result$stat.test )
        stat.test[,11] <- gsub(',','vs', stat.test[,11])
        # 输出
        write.csv(stat.test , file, row.names = F, fileEncoding = 'GB18030')
        
      }  ) 
    
    
    # 4、下载图形
    # 下载图片 .pdf
    output$downloadplot_pdf <- downloadHandler(
      filename = function() {
        paste('analysis.pdf')
      },
      content = function(file) {
        
        result <- result()
        plot <- result$plot
        
        pdf(file,width = input$w, height = input$h)
        
        print(plot)
        dev.off()
        
      }  )
    
    # .png
    output$downloadplot_png <- downloadHandler(
      filename = function() {
        paste('analysis.png')
      },
      content = function(file) {
        
        result <- result()
        plot <- result$plot
        
        png(file,width= input$w, height= input$h , unit="in",res=150)
        
        print(plot)
        dev.off()
        
      }  )
    
  })
  
  # 1、参考数据
  output$sample_table <- renderDataTable({
    sample_table <- data()
    
  } )
  
  # 2、下载参考数据
  output$downloadSampleData <- downloadHandler(
    
    filename = function() {
      paste('方差分析参考数据.csv')
    } ,
    
    content = function(file) {
      
      sample_table <- data()
      write.csv(sample_table , file, row.names = F, fileEncoding = "GB18030")
      
    }  ) 
  
  
}

# Run the app
shinyApp(ui = ui, server = server)
