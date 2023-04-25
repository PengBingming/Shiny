
library(shiny)
library(shinydashboard)

library(readxl)
library(reshape)

library(ggplot2)
library(ggpubr)

load("../file/qpcr_df0.Rdata")
 # 编写
 myfun1 <- function(df1, input){
  
  # 去重复

  colnames(df1) <- toupper(colnames(df1))
  
  
  dfs <- df1[which(df1$CQ.STD.DEV < input$Cq.Std.Dev), -which(colnames(df1)=='CQ.STD.DEV')]
  
  df1 <- melt( dfs, id.vars = c('GROUP','SAMPLE','TARGET') )

  m = t( reshape::cast(df1,TARGET~GROUP+SAMPLE,mean) )
  m <- as.data.frame.array(m)
  
  
  # 计算过程变量
  dat1 <- data.frame( matrix(0, ncol = 4*(ncol(m)-1) , nrow = nrow(m) )   )
  rownames(dat1) <- rownames(m)
  
  gene <- setdiff(  unique(df1$TARGET) , input$gene )
  
  n1 <- vector()
  for (i in gene ) {
    n1 <- c(n1 , paste0(  i ,"_",c("▲ct","2^(-▲ct)","mean","2^(-▲▲ct)" )))
  }
  
  colnames(dat1 ) <- n1
  
  m <- cbind(Sample=do.call( rbind, strsplit( rownames(m),'_') )[, 2],m)
  m <- cbind(Group=do.call( rbind, strsplit( rownames(m),'_') )[, 1],m)
  
  m <- cbind(m, dat1)
  
  for (i in 1:length(gene) ) {
    
    n <- i*4+ length(gene) 
    
    m[ , n] <-      m[, gene[i] ] - m[, input$gene]  # ▲ct
    
    m[ , n+1] <- 2^(-(m[, gene[i] ] - m[, input$gene]) ) # '2^(-▲ct)'
    

    mean <- mean(na.omit( m[ which(m$Group == input$group),n+1 ]) ) # 对照组 '2^(-▲ct)' 平均值
    
    m[ which(m$Group == input$group), n+2 ] <- mean # 赋值
     
    m[ , n+3 ] <- m[ , n+1] / mean    # 所有组除以对照组平均（'2^(-▲ct)'）
    
  }
  
  result <- list()
  result$dat1 <- m
  result$dat2 <- m[, c( 1:(3+length(gene) ), 1:length(gene)*4+(3+length(gene) ) ) ]
                 
  return(result)
  
  }

 
 # 设置shiny界面参数
 header <- dashboardHeader(title = "Q-PCR 数据分析") # 页面标题
 
 # 可输入、输出的按钮 
 sidebar <- dashboardSidebar(
   fileInput("file", "输入 xlsx 文件",
             multiple = TRUE),

   actionButton("submit", "展示示例 / 开始分析"),
   tags$hr(),
   selectInput("gene",  "内参", c("") ),
   selectInput("group", "对照", c("") ), 

   numericInput(inputId = "Cq.Std.Dev",
                label = "Cq.Std.Dev < ",
                value = 0.5) ,
   radioButtons("type", "结果选择",
                choices = c('计算结果' = "dat2",
                            '包含过程' = "dat1"),
                selected = "dat2"),
   tags$hr(),
   tags$br(),
   downloadButton("downloadSampleData", "下载参考数据"),
   hr()
   
 )
 
 # 网页呈现的内容 
 body <- dashboardBody(
      tabsetPanel(

      tabPanel('参考数据',
         fluidRow( 
           column(width = 6, 
                  box(dataTableOutput("sample"), 
                      title = "参考数据：Group为分组，Target为基因，Sample为样本",
                      width = 12, 
                      height = 10, 
                      status = "primary",
                      solidHeader = T))
           )),
      tabPanel('计算结果',
               fluidRow( 
                 column(width = 6, 
                        downloadButton("downloadData", "下载计算结果"),
                        box( title = "处理数据",
                             width = NULL, 
                             height = 10, 
                             status = "primary",
                             solidHeader = T)),
               column(width = 6, 
                      downloadButton("downloadplot","下载表达图形"),
                      box(selectInput("select", "表达情况", c("") ), 
                          title = "图形",
                          width = 12, 
                          height = 10, 
                          status = "primary",
                          solidHeader = T)
                       ) ),
               tags$br(),
               tags$br(),
               tags$br(),
               tags$br(),
               splitLayout(cellWidths = c("50%","50%"),
                           dataTableOutput("results"),
                           plotOutput("plot") ) 
      )

 )
 )
 
 # 设置ui
 ui <- dashboardPage(header, sidebar, body)
 
 
 
 server <- function(input,output,session){
   
   
   # 读取 qpcr数据
   df1 <- reactive({
     
     file1 <- input$file
     if(is.null(file1)){
       df1 <- df0
     }
     else{
       
       d <- tail( unlist(strsplit(file1$datapath,'[.]') ), 1)
       
       if(d=='csv'){
         df1 <- data.frame( read.csv(file1$datapath,1) )
       } else{
         df1 <- data.frame( read_excel(file1$datapath,1) )
       } 
     } 
     
     return(df1)
   })

   
   # 1、展示参考数据
   output$sample <- renderDataTable({
     
     df1 <- df1()
     
   }) 
   
   
 observe({
     
   df1 <- df1()
   colnames(df1) <- toupper(colnames(df1))

     # 实验组与对照组情况
     updateSelectInput(session, "gene",  label = '内参', choices = unique(df1$TARGET) , selected = unique(df1$TARGET)[1] )
     updateSelectInput(session, "group", label = '对照', choices = unique(df1$GROUP)  , selected = unique(df1$GROUP )[1] )
 })
 

     


   
 observeEvent(input$submit, {
   # 3、展示计算结果
   output$results <- renderDataTable({
     
     if ( is.null( df1() ) ) { return() }
     
     result <- myfun1(df1=df1(), input = input)
     
       if(input$type=='dat1'){
           dat <- result$dat1
       } 
       else if(input$type=='dat2'){
           dat <- result$dat2
          } 
     return(dat)
   })
   
 
     
   observe({

     if ( is.null( df1() ) ) { return() }
     df1 <- df1()
     colnames(df1) <- toupper(colnames(df1))

     gene <- setdiff( unique(df1$TARGET), input$gene )
     updateSelectInput(session, "select", label = '图形', choices =  gene , selected = gene[1] )

    output$plot <- renderPlot({
     result <- myfun1(df1=df1(), input = input)

     dat2 <-result$dat2

     n1 <- which(gene==input$select)

     p <- ggplot(dat2, aes(x= Group, y= dat2[,n1+3+length(gene)], color= Group ) ) +
           geom_boxplot(alpha = 1, size = .5 ) +
           geom_jitter(size = .5, alpha= .5 , position =  position_jitter(width = .1) ) +
           ggtitle("") +
           labs( y= paste0( input$select,'/', input$gene) )+
           stat_compare_means( # method = "anova"
                              ) +
     theme_classic(base_size =  9)
       
     return( p )
   })
    
    # 下载图片
    output$downloadplot <- downloadHandler(
      filename = function() {
        paste('exp.png')
      },
      content = function(file) {
        
        png(file,width=5,height=6,unit="in",res=150)
        
        result <- myfun1(df1=df1(), input = input)
        
        dat2 <-result$dat2
        
        n1 <- which(gene==input$select)
        
        plot <- ggplot(dat2, aes(x= Group, y= dat2[,n1+3+length(gene)], color= Group ) ) +
          geom_boxplot(alpha = 1, size = .5 ) +
          geom_jitter(size = .5, alpha= .5 , position =  position_jitter(width = .1) ) +
          ggtitle("") +
          labs( y= paste0( input$select,'/', input$gene) )+
          stat_compare_means( # method = "anova"
          ) +
          theme_classic(base_size =  9)
        
        print( plot )
        
        dev.off()
        
      }  )
    
   } )
   

 } )
   
 # 2、下载参考数据
 output$downloadSampleData <- downloadHandler(
   filename = function() {
     paste('q-pcr参考数据.csv')
   },
   content = function(file) {
     
     df1<- df1()
     write.csv(df1 , file, row.names = F,fileEncoding = 'GB18030')
   }
 ) 
 
 # 3、输出计算结果
 output$downloadData <- downloadHandler(
   
   filename = function() {
     paste("计算结果.csv")
   },
   
   content = function(file) {
     result <- myfun1(df1=df1(), input)
     dat2 <- result$dat2
     write.csv(dat2 , file, row.names = F,fileEncoding = 'GB18030')
   }
 )
 

 }
 
 
 shinyApp(ui,server)
