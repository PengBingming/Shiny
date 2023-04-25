
# 加载包
library(shiny)
library(shinydashboard)
library(readxl)
library(xlsx)

library(showtext)
showtext_auto() 

myfun <- function(df0){
  
  df1 <- df0 ; df1 <- as.data.frame(df1)
  
  colnames(df1) <- c('parameter', 'Penh' )
  
  # 取各组小鼠名称与计数
  a <- names( which( table(df1$parameter) >= 2  ) ) 
  num <- table(df1$parameter) 
  
  # 将组名加到数据列
  for (i in 1:length(a)) { 
    df1[which(df1== a[i]),2]<-a[i]
  }
  
  # 去除缺失行
  df1 <- df1[which(!is.na(df1$Penh)),]
  
  
  v0 <- df1[,2] #取数据列赋值到新建向量v0
  
  # 循环参数准备
  l <- length(a)
  df2 <- data.frame(matrix(NA,ncol = l))
  colnames(df2) <- c(a)
  
  # 提取数据循环
  for(i in 1:length(a)) {
    
    n1 <- num[a[i]]
    
    v1 <- v0[(which(v0 == a[i])[1]+1):(which(v0 == a[i])[n1]-1)]
    df2[1:length(v1),i] <- v1
    
    # v2 <- v0[(which(v0 == a[i])[1]):(which(v0 == a[i])[n1] )]
    # df2[1:length(v2),i+l+1] <- v2
    
  }
  
  for(i in 1:dim(df2)[2]) {
    
    df2[which(is.na(df2[,i])),i] <- ''
    
  }
  
  return(df2)
  
}

# 设置shiny界面参数
header <- dashboardHeader(title = "肺功能数据分析") # 页面标题

# 可输入、输出的按钮 
sidebar <- dashboardSidebar(
  fileInput("file", "输入 excel 文件",
            multiple = TRUE,),
  tags$hr(),
  
  # Horizontal line ----
  downloadButton("downloadData", "下载处理结果"),
  tags$hr(),
  downloadButton("downloadSampleData", "下载参考数据")
)

# 网页呈现的内容 
body <- dashboardBody( 
  fluidRow( column(width =12, 
                   box(dataTableOutput("results"),
                       title = "处理结果",
                       width = NULL, 
                       height = 10, 
                       status = "primary",
                       solidHeader = T ) ) 
  ),
  tags$br(),
  tags$br(),
  
  fluidRow( column(width = 12, 
                   box(
                       dataTableOutput("sample"), 
                       title = "参考格式",
                       width = NULL, 
                       height = 10, 
                       status = "primary",
                       solidHeader = T) 
  )
  )
)

# 设置 ui
ui <- dashboardPage(header, sidebar, body)


server <- function(input, output) {
  
  df1 <- reactive({
    
    file1 <- input$file
    if(is.null(file1)){return(NULL)}
    
    library(readxl)
    data.frame(read_excel(file1$datapath,1))
    
  })
  
  # 1、展示参考数据
  output$sample <- renderDataTable({
    if (!is.null(df1() ) ) { return(df1() ) }
    
    library(readxl)
    sample1 <- read_excel('/home/shiny/shiny-server/lung/2023.3.23.xlsx',1)
    
  }) 
  
  # 2、下载参考数据
  output$downloadSampleData <- downloadHandler(
    filename = function() {
      paste('肺功能参考数据.xlsx')
    },
    content = function(file) {
      
      sample2 <- read_excel('/home/shiny/shiny-server/lung/2023.3.23.xlsx',1)
      write.xlsx(as.data.frame(sample2), file, sheetName = '数据',row.names = F, showNA = F)
      
    }
  ) 
  
  # 3、分析情况展示
  output$results <- renderDataTable({
    
    if (is.null(df1() ) ) { return() }
    df2 <- myfun(df0 = df1() )
    
  })
  
  
  # 5、依据函数计算结果，并设置可下载
  output$downloadData <- downloadHandler(
    
    filename = function() {
      paste("肺功能处理结果.xlsx")
    },
    content = function(file) {
      
      df2 <- myfun(df0 = df1() )
      write.xlsx(as.data.frame(df2), file, sheetName = 'data',row.names = F, showNA = T)
      
    } )
  
}

shinyApp(ui, server)
