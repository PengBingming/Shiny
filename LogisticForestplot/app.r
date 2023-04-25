
library(shiny)
library(shinydashboard)

library(showtext) # 解决画图中文乱码
showtext_auto()

library(AER)
library(readxl) # 读取 Excel

library(ggplot2)

# 编写函数
myfun_fp <- function(data){
  fit <- data
  df1 <- cbind( "Pvalue"= round( summary(fit)$coefficients[,4], 3), # P
                'OR' = round(exp( coef(fit) ), 2) ,  # OR
                "Lower"= round(exp( confint(fit) )[,1], 2 ), # 置信区间
                "Upper"= round(exp( confint(fit) )[,2] , 2)
  ) 
  
  df1 <- data.frame(df1)[-1,]
  df1 <- cbind('Var'=rownames(df1),df1)
  
  df1$Factor <- ifelse(df1$Lower>1,'Risk',ifelse(df1$Upper<1,'Protective','Not sig.'))
  
  df1$Pvalue <- ifelse(df1$Pvalue >= 0.001,df1$Pvalue,'<0.001')
  
  df1 <- df1[order(df1$OR,decreasing = T),]
  df1$`OR (95% CI)` <- paste0(df1$OR,'(',df1$Lower,'-',df1$Upper,')')
  
  # # 整理数据
  # fit <- data
  # 
  #   df1 <- cbind( exp( coef(fit) ) ,  # OR
  #                 summary(fit)$coefficients[,4], # P
  #                 exp( confint(fit) ) # 置信区间
  #   ) 
  # df1 <- data.frame(df1)[-1,]
  # df1 <- cbind('Var'=rownames(df1),df1)
  # colnames(df1)[-1] <- c("OR","Pvalue","OR_lower","OR_upper")
  # 
  # df2 <- df1
  # df2$OR_mean <- df2$OR
  # df2$OR <- paste0(   round(df2$OR,2),   # OR
  #                     "(", round(df2$OR_lower,2), # OR_lower
  #                     "~", round(df2$OR_upper,2), # OR_upper
  #                     ")")
  # df2$Pvalue <- ifelse( df2$Pvalue>=0.001 , round(df2$Pvalue,3), "<0.001")
  # df2
  # 
  # fp <- rbind("lables"=NA,df2)
  # fp[1, 1:3]  <-  c('Variable', 'OR(95% CI)', 'Pvalue')
  
  return(df1)
}

# 画图
myfun_plot <-  function(data,input){

  df1$Var1 <- factor(1:nrow(df1))

  annotation <- data.frame(matrix("",ncol = 3,nrow = c(3*nrow(df1)+3) ))
  colnames(annotation) <- c('x','y','label')
  annotation$label <- c('OR (95% CI)','Odds  Ratio','P Value',paste0(df1$OR,'(',df1$Lower,'-',df1$Upper,')'),df1$Pvalue, df1$Var)
  annotation$x <- c( c(-0.3,1,-0.75),rep(-0.3, nrow(df1) ),rep(-0.75,nrow(df1) ), rep(-1.1, nrow(df1)) )
  annotation$y <- c(rep(nrow(df1)+0.47,3),seq(1, nrow(df1), 1),seq(1, nrow(df1), 1) ,seq(1, nrow(df1),1 ) )

    plot <-  ggplot(df1, aes(OR, Var1)) +
    geom_point(size=3.6, aes(col=Factor)) +
    geom_errorbarh(aes(xmax =Upper, xmin = Lower), col=input$lines, height = 0.4)+
    geom_vline(aes(xintercept=1),col=input$zero)+
    scale_x_continuous(limits=c(-1.2, max(df1$Upper)), breaks=seq(0,max(df1$Upper) , 0.5)) +
    theme_bw() + 
    theme(legend.position ="top") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.background = element_rect(fill =input$background, color ='red'),
          axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=10, face = "bold"),
          legend.text=element_text(size=11)) +
    labs(x='',y='',title = input$title)+
    geom_text(data=annotation,aes(x=x,y=y,label=label))
  
  return(plot)
}


# 参考数据
data(Affairs,package = "AER")
df <- Affairs

df <- cbind("group"= ifelse(df$affairs>0,1,0), df[,-1])
df$group <- factor(df$group, levels = c(0,1) )

# 设置shiny界面参数
header <- dashboardHeader(title = "Logistic回归森林图") # 页面标题

sidebar <- dashboardSidebar(
  fileInput("file", "输入文件：",
            multiple = TRUE,),
  h5('(支持 .csv .xlsx .xls文件)'),
  actionButton("submit", "展示示例 / 开始分析"),
  tags$hr(),
  # 筛选变量
  selectInput("type", "变量筛选：",
              c('全部变量'='full',
                '逐步回归'='step'), 
              selected = "full" ),
  
  # 确定效应变量
  
  tags$hr(),
  
  selectInput("background",  "背景颜色", colors() , selected = "skyblue" ),
  selectInput("zero",        "竖线颜色", colors() , selected = "black" ),
  selectInput("lines",       "横线颜色", colors() , selected = "black" ),
  
  tags$hr()
)

# 网页呈现的内容 
body <- dashboardBody(
  tabsetPanel(
    tabPanel(h4("Logistics回归"),
             fluidRow( 
               column(width = 7,
                      downloadButton("downloadtable","下载参考数据"),
                      box( title = "输入数据/参考数据",
                           width = NULL, 
                           height = 20, 
                           status = "primary",
                           solidHeader = T ),
                      box( selectInput("group", "响应变量：二分类", c("") ),
                           width = 4, 
                           height = 0,
                           solidHeader = T ),
                      box( selectInput("num",    "连续变量：数值",  multiple = T,c("") ), # 连续变量
                           width = 3, 
                           height = 0,
                           solidHeader = T ),
                      box( selectInput("factor", "分类变量：因子",  multiple = T,c("") ), # 分类变量
                           width = 3, 
                           height = 0,
                           solidHeader = T )),
               
               column(width = 5,
                      downloadButton("downloadtable1","下载回归结果"),
                      box( title = "Logistic 回归结果",
                           width = NULL, 
                           height = 20, 
                           status = "primary",
                           solidHeader = T ) ) ),
             br(),
             br(),
             br(),
             splitLayout(cellWidths = c("58%",'2',"40%"),
                         dataTableOutput("contents"),"",
                         dataTableOutput("contents1"))
    ),
    
    tabPanel(h4("森林图"),
             fluidRow( 
               column("下载图形：",
                      downloadButton("downloadplot1",".png"),
                      downloadButton("downloadplot2",".pdf"),
                      width = 7,
                      box( title = "森林图/Forestplot",
                           width = NULL, 
                           height = 20, 
                           status = "primary",
                           solidHeader = T ),
                      box( numericInput(inputId = 'w',label = '下载图形长：',value = 10),
                           width = 4, 
                           height = 0,
                           solidHeader = T ),
                      box( numericInput(inputId = 'h',label = '下载图形高：',value = 5),
                           width = 4, 
                           height = 0,
                           solidHeader = T ),
                      box( textInput(inputId = 'title',label = '图形标题',value="Forestplot"),
                           width = 4, 
                           height = 0,
                           solidHeader = T )
               ),
               
               column(width = 5,
                      downloadButton("downloadtable2","下载表格"),
                      box( title = "表格/Table",
                           width = NULL, 
                           height = 20, 
                           status = "primary",
                           solidHeader = T ) ) ),
             br(),
             br(),
             br(),
             splitLayout(cellWidths = c("53%",'5',"42%"),
                         plotOutput("plot2"),"",
                         dataTableOutput("contents2"))
    )
    
  ) # tabsetPanel
) # dashboardBody

# 设置ui
ui <- dashboardPage(header, sidebar, body)

# 设置输出
server <- function(input,output,session){
  
  # 数据：读取输入文件或参考数据
  df1 <- reactive({
    
    file1 <- input$file
    
    if(is.null(file1)){
      df1 <- df
    } 
    else{
      
      d <- tail( unlist(strsplit(file1$datapath,'[.]') ), 1)
      
      if(d=='csv'){
        df1 <- data.frame( read.csv(file1$datapath,fileEncoding = "GB18030") )
      } else{
        df1 <- data.frame( read_excel(file1$datapath,1) )
      } 
      
    } # else
    
    return(df1)
    
  })
  
  
  # 输入数据/参考数据
  output$contents <- renderDataTable({
    
    df <- df1()
    
  })  
  
  observe({
    
    df <- df1()
    
    updateSelectInput(session, 'group',  label = "响应变量：二分类",     choices = c(colnames(df) ),selected = colnames(df)[1] )
    updateSelectInput(session, 'num',    label = "连续变量：数值", choices = c(colnames(df) ) )
    updateSelectInput(session, 'factor', label = "分类变量：因子", choices = c(colnames(df) ) )
    
  })
  
  
  
  observeEvent(input$submit, {
    
    
    
    myfun <- reactive({
      
      df <- df1()
      
      n  <-  which(colnames(df) == input$group ) 
      
      # logistics回归
      
      colnames(df)[n] <- 'group_yn'
      
      df$group_yn <- factor(df$group_yn )
      
      # 数值型数据
      if(length(input$num)==1){
        df[,input$num ] <- as.numeric(df[,input$num ]  )
      }
      else if(length(input$num)>1){
        df[,input$num ] <- apply( df[,input$num ] ,2, as.numeric)
      }
      
      # 因子型数据
      if( length(input$factor)==1 ){
        df[,input$factor ] <- factor(df[,input$factor ]  )
      }
      else if( length(input$factor)>1 ){
        df[,input$factor ] <- apply( df[,input$factor ] ,2, factor)
      }
      
      fit.full <- glm(group_yn~ . , # 所有变量
                      data = df, # 数据集
                      family = binomial(link='logit') # 拟合方式
      )
      
      fit.step <- step(object = fit.full,trace = 0)
      
      if(input$type == 'full'){
        
        fp <- myfun_fp( data = fit.full )
        plot <- myfun_plot(fp, input)
        fit <- as.data.frame(summary( fit.full )$coefficients)
        
        
      }
      else if(input$type == 'step'){
        
        fp <- myfun_fp( data = fit.step )
        plot <- myfun_plot(fp,input)
        fit <- as.data.frame(summary( fit.step )$coefficients)
      }
      
 
      fit <- cbind( 'Var'=rownames(fit), fit )
      
      results <-list()
      
      results$fp <- fp
      results$fit <- fit 
      results$plot <- plot
      
      return( results )
      
    })
    
    # 森林图图形
    output$plot2 <- renderPlot({
      
      results <- myfun()
      plot <- results$plot
      
      return(plot)
      
    }  )
    
    
    # 展示 logistic 结果表格
    output$contents1 <- renderDataTable({
      
      results <- myfun()
      fit <- results$fit
      
      return(fit)
      
    })  
    
    # 展示森林图表格
    output$contents2 <- renderDataTable({
      
      results <- myfun()
      fp <- results$fp
      
      return(fp)
      
    })  
    
    
    
    # 下载图片 png
    output$downloadplot1 <- downloadHandler(
      filename = function() {
        paste('Forestplot.png')
      },
      content = function(file) {
        
        png(file,width= input$w, height= input$h , unit="in",res=150)
        
        results <- myfun()
        plot <- results$plot
        
        print(plot)
        
        dev.off()
        
      }  )
    
    # 下载图片 pdf
    output$downloadplot2 <- downloadHandler(
      filename = function() {
        paste('Forestplot.pdf')
      },
      content = function(file) {
        
        pdf(file,width=input$w, height=input$h)
        
        results <- myfun()
        plot <- results$plot
        
        print(plot)
        
        dev.off()
        
      }  )
    
    
    
    # 1、下载logistic回归结果
    output$downloadtable1 <- downloadHandler(
      filename = function() {
        paste('logistic.csv')
      },
      content = function(file) {
        
        results <- myfun()
        fit <- results$fit
        write.csv(fit , file,  row.names = F, fileEncoding = 'GB18030')
        
      } )
    
    # 2、下载森林图结果
    output$downloadtable2 <- downloadHandler(
      filename = function() {
        paste('forestplot.csv')
      },
      content = function(file) {
        
        results <- myfun()
        fp <- results$fp
        write.csv(fp , file,  row.names = F, fileEncoding = 'GB18030')
        
      } )
    
    
  })
  
  # 0、下载参考数据
  output$downloadtable <- downloadHandler(
    filename = function() {
      paste('sample.csv')
    },
    content = function(file) {
      df <- df
      write.csv(df , file,  row.names = F, fileEncoding = 'GB18030')
    } )
  
}

shinyApp(ui=ui, server = server)

