
library(shiny)
library(shinydashboard)

library(showtext) # 解决画图中文乱码
showtext_auto()

library(AER)
library(readxl) # 读取 Excel

library(rhandsontable)
library(ggplot2)

data(Affairs, package = "AER")

myfun_table <- function(data ){
  
  # 参考数据
  df <- data
  
  df <- cbind("group"= ifelse(df$affairs>0,1,0), df[,-1])
  df$group <- factor(df$group, levels = c(0,1) )
  
  fit.full <- glm(group~ . , # 所有变量
                  data = df, # 数据集
                  family = binomial(link='logit')) # 拟合方式
  fit <- fit.full
  df1 <- cbind( "Pvalue"= round( summary(fit)$coefficients[,4], 3), # P
                'OR' = round(exp( coef(fit) ), 2) ,  # OR
                "Lower"= round(exp( confint(fit) )[,1], 2 ), # 置信区间
                "Upper"= round(exp( confint(fit) )[,2] , 2)
  ) 
  
  df1 <- data.frame(df1)[-1,]
  df1 <- cbind('Var'=rownames(df1),df1)
  
  # df1$Factor <- ifelse(df1$Lower>1,'Risk',ifelse(df1$Upper<1,'Protective','Not sig.'))
  
  # df1$Pvalue <- ifelse(df1$Pvalue >= 0.001,df1$Pvalue,'<0.001')
  
  df1 <- df1[order(df1$OR,decreasing = T),]
  # df1$`OR (95% CI)` <- paste0(df1$OR,'(',df1$Lower,'-',df1$Upper,')')
  rownames(df1) <- 1:nrow(df1)
  
  return(df1)
}

sample <- myfun_table(Affairs)

myfun_plot <- function(data,input){
  
  df1 <- data

  df1$Factor <- ifelse(df1$Lower > 1,'Risk',ifelse(df1$Upper < 1,'Protective','Not sig.'))
  df1$`OR (95% CI)` <- paste0(df1$OR,'(',df1$Lower,'-',df1$Upper,')')
  
  df1$Pvalue.sig <- ifelse(df1$Pvalue >= 0.001,df1$Pvalue,'<0.001')
  
  df1 <- df1[order(df1$OR,decreasing = T), ]
  df1$Var1 <- factor( 1:nrow(df1) )
  
  annotation <- data.frame(matrix("",ncol = 3,nrow = 3*nrow(df1)+3 ))
  colnames(annotation) <- c('x','y','label')

  annotation$label <- c('OR (95% CI)','Odds Ratio','P Value',paste0(df1$OR,'(',df1$Lower,'-',df1$Upper,')'),df1$Pvalue.sig, df1$Var)

  annotation$x <- c( c(-0.3,1,-0.75),rep(-0.3, nrow(df1) ),rep(-0.75,nrow(df1) ), rep(-1.1,nrow(df1)) )

  annotation$y <- c(rep(nrow(df1)+ 0.45,3),seq(1, nrow(df1), 1),seq(1, nrow(df1), 1) ,seq(1, nrow(df1),1 ) )
  
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

# 设置shiny界面参数
header <- dashboardHeader(title = "森林图") # 页面标题

sidebar <- dashboardSidebar(
  fileInput("file", "输入文件：",
            multiple = TRUE,),
  h5('(支持 .csv .xlsx .xls文件)'),
  actionButton("submit", "开始分析/Analysis"),
  tags$hr(),

  selectInput("background",  "背景颜色", colors() , selected = "skyblue" ),
  selectInput("zero",        "竖线颜色", colors() , selected = "black" ),
  selectInput("lines",       "横线颜色", colors() , selected = "black" ),
  
  tags$hr()
)

# 网页呈现的内容 
body <- dashboardBody(
  tabsetPanel(
    
    tabPanel(h4("森林图/Forestplot"),
             fluidRow( 
               column(width = 8,
                      box( title = "图形/Plot",
                           br(),br(),br(),
                           h6("运行后可查看 数据/图形。可直接编辑表格，鼠标位于表格右下单元格时，
                              右键即可插入或删除行；也可下载参考数据表格，编辑后上传"),
                           width = NULL, 
                           height = 20, 
                           status = "primary",
                           solidHeader = T ) ),
               column(width = 4,
                      box( title = "表格/Table",
                           width = NULL, 
                           height = 20, 
                           status = "primary",
                           solidHeader = T ) )
                      ),
             fluidRow(  
              column(width = 1, br(), downloadButton("downloadplot2",".pdf")
                      ),
               column(width = 1, br(),downloadButton("downloadplot1",".png")
               ) ,
              column(width = 2,numericInput(inputId = 'w',label = '下载图形长：',value = 10)
              ) ,
              column(width = 2,numericInput(inputId = 'h',label = '下载图形高：',value = 5)
              ),
              column(width = 2,textInput(inputId = 'title',label = '图形标题',value="Forestplot")
              ),
              column(width = 1, br(), downloadButton("downloadtable","表格/Table")
              )
              ),
             br(),br(),
             splitLayout(cellWidths = c("66%",'4',"30%"),
                         plotOutput("plot"),"",
                         rHandsontableOutput('table') )
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
      df1 <- sample
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

    
    observeEvent(input$submit, {
      
      df <- df1()

      output$table <- renderRHandsontable(
        
        df %>% rhandsontable( rowHeaderWidth = 22, height = 500 )%>% 
          hot_cols(columnSorting = TRUE) %>% 
          hot_col("Lower",type = 'numeric') %>%
          hot_col("Upper", type = 'numeric') %>%
          hot_col('OR', type = 'numeric')
      )

      observeEvent(input$table, {
        
        myfun <- reactive({
          
          mydata <- hot_to_r( input$table )
          plot <- myfun_plot(mydata,input)
          
          return( plot )
        })
        
        output$plot <- renderPlot({
          
          plot <- myfun()
          return(plot)
          
        } ) 

        # 下载图片 png
        output$downloadplot1 <- downloadHandler(
          filename = function() {
            paste('Forestplot.png')
          },
          content = function(file) {
            
            png(file,width= input$w, height= input$h , unit="in",res=150)
            
            plot <- myfun()
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
            
            plot <- myfun()
            print(plot)
            
            dev.off()
            
          }  )
        
        # 1、下载logistic回归结果
        output$downloadtable <- downloadHandler(
          filename = function() {
            paste('Forestplot.csv')
          },
          content = function(file) {

            mydata <- hot_to_r( input$table )
            write.csv(mydata, file,  row.names = F, fileEncoding = 'GB18030')
            
          } )
  })
  
 })
  
}


shinyApp(ui=ui, server = server)
