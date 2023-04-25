
options(digits = 6)
# 加载包
library(shiny)
library(shinydashboard)

library(readxl)
library(rhandsontable)
# library(xlsx)

library(drc)

library(showtext)
showtext_auto() 

library(ggplot2)
library(patchwork) # 合并图片

# 设置shiny界面参数
header <- dashboardHeader(title = "ELISA 数据分析") # 页面标题

# 可输入、输出的按钮 
sidebar <- dashboardSidebar(
  fileInput("file", "输入文件：",
            multiple = TRUE,),
  h5('(支持 .csv .xlsx .xls文件)'),
  actionButton("submit", "展示示例 / 开始分析"),
  tags$hr(),
  selectInput("type", "拟合类型",
              choices = c('直线' = "line",
                          '二项式' = "binomial",
                          '4参数logistic' = "logistic"),
              selected = "line"),

  tags$hr(),

  selectInput("rep", "有无副孔",
              choices = c(无副孔 = "wu",
                          有副孔 = "you"),
              selected = "you"),
  tags$hr(),
  numericInput("round", "保留结果小数位数", min=0,max = 10,value = 3 ),
  downloadButton("downloadData", "下载计算结果"),
  tags$br(),
  tags$br(),
  downloadButton("downloadSampleData", "下载参考数据"),
  tags$hr()
)

# 网页呈现的内容 
body <- dashboardBody( 
  
  tabsetPanel(
    tabPanel('计算结果',
             fluidRow( column(width = 6, 
                              box( 
                                title = "数据/Data：",
                                h6('支持复制粘贴'),
                                width = NULL, 
                                height = 10, 
                                status = "primary",
                                solidHeader = T ) ),
                       column(width = 6, 
                              box( 
                                title = "计算结果/Result：",
                                width = NULL, 
                                height = 10, 
                                status = "primary",
                                solidHeader = T ) )
             ),
             tags$br(),
             tags$br(),
             tags$br(),
             splitLayout(cellWidths = c("48%","2%","50%"),
                         rHandsontableOutput("table"),"",
                         rHandsontableOutput('results')),
             fluidRow( column(width =12, 
                              box( 
                                title = "拟合图形/Plot：",
                                width = NULL, 
                                height = 10, 
                                status = "primary",
                                solidHeader = T ) )
             ),
             tags$br(),
             splitLayout(cellWidths = c("90%"),
                         plotOutput("standard_plot") ) )

  )
  
)

# 设置 ui
ui <- dashboardPage(header, sidebar, body)

server <- function(input, output) {
  
  # 读取 ELISA 数据
  df1 <- reactive({
    
    file1 <- input$file
    if(is.null(file1)){
      if(input$rep == "you"){
        df1 <- data.frame(read_excel('../file/IL-13_1.xlsx',1) )
      }
      else if(input$rep == "wu"){
        df1 <- data.frame(read_excel('../file/IL-13.xlsx',1) )
      }
    }
    else if(!is.null(file1)){
      
      d <- tail( unlist(strsplit(file1$datapath,'[.]') ), 1)
      
      if(d=='csv'){
        df1 <- data.frame( read.csv(file1$datapath, fileEncoding = "GB18030") )
      } else{
        df1 <- data.frame( read_excel(file1$datapath,1) ) 
      } 
    } 
    return(df1)
  })
  

  
  observeEvent(input$submit, {
    
    df0 <- df1()
    
    output$table <- renderRHandsontable(
      rhandsontable(df0) %>% hot_col(colnames(df0),type = 'numeric')
    )
    
   observeEvent(input$table, {
    # 分析处理 ELISA数据
    myfun1 <- reactive({
      
      # 1、读取数据
      if(is.null(df1 ) ){return(NULL)} 
      
      df1 <- hot_to_r( input$table )
      df1 <- as.data.frame( df1 )
      colnames(df1) <- tolower(colnames(df1))
      
      # 2、判断有无副孔
      if(input$rep == "you"){
        
        nm2 <- vector() 
        for (i in 1:( ncol(df1)/2-1) ) {
          nm1 <- c(paste0('孔',i) , paste0('副孔',i) )
          nm2 <-c(nm2,nm1)
        }
        
        colnames(df1)[3:ncol(df1)] <- nm2
        
        k1 <- 2*(2:(dim(df1)[2]/2))-1 ; k2 <- 2*(2:(dim(df1)[2]/2))
        df1 <- cbind(df1[,1:2],(df1[,k1] + df1[,k2])/2)
        df2 <- df1
        
      }  else  (df2 <- df1 )
      
      # 3、直线回归
      if(input$type == "line"){
        
        fit1 <-lm(y ~ x, data=df1) 
        a <- fit1$coefficients[2]; b <- fit1$coefficients[1] # 参数赋值
        
        f <- paste0('y  ≈  ',round(a, 6),'x +  ',round(b, 6) ) # 函数
        r <- paste0('Adjusted R-squared ≈ ', round(as.numeric(summary(fit1)[9]), 3) ) # R方
        
        p1 <- ggplot( data = df1, aes(x = x, y = y) ) +
          geom_point(shape = 1 ) +
          geom_smooth(data = df1,  aes(x = x, y = y), color= 'gray60', method = "gam", formula = y ~ poly(x,1) ) +
          annotate("text", x = (max(df1$x)-min(df1$x))/2 , y = max(df1$y)*0.98, label = paste0(f , '\n' , r) , colour="black") +
          labs(title = 'ELISA : binomial',
               x = "x : conc \n \n 拟合曲线", 
               y = "y : OD") 
        
        # 构建由 y 求 x 函数
        myfun <- function(y,a,b){
          
          x <- round((y-b)/a , input$round) # 求 k 值
            
          return(x) # 输出计算值
        }
        
        x1 <- myfun(y=df1$y,a=a,b=b)
        
        y0 <- as.matrix(df1[,3:dim(df1)[2] ]) # 取96孔板其余列 OD值（y）
        x0 <- myfun(y=y0,a=a,b=b) # 求浓度（x）,如果报错可能为OD值超过函数界值。

        df2[, 3:dim(df2)[2] ]<-matrix(x0, ncol = dim(df2)[2]-2, nrow = 8,byrow = F)
        df2[, dim(df2)[2]+1 ]<-x1 ; colnames(df2)[(dim(df2)[2])]<-'x1'
        
        df2 <- as.data.frame(df2)
        
        df_f <- data.frame('od' = runif(1000, min = min(df1$y), max = max(df1$y) )) 
        df_f$conc <- as.numeric(myfun(y= df_f$od,a=a,b=b))
        
        p2 <- ggplot(data = df_f , aes(x = od , y = conc) ) + 
          geom_line( color= 'pink3' ) +
          geom_point( aes(x= y, y = x ) , shape= 1 ,data = df1 ) +
          ggtitle("ELISA : binomial") + xlab("OD \n \n 反函数曲线") + ylab("Conc")
        
        standard_plot <- p1 + p2
        
      }
      
      # 4、二项式回归
      else if(input$type == "binomial"){
        
        fit2<-lm(y ~ x+ I(x^2), data=df1) # 拟合曲线，x为浓度，y为OD值
        
        a <- fit2$coefficients[3]; b <- fit2$coefficients[2]; c <- fit2$coefficients[1] # 参数赋值
        
        f <- paste0('y  ≈  ',round(a,8),'x^2  +  ',round(b, 4),'x +  ',round(c, 4) ) # 函数
        r <- paste0('Adjusted R-squared ≈ ', round(as.numeric(summary(fit2)[9]), 3) ) # R方
        
        p1 <- ggplot( data = df1, aes(x = x, y = y) ) +
          geom_point(shape = 1 ) +
          geom_smooth(data = df1,  aes(x = x, y = y), color= 'gray60', method = "gam", formula = y ~ poly(x,2) ) +
          annotate("text", x = (max(df1$x)-min(df1$x))/2 , y = max(df1$y)*0.98, label = paste0(f , '\n' , r) , colour="black") +
          labs(title = 'ELISA : binomial',
               x = "x : conc \n \n 拟合曲线", 
               y = "y : OD") 
        
        # 构建由 y 求 x 函数
        myfun <- function(y,a,b,c){
          k<- 1/a # 求 k 值
          n<- -b/(2*a) # 求 n 值
          m<-(n)^2-(c/a) # 求 m 值
          
          # sqrt() 为开方，前面符号为±，按实际情况(是否大于0)选择 +（x1） 或 -（x2）
          if(a<0){
            x <- round(-sqrt(k*y+m)+n, input$round) # a大于0时的解
          }
          else if(a>0){
            x <- round(sqrt(k*y+m)+n, input$round) # a小于0时的解
          }
          
          return( x ) # 输出计算值
        }
        
        x1<-myfun(y=df1$y,a=a,b=b,c=c)
        
        y0 <- as.matrix(df1[,3:dim(df1)[2] ]) # 取96孔板其余列 OD值（y）
        x0 <- myfun(y=y0,a=a,b=b,c=c) # 求浓度（x）,如果报错可能为OD值超过函数界值。
        
        df2[, 3:dim(df2)[2] ]<-matrix(x0, ncol = dim(df2)[2]-2, nrow = 8,byrow = F)
        df2[, dim(df2)[2]+1 ]<-x1 ; colnames(df2)[(dim(df2)[2])]<-'x1'
        
        df_f <- data.frame('od' = runif(1000, min = min(df1$y), max = max(df1$y) )) 
        df_f$conc <- as.numeric(myfun(y= df_f$od,a=a,b=b,c=c))
        
        p2 <- ggplot(data = df_f , aes(x = od , y = conc) ) + 
          geom_line( color= 'pink3' ) +
          geom_point( aes(x= y, y = x ) , shape= 1 ,data = df1 ) +
          # xlim( min(df1$y) , max(df1$y) )+
          # ylim( min(df1$x) , max(df1$x) )+
          ggtitle("ELISA : binomial") + xlab("OD \n \n 反函数曲线") + ylab("Conc")
        
        standard_plot <- p1 + p2
        
      }
      # 5、logistic 回归
      else if(input$type == "logistic"){
        
        library(drc)
        
        # 拟合
        pl4 <- drm( y~x , fct=LL.4( names=c("Slope", "Lower", "Upper", "ED50") ), data= df1)
        RSD <- abs(sqrt(summary(pl4)$"resVar") / mean(fitted(pl4))) # 相对标准偏差RSD
        
        # 提取参数
        var <- as.numeric(pl4[["coefficients"]])
        
        UpperLimit <- var[3] ; EC50 <- var[4] 
        Slope <- -var[1] ; LowerLimit <- var[2] 
        
        # 设置函数
        myfun<- function(y,a,b,c,d){
          x <- ((a-d)/(y-d)-1)^(1/b)*c
          return( round(x, input$round) )
        }
        
        # 查看标曲拟合度
        x1 <- myfun( y=df1$y, a=LowerLimit , b= Slope, c=EC50, d=UpperLimit)
        
        # 计算后赋值给df2
        df_r <- dim(df1)[1]; df_n <- dim(df1)[2]
        
        y0 <- as.matrix(df1[,3:df_n]) # 取96孔板其余列 OD值（y）
        
        x0 <- myfun(y0, a = LowerLimit , b = Slope , c = EC50, d = UpperLimit) # 求浓度（x）,如果报错可能为OD值超过函数界值。
        
        df2[,3:df_n] <- matrix(x0,ncol = df_n-2,nrow = df_r,byrow = F)
        
        df2[,df_n+1] <- x1 ; colnames(df2)[df_n+1] <- "x1"
        
        # 四参数logistic拟合图
        
        df_f <- data.frame('conc' = runif(1000, min = min(df1$x), max = max(df1$x) ) )
        
        pm <- predict(pl4, newdata=expand.grid(df_f$conc), interval="confidence")
        
        df_f$od <- pm[,1]
        df_f$pmin <- pm[,2]
        df_f$pmax <- pm[,3]
        
        f <- paste0("y==", UpperLimit , "+", "frac(", LowerLimit ,"-",UpperLimit,",",
                    "1+(","frac(","x",",",EC50,")" ,")^",Slope," )")
        
        rsd <- paste0('RSD ≈ ', round(as.numeric(RSD ), input$round) ) 
        
        p1 <- ggplot(data = df_f , aes(x = conc , y = od) ) +
          geom_line( color ='black' ) +
          geom_point( aes(x, y ) , shape= 1 ,data = df1 ) +
          geom_ribbon(data=df_f, aes(x= conc, y= od, ymin=pmin, ymax=pmax), alpha=0.2)+
          annotate("text", x = (max(df1$x)-min(df1$x))/2 , y = max(df1$y)*0.95, parse=T, 
                   label= f, size = 4)+
          annotate("text", x = (max(df1$x)-min(df1$x))/2 , y = max(df1$y)*0.85, 
                   label= rsd, size = 4)+
          ggtitle("ELISA : 4pl_logistic") + xlab("x : Conc \n \n 拟合曲线\n(RSD：相对标准偏差)") + ylab("y : OD")
        
        # 反函数图
        df_f <- data.frame('od' = runif(1000, min = min(df1$y), max = max(df1$y) ) )
        df_f$conc <- myfun( y= df_f$od, a=LowerLimit , b= Slope, c=EC50, d=UpperLimit)
        
        p2 <- ggplot(data = df_f , aes(x = od , y = conc) ) +
          geom_line( color ='pink3' ) +
          geom_point( aes(x= y, y= x) , df1, shape= 1) +
          ggtitle("ELISA : 4pl_logistic") + xlab("OD \n \n 反函数曲线") + ylab("Conc")
        
        standard_plot <- p1 + p2
        
      }
      
      results               <- list()
      results$df1           <- df1
      results$df2           <- df2
      results$standard_plot <- standard_plot
      
      return(results)
    })
    
    
    # 3、分析情况展示
    
    output$results <- renderRHandsontable(
     
         rhandsontable( myfun1()$df2 ) 
    )
    
    
    # 4、画图ELISA的标曲图像，并呈现在网页上
    output$standard_plot <-  renderPlot({
      
      if ( is.null(myfun1() )) { return() }
      
      results <- myfun1()
      standard_plot <- results$standard_plot
      return(standard_plot)
      
    })
    
    # 5、依据函数计算结果，并设置可下载
    output$downloadData <- downloadHandler(
      
      filename = function() {
        paste("计算结果.csv")
      },
      content = function(file) {
        
        results <- myfun1()
        df2 <- as.data.frame( results$df2 )
        
        write.csv( df2 , file , row.names = F, fileEncoding = 'GB18030') 
      } )
  })
  
  })
  
  # 2、下载参考数据
  output$downloadSampleData <- downloadHandler(
    
    filename = function() {
      paste('Elisa_data.csv')
    },
    content = function(file) {
      sample  <- hot_to_r( input$table )
      write.csv( sample , file , row.names = F, fileEncoding = 'GB18030') 
      
    } )  
  
}

shinyApp(ui, server)
