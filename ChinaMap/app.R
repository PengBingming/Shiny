
library(shiny)
library(shinydashboard)

library(showtext) # 解决画图中文乱码
showtext_auto()

library(plyr) # 合并地图数据与输入数据
library(readxl)
library(rhandsontable)

library(ggplot2)

load("../file/china_map.Rdata") # 地图数据

# 设置shiny界面参数
header <- dashboardHeader(title = "中国地图（省级）") # 页面标题

# 可输入、输出的按钮 
sidebar <- dashboardSidebar(
  fileInput("file", "输入文件：",
            multiple = TRUE,),
  h5('(支持 .csv .xlsx .xls文件)'),
  actionButton("submit", "开始画图/Analysis"),
  tags$hr(),
  selectInput("type", "标签 Label",
               c('中文（省）'='Chinese',
                 '省会城市'='City',
                 'English'='English',
                 '无标签'='wu',
                 '自定义/Customization'='Customization'
               ), 
               selected = "Chinese" ),
  radioButtons("netline", "经纬度线   Lat & long", 
               c('有'='you',
                 '无'='wu'),
               selected = "wu" ),
  tags$hr(),
  selectInput("low",  "低值颜色 LowColor", colors() , selected = "white" ),
  selectInput("high", "高值颜色 HighColor", colors() , selected = "red" ),
  selectInput("line", "界限颜色 LineColor", colors() , selected = "grey60" ),

  tags$hr()
)

# 网页呈现的内容 
body <- dashboardBody(
  tabsetPanel(
    tabPanel(h4("中国地图/ChinaMap"),
             fluidRow( 
               column(width = 8,
                      box( title = "图形/Plot",
                           width = NULL, 
                           height = 20, 
                           status = "primary",
                           solidHeader = T )),
               column(width = 4,
                      box( title = "数据/Data",
                           br(),
                           br(),
                           h5('Customization 列为自定义标签，value为数据（数值型）'),
                           width = NULL, 
                           height = 20, 
                           status = "primary",
                           solidHeader = T ) ),
               column(width = 1,
                      br(),
                      downloadButton("downloadplot",".Rdata")
               ),
               column(width = 1,
                      br(),
                      downloadButton("downloadplot1",".png")
               ),
               column(width = 1,
                      br(),
                      downloadButton("downloadplot2",".pdf") 
               ),
                 column( width = 2,
                         br(),
                         numericInput(inputId = 'w',label = '下载图形长：Weight',value = 20)
                         ),
                 column( width = 2, 
                         br(),
                         numericInput(inputId = 'h',label = '下载图形高：High',value = 10)
                         ),
               column(width = 1),
                 column(width = 2,
                        br(),
                        downloadButton("downloadtable","下载画图数据")
                      )
             ),
             splitLayout(cellWidths = c("68%","32%"),
                         plotOutput("plot1"),
                         rHandsontableOutput('table') )
    )
    
  ) # tabsetPanel
) # dashboardBody

# 设置ui
ui <- dashboardPage(header, sidebar, body)

# 设置输出
server <- function(input,output){

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
      rhandsontable(df,rowHeaderWidth = 22,width = 360, height = 400) %>% 
        hot_cols(columnSorting = TRUE) %>% 
        ot_col("NAME", readOnly = TRUE) %>% 
        hot_col('value',type = 'numeric')
    )
    
  observeEvent(input$table, {
   
    # 编写函数
    myfun <- reactive({
  
      mydata <- hot_to_r( input$table )
      
      china_data <- join(china_map_data, mydata, type="full")  #合并 输入数据 与 坐标数据
      china_province <- join(province_city, mydata, type="full")  #合并 输入数据 与 标签数据
      china_data$value <- as.numeric(china_data$value)
      
      # library(ggplot2)
      
      # 一、无经纬度网格线
      if(input$netline=='wu'){ 
        
        # 1 添加各省地名（中文）标签
        plot1 <- ggplot( china_data,aes(long,lat) )+
          geom_polygon( aes(group=group, fill=value), colour = input$line)+
          scale_fill_gradient( low= input$low, high= input$high ) +
          geom_line(data=l9,aes(x=long,y=lat,group=group),color= input$line,size=0.5)+ # 9段线
          # coord_cartesian(xlim=c(70,140),ylim=c(4,52))+ #缩小显示范围在南部区域；不用
          coord_map("polyconic") +
          theme(
            aspect.ratio = .8, #调节长宽比
            panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank()
          ) +
          # 添加各省地名（中文）标签
          geom_text(aes(x = jd,y = wd,label = province), data =province_city)
        
        # 2 添加各省地名（英文）标签
        plot2 <- ggplot( china_data,aes(long,lat) )+
          geom_polygon( aes(group=group, fill = value), colour = input$line)+
          scale_fill_gradient( low= input$low, high= input$high ) +
          geom_line(data=l9,aes(x=long,y=lat,group=group),color= input$line, size=0.5)+ # 9段线
          coord_map("polyconic") +
          theme(
            aspect.ratio = .8, #调节长宽比
            panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank()
          ) +
          # 添加各省地名（英文）标签
          geom_text(aes(x = jd,y = wd,label = city), data = province_city)
        
        # 3 添加各省（省会）标签
        plot3 <- ggplot( china_data,aes(long,lat) )+
          geom_polygon( aes(group=group, fill=value), colour = input$line)+
          scale_fill_gradient( low= input$low, high= input$high ) +
          geom_line(data=l9,aes(x=long,y=lat,group=group),color= input$line,size=0.5)+ # 9段线
          coord_map("polyconic") +
          theme(
            aspect.ratio = .8, #调节长宽比
            panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank()
          ) +
          # 添加各省（省会）标签
          geom_text(aes(x = jd,y = wd,label = province_english), data = province_city)
        
        
        # 4 添加各省省会（自定义）标签
        plot4 <- ggplot( china_data,aes(long,lat) )+
          geom_polygon( aes(group=group, fill=value), colour = input$line)+
          scale_fill_gradient( low= input$low, high= input$high ) +
          geom_line(data=l9,aes(x=long,y=lat,group=group),color= input$line,size=0.5)+ #9段线
          coord_map("polyconic") +
          theme(
            aspect.ratio = .8, #调节长宽比
            panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank()
          ) +
          # 添加各省省会（自定义）标签
          geom_text(aes(x = jd,y = wd,label = Customization), data = china_province)
        
        # 5 无标签
        plot5 <- ggplot( china_data,aes(long,lat) )+
          geom_polygon( aes(group=group, fill=value), colour = input$line)+
          scale_fill_gradient( low= input$low, high= input$high ) +
          geom_line(data=l9,aes(x=long,y=lat,group=group),color= input$line,size=0.5)+ #9段线
          coord_map("polyconic") +
          theme(
            aspect.ratio = .8, #调节长宽比
            panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank()
          )
        
      }
      # 二、添加经纬度网格线
      else if(input$netline=='you'){
        # 1 添加各省地名（中文）标签
        plot1 <- ggplot( china_data,aes(long,lat) )+
          geom_polygon( aes(group=group, fill=value), colour = input$line)+
          scale_fill_gradient( low= input$low, high= input$high ) +
          geom_line(data=l9,aes(x=long,y=lat,group=group),color= input$line,size=0.5)+ # 9段线
          coord_map("polyconic") +
          # 添加各省地名（中文）标签
          geom_text(aes(x = jd,y = wd,label = province), data =province_city)
        
        # 2 添加各省地名（英文）标签
        plot2 <- ggplot( china_data,aes(long,lat) )+
          geom_polygon( aes(group=group, fill = value), colour = input$line)+
          scale_fill_gradient( low= input$low, high= input$high ) +
          geom_line(data=l9,aes(x=long,y=lat,group=group),color= input$line, size=0.5)+ # 9段线
          coord_map("polyconic") +
          # 添加各省地名（英文）标签
          geom_text(aes(x = jd,y = wd,label = city), data = province_city)
        
        # 3 添加各省（省会）标签
        plot3 <- ggplot( china_data,aes(long,lat) )+
          geom_polygon( aes(group=group, fill=value), colour = input$line)+
          scale_fill_gradient( low= input$low, high= input$high ) +
          geom_line(data=l9,aes(x=long,y=lat,group=group),color= input$line,size=0.5)+ # 9段线
          coord_map("polyconic") +
          # 添加各省（省会）标签
          geom_text(aes(x = jd,y = wd,label = province_english), data = province_city)
        
        
        # 4 添加各省省会（自定义）标签
        plot4 <- ggplot( china_data,aes(long,lat) )+
          geom_polygon( aes(group=group, fill=value), colour = input$line)+
          scale_fill_gradient( low= input$low, high= input$high ) +
          geom_line(data=l9,aes(x=long,y=lat,group=group),color= input$line,size=0.5)+ #9段线
          coord_map("polyconic") +
          # 添加各省省会（自定义）标签
          geom_text(aes(x = jd,y = wd,label = Customization), data = china_province)
        
        # 5 无标签
        plot5 <- ggplot( china_data,aes(long,lat) )+
          geom_polygon( aes(group=group, fill=value), colour = input$line)+
          scale_fill_gradient( low= input$low, high= input$high ) +
          geom_line(data=l9,aes(x=long,y=lat,group=group),color= input$line,size=0.5)+ #9段线
          coord_map("polyconic")
      }
      
      plots <-list()
      plots$plot1 <- plot1
      plots$plot2 <- plot2
      plots$plot3 <- plot3
      plots$plot4 <- plot4
      plots$plot5 <- plot5
      
      return( plots )
      
    })
    
    # 中文画图
    output$plot1 <- renderPlot({
      
      plots <- myfun()
      
      if( input$type=='Chinese' ){
        plot <- plots$plot1
      }
      else if( input$type=='City' ){
        plot <- plots$plot2
      }
      else if( input$type=='English' ){
        plot <- plots$plot3
      }
      else if( input$type=='Customization' ){
        plot <- plots$plot4
      }
      else if( input$type=='wu' ){
        plot <- plots$plot5
      }
      return(plot)
      
    }  )
    
    # 下载图片 png
    output$downloadplot1 <- downloadHandler(
      filename = function() {
        paste('ChinaMap.png')
      },
      content = function(file) {
        
        png(file,width= input$w, height= input$h , unit="in",res=150)
        
        plots <- myfun()
        
        if( input$type=='Chinese' ){
          plot <- plots$plot1
        }
        else if( input$type=='City' ){
          plot <- plots$plot2
        }
        else if( input$type=='English' ){
          plot <- plots$plot3
        }
        else if( input$type=='Customization' ){
          plot <- plots$plot4
        } 
        else if( input$type=='wu' ){
          plot <- plots$plot5
        }
        
        print(plot)
        
        dev.off()
        
      }  )
    
    # 下载图片 pdf
    output$downloadplot2 <- downloadHandler(
      filename = function() {
        paste('ChinaMap.pdf')
      },
      content = function(file) {
        
        pdf(file,width=input$w, height=input$h)
        
        plots <- myfun()
        
        if( input$type=='Chinese' ){
          plot <- plots$plot1
        }
        else if( input$type=='City' ){
          plot <- plots$plot2
        }
        else if( input$type=='English' ){
          plot <- plots$plot3
        }
        else if( input$type=='Customization' ){
          plot <- plots$plot4
        } 
        else if( input$type=='wu' ){
          plot <- plots$plot5
        }
        
        print(plot)
        
        dev.off()
        
      }  )
    
    # # 下载图形 .Rdata
    output$downloadplot <- downloadHandler(
      
      filename = function() {
        paste("plots.Rdata")
      },
      
      content = function(file) {
        
        plots <- myfun()
        
        if( input$type=='Chinese' ){
          plot <- plots$plot1
        }
        else if( input$type=='City' ){
          plot <- plots$plot2
        }
        else if( input$type=='English' ){
          plot <- plots$plot3
        }
        else if( input$type=='Customization' ){
          plot <- plots$plot4
        } 
        else if( input$type=='wu' ){
          plot <- plots$plot5
        }
        
        save( plot, file = file)
        
      }  )
    
  } )
    
  } )
  # 下载参考数据
  output$downloadtable <- downloadHandler(
    filename = function() {
      paste('ChinaMap_data.csv')
    },
    content = function(file) {
      sample <- hot_to_r( input$table )
      write.csv(sample, file,  row.names = F, fileEncoding = 'GB18030')
    } )
  
}
shinyApp(ui=ui, server = server)
