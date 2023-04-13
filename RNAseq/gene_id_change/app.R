
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

library(airway)

header <- dashboardHeader(title = "基因 ID 转换")

geneName <- intersect(columns(org.Mm.eg.db),columns(org.Hs.eg.db) )
# # ②.2 侧边界面 ------------------------------------------------------------------
sidebar <- dashboardSidebar(
  
  fileInput(inputId = "file1", "选择文件：",
            multiple = F),
  h5('支持：.csv .xlsx .xls文件'),
  h5('转换的 基因列名 需为：ID'),
  tags$hr(),
  radioButtons("species", "物种：",
               choices = c('人' = "hsa",
                           '鼠' = "mmu"),
               selected = "hsa"),
  tags$hr(),
  selectInput("from_geneName","基因名：form", geneName, selected = "ENTREZID"),
  selectInput("to_geneName","基因名：to", geneName , selected = "ENTREZID"),
  actionButton("submit1", "开始转换"),
  tags$hr()
)


# # ②.3 网页呈现的内容  --------------------------------------------------------------
body <- dashboardBody( 
  
  # navbarPage 1 DEseq2
  tabsetPanel(
    tabPanel('转换前表格',
             fluidRow(  
               column(width = 12, 
                      box(downloadButton("downloadSampleData", "参考数据"),
                          title = '原始表格',
                          width = NULL, 
                          height = 10,
                          status = "primary",
                          solidHeader = T)
               )
             ),
             tags$br(),
             tags$br(),
             tags$br(),
             tags$br(),
             splitLayout(cellWidths = c("95%"),
                         dataTableOutput("exp_orgin") ) 
    ), # tabPanel 1.1
    tabPanel('转换后表格',
             fluidRow(
               column(width = 12, 
                      box(downloadButton("downloadResultData", "转换结果"),
                        title = "转换后表格",
                        width = NULL, 
                        height = 10, 
                        status = "primary",
                        solidHeader = T) 
               )
             ),
             tags$br(),
             tags$br(),
             tags$br(),
             tags$br(),
             splitLayout(cellWidths = c("95%"),
                         dataTableOutput("exp_change")) )
  ) # navbarPage
  
) # dashboardBody

# # ②.4 ui -----------------------------------------------------------------
ui <- dashboardPage(header, sidebar, body)


server <- function(input, output) {
  
  
  
 observeEvent(input$submit1, {
    id_change <- reactive({
    
    file1 <- input$file1
    if(is.null(file1 )){
      
      data("airway",package = "airway")
      dat <- data.frame( assay(airway) )
      dat<- cbind('ID'=rownames(dat), dat)
      
      }
    else{

      d <- tail( unlist(strsplit(file1$datapath,'[.]') ), 1)

      if(d=='csv'){
        dat <- data.frame( read.csv(file1$datapath,1) )
      } else{
        dat <- data.frame( read_excel(file1$datapath,1) )
      }

      }

    if(input$from_geneName==input$to_geneName){
      dat1 <- dat
      dat1$m <- "基因 ID 相同，无需转换"
      colnames(dat1)[ncol(dat1)] <- input$to_geneName

    } else if(!input$from_geneName==input$to_geneName){

      # dat$ID = rownames(dat)
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

    df_exp$dat  <- dat
    df_exp$dat1 <- dat1

    return(df_exp)

  } )
    output$exp_change <- renderDataTable({
      
      if (is.null(id_change() ) ) { return() }
      
      df_exp <- id_change()
      
      dat1   <- df_exp$dat1
      
      return(dat1)
    
      })
    
    # 下载结果
    output$downloadResultData <- downloadHandler(
      filename = function() {
        paste('转换结果数据.csv')
      },
      content = function(file) {
        df_exp <- id_change()
        dat1   <- df_exp$dat1
        write.csv(dat1, file,row.names = F, fileEncoding = "GB18030") 
      } )
    
 } )
    
  output$exp_orgin <- renderDataTable({
    
    file1 <- input$file1
    if(is.null(file1 )){
     
      data("airway",package = "airway")
      dat <- data.frame( assay(airway) )
      dat<- cbind('ID'=rownames(dat),dat)
    }
    else{
      
      d <- tail( unlist(strsplit(file1$datapath,'[.]') ), 1)
      if(d=='csv'){
        dat <- data.frame( read.csv(file1$datapath) )
      } else{
        dat <- data.frame( read_excel(file1$datapath,1) )
      }

    }

    return(dat)
    
 })
  # 3.2.1 下载参考数据  DEseq2
  output$downloadSampleData <- downloadHandler(
    filename = function() {
      paste('参考数据.csv')
    },
    content = function(file) {
      
      data("airway",package = "airway")
      dat <- data.frame( assay(airway) )
      dat <- cbind('ID'=rownames(dat),dat)
      write.csv(dat, file,row.names = F, fileEncoding = "GB18030") 
    }
  )
  

}

shinyApp(ui, server)
