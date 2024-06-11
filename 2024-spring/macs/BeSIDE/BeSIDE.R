library(tidyverse)
library(pliman) # requires Bioconductor package - EBImage
library(EBImage)
# upload file by using only one FOLDER!!!


microsizer <- function(a){
  
  input1 <- a
  
  img.files1 = list.files(path= a,
                          pattern="*.jpg",
                          full.names = TRUE)
  print(img.files1)
  
  img1 = readImage(img.files1) # each images in an array 2040 1528 3 5
  dim(img1) # each image is in an array of 2040 x1528 x 3 x 5
  
  # resize images
  img1 = resize(img1, w=1020, h=764)
  dim(img1)
  
  ##
  ## save resized images in ./result folder (THIS part is NOT necessary!)
  ##
  path_result <- paste(input1, "result", sep = "/")
  print(path_result)
  display(img1)
  
  # naming files
  i <- (length(img.files1))
  
  files = paste(rep("Img",i),1:i,rep("resized.jpg",i), sep="")
  paste(rep(path_result,i),files, sep="/")
  
  writeImage(img1,
             files=paste(rep(path_result,i),files,sep="/")) 
  # result/ folder should be created before saving...
  ##
  ##
  
  ##
  ## for calculating
  ##  
  
  # adjust threshold and analyze picture
  str(img1)
  dim(img1)
  colorMode(img1) = "Color"
  colorMode(img1)
  numberOfFrames(img1, type = "total") # 15 frames (RGB 3 frames x 5 images)
  numberOfFrames(img1, type = "render") # 5 images
  getFrame(img1,1, type="render") # "render" gets a single image
  getFrames(img1, type="render")
  ##
  
  ##
  ## each images in 'img1' object is called, analyzed...
  ##
  
  npic = numberOfFrames(img1, type = "render") # total number of images
  combine_sum = NULL # your final data?
  
  for (i in 1:npic) {
    # 'i'th image in 'img1' object
    #i = 1
    x = analyze_objects(getFrame(img1,i,type="render"),
                        marker = "id", 
                        index = "I",          # <--- BI, I
                        object_size = "large", # <--- small, medium, large, elarge
                        lower_noise = 0.1,     # <--- change this?
                        lower_size = NULL,     # <--- change this?
                        upper_size = NULL,
                        lower_eccent = NULL,
                        lower_circ = NULL)
    # index https://tiagoolivoto.github.io/pliman/articles/indexes.html
    # 'BI' using brightness for detecting objects??
    # 'x' is the result of analyze_objects
    plot_measures(x, measure = "id")
    names(x)
    #[1] "results"          "statistics"       "object_rgb"       "object_index"     "efourier"        
    #[6] "efourier_norm"    "efourier_error"   "efourier_power"   "efourier_minharm" "veins"           
    #[11] "angles"           "width_at"         "mask"             "pcv"              "contours"        
    #[16] "parms
    x$results %>% 
      filter(circularity_norm > 0.6 & circularity > 2) %>%
      select(id, diam_mean, area, circularity, circularity_norm) %>% 
      arrange(circularity_norm)
    #
    y = get_measures(x, dpi=7.9016)
    names(y)
    # what does this number (79016) determine?
    # what for 'micrometer()' function?
    class(y)
    class(y) = "data.frame"
    y %>% 
      filter(circularity_norm > 0.6 & circularity > 2) %>%
      select(id, diam_mean, area, circularity, circularity_norm) %>% 
      arrange(circularity_norm)
    # collecting/combining results
    combine_sum = bind_rows(combine_sum, y, .id="img")
  }
  combinebead <- combine_sum %>% select(diam_mean)
  write_csv(combinebead, paste(path_result, "RAWDATA.csv", sep = "/"))
  
}
#==========================================================================

library(shiny)
library(shinyFiles) 
library(bslib)
library(DT)
library(ggplot2)
library(shinyDirectoryInput)
library(colourpicker)
library(shinyalert)
library(readxl)


ui <- page_sidebar(
  title = h1(strong("BeSIDE")),
  
  mainPanel( h3(strong("BeSIDE - Bead Size Distribution Estimator")),
  "This application is the useful for estimating the size and number of microparticles by using unlimited microscopic pictures and automatically save in the folder",
  textOutput("progress"),
  tags$head(tags$style("#progress{color: green;
                        font-size: 18px;
                        font-style: bold;}"))),
  
 
  sidebar = sidebar(
    bg = "white",
    accordion(
              
      accordion_panel(
        
        "Select color",
        strong("Histogram:", textOutput("value_H", inline = TRUE)),
        colourInput("col_H", "Choose colour", "#666699"),
  
        
        
        hr(),
        
        strong("Density histogram:", textOutput("value_D", inline = TRUE)),
        colourInput("col_D", "Choose colour", "#666699"),
       
        
        
        hr(),
        
        strong("Size histogram:", textOutput("value_S", inline = TRUE)),
        colourInput("col_S", "Choose colour", "#666699"),
     
        
        br(), 
        fluidRow(
          column(3,
                 actionButton("btn", "SET")
          ),
          column(3),
          column(3, actionButton("reset_C", "RESET"),
          )),
        br(),
        textOutput("complete"),
        tags$head(tags$style("#complete{color: green;
                        font-size: 15px;
                        font-style: bold;}")),
        hr(),
        
      ),
      
      accordion_panel(
        "Select Folder",
                      br(),
                      directoryInput('directory', label = 'select a directory', value = '~'),
                      verbatimTextOutput("value"),
                       actionButton("start", label = "RUN",
                                               style="color:white; background-color: black; border-color: white"),
                      
                      
        br()
        
   
    ),
    
      )
    ),
  
  accordion(
 
    accordion_panel(
      "Overall",
      DTOutput(outputId = "SUM")
    ),
    accordion_panel(
      "Histogram",
      plotOutput(outputId = "HTSUM")
    ),
    accordion_panel(
      "Density",
      plotOutput(outputId = "HTDEN")
    ),
    accordion_panel(
      "Size",
      plotOutput(outputId = "HTSIZE")
  )
  ),
  
  mainPanel(
    fluidRow(
      column(3, actionButton("modify", label = "UPDATE")),
      column(3),
      column(3, actionButton("exit", label = "CLOSE",
                             style="color:white; background-color: coral; border-color: white"))
      
    )
  )
)

server <- function(input, output, session) {
  
  observe({
    if (input$exit > 0) stopApp()  # stop shiny
    
  })
  
  observeEvent(input$btn,
               {
          
                 output$complete <- renderText("Recent update is completed.")
               })
  output$value_S <- renderText(input$col_S)
  output$value_D <- renderText(input$col_D)
  output$value_H <- renderText(input$col_H)
  
  observeEvent(input$reset_C,
               {
                 updateColourInput(session, "col_S",
                                   value = "#666699")
                 updateColourInput(session, "col_D",
                                   value = "#666699")
                 updateColourInput(session, "col_H",
                                   value = "#666699")
                 output$value_S <- renderText("#666699")
                 output$value_D <- renderText("#666699")
                 output$value_H <- renderText("#666699")
                 
                 output$complete <- renderText("Reset is done")
                 
               })
  
  
  observeEvent(
    input$start,
    
    isolate ({
      ignoreNULL = F
      eventExpr = {
        input$directory
      }
      handlerExpr = {
        if (input$directory > 0) {
          # condition prevents handler execution on initial app launch
          
          # launch the directory selection dialog with initial path read from the widget
          pp = choose.dir(default = readDirectoryInput(session, 'directory'))
          
          # update the widget value
          updateDirectoryInput(session, 'directory', value = pp)
        }}
      
    })
  )
  
  
  
  observeEvent(
    
    input$start,
    
    isolate ({
      isolate({
        pp = choose.dir(default = readDirectoryInput(session, 'directory'))
        output$value = renderText(pp)
        EEE <- microsizer(pp)
        output$HTSUM <- renderPlot({
          isolate({
            BuildHT_SUM <- ggplot(data = EEE)+
              aes(x = EEE$diam_mean, y = ..density..)+
              geom_histogram(alpha =0.4, fill = "azure2", col = "azure3")+
              geom_density(alpha = 0.3, 
                           fill = input$col_H)+
              scale_x_continuous(expand=c(0, 0)) +
              scale_y_continuous(expand=c(0, 0)) +
              labs(x="Bead size (um)", y="Density") +
              theme_minimal() +
              theme(axis.text.x = element_text(size=10),
                    axis.text.y = element_text(size=10),
                    legend.position = "bottom",
                    plot.title = element_text(hjust = 0.5),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "black"),
                    legend.text = element_text(size=15),
                    title = element_text(size=12))
            plot(BuildHT_SUM)
            
          })
          
          path_result <- paste(pp, "result", sep = "/")
          png(paste(path_result, "Histogram summary.png", sep = "/"))
          plot(BuildHT_SUM)
          dev.off()
          dev.set()
          
          
          
        })
        
        ###
        output$HTDEN <- renderPlot({
          isolate({
            BuildHT_DEN <- ggplot(data = EEE)+
              aes(x = EEE$diam_mean, y = ..density..)+
              geom_density(alpha = 0.4, fill = input$col_D) +
              scale_x_continuous(expand=c(0, 0)) +
              scale_y_continuous(expand=c(0, 0)) +
              labs(x="Bead size (um)", y="Density") +
              theme_minimal() +
              theme(axis.text.x = element_text(size=10),
                    axis.text.y = element_text(size=10),
                    legend.position = "bottom",
                    plot.title = element_text(hjust = 0.5),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "black"),
                    legend.text = element_text(size=15),
                    title = element_text(size=12))
            
            plot(BuildHT_DEN)
            
          })
          
          path_result <- paste(pp, "result", sep = "/")
          png(paste(path_result, "Histogram density.png", sep = "/"))
          plot(BuildHT_DEN)
          dev.off()
          dev.set()
        })
        
        
        output$HTSIZE <- renderPlot({
          isolate({
            BuildHT_HT <- ggplot(data = EEE)+
              aes(x = EEE$diam_mean, y = ..density..)+
              geom_histogram(alpha =0.4, fill = input$col_S, col = "black") +
              scale_x_continuous(expand=c(0, 0)) +
              scale_y_continuous(expand=c(0, 0)) +
              labs(x="Bead size (um)", y="Density") +
              theme_minimal() +
              theme(axis.text.x = element_text(size=10),
                    axis.text.y = element_text(size=10),
                    legend.position = "bottom",
                    plot.title = element_text(hjust = 0.5),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "black"),
                    legend.text = element_text(size=15),
                    title = element_text(size=12))
            
            plot(BuildHT_HT)
          })
          path_result <- paste(pp, "result", sep = "/")
          png(paste(path_result, "Histogram.png", sep = "/"))
          plot(BuildHT_HT)
          dev.off()
          dev.set()
        })
        
        output$SUM <- renderDT(EEE, options = list(pageLength =5))
        output$progress = renderText("Analysis is done")
        

        ################################## Update again   
        observe({
          if (input$modify >0) {
            path_result_csv <- paste(pp, "result","RAWDATA.csv", sep = "/")
            SSS <- read_csv(path_result_csv)
            output$HTSUM <- renderPlot({
              isolate({
                BuildHT_SUM <- ggplot(data = SSS)+
                  aes(x =SSS$diam_mean, y = ..density..)+
                  geom_histogram(alpha =0.4, fill = "azure2", col = "azure3")+
                  geom_density(alpha = 0.3, 
                               fill = input$col_H)+
                  scale_x_continuous(expand=c(0, 0)) +
                  scale_y_continuous(expand=c(0, 0)) +
                  labs(x="Bead size (um)", y="Density") +
                  theme_minimal() +
                  theme(axis.text.x = element_text(size=10),
                        axis.text.y = element_text(size=10),
                        legend.position = "bottom",
                        plot.title = element_text(hjust = 0.5),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(colour = "black"),
                        legend.text = element_text(size=15),
                        title = element_text(size=12))
                plot(BuildHT_SUM)
                
              })
              
              path_result <- paste(pp, "result", sep = "/")
              png(paste(path_result, "Histogram summary.png", sep = "/"))
              plot(BuildHT_SUM)
              dev.off()
              dev.set()
             
            })
            
            ###
            output$HTDEN <- renderPlot({
              isolate({
                BuildHT_DEN <- ggplot(data = SSS)+
                  aes(x = SSS$diam_mean, y = ..density..)+
                  geom_density(alpha = 0.4, fill = input$col_D) +
                  scale_x_continuous(expand=c(0, 0)) +
                  scale_y_continuous(expand=c(0, 0)) +
                  labs(x="Bead size (um)", y="Density") +
                  theme_minimal() +
                  theme(axis.text.x = element_text(size=10),
                        axis.text.y = element_text(size=10),
                        legend.position = "bottom",
                        plot.title = element_text(hjust = 0.5),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(colour = "black"),
                        legend.text = element_text(size=15),
                        title = element_text(size=12))
                
                plot(BuildHT_DEN)
                
              })
              
              path_result <- paste(pp, "result", sep = "/")
              png(paste(path_result, "Histogram density.png", sep = "/"))
              plot(BuildHT_DEN)
              dev.off()
              dev.set()
            })
            
            
            output$HTSIZE <- renderPlot({
              isolate({
                BuildHT_HT <- ggplot(data = SSS)+
                  aes(x = SSS$diam_mean, y = ..density..)+
                  geom_histogram(alpha =0.4, fill = input$col_S, col = "black") +
                  scale_x_continuous(expand=c(0, 0)) +
                  scale_y_continuous(expand=c(0, 0)) +
                  labs(x="Bead size (um)", y="Density") +
                  theme_minimal() +
                  theme(axis.text.x = element_text(size=10),
                        axis.text.y = element_text(size=10),
                        legend.position = "bottom",
                        plot.title = element_text(hjust = 0.5),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(colour = "black"),
                        legend.text = element_text(size=15),
                        title = element_text(size=12))
                
                plot(BuildHT_HT)
              })
              path_result <- paste(pp, "result", sep = "/")
              png(paste(path_result, "Histogram.png", sep = "/"))
              plot(BuildHT_HT)
              dev.off()
              dev.set()
            })
          }
        })
        


      })
    })
  )
  
}

shinyApp(ui = ui, server = server)





