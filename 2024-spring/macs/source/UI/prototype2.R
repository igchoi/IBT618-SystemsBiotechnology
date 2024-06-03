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
  files = paste(rep("Img",5),1:5,rep("resized.jpg",5), sep="")
  paste(rep(path_result,5),files, sep="/")
  
  writeImage(img1,
             files=paste(rep(path_result,5),files,sep="/")) 
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
##loading tab function
Loading_tab = function(n){
  dat <- data.frame(x = numeric(0), y = numeric(0))
  withProgress(message = 'Making plot', value = 0, {
    # Number of times we'll go through the loop
    n 
    
    for (i in 1:n) {
      # Each time through the loop, add another row of data. This is
      # a stand-in for a long-running computation.
      dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
      
      # Increment the progress bar, and update the detail text.
      incProgress(1/n, detail = paste("Doing part", i))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.3)
    }
  })
  
  plot(dat$x, dat$y)
}


#==========================================================================
library(shiny)
library(shinyFiles) 
library(bslib)
library(DT)
library(ggplot2)
library(shinyDirectoryInput)
library(colourpicker)


ui <- page_fillable(
  titlePanel(h1(strong("Microbead sizer"
                ))),
  sidebarLayout(
    sidebarPanel(
      strong(h2("Select color")),
      strong("Size histogram:", textOutput("value_S", inline = TRUE)),
      colourInput("col_S", "Choose colour", "#666699"),
      textInput("text_S", "New colour:HEX value)"),
      

    br(),
    
    strong("Density histogram:", textOutput("value_D", inline = TRUE)),
    colourInput("col_D", "Choose colour", "#666699"),
    textInput("text_D", "New colour:HEX value)"),
   
    
    br(),
    
    strong("histogram:", textOutput("value_H", inline = TRUE)),
    colourInput("col_H", "Choose colour", "#666699"),
    textInput("text_H", "New colour:HEX value)"),
    
    br(), 
    actionButton("btn", "Update"),
    br(),
    br(),
  
    strong(h2("Select Folder")),
    br(),
      directoryInput('directory', label = 'select a directory', value = '~'),
      verbatimTextOutput("value"),
      actionButton("start", label = "START"),
     
      
      ),
  
    mainPanel(
      card(card_header(h6(strong("OVERALL"))),
             DTOutput(outputId = "SUM")
           ),
      
      card(card_header(h6(strong("SIZE"))),
             plotOutput(outputId = "HTSUM")
           ),
    
      card(card_header(h6(strong("DENSITY"))),
             plotOutput(outputId = "HTDEN")
           ),
      card(card_header(h6(strong("HISTOGRAM"))),
             plotOutput(outputId = "HT")
           )
      
    )
     ))
    
    
    
server <- function(input, output, session) {  
  
  
  observeEvent(input$btn, 
               {
    updateColourInput(session, "col_S",
            value = input$text_S)
    updateColourInput(session, "col_D",
                      value = input$text_D)
    updateColourInput(session, "col_H",
                      value = input$text_H)
  
    })
        output$value_S <- renderText(input$col_S)
        output$value_D <- renderText(input$col_D)
        output$value_H <- renderText(input$col_H)
  
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
                         fill = input$col_S)+
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
      
      
      output$HT <- renderPlot({
        isolate({
          BuildHT_HT <- ggplot(data = EEE)+
            aes(x = EEE$diam_mean, y = ..density..)+
            geom_histogram(alpha =0.4, fill = input$col_H, col = "black") +
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
      
    })
  })
     )
  
  
  
}


shinyApp(ui = ui, server = server)

