library(ggplot2)
library(lazyeval)
library(plotly)


plotTheme <- function(base_size = 12) {
    theme(
        text = element_text( color = "black"),
        plot.title = element_text(size = 12,colour = "black",hjust=0.5),
        plot.subtitle = element_text(face="italic"),
        plot.caption = element_text(hjust=0),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line("grey80", size = 0.1),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey80", color = "white"),
        strip.text = element_text(size=8),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        axis.title.x = element_text(hjust=1,size=15),
        axis.title.y = element_text(hjust=1,size=15),
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_text(colour = "black", face = "bold"),
        legend.text = element_text(colour = "black", face = "bold"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(vjust=-1,angle=90,size=15))
}

shinyServer(function(input, output, session) {
  
  
  output$contents <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    
    df <- read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote)
    
    if(input$disp == "head") {
      return(head(df))
    }
    else {
      return(df)
    }
    
  })
  
})
