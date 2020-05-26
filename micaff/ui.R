library(shiny)
library(shinyjs)

jsfile <- "getdata.js"
cssfile <- "style.css"
js_scroll_file <-'scroll.js'

shinyUI(
    fluidPage(
      titlePanel("MicAff Application!"),
      sidebarLayout(
        sidebarPanel(
          fileInput("file1", "Choose Files for analysis",
                    multiple = TRUE,
                    accept = c(".CEL", ".cel"))
        ),
        mainPanel(
        )
      )
    )
)