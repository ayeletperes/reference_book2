#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(rdrop2)
library(rhandsontable)
# Authenticate and save token for later use2
token <- drop_auth(rdstoken = "dropbox_token.rds")

# Retrieveing your file is as simple as

ui <- fluidPage(
    
    # Application title
    mainPanel(
      rHandsontableOutput('hot')
    )
)

server <- function(input, output, session) {
    
    values = reactiveValues(data=NULL, update = NULL)
    
    observe({
      tab <- paste0("public/conclusions/check_table.tsv")
      drop_download(tab, local_path = "check_table.tsv",
                    overwrite = TRUE, verbose = FALSE, progress = FALSE, dtoken = token)
      df <- read.delim("check_table.tsv")
      values$data <- df
    })
    
    observe({
      if(!is.null(input$hot))
        values$data <- hot_to_r(input$hot)
    })
    
    output$hot <- renderRHandsontable({
      rhandsontable(values$data)
    })
    
    #df <- eventReactive(df = hot_to_r(input$markdown))
    
    observe({
      req(input$hot)
      values$update <- hot_to_r(input$hot)
    })

    
    
    observe({
      req(values$update)
      write.table(values$update,
                  paste0("check_table.tsv"), row.names = F, sep = "\t")
      print(class(values$update))
      rdrop2::drop_upload(paste0("check_table.tsv"),
                          path = "public/conclusions/",
                          mode = "overwrite", 
                          dtoken = token, verbose = T)
    })
}

shinyApp(ui = ui, server = server)
