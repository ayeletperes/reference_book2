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
# Authenticate and save token for later use2
token <- drop_auth(rdstoken = "dropbox_token.rds")

# Retrieveing your file is as simple as

edit_links <- readLines("edit_links.txt")
share_links <- readLines("share_links.txt")
ui <- fluidPage(
    
    # Application title
    mainPanel(
        uiOutput('ab1_out'),
        shiny::actionButton(inputId='ab2', label="Reload text", 
                            icon = icon("sync-alt")),
        actionButton("refresh", "Refresh app"),
        #htmlOutput("df_output"),
        uiOutput('markdown')
    )
)

server <- function(input, output, session) {
    
    observeEvent(input$refresh, {
      session$reload();
    })
    
    observe({
        query <- parseQueryString(session$clientData$url_search)
        if(!is.null(query$group)) {
            group <- strsplit(query$group, "\"")[[1]][1]
            onclick <- edit_links[grep(group,edit_links)]
            output$ab1_out <- renderUI({shiny::actionButton(inputId='ab1', label="Edit text", 
                                         icon = icon("edit"), 
                                         onclick =paste0("window.open('",onclick,"')"))})
        }else{
            group <- "IGHV1-2G2"
            onclick <- edit_links[grep(group,edit_links)]
            output$ab1_out <- renderUI({shiny::actionButton(inputId='ab1', label="Edit text", 
                                                            icon = icon("edit"), 
                                                            onclick =paste0("window.open('",onclick,"')"))})
        }
    
        observeEvent(input$ab2, {
            drop_download(paste0("public/conclusions/",group,".docx"), local_path = paste0(group,".docx"),
                          overwrite = TRUE, verbose = FALSE, progress = FALSE, dtoken =token)

            invisible(rmarkdown::pandoc_convert(input =  paste0(group,".docx"), to = "markdown", output = paste0(group,"_docx.md"), options = c("--wrap=none","--reference-links","--metadata-file=metadata.yaml","--standalone"), verbose = FALSE))
            output$markdown <- renderUI({
                HTML(markdown::markdownToHTML(knitr::knit( paste0(group,"_docx.md"), quiet = TRUE)))
            })
            
        })
        
        output$markdown <- renderUI({
            drop_download(paste0("public/conclusions/",group,".docx"), local_path = paste0(group,".docx"),
                          overwrite = TRUE, verbose = FALSE, progress = FALSE, dtoken = token)
            
            invisible(rmarkdown::pandoc_convert(input =  paste0(group,".docx"), to = "markdown", output = paste0(group,"_docx.md"), options = c("--wrap=none","--reference-links","--metadata-file=metadata.yaml","--standalone"), verbose = FALSE))
            
            HTML(markdown::markdownToHTML(knitr::knit( paste0(group,"_docx.md"), quiet = TRUE)))
        })
    })
}

shinyApp(ui = ui, server = server)
