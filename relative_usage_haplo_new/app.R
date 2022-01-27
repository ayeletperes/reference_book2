options(repos = BiocManager::repositories())

require(dplyr)
require(tidyr)
require(htmltools)
require(bbplot)
require(ggplot2)
require(shiny)
require(shinyWidgets)
require(data.table)
require(alakazam)
require(unikn)
require(plotly)
require(RColorBrewer)
require(unikn)
require(bbplot)
require(plotly)
require(ggplot2)
require(shinythemes)
library(htmlwidgets)


hline <- function(y = 0,
                  color = "red",
                  x0 = 0,
                  x1 = 1) {
  list(
    type = "line",
    x0 = x0,
    x1 = x1,
    xref = "x",
    y0 = y,
    y1 = y,
    line = list(
      color = color,
      dash = "dot",
      width = 4
    )
  )
}



load("alleles_dbs.rda")
allele_db <-
  alleles_dbs$IGH$functional$nonsingle$all$`318`$complete$`95`
allele_db <-
  allele_db %>% dplyr::rowwise() %>% dplyr::mutate(
    gene = alakazam::getGene(or_allele, strip_d = F, omit_nl = F),
    group = strsplit(gsub(gene, "", new_allele), "[*]")[[1]][1],
    gene_group = alakazam::getGene(new_allele, strip_d = F, omit_nl = F)
  )
load("functional_groups.rda")
func_groups <-
  functional_groups$IGH$functional$nonsingle$all$`318`$complete$`95`
cols <- c("#FAAB18", "#1380A1", "#990000", "#588300")
pal <- cols %>%
  unikn::newpal(names = c("orangy", "bluish", "redish", "greeny"))
load("data_frac_new2.rda")
chain <- "IGH"

data <-
  setDT(data_frac$IGH$functional$nonsingle$all$`318`$complete$`95`)
data[, v_call := paste0(v_gene, "*", v_allele)]
data$group_plot <- ifelse(is.na(data$j_call), 1, 2)

ui <- fluidPage(
  title = "Groups relative usage",

  tags$head(
    tags$style(
      "[type = 'number'] {font-size:10px;height:5px;}",
      ".butt{background-color:#FFFAFA;} .butt{color: #FFD700;} .butt{height:40px;
    width:57px;
    padding-top: 3 !important;}",
      ".butt:hover {
  background-color: #FFD700;
             color:#FFFAFA;}",
      HTML(
        '
      .same-row {
        max-width: 200px;
        display: table-cell;
        padding-right: 10px;
      }
    '
      )
    ),

    tags$script(
      "$(document).on('change', '.dynamicSI input', function () {
                               Shiny.setInputValue('lastSelectId', this.id, {priority: 'event'});
                              });"
    )
  ),
  theme = shinytheme("sandstone"),
  wellPanel(
    width = 12,
    style = "background: #FFFAFA",
    fluidRow(
      column(
        12,
        switchInput(
          inputId = paste0("hetro"),
          label = "J6 heterozygous",
          value = F,
          width = "150%",
          size = "small",
          labelWidth = "100px"

        ),
        # switchInput(
        #   inputId = paste0("naive"),
        #   label = "Naive repertoires",
        #   value = T, size = "small", inline = T, labelWidth = "100px"
        #
        # ),
        shinyWidgets::awesomeRadio(
          inputId = paste0("mut"),
          label = "# max mutations allowed",
          choices = c(
            "0=" = 0,
            "1=<" = 1,
            "2=<" = 2,
            "3=<" = 3
          ),
          selected = 0,
          inline = T

        ),
        switchInput(
          inputId = paste0("abs"),
          label = "absolute thresh",
          value = T,
          width = "150%",
          size = "small",
          labelWidth = "100px"

        )
      )
    ),
    fluidRow(
      #column(8,
      conditionalPanel(
        condition = "input.abs == false",
        # shinyWidgets::dropdownButton(
        #   inputId = "drop_allele",
        column(
          width = 2,
          actionBttn(
            inputId = "reset_input",
            label = "reset thresholds",
            style = "minimal",
            color = "success",
            size = "xs",
            icon = icon("undo")
          )
        ),
        br(),
        column(width = 12, uiOutput("thresh")),
        # icon = icon("sliders-h"),
        # label = "alleles thresholds",
        # circle = FALSE,
        # status = "btn-dark",
        # size = "sm"
        #)

      )
    ),
    conditionalPanel(
      condition = "input.abs == true",
      # shinyWidgets::dropdownButton(
      #   inputId = "drop_allele_abs",
      column(
        width = 2,
        actionBttn(
          inputId = "reset_input2",
          label = "reset thresholds",
          style = "minimal",
          color = "success",
          size = "xs",
          icon = icon("undo")
        )
      ),
      br(),
      column(width = 12, uiOutput("thresh_abs"))
      #   icon = icon("sliders-h"),
      #   label = "absolute alleles thresholds",
      #   circle = FALSE,
      #   status = "btn-dark",
      #   size = "sm"
      # )


      #)

    )
  ),
  mainPanel(tabsetPanel(id = "tabs"))
)


absolute_thresholds_dict <- list(
  "IGHV1-18G1" = c(
    "V1-18*01" = 0.001,
    "V1-18*04" = 0.001,
    "V1-18*03" = 0.0001,
    "V1-18*03" = 0.0001,
    "V1-18*01_T189G" = 1e-04,
    "V1-18*01_A190G" = 1e-04,
    'V1-18*01_A196G' = 1e-04
  ),
  "IGHV1-2G2" = c(
    "V1-2*02" = 0.001,
    "V1-2*04" = 0.0001,
    "V1-2*05" = 0.0001,
    "V1-2*06" = 0.001,
    "V1-2*07" = 0.0001,
    'V1-2*01' = 1e-04,
    'V1-2*03' = 1e-04,
    'V1-2*04_A119C' = 1e-04,
    'V1-2*02_A85C' = 1e-04,
    'V1-2*04_A85C' = 1e-04,
    'V1-2*02_T174C' = 1e-04,
    'V1-2*06_T174C' = 1e-04
  ),
  "IGHV1-24G3" = c("V1-24*01" = 0.0001),
  "IGHV1-3G4" = c(
    "V1-3*01" = 0.0001,
    "V1-3*02" = 0.0001,
    "V1-3*04"  = 0.001,
    "V1-3*05" = 0.0001,
    'V1-3*03' = 1e-04
  ),
  "IGHV1-45G5" = c(
    "V1-45*02" = 0.0001,
    "V1-45*03" = 0.0001,
    'V1-45*01' = 1e-04
  ),
  "IGHV1-46G6" = c(
    "V1-46*01" = 0.001,
    "V1-46*02" = 0.0001,
    "V1-46*03" = 0.001,
    "V1-46*04" = 0.001,
    'V1-46*01_A103C' = 1e-04,
    'V1-46*01_A119C' = 1e-04
  ),
  "IGHV1-58G7" = c(
    "V1-58*01" = 0.0001,
    "V1-58*02" = 0.0001,
    "V1-58*03" = 0.001
  ),
  "IGHV1-69G8" = c(
    "V1-69*01/V1-69D*01" = 0.001,
    "V1-69*06" = 0.001,
    "V1-69*09" = 0.001,
    "V1-69*04" = 0.001,
    "V1-69*02" = 0.001,
    "V1-69*04_C184T" = 0.001,
    "V1-69*18" = 0.001,
    "V1-69*19" = 0.0001,
    "V1-69*14" = 0.001,
    "V1-69*12" = 0.001,
    "V1-69*05" = 0.001,
    "V1-69*13" = 0.001,
    "V1-69*10" = 0.001,
    "V1-69*17" = 0.001,
    "V1-69*08" = 0.001,
    "V1-69*15" = 0.001,
    'V1-69*10_A54G' = 1e-04,
    'V1-69*02_A18G_A244G' = 1e-04,
    'V1-69*04_G48A_A163G_A244G' = 1e-04,
    'V1-69*10_A225G' = 1e-04,
    'V1-69*04_A225G' = 1e-04,
    'V1-69*01_G54A' = 1e-04,
    'V1-69*05_G54A' = 1e-04,
    'V1-69*01_A316G' = 1e-04,
    'V1-69*04_T191C' = 1e-04,
    'V1-69*01_C26T' = 1e-04,
    'V1-69*04_A163G' = 1e-04,
    'V1-69*14_G240A' = 1e-04,
    'V1-69*06_G240A' = 1e-04,
    'V1-69*04_A85C' = 1e-04,
    'V1-69*16' = 1e-04,
    'V1-69*11' = 1e-04,
    'V1-69*03' = 1e-04
  ),
  "IGHV1-69-2G9" = c("V1-69-2*01" = 0.0001),
  "IGHV1-8G10" = c(
    "V1-8*01" = 0.0001,
    "V1-8*02" = 0.001,
    "V1-8*03" = 0.001
  ),
  "IGHV2-26G11" = c(
    "V2-26*01" = 0.0001,
    "V2-26*02" = 0.0001,
    "V2-26*03" = 0.0001,
    "V2-26*04" = 0.0001
  ),
  "IGHV2-5G12" = c(
    "V2-5*01" = 0.001,
    "V2-5*02" = 0.0005,
    "V2-5*04" = 0.001,
    "V2-5*05" = 0.001,
    "V2-5*08" = 0.001,
    'V2-5*06' = 1e-04,
    'V2-5*09' = 1e-04
  ),
  "IGHV2-70G13" = c(
    "V2-70*01" = 0.0001,
    "V2-70*13" = 0.001,
    "V2-70*15" = 0.0001,
    "V2-70*17" = 0.0001,
    "V2-70*20" = 0.0001,
    "V2-70*04/V2-70D*04" = 0.007,
    'V2-70D*14' = 1e-04,
    'V2-70*19' = 1e-04,
    'V2-70*18' = 1e-04,
    'V2-70*16' = 1e-04,
    'V2-70*12' = 1e-04,
    'V2-70*11' = 1e-04,
    'V2-70*10' = 1e-04,
    'V2-70*08' = 1e-04,
    'V2-70*07' = 1e-04,
    'V2-70*06' = 1e-04,
    'V2-70*03' = 1e-04,
    'V2-70*02' = 1e-04
  ),
  "IGHV3-11G14" = c(
    "V3-11*01" = 0.001,
    "V3-11*04" = 0.001,
    "V3-11*05" = 0.001,
    "V3-11*06" = 0.0001,
    'V3-11*03' = 1e-04,
    'V3-11*01_A152G' = 1e-04,
    'V3-11*05_A152G' = 1e-04,
    'V3-11*01_A85C' = 1e-04,
    'V3-11*06_A152G' = 1e-04
  ),
  "IGHV3-13G15" = c(
    "V3-13*01" = 0.001,
    "V3-13*04" = 0.0001,
    "V3-13*05" = 0.0001,
    "V3-13*03" = 0.001,
    'V3-13*02' = 1e-04,
    'V3-13*01_G290A_T300C' = 1e-04,
    'V3-13*04_A85C' = 1e-04
  ),
  "IGHV3-15G16" = c(
    "V3-15*01" = 0.0001,
    "V3-15*01_A313T" = 0.001,
    "V3-15*02" = 0.0001,
    "V3-15*04" = 0.001,
    "V3-15*05" = 0.001,
    "V3-15*07" = 0.0001,
    'V3-15*03' = 1e-04,
    'V3-15*06' = 1e-04,
    'V3-15*08' = 1e-04,
    'V3-15*01_A152G' = 1e-04,
    'V3-15*01_A85C' = 1e-04,
    'V3-15*07_A152G' = 1e-04
  ),
  "IGHV3-20G17" = c(
    "V3-20*01" = 0.0001,
    "V3-20*04" = 0.0001,
    'V3-20*04_A152G' = 1e-4
  ),
  "IGHV3-21G18" = c(
    "V3-21*01" = 0.001,
    "V3-21*02" = 0.001,
    "V3-21*03" = 0.001,
    "V3-21*04" = 0.001,
    "V3-21*06" = 0.001,
    "V3-21*01_A184G_T190A_A191C" = 0.001,
    "V3-21*05" = 0.001,
    'V3-21*01_A85C' = 1e-04,
    'V3-21*01_A152G' = 1e-04,
    'V3-21*01_A152G_G278C' = 1e-04,
    'V3-21*01_G278C' = 1e-04,
    'V3-21*01_C34T_A40C_C90T_A112T_C114G_A119G_A210C' = 1e-04,
    'V3-21*01_A131C' = 1e-04
  ),
  "IGHV3-23G19" = c(
    "V3-23D*01/V3-23*01" = 0.001,
    "V3-23*04" = 0.001,
    'V3-23*01_A131C' = 1e-04,
    'V3-23*01_A303G' = 1e-04,
    'V3-23*01_T158G' = 1e-04,
    'V3-23*01_G239T' = 1e-04,
    'V3-23*01_T154G' = 1e-04,
    'V3-23*01_A118C' = 1e-04,
    'V3-23*01_A118G' = 1e-04,
    'V3-23*04_A152G' = 1e-04,
    'V3-23*01_A152G' = 1e-04,
    'V3-23*01_A85C' = 1e-04,
    'V3-23*05' = 1e-04,
    'V3-23*03' = 1e-04,
    'V3-23*02' = 1e-04
  ),
  "IGHV3-30G20" = c(
    "V3-30-3*01" = 0.0001,
    "V3-33*01" = 0.0001,
    "V3-30*03" = 0.0001,
    "V3-30-3*03/V3-30*04" = 0.001,
    "V3-33*06" = 0.001,
    "V3-30*01" = 0.0001,
    "V3-30*10" = 0.001,
    "V3-30*11" = 0.0001,
    "V3-33*08" = 0.001,
    "V3-30*02_G49A" = 0.001,
    "V3-30*19_T189C" = 0.0001,
    "V3-30*16" = 0.001,
    "V3-30*14" = 0.0001,
    "V3-33*05" = 0.001,
    "V3-30*07" = 0.001,
    "V3-30-3*02" = 0.0001,
    "V3-30*09" = 0.0001,
    "V3-33*03" = 0.0001,
    "V3-30*15" = 0.0001,
    "V3-30*19" = 0.0001,
    'V3-30-5*01' = 1e-04,
    'V3-30*18_C75G' = 1e-04,
    'V3-30-3*01_A131C' = 1e-04,
    'V3-30*18_A131C' = 1e-04,
    'V3-30*04_T158G' = 1e-04,
    'V3-33*01_T158G' = 1e-04,
    'V3-30-3*01_T158G' = 1e-04,
    'V3-30*18_T158G' = 1e-04,
    'V3-30*18_T167C' = 1e-04,
    'V3-30*04_A119C' = 1e-04,
    'V3-30-3*01_T154G' = 1e-04,
    'V3-30-3*01_A119C' = 1e-04,
    'V3-33*01_A119C' = 1e-04,
    'V3-30*01_A119C' = 1e-04,
    'V3-30*02_A275G' = 1e-04,
    'V3-30*18_G113C_C114T' = 1e-04,
    'V3-30*04_C201T' = 1e-04,
    'V3-30*04_A152G' = 1e-04,
    'V3-30-3*01_A152G' = 1e-04,
    'V3-30*04_T154G' = 1e-04,
    'V3-33*08_A152G' = 1e-04,
    'V3-33*01_T154G' = 1e-04,
    'V3-33*01_A85C' = 1e-04,
    'V3-30*02_A152G' = 1e-04,
    'V3-30*18_A152G' = 1e-04,
    'V3-30*18_T154G' = 1e-04,
    'V3-33*01_A152G' = 1e-04,
    'V3-30*18_A85C' = 1e-04,
    'V3-33*07' = 1e-04,
    'V3-33*04' = 1e-04,
    'V3-33*02' = 1e-04,
    'V3-30*18/V3-30-5*02' = 1e-04,
    'V3-30*17' = 1e-04,
    'V3-30*13' = 1e-04,
    'V3-30*12' = 1e-04,
    'V3-30*08' = 1e-04,
    'V3-30*06' = 1e-04,
    'V3-30*05' = 1e-04,
    'V3-30*02' = 1e-04
  ),
  "IGHV3-35G21" = c("V3-35*02" = 0.0001),
  "IGHV3-43G22" = c(
    "V3-43*02" = 0.0001,
    "V3-43D*04_G4A" = 0.0001,
    "V3-43*01" = 0.0001,
    "V3-43D*03" = 0.0001,
    "V3-43D*04" = 0.0001,
    'V3-43*01_T177A' = 1e-04,
    'V3-43*01_A85C' = 1e-04
  ),
  "IGHV3-48G23" = c(
    "V3-48*03" = 0.001,
    "V3-48*02" = 0.001,
    "V3-48*01" = 0.001,
    "V3-48*04" = 0.001,
    'V3-48*02_A85C' = 1e-04,
    'V3-48*03_A85C' = 1e-04,
    'V3-48*03_T154G' = 1e-04,
    'V3-48*01_A85C' = 1e-04,
    'V3-48*02_T154G' = 1e-04,
    'V3-48*03_A152G' = 1e-04,
    'V3-48*01_T154G' = 1e-04,
    'V3-48*02_A152G' = 1e-04,
    'V3-48*04_A152G' = 1e-04,
    'V3-48*01_A152G' = 1e-04,
    'V3-48*01_A39C' = 1e-04
  ),
  "IGHV3-49G24" = c(
    "V3-49*05" = 0.0001,
    "V3-49*04" = 0.0001,
    "V3-49*03" = 0.0001,
    'V3-49*01' = 1e-04,
    'V3-49*02' = 1e-04,
    'V3-49*04_A152G' = 1e-04,
    'V3-49*05_A152G' = 1e-04,
    'V3-49*03_A152G' = 1e-04
  ),
  "IGHV3-53G25" = c(
    "V3-53*01" = 0.0001,
    "V3-53*04" = 0.0001,
    "V3-66*0/V3-66*041" = 0.0001,
    "V3-66*02" = 0.0001,
    "V3-53*02" = 0.001,
    "V3-66*03" = 0.0001,
    "V3-66*02_G303A" = 0.0001,
    'V3-53*04_C241G_G273A_C300T' = 1e-04,
    'V3-66*02_T158G' = 1e-04,
    'V3-66*01_A193G' = 1e-04,
    'V3-53*04_G81A' = 1e-04,
    'V3-66*02_G51T' = 1e-04,
    'V3-53*04_G37A' = 1e-04,
    'V3-53*02_G88A' = 1e-04,
    'V3-66*02_A152G' = 1e-04,
    'V3-66*01_A152G' = 1e-04,
    'V3-53*01_A118C' = 1e-04,
    'V3-53*05' = 1e-04,
    'V3-53*03' = 1e-04
  ),
  "IGHV3-62G26" = c('V3-62*04'=1e-04),
  "IGHV3-64G27" = c(
    "V3-64D*06" = 0.0001,
    "V3-64*02" = 0.0001,
    "V3-64*01" = 0.0001,
    "V3-64D*09" = 0.0001,
    "V3-64D*08" = 0.0001,
    "V3-64*07" = 0.001,
    'V3-64*05' = 1e-04,
    'V3-64*04' = 1e-04,
    'V3-64*03' = 1e-04
  ),
  "IGHV3-7G28" = c(
    "V3-7*03" = 0.001,
    "V3-7*01" = 0.001,
    "V3-7*04" = 0.0001,
    "V3-7*02" = 0.001,
    "V3-7*05" = 0.0001,
    'V3-7*01_A152G' = 1e-04,
    'V3-7*01_A85C' = 1e-04,
    'V3-7*01_T154G' = 1e-04,
    'V3-7*03_A152G' = 1e-04,
    'V3-7*01_T158G' = 1e-04,
    'V3-7*03_A131C' = 1e-04
  ),
  "IGHV3-72G29" = c("V3-72*01" = 0.0001),
  "IGHV3-73G30" = c(
    "V3-73*01" = 0.0001,
    "V3-73*02" = 0.0001,
    'V3-73*02_A152G' = 1e-04
  ),
  "IGHV3-74G31" = c(
    "V3-74*01" = 0.0001,
    "V3-74*03" = 0.001,
    'V3-74*02' = 1e-04
  ),
  "IGHV3-9G32" = c(
    "V3-9*01" = 0.0001,
    "V3-9*02"  = 0.001,
    "V3-9*01_T307C" = 0.0001,
    "V3-9*03" = 0.0001,
    'V3-9*01_A152G' = 1e-04,
    'V3-9*01_A85C' = 1e-04,
    'V3-9*01_A119C' = 1e-04
  ),
  "IGHV4-28G33" = c(
    "V4-28*05"  = 0.0001,
    "V4-28*07" = 0.0001,
    "V4-28*02" = 0.0001,
    "V4-28*06" = 0.0001,
    'V4-28*01' = 1e-04,
    'V4-28*03' = 1e-04,
    'V4-28*04' = 1e-04
  ),
  "IGHV4-30-2G34" = c(
    "V4-30-2*01" = 0.0001,
    "V4-30-4*07" = 0.0001,
    "V4-30-2*06" = 0.0001,
    "V4-30-2*05" = 0.0001,
    'V4-30-2*01_G139C' = 1e-04,
    'V4-30-2*01_C285T' = 1e-04,
    'V4-30-2*03' = 1e-04,
    'V4-30-2*02' = 1e-04
  ),
  "IGHV4-30-4G35" = c(
    "V4-30-4*01" = 0.001,
    "V4-30-4*01_A70G_A107G" = 0.0001,
    "V4-30-4*08" = 0.001,
    'V4-30-4*02' = 1e-04,
    'V4-30-4*03' = 1e-04,
    'V4-30-4*04' = 1e-04
  ),
  "IGHV4-31G36" = c(
    "V4-31*03" = 0.0001,
    "V4-31*01" = 0.0001,
    "V4-31*11" = 0.0001,
    "V4-30-4*08" = 0.001,
    "V4-31*02" = 0.001,
    'V4-31*04' = 1e-04,
    'V4-31*05' = 1e-04,
    'V4-31*06' = 1e-04,
    'V4-31*07' = 1e-04,
    'V4-31*08' = 1e-04,
    'V4-31*09' = 1e-04,
    'V4-31*10' = 1e-04,
    'V4-31*03_T76C' = 1e-04,
    'V4-31*11_G4C_G21C_C25T_A113C' = 1e-04,
    'V4-31*11_G70C' = 1e-04,
    'V4-31*03_T285C' = 1e-04
  ),
  "IGHV4-34G37" = c(
    "V4-34*01"  = 0.0001,
    "V4-34*12" = 0.001,
    "V4-34*02" = 0.002,
    'V4-34*03' = 1e-04,
    'V4-34*04' = 1e-04,
    'V4-34*05' = 1e-04,
    'V4-34*06' = 1e-04,
    'V4-34*07' = 1e-04,
    'V4-34*08' = 1e-04,
    'V4-34*01_T208C' = 1e-04,
    'V4-34*01_A220G_A225G' = 1e-04,
    'V4-34*01_A225G' = 1e-04,
    'V4-34*01_T300C' = 1e-04,
    'V4-34*01_A220G' = 1e-04
  ),
  "IGHV4-34G38" = c(
    'V4-59*12_G129C_A147G_T163G_T165A_T169A_T172C_C174T_G189A_C207G_C300T'=1e-04,'V4-34*10'=1e-04,'V4-34*09'=1e-04
  ),
  "IGHV4-34G39" = c('V4-34*11'=1e-04),
  "IGHV4-38-2G40" = c(
    "V4-38-2*01" = 0.0001,
    "V4-38-2*02" = 0.001,
    'V4-38-2*02_G246A' = 1e-04,
    'V4-38-2*02_A291G_C300T' = 1e-04,
    'V4-38-2*02_A220G' = 1e-04,
    'V4-38-2*02_A291G' = 1e-04
  ),
  "IGHV4-39G41" = c(
    "V4-39*01/V4-39*02_C258G" = 0.0001,
    "V4-39*07" = 0.001,
    "V4-39*01_C66G" = 0.0001,
    "V4-39*07_C288A" = 0.001,
    "V4-39*05" = 0.001,
    "V4-39*02" = 0.0001,
    'V4-39*03' = 1e-04,
    'V4-39*06' = 1e-04,
    'V4-39*01_C66G_A142G_G144C' = 1e-04,
    'V4-39*01_G315A' = 1e-04,
    'V4-39*01_A200C' = 1e-04,
    'V4-39*01_A202C' = 1e-04,
    'V4-39*01_A220G' = 1e-04,
    'V4-39*01_A291G' = 1e-04,
    'V4-39*07_A220G_A225G' = 1e-04,
    'V4-39*07_C300T' = 1e-04,
    'V4-39*01_C66G_G315A' = 1e-04
  ),
  "IGHV4-4G42" = c(
    "V4-4*02" = 0.0001,
    "V4-4*03" = 0.001,
    "V4-4*02_A106G" = 0.001,
    "V4-4*01" = 0.001,
    'V4-4*04' = 1e-04,
    'V4-4*05' = 1e-04,
    'V4-4*02_C20T' = 1e-04
  ),
  "IGHV4-4G43" = c(
    "V4-4*07" = 0.0001,
    'V4-4*07_T208C' = 1e-04,
    'V4-4*07_A70G' = 1e-04,
    'V4-59*10' = 1e-04
  ),
  "IGHV4-4G44" = c(
    "V4-59*01" = 0.0001,
    "V4-59*08" = 0.0001,
    "V4-59*11" = 0.001,
    "V4-59*07" = 0.0001,
    "V4-59*12" = 0.001,
    "V4-59*01_G267A" = 0.0001,
    "V4-59*02" = 0.005,
    "V4-59*13" = 0.0001,
    "V4-59*04" = 0.0001,
    'V4-59*01_A220G_A225G' = 1e-04,
    'V4-59*12_C300T' = 1e-04,
    'V4-59*08_T208C' = 1e-04,
    'V4-59*01_T208C' = 1e-04,
    'V4-59*01_T76C' = 1e-04,
    'V4-59*05' = 1e-04,
    'V4-59*03' = 1e-04,
    'V4-4*09' = 1e-04,
    'V4-4*08' = 1e-04
  ),
  "IGHV4-59G45" = c('V4-59*06'=1e-04),
  "IGHV4-61G46" = c(
    "V4-61*01" = 0.0001,
    "V4-61*08" = 0.001,
    "V4-61*01_A41G" = 0.0001,
    "V4-61*03" = 0.0001,
    'V4-61*04' = 1e-04,
    'V4-61*05' = 1e-04,
    'V4-61*10' = 1e-04
  ),
  "IGHV4-61G47" = c(
    "V4-61*02" = 0.001,
    "V4-61*02_A234G" = 0.001,
    'V4-61*09' = 1e-04,
    'V4-61*02_A291G' = 1e-04
  ),
  "IGHV5-10-1G48" = c(
    "V5-10-1*03" = 0.001,
    "V5-10-1*01" = 0.001,
    "V5-10-1*02" = 0.0001,
    'V5-10-1*04' = 1e-04
  ),
  "IGHV5-51G49" = c(
    "V5-51*01" = 0.0001,
    "V5-51*03" = 0.0001,
    "V5-51*04" = 0.0001,
    "V5-51*06" = 0.0001,
    "V5-51*07" = 0.0001,
    'V5-51*02' = 1e-04
  ),
  "IGHV6-1G50" = c(
    "V6-1*01" =  0.0001,
    "V6-1*01_T91C" =  0.0001,
    'V6-1*02' = 1e-04
  ),
  "IGHV7-4-1G51" = c(
    "V7-4-1*01" =  0.0001,
    "V7-4-1*02" =  0.0001,
    'V7-4-1*03' = 1e-04,
    'V7-4-1*04' = 1e-04,
    'V7-4-1*05' = 1e-04,
    'V7-4-1*02_A200T' = 1e-04
  )
)

absolute_thresholds_dict <- read.delim("alleles_db.tsv", stringsAsFactors = F)

absolute_thresholds_dict <- sapply(unique(absolute_thresholds_dict$func_group), function(x){
  tmp <- absolute_thresholds_dict[absolute_thresholds_dict$func_group==x,]
  setNames(tmp$thresh,gsub("IGH","",tmp$or_allele))
})

server <- function(input, output, session) {
  input_vals <-
    reactiveValues(
      tabs_count = 0,
      g_group = "IGHV3-30-3G20",
      g = allele_db %>% dplyr::filter(gene_group == "IGHV3-30-3G20") %>% pull(gene) %>% unique(),
      v_gene_cut = "IGHV3-30-3G20",
      allele_thresh = NULL,
      allele_thresh_names = NULL,
      allele_thresh_abs = NULL,
      allele_thresh_names_abs = NULL,
      tabs = NULL,
      tabs_left = NULL
    )
  allele_num <- reactiveValues(len = 0)
  a_thresh <- reactiveValues(a_thresh = NULL, a_thresh_abs = NULL)

  query <- eventReactive(session$clientData$url_search, {
    parseQueryString(session$clientData$url_search)
  })

  observeEvent(query()$g_group, {
    if (!is.null(query()$g_group)) {
      g_group <- strsplit(query()$g_group, "\"")[[1]][2]
      g <-
        allele_db %>% dplyr::filter(gene_group == g_group) %>% pull(gene) %>% unique()
      v_gene_cut <-
        ifelse(grepl("G", g_group), g_group, func_groups[as.character(g_group)])
      input_vals$g_group <- g_group
      input_vals$g <- g
      input_vals$v_gene_cut <- v_gene_cut
    }
  })


  tmp_allele_db =
    reactive({
      tmp_allele_db = allele_db %>% dplyr::filter(grepl(as.character(input_vals$g_group), new_allele)) %>%
        dplyr::group_by(new_allele) %>% dplyr::summarise(or_allele = paste0(sort(or_allele), collapse = "/"))
    })

  or_allele = reactive({
    or_allele <-
      setNames(gsub(chain, "", as.character(tmp_allele_db()$or_allele)), as.character(gsub(
        paste0(input_vals$g_group, "[*]"),
        "",
        tmp_allele_db()$new_allele
      )))

  })


  data_cut <- reactive({
    tmp <- data %>%
      dplyr::filter(v_gene == input_vals$v_gene_cut,!is.na(v_allele),
                    group_plot == 1) %>%
      dplyr::ungroup() %>% dplyr::mutate(v_allele_axis = or_allele()[v_allele])

    input_vals$allele_thresh <-
      setNames(rep(0.5, length(unique(
        tmp$v_allele_axis
      ))), unique(tmp$v_allele_axis))

    input_vals$allele_thresh_abs <-
      absolute_thresholds_dict[[input_vals$g_group]][unique(tmp$v_allele_axis)]

    if(length(input_vals$allele_thresh_abs)<length(unique(tmp$v_allele_axis))){
      tmp_names <- unique(tmp$v_allele_axis)
      tmp_names <-   tmp_names[ !tmp_names %in% names(absolute_thresholds_dict[[input_vals$g_group]])]
      input_vals$allele_thresh_abs <- c(input_vals$allele_thresh_abs, setNames(rep(1e-4, length(tmp_names)), tmp_names))
    }


    return(tmp)

  })

  # observeEvent(data_cut(), {
  #
  #   input_vals$allele_thresh <-
  #     setNames(rep(0.5, length(unique(
  #       data_cut()$v_allele_axis
  #     ))), unique(data_cut()$v_allele_axis))
  #   input_vals$allele_thresh_abs <-
  #     absolute_thresholds_dict[[input_vals$g_group]]
  #
  #   if(length(input_vals$allele_thresh_abs)<length(unique(data_cut()$v_allele_axis))){
  #     tmp_names <- unique(data_cut()$v_allele_axis)
  #     tmp_names <-   tmp_names[ !tmp_names %in% names(absolute_thresholds_dict[[input_vals$g_group]])]
  #     input_vals$allele_thresh_abs <- c(input_vals$allele_thresh_abs, setNames(rep(1e-4, length(tmp_names)), tmp_names))
  #   }else{
  #
  #     tmp_names <- unique(data_cut()$v_allele_axis)
  #     keep <-   names(input_vals$allele_thresh_abs)[names(input_vals$allele_thresh_abs) %in% tmp_names]
  #     input_vals$allele_thresh_abs <- input_vals$allele_thresh_abs[keep]
  #   }
  #
  # })


  output$thresh <-
    renderUI({
      a_thresh$a_thresh <- input_vals$allele_thresh
      input_vals$allele_thresh_names <-
        setNames(paste0("allele", 1:length(unique(
          data_cut()$v_allele_axis
        ))),
        names(reactiveValuesToList(input_vals)$allele_thresh))
      allele_num$len <- length(unique(data_cut()$v_allele_axis))
      print(allele_num$len)
      div(class = "dynamicSI",
          lapply(1:length(unique(
            data_cut()$v_allele_axis
          )), function(i) {
            lab <- unique(data_cut()$v_allele_axis)[i]
            val <- as.numeric(input_vals$allele_thresh[lab])
            div(
              style = "display: inline-block;font-size:80%;",
              numericInput(
                inputId = paste0("allele", i),
                label = lab,
                value = val,
                min = 0.5,
                max = 20,
                step = 0.1,
                width = "100px"
              )
            )
          }))

      #lab <- unique(data_cut()$v_allele_axis)[1]
      #val <- as.numeric(a_thresh[lab])
      #numericInput(inputId = "dat", label = "dfkjs", value = 0.5, min = 0.5, max = 20, step = 0.1)
    })

  outputOptions(output, "thresh", suspendWhenHidden = FALSE)#, priority = 20)


  output$thresh_abs <-
    renderUI({
      a_thresh$a_thresh_abs <- input_vals$allele_thresh_abs
      input_vals$allele_thresh_names_abs <-
        setNames(
          paste0("allele_abs", 1:length(unique(
            data_cut()$v_allele_axis
          ))),
          names(reactiveValuesToList(input_vals)$allele_thresh_abs)
        )
      allele_num$len <- length(unique(data_cut()$v_allele_axis))
      div(class = "dynamicSI",
          lapply(1:length(unique(
            data_cut()$v_allele_axis
          )), function(i) {
            lab <- unique(data_cut()$v_allele_axis)[i]
            val <- as.numeric(input_vals$allele_thresh_abs[lab])
            div(
              style = "display: inline-block;font-size:80%;",
              numericInput(
                inputId = paste0("allele_abs", i),
                label = lab,
                value = val,
                min = 1e-6,
                max = 0.01,
                step = 0.0001
              )
            )
          }))

      #lab <- unique(data_cut()$v_allele_axis)[1]
      #val <- as.numeric(a_thresh[lab])
      #numericInput(inputId = "dat", label = "dfkjs", value = 0.5, min = 0.5, max = 20, step = 0.1)
    })

  outputOptions(output, "thresh_abs", suspendWhenHidden = FALSE)#, priority = 20)


  observeEvent(input$lastSelectId, once = F, {
    vals <- reactiveValuesToList(input)
    vals <- vals[grep("allele[0-9]", names(vals))]
    if (length(vals) != 0) {
      names_ <- reactiveValuesToList(input_vals)$allele_thresh
      names_input <-
        reactiveValuesToList(input_vals)$allele_thresh_names

      for (x in 1:length(vals)) {
        a_names = names(names_)[x]

        a_thresh$a_thresh[a_names] <-
          as.numeric(vals[[names_input[a_names]]])
      }
    }
  })


  observeEvent(input$lastSelectId, once = F, {
    vals <- reactiveValuesToList(input)
    vals <- vals[grep("allele_abs[0-9]", names(vals))]
    if (length(vals) != 0) {
      names_ <- reactiveValuesToList(input_vals)$allele_thresh_abs
      names_input <-
        reactiveValuesToList(input_vals)$allele_thresh_names_abs

      for (x in 1:length(vals)) {
        a_names = names(names_)[x]
        a_thresh$a_thresh_abs[a_names] <-
          as.numeric(vals[[names_input[a_names]]])
      }
    }
  })

  observe({
    input$lastSelect

    if (!is.null(input$lastSelectId)) {
      cat("lastSelectId:", input$lastSelectId, "\n")
      cat("Selection:", input[[input$lastSelectId]], "\n\n")
    }
  })

  #observeEvent(, once = F,{

  data_gene <- eventReactive(list(#input$lastSelect
    a_thresh$a_thresh,
    a_thresh$a_thresh_abs,
    input$mut),
    ignoreInit = TRUE,
    {
      tmp <- isolate(data_cut())
      print(tmp)
      if (!input$abs) {
        tmp <-
          tmp %>% dplyr::group_by(v_allele_axis) %>% dplyr::filter(freq >=  a_thresh$a_thresh[v_allele_axis] /
                                                                     100, mut == input$mut)
      } else{
        #print("ok")
        #print(a_thresh$a_thresh_abs)
        tmp <-
          tmp %>% dplyr::group_by(v_allele_axis) %>% dplyr::filter(freq2 >=  a_thresh$a_thresh_abs[v_allele_axis], mut == input$mut)
      }
      print(tmp %>% filter(subject == "P11_I39_S1") %>% as.data.frame())
      tmp <- tmp %>% dplyr::arrange(desc(freq)) %>%
        dplyr::group_by(subject, v_gene) %>% dplyr::mutate(
          zygousity_state = as.numeric(length(unique(v_allele_axis))),
          v_alleles = paste0(1:unique(zygousity_state), " - ", or_allele()[v_allele[1:unique(zygousity_state)]], collapse = ";"),
          v_alleles_abc = paste0(sort(or_allele()[v_allele[1:unique(zygousity_state)]]), collapse = ";")
        ) %>% arrange(subject)
      #tmp <- tmp %>% dplyr::group_by(subject, zygousity_state) %>%
      #  dplyr::mutate(loc_state = loc <= zygousity_state) %>% dplyr::filter(loc_state) %>% ungroup()
      print(a_thresh$a_thresh_abs)
      print(tmp %>% filter(subject == "P11_I39_S1") %>% as.data.frame())
      # if (input$naive) {
      #   tmp <- tmp %>% dplyr::filter(!grepl("P4",subject))
      # } else{
      #   tmp <- tmp
      # }

      return(tmp)
    })

  #})
  plt_data_het = eventReactive(input$hetro, {
    #req(data_gene())
    tmp_plot <- isolate(data_gene())
    if (input$hetro) {
      tmp_plot <- tmp_plot %>% dplyr::filter(J6 == 1)
    } else{
      tmp_plot <- tmp_plot
    }
    return(tmp_plot)
  })

  observe({
    req(data_gene())
    states <- sort(as.numeric(unique(data_gene()$zygousity_state)))

    if (input_vals$tabs_count < 1) {
      input_vals$tabs_count <- 1
      lapply(states, function(i) {
        insertTab(
          session = session,
          inputId = "tabs",
          tabPanel(
            title = uiOutput(paste0("State_", i)),
            value = paste0("tab_", i),
            downloadButton(
              outputId = paste0("downloadData_", i),
              label = "Download frequency table",
              style = "width:50%;",
              class = "butt"
            ),
            br(),
            br(),
            #tags$head(tags$style(type = "text/css", "#",paste0('scatter', i)," {width:95vh !important;}")),
            plotlyOutput(paste0('scatter', i), width = "150%"),
            br(),
            br(),
            #tags$head(tags$style(type = "text/css", "#",paste0('scatter', i)," {height:95vh !important;}")),
            plotlyOutput(paste0('hover', i))

          )
        )
        output[[paste0("State_", i)]] = renderText({
          paste0("<font color=\"#000000\">",
                 paste0("State ", i),
                 "</font>")
        })
      })

      input_vals$tabs <- paste0("State ", states)
      input_vals$tabs_left <- paste0("State ", states)
    } else{
      lapply(input_vals$tabs, function(i) {
        if (i %in% paste0("State ", states)) {
          output[[gsub(" ", "_", i)]] = renderText({
            paste0("<font color=\"#000000\">", i, "</font>")
          })
        } else{
          output[[gsub(" ", "_", i)]] = renderText({
            paste0("<font color=\"#eeeeee\">", i, "</font>")
          })

        }

      })

    }

    updateTabsetPanel(
      session = session,
      inputId = "tabs",
      selected = paste0("tab_", 1)
    )


    Map(function(state) {
      state_val <- reactiveValues(i = as.numeric(state))


      plt_data = eventReactive(c(input$hetro, input$abs), {
        req(plt_data_het())
        req((input$hetro == T | input$hetro == F))

        tmp_plot <-  isolate(plt_data_het()) %>%
          dplyr::filter(zygousity_state == as.numeric(state_val$i),
                        is.na(j_call)) %>% ungroup() %>% dplyr::rowwise() %>% dplyr::mutate(
                          v_alleles_p = v_alleles_abc,
                          v_alleles_p = gsub(";", "\n", v_alleles_p),
                          text = HTML(
                            paste(
                              '</br>Project: ',
                              project,
                              '</br>Subject: ',
                              subject,
                              '</br>Alleles: ',
                              v_alleles,
                              '</br># assignments: ',
                              count,
                              '</br>Group normalization freq.: ',
                              freq,
                              '</br>Rep. normalization freq.: ',
                              as.character(freq2),
                              '</br>',
                              J6_TAG
                            )
                          )
                        ) %>% ungroup()




        return(tmp_plot)
      })
      # plt_data_haplo <- reactive({
      #   tmp_haplo <- data %>%
      #     dplyr::filter(v_gene == input_vals$v_gene_cut,
      #                   !is.na(v_allele),
      #                   group_plot == 2) %>%
      #     ungroup()
      #
      #   tmp_haplo <- tmp_haplo %>% dplyr::group_by(subject) %>%
      #     dplyr::filter(v_allele %in% plt_data()$v_allele[plt_data()$subject ==
      #                                                       unique(subject)]) %>%
      #     dplyr::mutate(
      #       v_alleles_p = unique(plt_data()$v_alleles_p[plt_data()$subject ==
      #                                                     unique(subject)]),
      #       v_allele_axis = or_allele()[v_allele]
      #     )
      #   return(tmp_haplo)
      # })

      colors <- eventReactive(plt_data(), {
        n_alleles <- length(unique(plt_data()$v_alleles_p))
        if (n_alleles <= 4)
          setNames(cols[1:n_alleles], unique(plt_data()$v_alleles_p))
        else
          setNames(pal %>% usecol(n = n_alleles),
                   unique(plt_data()$v_alleles_p))
      })

      plt_data_haplo <- reactiveValues(db = NULL)

      filteredScores <- eventReactive(plt_data(), {
        req((input$hetro == T | input$hetro == F))
        req((input$abs == T | input$abs == F))
        tmp <- isolate(data_cut())
        if (!input$abs) {
          tmp <-
            tmp %>% dplyr::group_by(v_allele) %>% dplyr::filter(freq >=  a_thresh$a_thresh[v_allele_axis] /
                                                                  100, mut == input$mut)
        } else{
          tmp <-
            tmp %>% dplyr::group_by(v_allele) %>% dplyr::filter(freq2 >=  a_thresh$a_thresh_abs[v_allele_axis], mut == input$mut)
        }
        tmp <- tmp %>% dplyr::arrange(desc(freq)) %>%
          dplyr::group_by(subject, v_gene) %>% dplyr::mutate(
            zygousity_state = as.numeric(length(unique(v_allele))),
            v_alleles = paste0(
              1:unique(zygousity_state),
              " - ",
              v_allele_axis[1:unique(zygousity_state)],
              collapse = ";"
            ),
            v_alleles_abc = paste0(sort(v_allele_axis[1:unique(zygousity_state)]), collapse = ";")
          ) %>% arrange(subject)
        #tmp <- tmp %>% dplyr::group_by(subject, zygousity_state) %>%
        #  dplyr::mutate(loc_state = loc <= zygousity_state) %>% dplyr::filter(loc_state) %>% ungroup()

        if (input$hetro) {
          tmp <- tmp %>% dplyr::filter(J6 == 1)
        } else{
          tmp <- tmp
        }

        # if (input$naive) {
        #   tmp <- tmp %>% dplyr::filter(!grepl("P4",subject))
        # } else{
        #   tmp <- tmp
        # }

        plot_data <-  tmp %>%
          dplyr::filter(zygousity_state == as.numeric(state_val$i),
                        is.na(j_call)) %>% ungroup() %>% dplyr::rowwise() %>% dplyr::mutate(
                          v_alleles_p = v_alleles_abc,
                          v_alleles_p = gsub(";", "\n", v_alleles_p),
                          text = HTML(
                            paste(
                              '</br>Project: ',
                              project,
                              '</br>Subject: ',
                              subject,
                              '</br>Alleles: ',
                              gsub(";", "\n", v_alleles),
                              '</br># assignments: ',
                              count,
                              '</br>Group norm. freq.: ',
                              freq,
                              '</br>Rep. norm. freq.:',
                              as.character(freq2),
                              '</br>',
                              J6_TAG
                            )
                          )
                        ) %>% ungroup()

        print(plot_data %>% filter(subject == "P11_I39_S1") %>% as.data.frame())
        tmp_haplo <- data %>%
          dplyr::filter(
            v_gene == input_vals$v_gene_cut,!is.na(v_allele),
            group_plot == 2,
            mut == input$mut
          ) %>%
          ungroup()

        tmp_haplo <- tmp_haplo %>% dplyr::group_by(subject) %>%
          dplyr::filter(v_allele %in% plot_data$v_allele[plot_data$subject ==
                                                           unique(subject)]) %>%
          dplyr::mutate(
            v_alleles_p = unique(plot_data$v_alleles_p[plot_data$subject ==
                                                         unique(subject)]),
            v_allele_axis = or_allele()[v_allele]
          )

        plt_data_haplo$db <- tmp_haplo
        shiny:::flushReact()
        #plot_data <- plt_data()

        loc2 <-
          setNames(1:length(unique(plot_data$v_allele_axis)), sort(unique(plot_data$v_allele_axis)))

        plot_data$loc2 <- loc2[plot_data$v_allele_axis]

        validate(need(
          nrow(plot_data) > 0,
          paste0("No rearrangments for this state")
        ))

        if (state_val$i != 1 &
            length(unique(plot_data$v_alleles_p)) != 1) {
          loc_jitter <- list()
          for (ii in 1:length(unique(plot_data$loc2))) {
            loc_c <- as.character(unique(plot_data$loc2)[ii])
            loc_jitter[[loc_c]] <-
              seq(0, 0.5, length.out = length(unique(plot_data$v_alleles_p[plot_data$loc2 ==
                                                                             unique(plot_data$loc2)[ii]])))

            loc_jitter[[loc_c]]  <-
              setNames(loc_jitter[[loc_c]] , sort(unique(plot_data$v_alleles_p[plot_data$loc2 ==
                                                                                 unique(plot_data$loc2)[ii]])))
          }


          plot_data <-
            plot_data %>% dplyr::arrange(loc2, v_alleles_p) %>% dplyr::group_by(loc2) %>%
            dplyr::mutate(loc_plot = loc2 + loc_jitter[[as.character(unique(loc2))]][v_alleles_p],) %>% ungroup()

          if (length(plot_data$loc_plot))
            plot_data <-
            plot_data %>% dplyr::mutate(jitter_offset = jitter(loc_plot, factor = 1))
        } else{
          plot_data <-
            plot_data %>% dplyr::arrange(loc2, v_alleles_p) %>% dplyr::group_by(loc2) %>%
            dplyr::mutate(loc_plot = loc2,) %>% ungroup()
          if (length(plot_data$loc_plot))
            plot_data <-
              plot_data %>% dplyr::mutate(jitter_offset = jitter(loc_plot, factor = 1))
        }

        tickvals_tmp <-
          plot_data %>% dplyr::pull(loc_plot) %>% unique() %>% sort()

        tickvals <- c()

        for (i in 1:length(loc2)) {
          tickvals <-
            c(tickvals, mean(tickvals_tmp[floor(tickvals_tmp) == i]))
        }


        ticktext <-
          plot_data %>% dplyr::pull(v_allele_axis) %>% unique() %>% sort()

        validate(need(
          nrow(plot_data) > 0,
          paste0("No rearrangments for this state")
        ))

        plot_data_final <- plot_data %>%
          highlight_key(., ~ subject)

        #a_names = isolate(input_vals$allele_thresh)

        #a_thresh = sapply(1:length(input_vals$allele_thresh), function(x) input[[paste0("allele",x)]])
        #names(a_thresh) <- names(input_vals$allele_thresh)
        vals <- reactiveValuesToList(a_thresh)$a_thresh
        allele_thresh_state <- vals[names(vals) %in% names(loc2)]

        names_tmp <- names(loc2)[names(loc2) %in% names(vals)]

        allele_thresh_state <-
          allele_thresh_state[order(factor(names(allele_thresh_state), levels = names_tmp))]

        vals <- reactiveValuesToList(a_thresh)$a_thresh_abs
        allele_thresh_state_absolute <-
          vals[names(vals) %in% names(loc2)]

        names_tmp <- names(loc2)[names(loc2) %in% names(vals)]

        allele_thresh_state_absolute <-
          allele_thresh_state_absolute[order(factor(names(allele_thresh_state_absolute), levels =
                                                      names_tmp))]

        colors_line <-
          setNames(c("red", "green", "blue", "purple", "black"),
                   unique(allele_thresh_state))

        colors_line_abs <-
          setNames(
            c("red", "green", "blue", "purple", "black"),
            unique(allele_thresh_state_absolute)
          )
        print(table(plot_data$project))
        plotly1 <-
          plot_data_final %>%
          plotly::plot_ly(colors = colors()) %>%
          plotly::add_trace(
            type = "scatter",
            x = ~ (jitter_offset),
            y = ~ freq,
            text = ~ text,
            symbol = ~ project,
            mode = 'markers',
            marker = list(color = "grey", size = 12),
            showlegend = TRUE,
            opacity = 0.9,
            hoverinfo = 'text',
            legendgroup = ~ project
          ) %>%
          plotly::add_trace(
            type = "scatter",
            x = ~ (jitter_offset),
            y = ~ freq,
            text = ~ text,
            color = ~ v_alleles_p,
            colors = colors(),
            mode = 'markers',
            showlegend = FALSE,
            opacity = 0.8,
            hoverinfo = 'none',
            legendgroup = ~ v_alleles_p
          ) %>%
          plotly::add_trace(
            x = ~ as.numeric(loc_plot),
            y = ~ freq,
            color = ~ v_alleles_p,
            colors = colors(),
            type = "box",
            hoverinfo = "none",
            showlegend = FALSE,
            fillcolor = "transparent",
            name = ~ v_alleles_p,
            legendgroup = ~ v_alleles_p
          ) %>%
          plotly::layout(
            hovermode = 'closest',
            showlegend = T,
            shapes = lapply(1:length(names(
              allele_thresh_state
            )), function(ia) {
              a = names(allele_thresh_state)[ia]
              xx = tickvals
              hline(
                unname(allele_thresh_state[a]) / 100,
                x0 = xx[ia] -
                  0.25,
                x1 = xx[ia] + 0.25,
                color = colors_line[as.character(allele_thresh_state[a])]
              )
            }),
            legend = list(
              tracegroupgap = 20,
              title = list(text =
                             '<b>  </b>'),
              orientation = "V"
            ),
            xaxis = list(
              title = paste0(input_vals$g, " Alleles"),
              autotick = F,
              tickmode = "array",
              tickvals = tickvals,
              ticktext = ticktext
            ),
            yaxis = list(title = "Group\nnormalization",
                         range = c(0, 1.05))
          )


        plotly2 <-
          plot_data_final %>%
          plotly::plot_ly(colors = colors()) %>%
          plotly::add_trace(
            type = "scatter",
            x = ~ (jitter_offset),
            y = ~ freq2,
            text = ~ text,
            symbol = ~ project,
            mode = 'markers',
            marker = list(color = "grey", size = 12),
            showlegend = FALSE,
            opacity = 0.9,
            hoverinfo = 'text',
            legendgroup = ~ project
          ) %>%
          plotly::add_trace(
            type = "scatter",
            x = ~ (jitter_offset),
            y = ~ freq2,
            text = ~ text,
            color = ~ v_alleles_p,
            colors = colors(),
            mode = 'markers',
            showlegend = FALSE,
            opacity = 0.8,
            hoverinfo = 'none',
            legendgroup = ~ v_alleles_p
          ) %>%
          plotly::add_trace(
            x = ~ as.numeric(loc_plot),
            y = ~ freq2,
            color = ~ v_alleles_p,
            colors = colors(),
            type = "box",
            hoverinfo = "none",
            fillcolor = "transparent",
            showlegend = FALSE,
            name = ~ v_alleles_p,
            legendgroup = ~ v_alleles_p
          ) %>%
          plotly::layout(
            hovermode = 'closest',
            shapes = lapply(1:length(
              names(allele_thresh_state_absolute)
            ), function(ia) {
              a = names(allele_thresh_state_absolute)[ia]
              xx = tickvals
              hline(
                unname(allele_thresh_state_absolute[a]),
                x0 = xx[ia] -
                  0.25,
                x1 = xx[ia] + 0.25,
                color = colors_line_abs[as.character(allele_thresh_state_absolute[a])]
              )
            }),
            legend = list(
              tracegroupgap = 20,
              title = list(text =
                             '<b>  </b>'),
              orientation = "V"
            ),
            xaxis = list(
              title = paste0(input_vals$g, " Alleles"),
              autotick = F,
              tickmode = "array",
              tickvals = tickvals,
              ticktext = ticktext
            ),
            yaxis = list(title = "Rep.\nnormalization")
            #,range = c(0, 0.5))
          )

        s <- plotly::subplot(
          plotly1,
          plotly2,
          nrows = 2,
          shareY = T,
          titleX = F,
          titleY = T,
          shareX = T,
          margin = 0.1,
          which_layout = 1
        ) %>% plotly::layout(
          legend = list(
            orientation = 'h',
            y = -0.5 - max(nchar(
              names(allele_thresh_state_absolute)
            )) / 100,
            x = 0
          ),
          xaxis = list(
            titlefont = list(size = 14),
            tickfont = list(size = 14)
          ),
          yaxis = list(
            titlefont = list(size = 16),
            tickfont = list(size = 16)
          ),
          yaxis2 = list(
            titlefont = list(size = 16),
            tickfont = list(size = 16)
          ),
          legend = list(font = list(size = 12)),
          hoverlabel = list(bgcolor = 'rgba(255,255,255,0.75)',
                            font = list(color = 'black'))
        )  %>%
          plotly::highlight(
            on = "plotly_click",
            opacityDim = 0.3,
            #selectize = TRUE,
            selected = attrs_selected(showlegend = F),
            persistent = T
          )

        s$x$source <- paste0("scatter", as.numeric(state_val$i))

        return(
          #onRender(
          list(
            plot = s %>% plotly::layout(margin = list(t = 50)) %>%
              plotly::event_register('plotly_click') %>% plotly::config(displayModeBar = T, #edits = list(shapePosition = TRUE),
                                                                        scrollZoom = FALSE),
            table = plot_data
          )
        )
        #     "function (el, x) {
        #
        #   // Get selectized search widget
        #   s = $('select')[0].selectize;
        #
        #   // Update placeholder
        #   s.settings.placeholder = 'Select a subject';
        #   s.updatePlaceholder();
        #
        #   // Hide group name (random string) that would appear above placeholder
        #   $('#htmlwidget_container').find('label').css('display', 'none');
        #
        #
        # }"
        #   )

      })

      output[[paste0("downloadData_", as.numeric(state_val$i))]] <-
        downloadHandler(
          filename = paste0(
            input_vals$g_group,
            "_state",
            as.numeric(state_val$i),
            ".tsv"
          ),
          content = function(file) {
            tab = isolate(filteredScores()$table)
            tab$v_allele <- as.character(tab$v_allele)
            columns_names <-
              strsplit(
                "subject,v_gene,v_allele,count,freq,freq2,project,J6_TAG,mut,v_allele_axis,zygousity_state,v_alleles",
                ","
              )[[1]]
            columns_names2 <-
              strsplit(
                "subject,v_gene,v_allele,count,group_freq,rep_freq,project,j6_hetro,mut,original_allele,zygousity_state,state_alleles_combination",
                ","
              )[[1]]
            tab <- tab[, columns_names]
            names(tab) <- columns_names2
            data.table::fwrite(x = as.data.table(tab),
                               file = file,
                               sep = "\t")
          }
        )

      output[[paste0("scatter", as.numeric(state_val$i))]] <-
        renderPlotly({
          filteredScores()$plot
        })

      #.update x/y reactive values in response to changes in shape anchors
      # observe({
      #   ed <- event_data("plotly_relayout", source=paste0("scatter", as.numeric(isolate(state_val$i))))
      #   shape_anchors <- ed[grepl("^shapes", names(ed))]
      #   if (length(shape_anchors) != 4) return()
      #   row_index <- unique(readr::parse_number(names(shape_anchors)) + 1)
      #   pts <- as.numeric(shape_anchors)
      #   a_thresh$a_thresh[row_index] <- pts[3]
      # })
      observeEvent(input$hetro, once = F, {
        output[[paste0("hover", as.numeric(state_val$i))]] <- renderPlotly({
          eventdat <-
            event_data('plotly_click', source = paste0("scatter", as.numeric(state_val$i))) # get event data from source main

          validate(need(
            !is.null(eventdat),
            "Please choose a J6 heterozygouse sample"
          ))

          eventdat <- eventdat[nrow(eventdat), ]
          subject_key <-
            eventdat[['key']] # Index of the data point being charted


          plot_bar <- plt_data_haplo$db

          validate(need(
            !is.null(plot_bar),
            paste0(
              "Sample ",
              subject_key,
              " doesn't have J6 sequences with the selected gene for this sample.\nPlease choose a different sample"
            )
          ))

          plot_bar <-
            dplyr::filter(plot_bar, !is.na(j_call), subject == subject_key)

          validate(need(
            nrow(plot_bar) > 0,
            paste0(
              "Sample ",
              subject_key,
              " doesn't have J6 sequences with the selected gene for this sample.\nPlease choose a different sample"
            )
          ))

          # draw plot according to the point number on hover
          #shiny:::flushReact()
          plot_bar <- plot_bar  %>% rowwise() %>% mutate(text = HTML(
            paste(
              '</br>Project: ',
              project,
              '</br>Subject: ',
              subject,
              '</br>Allele: ',
              v_allele,
              '</br># assignments: ',
              count,
              '</br>Group normalization freq.: ',
              freq,
              '</br>Rep. normalization freq.: ',
              freq2
            )
          ))
          print(plot_bar %>% as.data.frame())
          plotly::ggplotly(
            ggplot(
              plot_bar,
              aes(
                v_allele_axis,
                freq,
                fill = v_alleles_p,
                text = text
              )
            ) +
              geom_col(width = 0.2) + facet_grid(j_call ~ .) + theme_minimal(base_size = 12) +
              labs(
                y = "Relative allele frequency",
                x = paste0(subject_key, " Alleles"),
                fill = ""
              ) + scale_y_continuous(limits = c(0, 1)) +
              scale_fill_manual(values = colors()[unique(plot_bar$v_alleles_p)]) + theme(legend.position = "none"),
            tooltip = "text"
          ) %>%
            plotly::layout(showlegend = F) %>% plotly::config(displayModeBar = F)
          #) %>% layout(
          #  legend = list(orientation = 'h') ) %>% config(displayModeBar = F)
        })
      })
    }, states)
  })

  # observeEvent(req(any(a_thresh$a_thresh!=input_vals$allele_thresh)),{
  #
  #   lapply(1:length(unique(data_cut()$v_allele_axis)), function(i){
  #     lab <- unique(data_cut()$v_allele_axis)[i]
  #     if(a_thresh$a_thresh[i]!=input_vals$allele_thresh[lab]) updateNumericInput(session, paste0("allele",i), value = as.numeric(a_thresh$a_thresh[i]))
  #
  #   }
  #
  #   )
  #
  # })


  observeEvent(input$reset_input, once = F, {
    lapply(1:length(input_vals$allele_thresh), function(i) {
      lab <- names(input_vals$allele_thresh)[i]

      val <- as.numeric(input_vals$allele_thresh[lab])

      updateNumericInput(session, paste0("allele", i), value = val)

    })

  })


  observeEvent(input$reset_input2, once = F, {
    lapply(1:length(input_vals$allele_thresh_abs), function(i) {
      lab <- names(input_vals$allele_thresh_abs)[i]

      val <- as.numeric(input_vals$allele_thresh_abs[lab])

      updateNumericInput(session, paste0("allele_abs", i), value = val)

    })

  })

}

shinyApp(ui = ui, server = server)
