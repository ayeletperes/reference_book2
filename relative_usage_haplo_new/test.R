pacman::p_load(
  'dplyr',
  'tidyr',
  'htmltools',
  'bbplot',
  'ggplot2',
  'shiny',
  'shinyWidgets',
  'dendextend',
  'data.table',
  'alakazam',
  "unikn",
  'plotly',
  'ggdendro',
  "RColorBrewer",
  install = F
)
library(shinyWidgets)
library(unikn)
library(bbplot)
library(plotly)
library(ggplot2)
options(repos = BiocManager::repositories())

hline <- function(y = 0,
                  color = "red",
                  x0 = 0,
                  x1 = 1) {
  list(
    type = "line",
    x0 = x0,
    x1 = x1,
    xref = "paper",
    y0 = y,
    y1 = y,
    line = list(color = color, dash = "dot")
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

ui <- fluidPage(mainPanel(tabsetPanel(id = "tabs")))

library(htmlwidgets)

server <- function(input, output, session) {
  input_vals <- reactiveValues(
    allele_thresh = 0.05,
    g_group = "IGHV2-26G11",
    g = allele_db %>% filter(gene_group == "IGHV2-26G11") %>% pull(gene) %>% unique(),
    v_gene_cut = ifelse(grepl("G", "IGHV2-26G11"), "IGHV2-26G11", func_groups[as.character("IGHV2-26G11")])
  )
  observe({
    query <- parseQueryString(session$clientData$url_search)

    if (!is.null(query$allele_thresh)) {
      allele_thresh <-
        as.numeric(strsplit(query$allele_thresh, "\"")[[1]][2])
      g_group <- strsplit(query$g_group, "\"")[[1]][2]
      g <-
        allele_db %>% filter(gene_group == g_group) %>% pull(gene) %>% unique()
      v_gene_cut <-
        ifelse(grepl("G", g_group), g_group, func_groups[as.character(g_group)])

      input_vals$allele_thresh <- allele_thresh
      input_vals$state <- state
      input_vals$g_group <- g_group
      input_vals$g <- g
      input_vals$v_gene_cut <- v_gene_cut
    }
  })

  data_gene <- reactive({
    tmp_allele_db =
      allele_db %>% dplyr::filter(grepl(as.character(input_vals$g_group), new_allele)) %>%
      dplyr::group_by(new_allele) %>% dplyr::summarise(or_allele = paste0(or_allele, collapse = "/"))

    or_allele =
      setNames(gsub(chain, "", as.character(tmp_allele_db$or_allele)), as.character(gsub(
        paste0(input_vals$g_group, "[*]"),
        "",
        tmp_allele_db$new_allele
      )))

    tmp <- data %>%
      dplyr::filter(v_gene == input_vals$v_gene_cut,
                    !is.na(v_allele),
                    group_plot == 1) %>%
      ungroup()

    tmp <-
      tmp %>% dplyr::group_by(subject) %>% dplyr::arrange(desc(freq)) %>%
      dplyr::group_by(subject, v_gene) %>% dplyr::mutate(
        zygousity_state = as.numeric(sum(
          freq > input_vals$allele_thresh / 100, na.rm = T
        )),
        v_alleles = paste0(1:unique(zygousity_state), " - ", or_allele[v_allele[1:unique(zygousity_state)]], collapse = ";"),
        v_alleles_abc = paste0(sort(or_allele[v_allele[1:unique(zygousity_state)]]), collapse = ";"),
        v_allele_axis = or_allele[v_allele]
      ) %>% arrange(subject)
    tmp <-
      tmp %>% dplyr::group_by(subject, zygousity_state, group_plot) %>%
      dplyr::mutate(loc_state = loc <= zygousity_state) %>% dplyr::filter(loc_state) %>% ungroup()
    return(tmp)
  })

  states <- sort(as.numeric(isolate(unique(data_gene()$zygousity_state))))


  lapply(states, function(i) {
    insertTab(inputId = "tabs",
              tabPanel(
                paste0("State ", i),
                column(
                  1,
                  shinyWidgets::materialSwitch(
                    inputId = paste0("hetro", i),
                    label = "J6 heterozygous",
                    status = "primary",
                    right = TRUE
                  )
                ),
                column(7, plotlyOutput(paste0('scatter', i))),

                column(4, plotOutput(paste0('hover', i)))
              ))
  })

  updateTabsetPanel(session = session, inputId = "tabs", selected = paste0("State ", states[1]))

  Map(function(state){
    state_val <- reactiveValues(i = as.numeric(state))

    plt_data = reactive({
      tmp_plot <-
        dplyr::filter(data_gene(), zygousity_state == as.numeric(isolate(state_val$i)), is.na(j_call)) %>% dplyr::rowwise() %>% dplyr::mutate(
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
              '</br>Relative freq.: ',
              round(freq, 4),
              '</br>Relative Rep. freq.: ',
              round(freq2, 4),
              '</br>',
              J6_TAG
            )
          )
        ) %>% ungroup()
      return(tmp_plot)
    })

    plt_data_haplo <- reactive({
      tmp_allele_db =
        allele_db %>% dplyr::filter(grepl(as.character(input_vals$g_group), new_allele)) %>%
        dplyr::group_by(new_allele) %>% dplyr::summarise(or_allele = paste0(or_allele, collapse = "/"))

      or_allele =
        setNames(gsub(chain, "", as.character(tmp_allele_db$or_allele)), as.character(gsub(
          paste0(input_vals$g_group, "[*]"),
          "",
          tmp_allele_db$new_allele
        )))


      tmp_haplo <- data %>%
        dplyr::filter(v_gene == input_vals$v_gene_cut,
                      !is.na(v_allele),
                      group_plot == 2) %>%
        ungroup()

      tmp_haplo <- tmp_haplo %>% dplyr::group_by(subject) %>%
        dplyr::filter(v_allele %in% plt_data()$v_allele[plt_data()$subject ==
                                                          unique(subject)]) %>%
        dplyr::mutate(v_alleles_p = unique(plt_data()$v_alleles_p[plt_data()$subject ==
                                                                    unique(subject)]),
                      v_allele_axis = or_allele[v_allele])
      return(tmp_haplo)
    })

    colors <- reactive({
      n_alleles <- length(unique(plt_data()$v_alleles_p))
      if (n_alleles <= 4)
        setNames(cols[1:n_alleles], unique(plt_data()$v_alleles_p))
      else
        setNames(pal %>% usecol(n = n_alleles),
                 unique(plt_data()$v_alleles_p))
    })

    output[[paste0("scatter",as.numeric(isolate(state_val$i)))]] <- renderPlotly({
      tmp_allele_db =
        allele_db %>% dplyr::filter(grepl(as.character(input_vals$g_group), new_allele)) %>%
        dplyr::group_by(new_allele) %>% dplyr::summarise(or_allele = paste0(or_allele, collapse = "/"))

      or_allele =
        setNames(gsub(chain, "", as.character(tmp_allele_db$or_allele)), as.character(gsub(
          paste0(input_vals$g_group, "[*]"),
          "",
          tmp_allele_db$new_allele
        )))

      plot_data <- isolate(plt_data())

      loc2 <-
        setNames(1:length(unique(plot_data$v_allele_axis)), sort(unique(plot_data$v_allele_axis)))

      plot_data$loc2 <- loc2[plot_data$v_allele_axis]

      if (isolate(state_val$i) != 1 &
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
          dplyr::mutate(loc_plot = loc2 + loc_jitter[[as.character(unique(loc2))]][v_alleles_p], ) %>% ungroup()

        if (length(plot_data$loc_plot))
          plot_data <-
          plot_data %>% dplyr::mutate(jitter_offset = jitter(loc_plot))
      } else{
        plot_data <-
          plot_data %>% dplyr::arrange(loc2, v_alleles_p) %>% dplyr::group_by(loc2) %>%
          dplyr::mutate(loc_plot = loc2, ) %>% ungroup()
        if (length(plot_data$loc_plot))
          plot_data <-
            plot_data %>% dplyr::mutate(jitter_offset = jitter(loc_plot))
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




      if (input[[paste0("hetro", as.numeric(isolate(state_val$i)))]]){
        plot_data <- plot_data %>% filter(J6 == 1)
      }

      plot_data <- plot_data %>%
        highlight_key(., ~ subject)



      plotly1 <-
        plot_data %>%
        plot_ly(colors = colors()) %>%
        add_trace(
          type = "scatter",
          x = ~ (jitter_offset),
          y = ~ freq,
          text = ~ text,
          symbol = ~ project,
          mode = 'markers',
          marker = list(color = "grey", size = 12),
          showlegend = TRUE,
          opacity = 0.9,
          hoverinfo = 'none',
          legendgroup = ~ project
        ) %>%
        add_trace(
          type = "scatter",
          x = ~ (jitter_offset),
          y = ~ freq,
          text = ~ text,
          color = ~ v_alleles_p,
          colors = colors(),
          mode = 'markers',
          showlegend = FALSE,
          opacity = 0.8,
          hoverinfo = 'text',
          legendgroup = ~ v_alleles_p
        ) %>%
        add_trace(
          x = ~ as.numeric(loc_plot),
          y = ~ freq,
          color = ~ v_alleles_p,
          colors = colors(),
          type = "box",
          hoverinfo = "none",
          fillcolor = "transparent",
          name = ~ v_alleles_p,
          legendgroup = ~ v_alleles_p
        ) %>%
        layout(
          hovermode = 'closest',
          shapes = list(hline(input_vals$allele_thresh / 100)),
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
          yaxis = list(title = "Relative\nallele frequency",
                       range = c(0, 1.05))
        )


      plotly2 <-
        plot_data %>%
        plot_ly(colors = colors()) %>%
        add_trace(
          type = "scatter",
          x = ~ (jitter_offset),
          y = ~ freq2,
          text = ~ text,
          symbol = ~ project,
          mode = 'markers',
          marker = list(color = "grey", size = 12),
          showlegend = FALSE,
          opacity = 0.9,
          hoverinfo = 'none',
          legendgroup = ~ project
        ) %>%
        add_trace(
          type = "scatter",
          x = ~ (jitter_offset),
          y = ~ freq2,
          text = ~ text,
          color = ~ v_alleles_p,
          colors = colors(),
          mode = 'markers',
          showlegend = FALSE,
          opacity = 0.8,
          hoverinfo = 'text',
          legendgroup = ~ v_alleles_p
        ) %>%
        add_trace(
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
        layout(
          hovermode = 'closest',
          shapes = list(hline(input_vals$allele_thresh / 100)),
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
          yaxis = list(title = "Relative\nrepertoire frequency",
                       range = c(0, 1.05))
        )

      s <- subplot(
        plotly1,
        plotly2,
        nrows = 2,
        shareY = T,
        titleX = F,
        titleY = T,
        shareX = T,
        margin = 0.05,
        which_layout = 1
      ) %>% layout(
        xaxis = list(
          titlefont = list(size = 22),
          tickfont = list(size = 22)
        ),
        yaxis = list(
          titlefont = list(size = 22),
          tickfont = list(size = 22)
        ),
        yaxis2 = list(
          titlefont = list(size = 22),
          tickfont = list(size = 22)
        ),
        legend = list(font = list(size = 16))
      )  %>%
        plotly::highlight(on = "plotly_click",
                          opacityDim = 0.3,
                          persistent = TRUE)%>% event_register('plotly_click')
      s$x$source <- paste0("scatter", as.numeric(isolate(state_val$i)))
      s

    })


    output[[paste0("hover",as.numeric(isolate(state_val$i)))]] <- renderPlot({
      eventdat <-
        event_data('plotly_click', source = paste0("scatter", as.numeric(isolate(state_val$i)))) # get event data from source main

      validate(need(
        !is.null(eventdat),
        "Please choose a J6 heterozygouse sample"
      ))

      eventdat <- eventdat[nrow(eventdat),]
      subject_key <-
        eventdat[['key']] # Index of the data point being charted

      plot_bar <- isolate(plt_data_haplo())

      plot_bar <-
        dplyr::filter(plot_bar,!is.na(j_call), subject == subject_key)

      validate(need(
        nrow(plot_bar) > 0,
        paste0(
          "Sample ",
          subject_key,
          ' is not haplotypable.\nPlease choose a different sample'
        )
      ))

      # draw plot according to the point number on hover


      tmp_allele_db =
        allele_db %>% dplyr::filter(grepl(as.character(input_vals$g_group), new_allele)) %>%
        dplyr::group_by(new_allele) %>% dplyr::summarise(or_allele = paste0(or_allele, collapse = "/"))

      or_allele =
        setNames(gsub(chain, "", as.character(tmp_allele_db$or_allele)), as.character(gsub(
          paste0(input_vals$g_group, "[*]"),
          "",
          tmp_allele_db$new_allele
        )))

      ggplot(plot_bar, aes(v_allele_axis, freq, fill = v_alleles_p)) +
        geom_col(width = 0.2) + facet_grid(j_call ~ .) + bbplot::bbc_style() +
        labs(y = "Relative\nallele frequency", x = "Alleles", fill = "") +
        scale_fill_manual(values = colors()[unique(plot_bar$v_alleles_p)])
    })
  },states)
}

shinyApp(ui = ui, server = server)
