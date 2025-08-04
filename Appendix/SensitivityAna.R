#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

##############################################################################
## Shinyapp – QALY‑weighted RMST               ##
##############################################################################
library(readxl)
library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)

df <- read_xlsx("df_opt.xlsx")
df <- df %>%
  mutate(
    across(starts_with("type"), ~ ifelse(.x < 1 | .x > 5 | is.na(.x), NA, .x))
  )

# ── UI -----------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("QALY Sensitivity Analysis: Different Weighting Principles"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("w2", "QALY Utility for Revascularization:",         min = 0.87, max = 1,   value = 0.90, step = 0.01),
      sliderInput("w3", "QALY Utility for Myocardio Infarction:",              min = 0.61, max = 0.9, value = 0.88, step = 0.01),
      sliderInput("w4", "QALY Utility for Hospitalized Angina:",    min = 0.69, max = 0.9, value = 0.69, step = 0.01),
      sliderInput("w5", "QALY Utility for Stroke:",          min = 0.15, max = 0.80,value = 0.15, step = 0.01)
    ),
    mainPanel(
      plotOutput("principle_plot", height = "550px")
    )
  )
)

# ── Server -------------------------------------------------------------------
server <- function(input, output) {
  
  ## Reactive block that does *all* the heavy lifting ------------------------
  principle_metrics <- reactive({
    ## 1.  Current utility vector from sliders
    w <- c(1, input$w2, input$w3, input$w4, input$w5, 0)
    
    # Print to console
    cat("Weights:", w, "\n")
    
    # Check for problems early
    if (any(is.na(w))) {
      stop("Some weights are NA")
    }
    
    #### Principle1 – Markov ----
    df_P1 <- df %>%
      mutate(
        p1 = time1 * w[1],
        p2 = ifelse(is.na(time2), 0, w[type1] * (time2 - time1)),
        p3 = ifelse(is.na(time3), 0, w[type2] * (time3 - time2)),
        p4 = ifelse(is.na(time4), 0, w[type3] * (time4 - time3)),
        p5 = ifelse(is.na(time5), 0, w[type4] * (time5 - time4)),
        p6 = ifelse(type5 == 0 | is.na(time5), 0, w[type5] * (36 - time5)),
        w_st = p1 + p2 + p3 + p4 + p5 + p6
      )
    
    #### Principle2 – Worst‑state update ----
    df_P2 <- df %>%
      rowwise() %>%
      mutate(
        p1 = time1 * w[1],
        worst1 = max(1, type1, na.rm = TRUE),
        p2 = ifelse(is.na(time2), 0, w[worst1] * (time2 - time1)),
        worst2 = max(c(1, type1, type2), na.rm = TRUE),
        p3 = ifelse(is.na(time3), 0, w[worst2] * (time3 - time2)),
        worst3 = max(c(1, type1, type2, type3), na.rm = TRUE),
        p4 = ifelse(is.na(time4), 0, w[worst3] * (time4 - time3)),
        worst4 = max(c(1, type1, type2, type3, type4), na.rm = TRUE),
        p5 = ifelse(is.na(time5), 0, w[worst4] * (time5 - time4)),
        worst5 = max(c(1, type1, type2, type3, type4, type5), na.rm = TRUE),
        p6 = ifelse(is.na(time5) | type5 == 0, 0, w[worst5] * (36 - time5)),
        w_st = p1 + p2 + p3 + p4 + p5 + p6
      ) %>%
      ungroup()
    
    #### Principle3 – Cumulative product ----
    df_P3 <- df %>%
      rowwise() %>%
      mutate(
        p1 = time1 * w[1],
        q2 = w[1] * w[type1],
        p2 = ifelse(is.na(time2), 0, q2 * (time2 - time1)),
        q3 = q2 * w[type2],
        p3 = ifelse(is.na(time3), 0, q3 * (time3 - time2)),
        q4 = q3 * w[type3],
        p4 = ifelse(is.na(time4), 0, q4 * (time4 - time3)),
        q5 = q4 * w[type4],
        p5 = ifelse(is.na(time5), 0, q5 * (time5 - time4)),
        q6 = q5 * w[type5],
        p6 = ifelse(is.na(time5) | type5 == 0, 0, q6 * (36 - time5)),
        w_st = p1 + p2 + p3 + p4 + p5 + p6
      ) %>%
      ungroup()
    
    #### Principle 4 – Naive RMST ----
    df_P4 <- df %>%
      mutate(
        w_st = time1 * w[1],
      )
    

    ## Summary metrics
    pr_labels <- c("Markov", "Worst‑state Update", "Cumulative Product", "Naive RMST")
    
    metrics <- bind_rows(
      tibble(principle = pr_labels[1],
             mean_evo  = mean(df_P1$w_st[df_P1$arm == 1], na.rm = TRUE),
             mean_pla  = mean(df_P1$w_st[df_P1$arm == 0], na.rm = TRUE)),
      tibble(principle = pr_labels[2],
             mean_evo  = mean(df_P2$w_st[df_P2$arm == 1], na.rm = TRUE),
             mean_pla  = mean(df_P2$w_st[df_P2$arm == 0], na.rm = TRUE)),
      tibble(principle = pr_labels[3],
             mean_evo  = mean(df_P3$w_st[df_P3$arm == 1], na.rm = TRUE),
             mean_pla  = mean(df_P3$w_st[df_P3$arm == 0], na.rm = TRUE)),
      tibble(principle = pr_labels[4],
             mean_evo  = mean(df_P4$w_st[df_P4$arm == 1], na.rm = TRUE),
             mean_pla  = mean(df_P4$w_st[df_P4$arm == 0], na.rm = TRUE))
    ) %>%
      mutate(
        diff_months = mean_evo - mean_pla,
        diff_days = diff_months * 30
      )
    
    metrics
  })
  
  output$principle_plot <- renderPlot({
    met <- principle_metrics()
    
    df_long <- bind_rows(
      met %>% transmute(principle, Treatment = "Evolocumab", RMST = mean_evo),
      met %>% transmute(principle, Treatment = "Placebo", RMST = mean_pla)
    ) %>%
      mutate(
        principle = factor(principle, levels = c("Markov", "Worst‑state Update", "Cumulative Product", "Naive RMST")),
        Treatment = factor(Treatment, levels = c("Evolocumab", "Placebo"))
      )
    
    ggplot(df_long, aes(x = principle, y = RMST, fill = Treatment)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
      geom_text(
        data = met,
        aes(x = principle,
            y = pmax(mean_evo, mean_pla) + 0.5,
            label = paste0("Δ = ", round(diff_days, 2), " days\n(",
                           round(diff_days / 30, 2), " months)")),
        inherit.aes = FALSE,
        color = "black",
        size = 4
      ) +
      scale_fill_manual(values = c("Evolocumab" = "#1f78b4", "Placebo" = "#e31a1c")) +
      coord_flip(ylim = c(30, 36)) +
      labs(
        title = "Q-RMST under Different Principles",
        x = "QALY Principle",
        y = "RMST (months)",
        fill = "Treatment"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 15, hjust = 1, size = 12),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11)
      )
  })
}


# ── Run the app --------------------------------------------------------------
shinyApp(ui = ui, server = server)
