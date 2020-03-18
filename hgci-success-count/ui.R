library(shiny)
options(shiny.sanitize.errors = FALSE)

shinyUI(fluidPage(
  tags$title("LCO Confidence Interval For Hygergeometric Success Count"),
  titlePanel("LCO Confidence Interval For Hygergeometric Success Count"),
  
  div(a(href="https://mfschilling.shinyapps.io/hgci-population-size/",target="_blank", "Click Here to Estimate Population Size Instead"),
      style = "font-size: 18pt;color:black"),br(),br(),
  
  sidebarPanel(
    
    sliderInput("level", 
                label = h5("Confidence Level %:"),
                min = 80, max = 99, value = 95, step=1),br(),
    
    sliderInput("sample_size", 
                label = h5("Sample Size:"),
                min = 1, max = 100, value = 25),br(),
    
    sliderInput("population_size", 
                label = h5("Population Size:"),
                min = 1, max = 100, value = 50),br(),
    
    
    div(submitButton("Submit"),align="right"), br(), br(), br(), br(), br(), 
    
    
    div("LCO-CI Generator Shiny app",align="right", style = "font-size: 8pt"), 
    
    div("maintained by", 
        a(href="http://www.csun.edu/~hcmth031/",target="_blank", 
          "Mark Schilling"),align="right", style = "font-size: 8pt"),
    
    div("Source Code:",
        a(href="https://github.com/mfschilling/HGCIs",
          target="_blank","GitHub Repo"),align="right", style = "font-size: 8pt"),
    
    div("For LCO", tags$i("binomial"), "confidence intervals click",
        a(href="http://shiny.calpoly.sh/LCO_CI_Generator/",
          target="_blank","here"),align="left", style = "font-size: 12pt"),
  ),
  
  mainPanel(
    p("Details on the Length/Coverage Optimal (LCO) confidence interval for 
      the hypergeometric success count can be found in the following journal 
      article:"),
    
    tags$blockquote(a(href="https://www.tandfonline.com/doi/full/10.1080/03610926.2020.1737879", 
                      target="_blank", 'Schilling, M.F. and Stanley, A., “A new approach to precise interval estimation for the parameters of the hypergeometric distribution",'),
                    a(href="https://www.tandfonline.com/doi/full/10.1080/03610926.2020.1737879", 
                      target="_blank", em("Communications in Statistics-Theory and Methods,")),
                    a(href="https://www.tandfonline.com/doi/full/10.1080/03610926.2020.1737879", 
                      target="_blank", "March 2020")),
    
    p("The ", a(href="https://github.com/mfschilling/HGCIs", target="_blank", "Github Repo"), 
      "contains the R code for the CI algorithm and Shiny app"),
    
    HTML("<hr style='height: 2px; color: #de7008; background-color: #df7109; border: none;'>"),
    
    textOutput("textlevel"),
    textOutput("textsample_size"),
    textOutput("textpopulation_size"),
    br(),
    verbatimTextOutput("LCOresults")
    )
  )
  )