library(shiny)

shinyUI(fluidPage(
  
  tags$title("Confidence Interval For Hypergeometric Population Size"),
  titlePanel("Confidence Interval For Hypergeometric Population Size"),
  
  div(a(href="https://mfschilling.shinyapps.io/hgci-success-count/",target="_blank", "Click Here to Estimate Success Count Instead"),
      style = "font-size: 18pt;color:black"),br(),br(),
  sidebarPanel(
    
    sliderInput("level", 
                label = h5("Confidence Level %:"),
                min = 80, max = 99, value = 95, step=1),br(),
    
    sliderInput("sample_size", 
                label = h5("Sample Size:"),
                min = 1, max = 100, value = 25),br(),
    
    sliderInput("num_successes", 
                label = h5("Number of Successes:"),
                min = 1, max = 100, value = 50),br(),
    
    
    div(submitButton("Submit"),align="right"), br(), br(), br(), br(), br(), 
    
    
    div("CI Generator Shiny app",align="right", style = "font-size: 8pt"), 
    
    div("maintained by", 
        a(href="http://www.csun.edu/~hcmth031/",target="_blank", 
          "Mark Schilling"),align="right", style = "font-size: 8pt"),
    
    div("Source Code:",
        a(href="https://github.com/mfschilling/HGCIs",
          target="_blank","GitHub Repo"),align="right", style = "font-size: 8pt")
  ),
  
  mainPanel(
    p("Details on this method of generating the confidence interval for 
      the hypergeometric population size can be found in the following journal 
      article:"),
    
    tags$blockquote("[The journal details will be inserted here after publication]",
                    em("example journal"),
                    a(href="[TODO: Insert journal link when published]", 
                      target="_blank", "(Online access)")),
    
    p("The ", a(href="https://github.com/mfschilling/HGCIs", target="_blank", "Github Repo"), 
      "contains the R code for the CI algorithm and Shiny app"),
    
    HTML("<hr style='height: 2px; color: #de7008; background-color: #df7109; border: none;'>"),
    
    textOutput("textlevel"),
    textOutput("textsample_size"),
    textOutput("textnum_successes"),
    br(),
    verbatimTextOutput("Results")
    )
  )
  )