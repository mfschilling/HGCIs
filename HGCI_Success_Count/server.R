library(shiny)
source("HyperM_Generator.R")

shinyServer(function(input, output, session) {
  
  output$textlevel <- renderText({ 
    paste0("Level = ", input$level,"%")
  })
  
  output$textsample_size <- renderText({ 
    paste("Sample size =", input$sample_size)
  })
  
  output$textpopulation_size <- renderText({ 
    paste("Population size =", input$population_size)
  })
  
  output$LCOresults <- renderPrint({
    
    level <- input$level/100
    sample_size <- input$sample_size
    population_size <- input$population_size
    
    LCO.CI.M(sample_size,population_size,level)
  })  
})