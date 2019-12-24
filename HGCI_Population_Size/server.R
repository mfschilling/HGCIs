library(shiny)
source("HyperN_Generator.R")

shinyServer(function(input, output, session) {
  
  output$textlevel <- renderText({ 
    paste0("Level = ", input$level,"%")
  })
  
  output$textsample_size <- renderText({ 
    paste("Sample size =", input$sample_size)
  })
  
  output$textnum_successes <- renderText({ 
    paste("Number of Successes =", input$num_successes)
  })
  
  output$Results <- renderPrint({
    
    level <- input$level/100
    sample_size <- input$sample_size
    num_successes <- input$num_successes
    
    CI.N(sample_size,num_successes,level)
    
  })  
})