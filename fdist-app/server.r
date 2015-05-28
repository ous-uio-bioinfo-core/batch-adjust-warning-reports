library(shiny)

#source("F-distribution_mod.R")
source("F-test.R")


shinyServer(
	function(input, output)
	{			
		#x = fdistcalc(input)
		
		x <- reactive({
			fdistcalc(input)
		})
		
		output$textoutput <- renderText({ 
			fdisttext(x())
		})
		
		output$fdistplot <- renderPlot({
	    fdistplot(x())
	  })
		
	}
)