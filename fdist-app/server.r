library(shiny)

source("F-distribution_mod.R")


shinyServer(
	function(input, output)
	{			
		#x = fdistcalc(input)
		
		x <- reactive({
			fdistcalc(input)
		})
		
 		output$fdistplot <- renderPlot({
 			#fdistplot(input)
 			fdistplot(x())
 		})		
		
		output$textoutput <- renderText({ 
			#x$a
			#fdisttext(x)
			fdisttext(x())
		})
		
	}
)