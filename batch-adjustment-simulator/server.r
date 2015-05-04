library(shiny)

source("batchadjustplot.r")

shinyServer(
	function(input, output)
	{			
		output$batchadjustplot <- renderPlot({
			batchadjustplot(input)			
		})		
		
		output$optionslink <- renderText({ 
			makeoptionslink(input)
		})
		
	}
)