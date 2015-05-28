library(shiny)

defaultdesign="a,b,c,d,e
10,10,0,0,0
0,10,10,0,0
0,0,10,10,0
0,0,0,10,10
10,0,0,0,10
"

shinyUI(
  fluidPage( 
  	headerPanel("Fdist app   - UNDER CONSTRUCTION -"),
    theme="bootstrap.css",
    includeCSS("www/styles.css"),
    tags$head(tags$script(src="helperfunctions.js")),
    #cat("getwd", getwd()),

    
    #titlePanel("title panel"),  
    sidebarLayout(
    	
      sidebarPanel(       
        width = 3,
        helpText(" helptext"),
                       
        # Experimental design
        wellPanel(
          h4("Trial design"),
          ("Group names and sample size (columns) per batch (rows)"),
          br(),
          tags$textarea(id="design", rows=10, cols=20, defaultdesign)
        ), # end wellPanel  
        
        # simulation input
        wellPanel(
          h4("Random samples"),
        	fluidRow(
        		column( width=4, "Samples" ),
        		column( width=8, numericInput( "samples", NA, 1000, min=1, max=10000, step=1) )
        	),
          fluidRow(
        		column( width=4, "Seed" ),
        		column( width=8, numericInput( "rngseed", NA, 139, min=0, step=1) )
          ),
          br()
        )

      ), # end sidebarPanel
      
      mainPanel(
      	div( h2("Summary"), htmlOutput("textoutput")),
        div( h2("Q-Q plot of F-statistic"), plotOutput("fdistplot",  height = "600px") ),
      	br(),br()
			)
      
    ) # end sidebarLayout
  ) # end fluidPage
) # end shinyUI