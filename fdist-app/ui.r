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
          strong("Number of samples from each group (column) in the batches (rows)"),          
          tags$textarea(id="design", rows=20, cols=20, defaultdesign)
        ), # end wellPanel  
        
        # simulation input
        wellPanel(
        	fluidRow(
        		#column( width=2, strong("Gene:"), align="right"),
        		#column( width=3, numericInput( "indexgene", NA, 1, min=1, max=1000, step=1) ) ,
        		
        		column( width=2, strong("seed:"), align="right" ),
        		column( width=3, numericInput( "rngseed", NA, 139, step=1) )
        	)
        )

      ), # end sidebarPanel
      
      mainPanel(
      	
      	plotOutput("fdistplot",  height = "1000px"),

      	div( htmlOutput("textoutput"), class="pull-right"),
      	
      	br(),br(),
        
    
          
      	br(),br()
			)
      
    ) # end sidebarLayout
  ) # end fluidPage
) # end shinyUI