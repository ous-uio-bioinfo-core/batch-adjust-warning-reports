library(shiny)


batches=1:20
groups=c("A", "B", "C")
baselinesignal=10
#groupsd=1

#initial settings
countsdefault = matrix(0, nrow=length(batches), ncol=length(groups), dimnames=list(batches, groups))
countsdefault[1,] = c(8,0,0)
countsdefault[2,] = c(2,2,0)
countsdefault[3,] = c(0,8,0)
groupeffectdefaults = rep(0, length(groups))
groupeffectdefaults[1:3] = c(0,0,0)
names(groupeffectdefaults) = groups
batcheffectdefaults = rep(1, length(batches))
#batcheffectdefaults[1:2] = c(0,2)


designheader = function(group)
{
  column( width=3, paste(group, "", sep="") , align="center")  
}

designrow = function(batch)
{
  fluidRow( 
    id = paste("batchrow", batch, sep=""),
    column( width=3, paste("Batch", batch, sep="") , align="right"),
    lapply( groups, groupcountselector, batch=batch)    
  )
}

groupcountselector = function(group, batch)
{
  thisid = paste("count", group, batch, sep="")
  column( 
    width=3,
    numericInput( thisid, NA, countsdefault[batch,group], min = 0, max = 100, step=1)
  )  
}

groupeffectselector = function(group)
{
  thisid = paste("groupeffect",group,  sep="")
  column( 
    width=4,
    numericInput( thisid, NA, groupeffectdefaults[group], min = -5, max = 5, step=1)
  )  
}

batcheffectselector = function(batch)
{
  thisid = paste("batcheffectvalue",batch, sep="")
  column( 
    width=4,
    numericInput( thisid, NA, batcheffectdefaults[batch], min = -5, max = 5, step=1)
  )  
}

#cat("adhoc.palette", adhoc.palette)

shinyUI(
  fluidPage( 
  	headerPanel("Batch Adjustment Simulator   - UNDER CONSTRUCTION -"),
    theme="bootstrap.css",
    includeCSS("www/styles.css"),
    tags$head(tags$script(src="helperfunctions.js")),
    #cat("getwd", getwd()),

    
    #titlePanel("title panel"),  
    sidebarLayout(
    	
      sidebarPanel(       
        width = 3,
        helpText("Set simulation parameters and plotting options.", " ", "[",a("Help", href="helptext.html"), "]"),
        
       
        
        
        # Experimental design
        wellPanel(          
          strong("Number of samples from each group (A,B,C) in the batches"),          
          #header
          fluidRow( column( width=3,""), lapply( groups, designheader)),          
          # batches          
          lapply( batches, designrow)          
        ), # end wellPanel
        
        # Group effects
        wellPanel(
          strong("Group effects, A,B and C"),
          # groupeffects          
          fluidRow( lapply( groups, groupeffectselector)),
          #helpText( "Added to the baseline signal (", baselinesignal,") for :", sep=""),
          br(),
          selectInput("groupeffectfraction", NA,
                      c("Added to all genes"=1,
                        "Added to 10%"=0.10,
                        "Added to 1%"=0.01,
                        "Added to 0.1%"=0.001,
                        "Added to 0.01%"=0.0001)),
          
          selectInput("groupeffecttype", NA,
                      c(
                        "As a per gene variable (mean=0, sd=abs(x))" = "sd",
                      	"As a constant (mean=x, sd=0)" = "mean"
                      	))
        

        ), # end wellPanel
        
        # Batch effects
        wellPanel(
          strong("Batch effects, batch 1, 2 and 3"),
          # groupeffects          
          fluidRow( lapply( batches, batcheffectselector)),
          
          selectInput("batcheffecttype", NA,
          						c(          							
          							"As a per gene variable (mean=0, sd=abs(x))" = "sd",
          							"As a constant (mean=x, sd=0)" = "mean"
          							))

        ), # end wellPanel
        
        
        # adjust methid
        wellPanel(
        selectInput("adjustmethod", "Adjust method:",
                    c("Mean-centring" = "Mean-centring",
                      "ComBat" = "ComBat",
                    	"ComBat, no covariates" = "ComBat no covariates",
                      "removeBatchEffect" = "removeBatchEffect"),
        						selected="removeBatchEffect")
        ),
        
        # simulation input
        wellPanel(
        fluidRow(
          #column( width=2, strong("Gene:"), align="right"),
          #column( width=3, numericInput( "indexgene", NA, 1, min=1, max=1000, step=1) ) ,
        
          column( width=2, strong("seed:"), align="right" ),
          column( width=3, numericInput( "rngseed", NA, 139, step=1) )
          )
        ),
        
        
        #Main plots
        wellPanel(          
          strong("Plot values from one example gene as:"),          
          checkboxInput("plottrue", "True values", TRUE),
          checkboxInput("plotbatchaffected", "Batch affected values", TRUE),          
          checkboxInput("plotbatchadjusted", "Batch adjusted values", TRUE),
          checkboxInput("plotlsmeans", "Batch as factor in 2-way ANOVA", TRUE)
        ), # end wellPanel
        
        # additional plot elements
        wellPanel(
          selectInput("plotCI", "Plot 95% Confidence Interval Boxes:",
                      c("Over values"="over",
                        "Right side of values"="right",
                        "Both"="both",
                        "Do not plot"="none")),
          checkboxInput("plotbatchbox", "Plot batch boxes", TRUE),
          checkboxInput("printbatchshift", "Print batch shift", FALSE),
          checkboxInput("zerocentre", "Zero-centre all values", FALSE),
          selectInput("plotdiagnosisplot", "Diagnostic plot using all simulated genes",
                      c("p-val histogram A vs. B" = "AB",
                        "p-val histogram A vs. C" = "AC",
                        "p-val histogram B vs. C" = "BC",
                      	"PCA,  SLOW!" = "pca",
                        "PVCA,  SLOW!" = "pvca",
                        "hclust" = "hclust",
                        "Do not plot" = "none"),
          						selected="AB")
        ),
        

        #Other
        wellPanel( 
          checkboxInput("blackbg", "Black background", FALSE)          
        ), # end wellPanel
        
        
        div( id="localhostonly",
          # parameters only visible when url is localhost or 127.0.0.1
          wellPanel(
            strong("Available as localhost"),
            fluidRow(            
              column( width=5, strong("Genes:"), align="right" ),
              column( width=3, numericInput( "ngenes", NA, value=1000) )
            ),
            fluidRow(            
              column( width=5, strong("Max p-value in plot:"), align="right" ),
              column( width=3,           selectInput("pvalueplotxlim", NA,
                                                     c(1,
                                                       0.25,
                                                       0.10,
                                                       0.05,
                                                       0.01)) )
            ),
            fluidRow(            
              column( width=5, strong("Number of batches:"), align="right" ),
              column( width=3,           tags$select( id="numberofbatches",
                                                      onchange="setbatchrows()",
                                                      tags$option(3),
                                                      tags$option(5),
                                                      tags$option(10),
                                                      tags$option(15),
                                                      tags$option(20)
                                                      )
              )
            )
            
            
          )

        ),
        
        
        
        

        

        br(),br(),br(),br(),br(),br(),br(),br(),br(),
        
        tags$script("initparams()")
      ), # end sidebarPanel
      
      mainPanel(
      	
      	plotOutput("batchadjustplot",  height = "1000px"),
      	
      	#a("test1", href="#", onclick="initparams();return false;"),      	      
      	"[", a("Help text, with examples", href="helptext.html"), "]",
      	
      	div( htmlOutput("optionslink"), class="pull-right"),
      	
      	br(),br(),
        
    
          
      	br(),br()
			)
      
    ) # end sidebarLayout
  ) # end fluidPage
) # end shinyUI