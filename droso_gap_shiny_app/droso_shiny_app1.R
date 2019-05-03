#------------------------------
# a test shiny app
#
#
#   author: Jianhong Chen
#   date: April 5th, 2019
#-------------------------------

library("shiny")
library("ggplot2")
library("tidyverse")
library("DT")
library("tools")
library("plotly")

# --------------------
# load the data set 
file_path = "/Users/jianhongchen/Desktop/UCSC_AMS_PhD/Github/Drosophila_gap_gene/droso_gap_shiny_app/droso_gap_bin.csv"
droso_bin = read_csv(file_path, col_names = TRUE)

# build UI of th app 
ui = fluidPage(
  # add title
  titlePanel("Drosophila Gap Gene Network", windowTitle =  "Drosophila"),
  
  # sidebar layout for both input and output
  sidebarLayout(
    
    # inputs
    sidebarPanel(
      # select variable for y-axis
      selectInput(inputId = "y", label = "Y-axis", 
                  choices = c("AP", "gt", "kr", "kni", "hb" ),
                  selected = "gt"),
      #select variable for x-axis
      selectInput(inputId = "x", label = "X-axis",
                  choices = c("AP", "gt", "kr", "kni", "hb"),
                  selected = "AP"),
      # change the transparency of the data point
      sliderInput(inputId = "alpha", label = "Alpha", 
                  min = 0, max = 1, value = 0.8),
      # show data table or not
      checkboxInput(inputId = "show_data", label = "Display Data Table", 
                    value = FALSE),
      # show time info. and plot with different timeframe
      sliderInput(inputId = "time", label = "Plot w.r.t. Time",
                         min = 1, max = 8, value = c(2, 4)),
      # plot title selection
      textInput(inputId = "plot_title", label = "Plot title",
                placeholder =  "Enter text for the plot title"),
      br(), br(),
      # display logos 
      h5("Built with",
         img(src = "https://www.rstudio.com/wp-content/uploads/2014/04/shiny.png", height = "30px"),
         "by",
         img(src = "https://www.rstudio.com/wp-content/uploads/2014/07/RStudio-Logo-Blue-Gray.png", height = "30px"),
         ".")
    ),
    # output
    mainPanel(
      plotOutput(outputId =  "scatterplot"),
      conditionalPanel("input.show_data == true", h3("Data Table")), # java script symtext here
      DT::dataTableOutput(outputId = "table")
    )
  )
)
# create server to create the outputs
server = function(input, output, session){
  
  # render scatter plot
  output$scatterplot = renderPlot({
    df = filter(droso_bin, time >= input$time[1] & time <= input$time[2])
    ggplot(data = df, aes_string(x = input$x, y = input$y)) + 
      geom_point(alpha = input$alpha) +
      labs(title = isolate(toTitleCase(input$plot_title))) + 
      theme(text = element_text(size = 20))
  })
  
  # print data table 
  output$table = DT::renderDataTable(
    if(input$show_data){
      DT::datatable(data = droso_bin, option = list(pageLength = 10), rownames = FALSE)
    }
  )
  
}
  

shinyApp(ui, server)
  
  
  
  
  
  
