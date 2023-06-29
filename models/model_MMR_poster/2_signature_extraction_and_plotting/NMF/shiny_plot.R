# To create an interactive Shiny app that allows users to change the "maxK" value, we need to modify your existing script to be compatible with Shiny's structure. This involves dividing the code into two parts: the User Interface (UI) and the Server.
# 
# The UI part is where you define the layout and appearance of the app. The server part is where the app's operations (like calculations, plot generation, etc.) are performed.
# 
# In the UI, we can use numericInput() to create an input field for "maxK". In the server function, we use renderPlot() to generate the plot, which will be updated every time the "maxK" input changes.
# 
# Here's how you can modify your existing script to turn it into a Shiny app:

# In this code, the numericInput() function creates a numeric input field with an id of "maxK". This id is used in the server part of the app to access the value entered by the user. The plotOutput() function in the UI corresponds to the renderPlot() function in the server. They both use the id "plot", which links the generated plot to the correct location in the UI.
# 
# When the user changes the value in the "maxK" input field, the renderPlot() function is automatically re-executed, which updates the plot with the new "maxK" value.
# 
# Please note that you need to replace the comment # Rest of your code goes here, replacing the hardcoded maxK with the input$maxK with your actual code. The print(res) line at the end of renderPlot() is used to display the plot in the Shiny app.


# Load necessary libraries
library(shiny)
library(tidyverse)
library(RcppML)
library(ggrepel)

# Define the UI
ui <- fluidPage(
  titlePanel("Interactive plot with maxK"),
  sidebarLayout(
    sidebarPanel(
      numericInput(inputId = "maxK", label = "Enter maxK:", value = 5)
    ),
    mainPanel(
      plotOutput(outputId = "plot")
    )
  )
)

# Define the server
server <- function(input, output) {
  output$plot <- renderPlot({
    maxK <- input$maxK
    topn = 5
    pos = position_jitter(w = 0.25, h = 0, seed = 1)
    # Rest of your code goes here, replacing the hardcoded maxK with the input$maxK
    # ...
    # Finally, print the plot
    print(res)
  })
}

# Run the app
shinyApp(ui = ui, server = server)
