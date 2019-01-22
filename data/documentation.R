observeEvent(input$outlier_about, {
  showModal(modalDialog(
    title = "Detecting and Removing Outliers",
    p("This button will calculate outliers by using Cook's distance. Once it is calculated and the plot appears, the option to remove outliers from further analyses will be available. Clicking the 'Remove Outliers' button will refresh all modules and remove outliers from the dataset. From this point onward, any calculations and plots will reflect the remaining data."), 
    p("If your uploaded data is fairly large and analyses are taking a long time to complete, removing outliers in this step should substantially reduce calculation times. You may also download the dataset before and after outliers are removed by using the 'Download Merged Data' button in the 'Merging Files' dropdown."),
    easyClose = TRUE,
    footer = NULL
    ))
})

observeEvent(input$trends_about, {
  showModal(modalDialog(
    title = "Shapes Trends Plots",
    p("This plot displays the changes over time in a selected 'shape' variable depending on the design elements selected under 'Color By' and 'Facet By'. If the imported data has less than 1000 observations, it will be fitted by local polynomial regression. If it has more than 1000 observations, it will be fitted by a generalized additive model."), 
    easyClose = TRUE,
    footer = NULL
  ))
})

observeEvent(input$heatmap_about, {
  showModal(modalDialog(
    title = "Shapes Heatmap Plots",
    p("The heatmap shows fluctuations over time in the chosen 'shape' variable under 'Color By' depending on the design elements selected for 'Group By' and 'Facet By'. The lighter the color, the higher the value for the chosen 'shape' variable."),
    easyClose = TRUE,
    footer = NULL
  ))
})

# observeEvent(input$boxplot_about, {
#   showModal(modalDialog(
#     title = "Shapes Boxplot Plots",
#     p("test"), 
#     easyClose = TRUE,
#     footer = NULL
#   ))
# })