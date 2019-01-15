observeEvent(input$outlier_about, {
  showModal(modalDialog(
    title = "Detecting and Removing Outliers",
    p("This button will calculate outliers by using Cook's distance. Once it is calculated and the plot appears, the option to remove outliers from further analyses will be available. Clicking the 'Remove Outliers' button will refresh all modules and remove outliers from the dataset. From this point onward, any calculations and plots will reflect the remaining data."), 
    p("If your uploaded data is fairly large and analyses are taking a long time to complete, removing outliers in this step should substantially reduce calculation times. You may also download the dataset before and after outliers are removed by using the 'Download Merged Data' button in the 'Merging Files' dropdown."),
    easyClose = TRUE,
    footer = NULL
    ))
})
