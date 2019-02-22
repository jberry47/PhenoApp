observeEvent(input$outlier_about, {
  showModal(modalDialog(
    title = "Detecting and Removing Outliers",
    p("This button will calculate outliers by using Cook's distance. Once it is calculated and the plot appears, the option to remove outliers from further analyses will be available. Clicking the 'Remove Outliers' button will refresh all modules and remove outliers from the dataset. From this point onward, any calculations and plots will reflect the remaining data.", strong("Removing outliers will cause ALL figures to then reflect the data without outliers. To undo this, refresh the application and re-load in the data.")), 
    p("If your uploaded data is fairly large and analyses are taking a long time to complete, removing outliers in this step should substantially reduce calculation times. You may also download the dataset before and after outliers are removed by using the 'Download Merged Data' button in the 'Merging Files' dropdown."),
    easyClose = TRUE,
    footer = NULL
    ))
})

observeEvent(input$shapes_anova_about, {
  showModal(modalDialog(
    title = "ANOVA Plot",
    p("This ANOVA attempts to estimate the percent variance explained in all 'shapes' elements due to the variables in the design file (disregarding Barcodes) on a specifiable DAP. It uses a fully random effect model and variances are calculated by Type3 sum of squares."),
    easyClose = TRUE,
    footer = NULL
  ))
})

observeEvent(input$anova_ts_about, {
  showModal(modalDialog(
    title = "Temporal ANOVA Plot",
    p("This ANOVA attempts to estimate the percent variance explained over time for a specifiable 'shape' element due to the variables in the design file. It uses a fully random effect model and variances are calculated by Type3 sum of squares."),
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

observeEvent(input$boxplot_about, {
  showModal(modalDialog(
    title = "Shapes Boxplots",
    p("The boxplots show fluctuations over time in the chosen 'shape' variable under 'Color By' depending on the design elements selected for 'Group By' and 'Facet By'. The analysis will only take the data from a particular day that can be specified under 'Which Day'. A histogram overlays the data in addition to the boxplot which is intended to show more detailed information about the spread of the data."),
    easyClose = TRUE,
    footer = NULL
  ))
})

observeEvent(input$vis_caps_about, {
  showModal(modalDialog(
    title = "Vis CAPS",
    p("This figure shows the constrained analysis of principal coordinates (CAP) which is intended to show variation amongst one design variable while also accounting for the remaining design variables. If there are only two varaiables in the design file, the CAPs calculation will be reduced to the formula that is used to calculate PCAs. "),
    p("Before the analysis begins, all rows and columns which contain only values of 0 are removed from the Color dataset."),
    p("To prevent the application from being overloaded by calculations, the 'Go' button will be disabled when a calculation begins. It can be re-activated by changing the variable in 'Main effect', 'Distance Type', or 'Which Day'. The option to download the plot will appear once the calculation is complete."),
    p("This analysis will likely take several minutes. It is possible to ",strong("reduce the calculation time by removing outliers from the data input"),". To do this, refer to the 'Outlier Detection and Removal' tab."),
    easyClose = T,
    footer = NULL
  ))
})

observeEvent(input$vis_joyplot_about, {
  showModal(modalDialog(
    title = "Vis Joyplot",
    p("This plot will display HUE channel data on a specified day for the first two elements in the design file (disregarding Barcodes). The HUE channel is shown as degrees on the x-axis and reflects the overall visually perceived color of the plants in each group."),
    easyClose = T,
    footer = NULL
  ))
})

observeEvent(input$nir_heatmap_withfacet_about, {
  showModal(modalDialog(
    title = "Nir Faceted Heatmap",
    p("The color-scale of these heatmaps represents the average grayscale density of entire plants. As the color becomes more gray, the average intensity increases. This heatmap specifically can be used to visualize average intensity fluctuations over time relative to the first two variables in the design file (disregarding Barcodes)."),
    easyClose = T,
    footer = NULL
  ))
})

observeEvent(input$nir_heatmap_nofacet_about, {
  showModal(modalDialog(
    title = "Nir Heatmap",
    p("The color-scale of these heatmaps represents the average grayscale density of entire plants. As the color becomes more gray, the average intensity increases. This heatmap specifically can be used to visualize average intensity fluctuations over time for the first variable in the design file (disregarding Barcodes)."),
    easyClose = T,
    footer = NULL
  ))
})

observeEvent(input$nir_caps_about, {
  showModal(modalDialog(
    title = "Nir CAPS",
    p("This figure shows the constrained analysis of principal coordinates (CAP) which is intended to show variation amongst one design variable while also accounting for the remaining design variables. If there are only two varaiables in the design file, the CAPs calculation will be reduced to the formula that is used to calculate PCAs."),
    p("Before the analysis begins, all rows and columns which contain only values of 0 are removed from the NIR dataset."),
    p("To prevent the application from being overloaded by calculations, the 'Go' button will be disabled when a calculation begins. It can be re-activated by changing the variable in 'Main effect', 'Distance Type', or 'Which Day'. The option to download the plot will appear once the calculation is complete."),
    p("This analysis will likely take several minutes. It is possible to ",strong("reduce the calculation time by removing outliers from the data input"),". To do this, refer to the 'Outlier Detection and Removal' tab."),
    easyClose = T,
    footer = NULL
  ))
})

 observeEvent(input$iqv_about, {
   showModal(modalDialog(
     title = "IQV Plot",
     p("The y-axis of this plot describes the whole color profile deviance of all images from a color-checker reference. This data is shown over a time period of 24 hours, displaying any fluctuations in image quality throughout the average day. More details about this scale can be", a("found here",target="_blank",href="https://peerj.com/articles/5727/"), "."),
     easyClose = T,
     footer = NULL
   ))
 })

observeEvent(input$water_about, {
  showModal(modalDialog(
    title = "Water Weight Plot",
    p("This plot can display fluctuations in weight(g) overtime either before watering or after watering. It can also show how much water (mL) was added each day for each pot to reach the target weight. Data will be displayed relative to two selected design variables in 'Facet By' and 'Color By'."),
    easyClose = T,
    footer = NULL
  ))
})

observeEvent(input$oof_about, {
  showModal(modalDialog(
    title = "Out of Frame Plot",
    p("This plot describes the approximate DAP in which plants begin to grow out of frame in the images. The plot does not take into account empty pots."),
    p("If the design file contains more than three variables (aside from Barcodes), the complexity will prevent the figure from loading. If this occurs, please subset the design file so it includes Barcodes and one to three design variables."),
    p("The plot will load in once a design variable is selected for all three input selections."),
    easyClose = T,
    footer = NULL
  ))
})

observeEvent(input$er_about, {
  showModal(modalDialog(
    title = "Emergence Rate Plot",
    p("This plot shows the approximate DAP in which plants begin to emerge from their pots and come into frame in the images. The plot does not take into account empty pots."),
    p("If the design file contains more than three variables (aside from Barcodes), the complexity will prevent the figure from loading. If this occurs, please subset the design file so it includes Barcodes and one to three design variables."),
    p("The plot will load in once a design variable is selected for all three input selections."),
    easyClose = T,
    footer = NULL
  ))
})