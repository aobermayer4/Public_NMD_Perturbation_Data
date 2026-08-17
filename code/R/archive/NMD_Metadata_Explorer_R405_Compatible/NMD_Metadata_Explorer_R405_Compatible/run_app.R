# Run this file from the project folder, or set app_dir explicitly.
app_dir <- if (file.exists("app.R")) "." else "NMD_Metadata_Explorer"
shiny::runApp(appDir = app_dir, launch.browser = TRUE)
