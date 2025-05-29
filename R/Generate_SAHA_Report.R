Generate_SAHA_Report <- function(SAHA, auto = NULL, output_file) {
   
   # Determine the full output path
   # Check if output_file is an absolute path (starts with / or C:/, D:/ etc.)
   # Note: The 'file.path()' function itself handles absolute paths correctly
   # when combined with a relative path. If 'output_file' is absolute,
   # 'file.path(getwd(), output_file)' will effectively return 'output_file'.
   # However, for clarity and explicit control, the check below is more robust
   # and prevents potential double slashes if getwd() was used naively.
   
   # This regex checks for paths starting with '/' (Unix-like) or 'X:/' (Windows)
   is_absolute_path <- grepl("^(/|[A-Za-z]:/)", output_file)
   
   if (is_absolute_path) {
      # If it's already an absolute path, use it directly
      final_output_path <- normalizePath(output_file, mustWork = FALSE)
   } else {
      # If it's a relative path, combine it with the current working directory
      final_output_path <- normalizePath(file.path(getwd(), output_file), mustWork = FALSE)
   }
   
   # Render the .Rmd file
   rmarkdown::render(input = system.file("extdata", "report_html.Rmd", package = "SAHAdata"),
                     output_file = final_output_path)
}
