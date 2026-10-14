# ============================================================================
# Shared helper: sourced by the analysis scripts, never run on its own.
# ============================================================================
#################
# Author: Ethan Tse
# Date: Sept 30th, 2024
#
# Utility Functions used throughout the project
#
#################

create_dir_if_needed <- function(directory_path) {
  # checks if directory exists. If not, recursively create.
  # no return value
  if (!dir.exists(directory_path)) {
    dir.create(directory_path, recursive = TRUE)
  }
}

ggsave_wrapper <- function(filename, plot, extension = "pdf", width = 7, height = 7) {
  # Wrapper for ggsave that saves to the correct directory
  # can be used to force saved images to be of a specific extension
  # extracting extension from the file: https://stackoverflow.com/questions/7779037/extract-file-extension-from-file-path
  require(tools)

  if (!is.null(extension)) {
    # remove extension in filename
    filepath <- tools::file_path_sans_ext(filename)
    filepath <- paste0(filepath, ".", extension)
  }

  if (extension == "png") {
    # need to use bg = "white" for png because black and white
    # are reversed in the background
    ggsave(
      filename = filepath,
      plot = plot,
      device = extension,
      width = width,
      height = height,
      create.dir = TRUE,
      bg = "white",
      dpi = 600
    )
  } else if (extension == "svg") {
    ggsave(
      filename = filepath,
      plot = plot + theme(text = element_text(family = "sans")),
      device = svglite::svglite,
      width = width,
      height = height,
      create.dir = TRUE,
      bg = "white",
      dpi = 600
    )
  } else {
    ggsave(
      filename = filepath,
      plot = plot,
      device = extension,
      width = width,
      height = height,
      create.dir = TRUE
    )
  }

  # save rds of the ggplot so we can edit later if needed
  saveRDS(plot, file = paste0(tools::file_path_sans_ext(filepath), ".rds"))
}
