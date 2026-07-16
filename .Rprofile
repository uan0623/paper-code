local({
  align_trailing_pipes <- function(text) {
    lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
    pipe_re <- "^(.*?)(\\s*)(%>%)(\\s*(#.*)?)$"

    align_block <- function(start, end) {
      pipe_lines <- start:end
      pipe_cols <- vapply(regmatches(lines[pipe_lines], regexec(pipe_re, lines[pipe_lines])), function(match) {
        nchar(match[2], type = "width") + nchar(match[3], type = "width") + 1L
      }, integer(1))
      target_col <- max(pipe_cols)

      for (idx in pipe_lines) {
        match <- regmatches(lines[idx], regexec(pipe_re, lines[idx]))[[1]]
        prefix <- trimws(match[2], which = "right")
        suffix <- match[5]
        pad <- strrep(" ", target_col - nchar(prefix, type = "width") - 1L)
        lines[idx] <<- paste0(prefix, pad, "%>%", suffix)
      }
    }

    block_start <- NULL
    for (idx in seq_along(lines)) {
      has_pipe <- grepl(pipe_re, lines[idx])
      if (has_pipe && is.null(block_start)) {
        block_start <- idx
      }
      if ((!has_pipe || idx == length(lines)) && !is.null(block_start)) {
        block_end <- if (has_pipe && idx == length(lines)) idx else idx - 1L
        if (block_end > block_start) {
          align_block(block_start, block_end)
        }
        block_start <- NULL
      }
    }

    paste(lines, collapse = "\n")
  }

  patch_languageserver_formatter <- function(...) {
    ns <- asNamespace("languageserver")
    original_style_text <- get("style_text", envir = ns)

    patched_style_text <- function(text, style, indention = 0L, trailing_empty_line = FALSE) {
      new_text <- original_style_text(text, style, indention, trailing_empty_line)
      if (is.null(new_text)) {
        return(NULL)
      }
      align_trailing_pipes(new_text)
    }

    unlockBinding("style_text", ns)
    assign("style_text", patched_style_text, envir = ns)
    lockBinding("style_text", ns)
  }

  setHook(packageEvent("languageserver", "onLoad"), patch_languageserver_formatter)
})
