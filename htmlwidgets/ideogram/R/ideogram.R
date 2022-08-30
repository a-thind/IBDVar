#' R binding for ideogram
#'
#' Plots ideogram.
#'
#' @import htmlwidgets
#' @import pals
#' @import DescTools
#' @param file IBIS IBD segments file
#' @export
ideogram <- function(file, width = NULL, height = NULL,
                     elementId = NULL) {
  ibd <- read.table(file, header=TRUE)

  # concatenate two sample names together
  ibd$name <- paste0(ibd$sample1, '-', ibd$sample2)
  # remove first two sample columns
  ibd <- ibd[, -c(1,2)]

  ibd$name <- as.factor(ibd$name)
  # create colours
  color <- pals::glasbey()[1:length(levels(ibd$name))]
  # add alpha
  color <- DescTools::SetAlpha(color, alpha=0.65)
  # map the colours to the samples
  ibd$color <- factor(ibd$name, labels=color)

  annots <- ibd[,c("name", "chrom", "phys_start_pos", "phys_end_pos", "color")]
  colnames(annots)[2:4] <- c("chr","start", "stop")
  annots$chr <- as.character(annots$chr)

  # forward options using x
  x = list(
    data = annots[,-1],
    samples = data.frame(name=levels(annots$name),
                         color=levels(factor(annots$color)),
                         shape=rep("rectangle",
                                   length(levels(annots$name))))
  )



  # create widget
  htmlwidgets::createWidget(
    name = 'ideogram',
    x,
    width = width,
    height = height,
    package = 'ideogram',
    elementId = elementId,
    sizingPolicy = htmlwidgets::sizingPolicy(
      browser.padding = 100
    )
  )
}

#' Shiny bindings for ideogram
#'
#' Output and render functions for using ideogram within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{'100\%'},
#'   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#'   string and have \code{'px'} appended.
#' @param expr An expression that generates a ideogram
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#'   is useful if you want to save an expression in a variable.
#'
#' @name ideogram-shiny
#'
#' @export
ideogramOutput <- function(outputId, width = '100%', height = '400px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'ideogram', width, height,
                                 package = 'ideogram')
}

#' @rdname ideogram-shiny
#' @export
renderIdeogram <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, ideogramOutput, env, quoted = TRUE)
}



