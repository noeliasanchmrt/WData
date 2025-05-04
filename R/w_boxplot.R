#' Box plot
#'
#' This function produces a box-and-whisker plot(s) of the given (grouped) values for biased data.
#'
#' @param formula A formula, such as `y ~ grp`, where `y` is a numeric vector of data values to be split into groups according to the grouping variable `grp` (usually a factor). Note that `~ g1 + g2` is equivalent to `g1:g2`.
#' @param data A `data.frame` (or list) from which the variables in `formula` should be taken.
#' @param subset An optional vector specifying a subset of observations to be used for plotting.
#' @param w A function representing the bias function applied to the data points. It must be evaluable and positive in each point of `x`. By default, it is set to the length-biased function.
#' @param qmethod Character string specifying the quantile estimator to be used. Options are `"qf.SBC"` ([`qf.SBC`][WData::qf.SBC()]) and `"qf.sen"` ([`qf.sen`][WData::qf.sen()]). Default is `"qf.SBC"`.
#' @param na.action A function indicating what should happen when the data contain `NA`s. The default is to ignore missing values in either the response or the group.
#' @param xlab,ylab X- and y-axis annotation. Can be suppressed by `ann = FALSE`.
#' @param ann Logical; if `TRUE`, axes are annotated (by `xlab` and `ylab`).
#' @param drop,sep,lex.order Passed to [`split.default`][base::split.default()].
#' @param x For specifying data from which the boxplots are to be produced. Either a numeric vector, or a single list containing such vectors. Additional unnamed arguments specify further data as separate vectors (each corresponding to a component boxplot). `NA`s are allowed in the data. Additional unnamed arguments specify further data as separate vectors.
#' @param range This determines how far the plot whiskers extend out from the box. If range is positive, the whiskers extend to the most extreme data point which is no more than range times the interquartile range from the box. A value of zero causes the whiskers to extend to the data extremes.
#' @param width A vector giving the relative widths of the boxes making up the plot.
#' @param varwidth If `varwidth` is `TRUE`, the boxes are drawn with widths proportional to the square-roots of the number of observations in the groups.
#' @param notch If notch is `TRUE`, a notch is drawn in each side of the boxes.
#' @param outline If outline is not true, the outliers are not drawn (as points whereas S+ uses lines).
#' @param names Character vector or expressions for group labels.
#' @param boxwex A scale factor to be applied to all boxes. When there are only a few groups, the appearance of the plot can be improved by making the boxes narrower.
#' @param staplewex Staple line width expansion, proportional to box width.
#' @param outwex Outlier line width expansion, proportional to box width.
#' @param plot If `TRUE` (the default) then a boxplot is produced. If not, the summaries which the boxplots are based on are returned.
#' @param border An optional vector of colors for the outlines of the boxplots. The values in border are recycled if the length of border is less than the number of plots.
#' @param col If col is non-null it is assumed to contain colors to be used to color the bodies of the box plots. By default they are in the background color.
#' @param log Character indicating if x or y or both coordinates should be plotted in log scale.
#' @param pars A list of (potentially many) more graphical parameters, e.g., `boxwex` or `outpch`; these are passed to `bxp` (if plot is `TRUE`); for details, see [`par`][graphics::par()].
#' @param horizontal Logical indicating if the boxplots should be horizontal; default `FALSE` means vertical boxes.
#' @param add Logical, if `TRUE` add boxplot to current plot.
#' @param at Numeric vector giving the locations where the boxplots should be drawn, particularly when `add = TRUE`; defaults to `1:n` where `n` is the number of boxes.
#' @param ... For the formula method, named arguments to be passed to the default method. For the default method, unnamed arguments are additional data vectors (unless `x` is a list when they are ignored), and named arguments are arguments and graphical parameters to be passed to `bxp` in addition to the ones given by argument pars (and override those in pars). Note that `bxp` may or may not make use of graphical parameters it is passed: see its documentation.
#' @return A list with components:
#' \item{stats}{A matrix, each column contains the extreme of the lower whisker, the lower hinge, the median, the upper hinge and the extreme of the upper whisker for one group/plot. If all the inputs have the same class attribute, so will this component.}
#' \item{n}{A vector with the number of (non-`NA`) observations in each group.}
#' \item{conf}{A matrix where each column contains the lower and upper extremes of the notch.}
#' \item{out}{The values of any data points which lie beyond the extremes of the whiskers.}
#' \item{group}{A vector of the same length as out whose elements indicate to which group the outlier belongs.}
#' \item{names}{A vector of names for the groups.}
#' @details
#' This function is a modification of the [`boxplot`][graphics::boxplot()] function, which is used to create boxplots for biased data by calculating the five-number summary of the data based on [`qf.sen`][WData::qf.sen()] quantile estimator or in [`qf.SBC`][WData::qf.SBC()] quantile estimator.
#' The generic function [`w_boxplot`][WData::w_boxplot()] has a default method ([`w_boxplot.default`][WData::w_boxplot.default()]) and a formula
#' interface ([`w_boxplot.formula`][WData::w_boxplot.formula()]). Parallel boxplots are plotted for multiple groups.
#' Missing values are ignored when forming boxplots.
#' @seealso [`w_boxstats`][WData::w_boxstats()]
#' @examples
#' w_boxplot(Width ~ Replica, data = shrub.data, qmethod = "qf.SBC", main = "Boxplot by Group")
#' w_boxplot(shrub.data$Width, qmethod = "qf.SBC", horizontal = TRUE, main = "Boxplot")
#' w_boxplot(Width ~ Replica, data = shrub.data, qmethod = "qf.sen", main = "Boxplot by Group")
#' w_boxplot(shrub.data$Width, qmethod = "qf.sen", horizontal = TRUE, main = "Boxplot")
#' @aliases w_boxplot w_boxplot.default w_boxplot.formula
#' @name w_boxplot

#' @rdname w_boxplot
#' @order 1
#' @export w_boxplot

w_boxplot <- function(x, ...) {
  UseMethod("w_boxplot")
}

#' @rdname w_boxplot
#' @order 2
#' @export w_boxplot.default
## Default S3 method:
w_boxplot.default <- function(x, w = function(y) {
                                ifelse(y >= 0, y, NA)
                              },
                              qmethod = "qf.SBC",
                              ..., range = 1.5, width = NULL, varwidth = FALSE,
                              notch = FALSE, outline = TRUE, names, plot = TRUE, border = par("fg"),
                              col = "lightgray", log = "", pars = list(
                                boxwex = 0.8, staplewex = 0.5,
                                outwex = 0.5
                              ), ann = !add, horizontal = FALSE, add = FALSE,
                              at = NULL) {
  args <- list(x, ...)
  namedargs <- if (!is.null(attributes(args)$names)) {
    attributes(args)$names != ""
  } else {
    rep_len(FALSE, length(args))
  }
  groups <- if (is.list(x)) {
    x
  } else {
    args[!namedargs]
  }
  if (0L == (n <- length(groups))) {
    stop("invalid first argument")
  }
  if (length(class(groups))) {
    groups <- unclass(groups)
  }
  if (!missing(names)) {
    attr(groups, "names") <- names
  } else {
    if (is.null(attr(groups, "names"))) {
      attr(groups, "names") <- 1L:n
    }
    names <- attr(groups, "names")
  }
  cls <- lapply(groups, class)
  cl <- NULL
  if (all(vapply(groups, function(e) {
    is.numeric(unclass(e)) && identical(
      names(attributes(e)),
      "class"
    )
  }, NA)) && (length(unique(cls)) == 1L)) {
    cl <- cls[[1L]]
  }
  for (i in 1L:n) {
    groups[i] <- list(w_boxstats(
      unclass(groups[[i]]),
      range,
      w = w,
      qmethod = qmethod,
    ))
  }

  stats <- matrix(0, nrow = 5L, ncol = n)
  conf <- matrix(0, nrow = 2L, ncol = n)
  ng <- out <- group <- numeric(0L)
  ct <- 1
  for (i in groups) {
    stats[, ct] <- i$stats
    conf[, ct] <- i$conf
    ng <- c(ng, i$n)
    if ((lo <- length(i$out))) {
      out <- c(out, i$out)
      group <- c(group, rep.int(ct, lo))
    }
    ct <- ct + 1
  }
  if (length(cl) == 1L && cl != "numeric") {
    oldClass(stats) <- oldClass(conf) <- oldClass(out) <- cl
  }
  z <- list(
    stats = stats, n = ng, conf = conf, out = out,
    group = group, names = names
  )
  if (plot) {
    if (is.null(pars$boxfill) && is.null(args$boxfill)) {
      pars$boxfill <- col
    }
    do.call(bxp, c(list(z,
      notch = notch, width = width,
      varwidth = varwidth, log = log, border = border,
      pars = pars, outline = outline, horizontal = horizontal,
      add = add, ann = ann, at = at
    ), args[namedargs]),
    quote = TRUE
    )
    invisible(z)
  } else {
    z
  }
}

#' @rdname w_boxplot
#' @order 3
#' @export w_boxplot.formula

## S3 method for class 'formula'
w_boxplot.formula <- function(formula, data = NULL, w = function(y) {
                                ifelse(y >= 0, y, NA)
                              },
                              qmethod = c("qf.SBC", "qf.sen"),
                              ..., subset, na.action = NULL,
                              xlab = mklab(y_var = horizontal), ylab = mklab(y_var = !horizontal),
                              add = FALSE, ann = !add, horizontal = FALSE, drop = FALSE,
                              sep = ".", lex.order = FALSE) {
  qmethod <- match.arg(qmethod)
  if (missing(formula) || (length(formula) != 3L)) {
    stop("'formula' missing or incorrect")
  }
  if (missing(xlab) || missing(ylab)) {
    mklab <- function(y_var) {
      if (y_var) {
        names(mf)[response]
      } else {
        paste(names(mf)[-response], collapse = " : ")
      }
    }
  }

  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame()))) {
    m$data <- as.data.frame(data)
  }
  m$... <- m$drop <- m$sep <- m$lex.order <- NULL
  m$xlab <- m$ylab <- m$add <- m$ann <- m$horizontal <- NULL
  m$na.action <- na.action
  m$qmethod <- NULL
  m[[1L]] <- quote(stats::model.frame.default)
  mf <- eval(m, parent.frame())
  response <- attr(attr(mf, "terms"), "response")

  w_boxplot(
    split(mf[[response]], mf[-response],
      drop = drop,
      sep = sep, lex.order = lex.order
    ),
    w = w, qmethod = qmethod,
    xlab = xlab, ylab = ylab,
    add = add, ann = ann, horizontal = horizontal, ...
  )
}
