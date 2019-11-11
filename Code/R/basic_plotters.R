########Common plotting helpers########
source('/Users/liz/Desktop/TactileSearch/basic_helper.R')

# constants
gr = 1.61803398875


###Ebar functions

# primary erorr bar function, makes guesses at y and sem if you don't specify them directly
# 3 use cases:  supply three vectors: x, y, sem
#               supply a matrix, treated as = {y, sem}, x is 1:length(y)
#               supply two vectors: treated as y and sem, x is 1:length(y)
#
ebars = function(x, y=NULL, sem=NULL, length = 0.05, type='n', col='black', pt.col=col, code=2, show_0=FALSE, ...) {
    if(is.null(y)) {
        if(is.matrix(x)) {
            y <- x[,1]
            sem <- x[,2]
        } else {
            y <- x
        }
        x <- seq_along(y)
    }

    if(is.matrix(y)) {
        sem <- y[,2]
        y <- y[,1]
    }

    if(is.null(sem)) {
        sem <- y
        y <- x
        x <- seq_along(y)
    }

    if(!show_0) {
        ind <- (sem > 0)
        ebars.y(x[ind], y[ind], sem[ind], length, code=code, col=col, ...)
    } else {
        ebars.y(x, y, sem, length, code=code, col=col, ...)
    }

    points(x, y, type=type, col=pt.col, ...)
}

ebars.x = function(x, y, sem, length = 0.05, code = 2, show0=FALSE, ...) {
    if(show0) {
        arrows(x - sem, y, x + sem, y, angle = 90, code = code, length = length, ...)
    }
    else {
        ind <- sem > 0
        arrows(x[ind] - sem[ind], y[ind], x[ind] + sem[ind], y[ind], angle = 90, code = code, length = length, ...)
    }
}

ebars.y = function(x, y, sem, length = 0.05, up = T, down = T, code = 2, ...) {
    if (up) {
        arrows(x0 = x, y0 = as.numeric(y), y1 = as.numeric(y + sem), angle = 90, code = code, length = length, ...)
    }
    if (down) {
        arrows(x0 = x, y0 = as.numeric(y), y1 = as.numeric(y - sem), angle = 90, code = code, length = length, ...)
    }
}

# the dots are sent only to the lines function
draw.line_and_error = function(x, y, sem, col ='black', type='o', pch=16, ebar.col=col, lwd=1, ebar.lwd=lwd, ebar.code=2, ...) {
    lines(x, y, col=col, lwd=lwd, type=type, pch=pch, ...)
    ebars.y(x,y,sem,col=ebar.col, lwd=ebar.lwd, code=ebar.code)
}


ebar_polygon = function(x, y, sem, col='black', alpha=100, border = NA, add_line=FALSE, lwd=1, ...) {
    # # have they already give a color with alpha?
    # if(startsWith(col, '#') & nchar(col)==9) {
    #     fill <- col
    #     stroke <- substr(col, 1, 7)
    # } else {
         stroke <- col
         fill <- getAlphaRGB(col, alpha)
    # }

    polygon(c(x, rev(x)), c(y + sem, rev(y - sem)), border = border, col = fill, ...)

    if(add_line) lines(x, y, col=stroke, lwd=lwd)
}



# clean barplot with no x axis
simpleBars = function(heights, col = c("blue", "darkgreen", "gray", "black"), ...) {
    barplot(heights, col = col, border = NA, axes = F, axisnames = F, ...)
    axis(2, labels = F, at = c(0, 0.5 * max(heights), max(heights)))
}

# the idea is to wrap any ole plot function with one a function that adds error bars
# we have to special case barplot
# right now this function only support y-axis error bars
with_error_bars = function(plot_func) {

    return(function(x, y = NULL, sem, col = "black", down = T, up = T, ...) {
        ## special case for barplot here
        if (body(plot_func) == body(barplot)) {
            x = plot_func(y, col = col, border = NA, ...)
            down = F
        } else {
            if (is.na(x)) {
                x = 1:length(y)
            }
            plot_func(x, y, col = col, ...)
        }

        ebars.y(x, y, sem, up = up, down = down, col = col)

        return(invisible(x))
    })
}

#Plots line with error bars, each column is point
ebar_plot = function(x=1:ncol(mat), mat, overlay = F, type = "o", lwd = 3, col = "black", ebar.lwd = (lwd - 1), code = 2, SUMM = m_sem,
    ...) {
    summ = SUMM(mat)

    if (overlay == F) {
        plot(x, summ[, 1], type = type, lwd = lwd, col = col, ...)
    } else {
        lines(x, summ[, 1], type = type, lwd = lwd, col = col, ...)
    }

    ebars.y(x, summ[, 1], summ[, 2], lwd = ebar.lwd, col = col, code = code)

    invisible(summ)
}


plotManyLines = function(x = 1:ncol(mat), mat, col = "black", ylim = NA, ...) {
    if (any(is.na(ylim))) {
        ylim = range(c(mat))
    }
    plot(x, mat[1, ], type = "n", ylim=ylim, ...)
    for (i in 1:nrow(mat)) lines(x, mat[i, ], col = col)
}

getAlphaRGB = function(colname, alpha) {
    c = col2rgb(colname)
    rgb(t(c), alpha = alpha, maxColorValue = 255)
}

getAlphaRGB_as_RGB = function(colname, alpha) {
    col = col2rgb(colname)
    col = col * (255-alpha) / 255
    rgb(t(col), maxColorValue = 255)
}

poly = function(x, lower, upper, col = "gray", ...) {
    lower = rep(lower, length.out = length(x))
    upper = rep(upper, length.out = length(x))
    polygon(c(x, rev(x)), c(upper, rev(lower)), col = col, ...)
}

bplot = function(x, heights, sems, width, col = rep("gray", length(heights)), ebar.col=col, ebar.code=3, ebar.lwd=3,
                 base = rep(0, length.out = length(heights))) {
    col <- rep_len(col, length.out = length(heights))
    base <- rep_len(base, length.out = length(heights))

    for (i in 1:length(heights)) {
        xp = c(x[i] - width, x[i] + width)
        poly(xp, base[i], heights[i], col[i], border = NA)
    }

    if(any(sems > 0))
        ebars.y(x, heights, sems, col = ebar.col, code=ebar.code, lwd = ebar.lwd)
}


#
# This is way too specifc, relying on ASYNCS to be defined. Move it out!!
barplot.group = function(y1, y2, col = c("orange", "blue")) {
    col = rep(col, length.out = ncol(y1) + ncol(y2))

    m1 = m_sem(y1)
    m2 = m_sem(y2)

    a = barplot(t(cbind(m1[, 1], m2[, 1])), beside = T, col = col, names.arg = ASYNCS, xlim = c(1, 40), space = c(0,
        0.75), border = F, ylim = c(0, 0.06))

    for (i in 1:nrow(m1)) {
        ebars.y(a[1, i], m1[i, 1], m1[i, 2], lwd = 2, col = col[1], down = F)
        ebars.y(a[2, i], m2[i, 1], m2[i, 2], lwd = 2, col = col[2], down = F)
    }
}

stack.plot = function(group, y, col, xAt = 10*0:4, ...) {
    groups = sort(unique(group))
    y.uni = sort(unique(y))
    counts = sapply(y.uni, function(g) max(sum(y[group == groups[1]] == g, na.rm = T), sum(y[group == groups[2]] ==
            g, na.rm = T)))
    plot(range(y.uni, na.rm = T), range(c(0, 5 + counts), na.rm = T), type = "n", axes = F, ylab = "", xlab = "", ...)
    axis(1, at = xAt, labels = T)
    for (i in 1:length(groups)) {
        yg = y[group == groups[i]]
        y.uni = sort(unique(yg))
        counts = sapply(y.uni, function(g) sum(yg == g, na.rm = T))

        putStack(y.uni, counts, col = col[i])

        abline(v = median(rep(y.uni, counts)), col = col[i], lwd = 3)
    }
}

putStack = function(y.uni, counts, col, pch = 20, horizontal = FALSE, offset=0, ...) {
    col = rep(col, length.out = length(counts))
    pch = rep(pch, length.out = length(counts))

    if (horizontal) {
        for (i in 1:length(counts)) {
            points(offset + 1:counts[i], rep(y.uni[i], abs(counts[i])), col = col[i], pch = pch[i], ...)
        }

    } else {
        for (i in 1:length(counts)) {
            points(rep(y.uni[i], abs(counts[i])), offset + 1:counts[i], col = col[i], pch = pch[i], ...)
        }

    }
}

plot.clean = function(xlim, ylim=xlim, x = 1, y = 1, type = "n", xlab="", ylab="", axes=FALSE, ...) {
    if(type!='n' & all(x==1) & all(y==1)) {
        x <- xlim
        y <- ylim
    }


    plot(x, y, type = type, axes = axes,
         ylab = ylab, xlab = xlab, xlim = range(xlim), ylim = range(ylim), ...)
}


plot.line = function(...) {
    plot(..., type='l')
}


# draw.axis = function(side, at, tcl=-0.3, labels=at, padj=0.5, adj=1, yline=0.65, xline=0.2, ...) {
#     if(length(side) > 1) {
#         return (invisible(sapply(side, draw.axis,
#             at=at, tcl=tcl, labels=labels, padj=padj, adj=adj, yline=yline, xline=xline, ...)
#         ))
#     }
#     axis(side, at=at, labels=F, tcl=tcl, ...)
#
#     if(side%%2 == 1)	mtext(labels, side=side, at=at, padj=padj, las=1, line=xline)
#     if(side%%2 == 0)	mtext(labels, side=side, at=at, adj=adj, las=1, line=yline)
#
#     invisible()
# }

draw.axis <- function(side, at, tcl=-0.3, labels=at, las=1, cex.axis=1, mgpy=c(3, .5, 0), mgpx=c(3, .4, 0), ...) {
  if(length(side) > 1) {
    return (invisible(sapply(side, draw.axis,
      at=at, tcl=tcl, labels=labels, cex.axis=cex.axis, las=las, ...)
    ))
  }
  mgp <- mgpy
  if(side %% 2) mgp <- mgpx

  axis(side, at=at, labels=labels, tcl=tcl, mgp=mgp, cex.axis=cex.axis, las=las, ...)
}

# for legacy reasons...
my_axis = draw.axis


# functions that add summaries to plots
# this is the base function, which draws the mean and optionally the SE
draw.FUN <- function(CENTER=mean, ERROR=se) {
    return (function(x,y,col='black', lwd=4, delta=0.1*range(x), add_error_poly=FALSE, error_alpha=150, ...) {
        lines(plus_minus(x,delta), rep(CENTER(y),2), lwd=lwd, col=col, ...)

        if(add_error_poly) {
            ebar_polygon(plus_minus(x, delta), CENTER(y) %>% rep(2), ERROR(y), col=getAlphaRGB(col, error_alpha))
        }
    })
}
draw.mean = draw.FUN()
draw.median = draw.FUN(median, mad)


# jitter points at each point
jpoints <- function(x, y, amt, jitter='x', col='black', ...) {
    x <- rep(x, length(y))
    if(jitter == 'x' || jitter == 'xy') {
        x = runif(length(y), x-amt, x+amt)
    }
    if(jitter == 'y' || jitter == 'xy') {
        y = runif(lengthy(y), y-amt, y+amt)
    }
    points(x,y, col=col, ...)
}

draw.points <- function(x, y, jx=0, jy=0, col='black', pch=16, ...) {
    if(missing(y)) {
        if(is.matrix(x)) {
            y <- x[,2]
            x <- x[,1]
        }
        else {
            y <- x
            x <- seq_along(y)
        }
    } else {
        # if length(x)==length(y) this does nothing
        x <- rep(x, length(y))
    }

    # jittering
    if(jx>0) x = runif(length(y), x-jx, x+jx)

    if(jy>0) y = runif(length(y), y-jy, y+jy)

    points(x,y, col=col, pch=pch, ...)

    invisible(cbind(x,y))
}




# confusion matrix
# the matrix is returned invisibly, allowing function chaining
plot.confusion_matrix <- function(mat, xlab='Predicted', ylab='Actual', cols=NA, zlim=NULL, max_col = 'dodgerblue3', color_steps=100) {
    stopifnot(diff(dim(mat))==0)

    if(is.na(cols))
        cols <- colorRampPalette(c('white', max_col), interpolate='linear')(color_steps)


    if(is.null(zlim))
        zlim = range(c(mat))

    sz = nrow(mat)

    rev_col <- function(mat) mat[,sz:1]

    mat %>% t %>% rev_col %>% image(col=cols, axes=F, x=seq(sz), y=seq(sz),
        ylab='', xlab='', xlim=c(0.5,sz+0.5), ylim=c(0.45,sz+0.55), zlim=zlim,
        asp=1, useRaster=TRUE)

    #hrz lines
    segments(x0=0.5, x1=sz+0.5, y0=0:sz + 0.5,y1=0:sz + 0.5)

    #vrt lines
    segments(x0=0:sz + 0.5,x1=0:sz + 0.5, y0=0.5, y1=sz+0.5)

    mtext(xlab, 2)
    mtext(ylab, 3)

    invisible(list('cmat'=mat, 'ramp' = cols, 'limits'=zlim, 'step_size'=color_steps))
}

# make a color bar using the information returned by plot.confusion_matrix
# oriented vertically unless orient == 'H'
plot.color_bar <- function(par_list, orient='V', ...) {
    mat <- matrix(seq(par_list$limits[1], par_list$limits[2], length=par_list$step_size)) %>% t

    if(substr(orient, 1, 1) == 'H')
        mat %<>% t

    image(mat, axes=F, col=par_list$ramp, zlim=par_list$limits, ...)
    box()
}

plot_pred = function(yhat, col = "black", ylim = c(-0.05, 1.05), names = NA) {
    x = barplot(yhat, col = col, ylim = ylim, axes = F, names = names, border = NA)
    axis(2, at = c(0, 0.5, 1), labels = F)
    abline(h = 0, lwd = 0.5)
    return(x)
}

plot.barfit = function(y, yhat, y.col = "black", yhat.col = "gray", error.col = "orangered", error_on_top = T, do_mean = T,
    names = NA) {
    par(mar = rep(1, 4)) #put this here for now

    #tack on the mean
    if (do_mean) {
        yhat = c(yhat, NA, mean(yhat))
        y = c(y, NA, mean(y))
    }

    x = plot_pred(yhat, col = yhat.col, names = names)

    if (error_on_top) {
        bar_stool(x, y, wid = 0.4, lwd = 3)
        error_line(x, y, yhat, col = error.col)
    } else {
        error_line(x, y, yhat, col = error.col, lwd = 1)
        bar_stool(x, y, wid = 0.4, lwd = 3)
    }
}

bar_stool = function(x, y, ylo = 0, wid = 0.2/length(x), col = "black", lwd = 2, legs = F, ...) {
    for (i in seq(x)) {
        if (legs) {
            sapply(c(-wid, +wid), function(q) lines(c(x[i] + q, x[i] + q), c(ylo, y[i]), col = "darkgray", lwd = 1))
        }
        lines(pm(x[i], wid), c(y[i], y[i]), col = col, lwd = lwd, ...)
    }
}

error_line = function(x, y, y2, col = "orangered", lwd = 2) {
    for (i in seq(x)) lines(c(x[i], x[i]), c(y[i], y2[i]), col = col, lwd = lwd)
}

# # # # #	Wrapper functions	# # # # #
to_pdf = function(PLOT, w, h, mar=rep(1,4)) {
    return(function(..., fname, width = w, height = h, margin=mar) {
        on.exit(dev.off())

        fname <- fix_pdf_name(fname)
        pdf(fname, width = width, height = height, useDingbats=FALSE)
        par(mar=margin)
        PLOT(...)
    })
}

fix_pdf_name <- function(fname) {
    if(!grepl("\\.pdf$", fname)) {
        fname = paste0(fname, ".pdf")
    }
    return(fname)
}

# a more flexible pdf wrapper that evaluates an arbitrary expression.
as_pdf = function(fname, w, h, expr, TEST=FALSE, bg='white') {
    if(! TEST) {
        on.exit(dev.off())

        fname <- fix_pdf_name(fname)
        pdf(fname, width=w, height=h, useDingbats=FALSE, bg=bg)
        res = eval(expr)
    } else {
        res = eval(expr)
    }
    return (invisible(res))
}
