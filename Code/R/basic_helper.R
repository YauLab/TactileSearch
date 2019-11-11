require(stringr)
require(magrittr)

########Common function helpers########

posneg = function(x) {
	c(x,-x) %>% unique %>% sort
}

# helpful when x and/or d are long expressions
plus_minus <- function(x,d) {
    if(missing(d) & is.matrix(x)) {
        d <- x[,2]
        x <- x[,1]
    }
    c(x-d,x+d)
}


# ensure low <= x <= hi
clip = function(x, low = min(x), high = max(x), limits=c(low, high)) {
	x[x < min(limits)] = min(limits)
	x[x > max(limits)] = max(limits)

	return(x)
}

`range<-` = function(x, values) {
    clip(x, limits=range(values))
}

# a range function that names it output
# helpful when used as part of a larger summary
range <- function(..., na.rm=FALSE) {
    res <- base::range(..., na.rm=na.rm)
    res %>% set_names(c('L', 'U'))
}

# a require function that installs packages that are not available
require <- function(package, lib.loc = NULL, quietly = FALSE, warn.conflicts = TRUE, character.only = FALSE) {
    # if a symbol is supplied, convert it to a character
    if(!character.only)
        package <- as.character(substitute(package))

    # we've already done the symbol conversion, so set character.only to TRUE to avoid
    # a second conversion
    if(!base::require(package, lib.loc, quietly, warn.conflicts, character.only=TRUE)) {
        install.packages(package)
        base::require(package, lib.loc, quietly, warn.conflicts, character.only=TRUE)
    }
}

get_column_from_list = function(ll, i) t(sapply(ll, function(m) m[,i]))

# convert a list of numeric vectors into a matrix whose rows
# are the vectors of the list. Names from the list are copied to the
# rows of the matrix, names from the (first) vector are copied to the columns
# NO checking is done to see if the names of the vectors are the same...
list_to_matrix = function(ll) {
    ll %>% unlist %>%
        matrix(nrow = length(ll[[1]])) %>%
        set_colnames(ll %>% names) %>%
        set_rownames(ll[[1]] %>% names) %>%
        t
}

#
# helper function for applying a list of functions to a matrix
# does some basic error checking and results formatting
mat_summary = function(mat, FUNS, NAMES='', IND=2, na.rm=TRUE) {
    if(is.null(mat)) return(rep(NULL, length=length(FUNS)))

    mat %<>% as.matrix

    # we're taking advantage of the scoping to
    # emphasize what is (and is not) changing each run through
    # sapply below
    summ <- function(FUN) {
        apply(mat, IND, FUN, na.rm=na.rm)
    }

    FUNS %>%
        sapply(summ) %>%
        set_colnames(rep(NAMES, length.out=ncol(.)))
}

# returns a function that will work on numeric vectors or
# the columns of a numeric matrix.
# Better than calling apply everywhere?
apply_to_vector_or_matrix = function(FUNS, NAMES){
    return(function(x, na.rm=TRUE, IND=2) {
        if(is.matrix(x))
            return (mat_summary(x, FUNS, NAMES, na.rm=na.rm, IND=IND))

        sapply(FUNS, function(FUN) FUN(x, na.rm=na.rm)) %>%
            set_names(NAMES)
    })

}

# mean +/- se
m_se = apply_to_vector_or_matrix(list(mean, se), c('mean', 'se'))
m_sem = m_se

#mean +/- sd
m_sd = apply_to_vector_or_matrix(list(mean, sd), c('mean', 'sd'))

# zscores along an axis
zmat = function(m, idx = 1) {
	apply(m, idx, scale) %>% t
}

# running mean smoother
# only accepts odd window sizes. If you give it an even, it will go +1
runmean <- function(x, k=3) {
    if(! (k %% 2)) k <- k + 1

    d <- (k-1) / 2

    n <- length(x)

    y <- rep(0,n)

    # after profiling, the for loop turns out to be no slower (and sometimes appreciably faster)
    # than a few other approaches using vectorization and/or par*pply and/or doing three separate loops without the if statements
    for(ii in 1:n) {
        #handle the early edge
        s <- ii-d
        if(s < 1) s <- 1

        #handle the late edge
        f <- ii+d
        if(f>n) f <- n

        y[ii] <- mean(x[s:f])
    }

    return(y)
}




# function helper to allow aggregate to be used in magrittr without the
# bracket notation
do_aggregate = function(data, formula, FUN, ...) {
    aggregate(formula, data=data, FUN, ...)
}

#wrapper for parLapply that puts the data FIRST
do_parL <- function(x, FUN, ..., CLUST) {
    if(missing(CLUST)) {
        CLUST <- makeForkCluster(detectCores())
        on.exit({stopCluster(CLUST)})
    }

    parLapply(CLUST, x, FUN, ...)
}

# decorator helpers., NO OPS
# Many functions work by side effect but fail to return their input.
# These wrappers help that... I'm thinking about axis deocrator etc or even plotting the output of a commadn
# but still returning what you want, e.g., m1 <- lm(...) %>% print_summmary %>% F_NOOP(plot)
#

# print function output but return input
print_summary = function(x, FUN=summary, ...) {
    x %>% FUN(...) %>% print
    invisible(x)
}

F_NOOP <- function(x, FUN, ...) {
    FUN(x, ...)
    invisible(x)
}

NOOP <- function(x, expr=NULL) {
    eval(expr); invisible(x)
}

# make (r|c)bind magrittr compatible
rbind_list <- function(ll) do.call(rbind, ll)
cbind_list <- function(ll) do.call(cbind, ll)

# easy column applier
column_apply = function(mat, FUN, ...) apply(mat, 2, FUN, ...)

# easy row applier
row_apply = function(mat, FUN, ...) apply(mat, 1, FUN, ...)

# row applier with an additional index variable
row_apply_ii <- function(mat, FUN_, ...) {
    ii <- 0
    row_apply(mat, function(...) {
        ii <<- ii + 1
        FUN_(..., ii=ii)
    })
}

# sapply with additional index variable
sapply_ii <- function(X, FUN_, simplify=TRUE, USE.NAMES=TRUE, ...) {
    ii <- 0
    sapply(X, function(...) {
        ii <<- ii + 1
        FUN_(..., ii=ii)
    }, simplify = simplify, USE.NAMES = USE.NAMES)
}


# helper to do row scaling
row_scale <- function(mat) apply(mat, 1, scl)

# make it easier to say not is.na in a pipe'd context
not_NA = function(x) !is.na(x)

# like which.min, but for equality
# useful when an expression for x or y is long
which.equal = function(x,y) which(x==y)

# NA aware standard error function
se = function(x, na.rm = TRUE) {
    # sqrt_length <- x %>% not_NA %>% sum %>% sqrt
    sqrt_length <- sqrt(sum(not_NA(x)))

    if(!na.rm) {
        sqrt_length <- sqrt(length(x))
    }

    sd(x, na.rm = na.rm) / sqrt_length
}

is_between = function(x, rng) {
    (x<=max(rng)) & (x>=min(rng))
}

is_within <- is_between

# this is related to %in% but checks for A contained_in the interval of B
# rather than A is_element_of B
`%within%` <- function(x, rng) {
  is_between(x, rng)
}

# quick way to arrange a variable by another one
sort_by = function(x, by) x[order(by)]

#sometimes we don't need the last item
# if you give me < 1 I will return the full vector with a warning
# this is helpful as it (I think) avoids error checking and
# should be reasonable for most events
remove_tail = function(x, k=1) {

    if(k<1) {
        warning('Tried removing less than 1 element, k = ', k)
        return(x)
    }

    stopifnot(k>=1 & k<length(x))
    #seems like it's faster to rewrite as selecting from the beginning
    #rather than using negative indexing
    # x[-(length(x):(length(x) - (k-1)))]
    x[1:(length(x)-k)]
}

# get the name from a file
get_fname = function(full_path, keep_extension = FALSE) {

    name <- full_path %>%
	    str_split("/") %>%
	    unlist %>%
	    tail(1)

	if (keep_extension)
		return(name)

    name %>%
        str_split("\\.") %>%
        unlist %>%
        remove_tail %>%
        # because we split on '.' above, collapse back
        paste(collapse='.')
}


# concatenate strings in sequence
`%&%` <- function(s1, s2) paste0(s1,s2)

# go into a directory, execute the code, then go back
do_in_directory = function(dirname, expr) {
    # dir.exists also catches mal-formed dirnames
    if(!dir.exists(dirname)) stop('Directory does not exist: ', dirname)

    # make sure we return to cwd (this works even if there is an error
    cwd = getwd()
    on.exit(setwd(cwd))

    setwd(dirname)
    retval = eval(expr)

    return(retval)
}


# Call FUN() NSIM_ times in parallel. all inputs to FUN must be named explicitly
# Wrapper for parSapply
parReplicate <- function(NSIM_, cl=NULL, FUN, ..., simplify=TRUE, USE.NAMES=TRUE, close_cluster=FALSE) {

    result <- parSapply(cl, seq(NSIM_), function(ii) {
        FUN(...)
    }, simplify = simplify, USE.NAMES = USE.NAMES)

    if(close_cluster) stopCluster(cl)

    return(result)
}


#
# add NAs to columns of a matrix so that you can output rectangular data
# we're relying on sapply to create a matrix for us
pad_matrix <- function(...) {
	len = max(sapply(list(...), length))
	sapply(list(...), function(x) c(x, rep(NA, len - length(x))))
}

#
# strip out NA elements in lists and vectors
#
# I think setting a list element to NULL also removes it...
remove_na_from_collection = function(collection) {
	if (is.list(collection))
		return(lapply(which(!is.na(collection)), getElement, object = collection))
	if (is.vector(collection))
		return(collection[!is.na(collection)])

	stop("Input must be one of {list|vector}")
}

#
# divide each value by the sum total
# although simple, it is often used to wrap results of other function calls,
# thus avoiding intermediate storage variables
# by default it removes NA values from the sum,
# but leaves them in the resulting vector
#
scl = function(x, na.rm = TRUE) {
	x/sum(x, na.rm = na.rm)
}


#
# this is really slow if you try to use pipes
rmse = function(y, yhat) sqrt(mean((y-yhat)^2))

# cat decorator that produces time stamps
tcat <- function() {
    t0 <- proc.time()
    t1 <- t0
    return (function(..., reset=FALSE){
        t2 <- proc.time()
        cat(...)
        if(reset) {
            t0 <<- t2
        } else {
            cat(', ', (t2-t1)[3], '| Total:', (t2-t0)[3])
        }
        cat('\n')
        t1 <<- t2
    })
}



# # # # #	Wrapper functions	# # # # #

# checks for the existence of aptly nm'd RDS file and
# returns its contents rather than repeating the expression evaluation
with_cache <- function(nm, expr, cache_dir='./', refresh=FALSE) {
    fname <- cache_dir %&% nm
    if(!endsWith(fname, '.RDS')) {
        fname <- fname %&% '.RDS'
    }

    if(file.exists(fname) & !refresh) {
        print('reading '%&% fname %&% ' from disk')

        res <- readRDS(fname)
    } else {
        res <- eval(expr)

        print('writing '%&% fname %&% ' to disk')
        saveRDS(res, file=fname)
    }

    return (res)
}

rm_cache <- function(var) {
    var = as.character(eval(substitute(quote(var))))
    fname <- var %&% '.RDS'
    if(!file.exists(fname)){
        warning(fname, ' does not exist')
    } else {
        unlink(fname)
    }
}

`%C%` <- function(var, expr, refresh=FALSE) {

    env = parent.frame()
    var = as.character(eval(substitute(quote(var))))
    expr = substitute(expr)

    fname <- var %&% '.RDS'

    # if the file exists and we aren't refreshing, read it in
    rds <- list(key='', value='')
    if(file.exists(fname) & !refresh) {
        rds <- readRDS(fname)
        res <- rds$value

        # this is needed for version 1 caches that didn't have keys
        if(!exists('key', rds)) {
            rds$key = 'none'
        }
        # print('old key: ' %&% rds$key)
    }

    key <- digest::digest(expr)
    # print('new key: ' %&% key)

    # if we're refreshing, or the file doesn't exist, or RHS has changed, then recalc
    # this should short-circuit after condition 2 so that we don't get an error on condition 3
    if(refresh | is.null(rds) | rds$key!=key ) {
        print('Calculating new value')
        res <- eval(expr)

        saveRDS(list('key'=key, 'value'=res), file=fname)
    } else {
        print('Loading old value')
    }

    env[[var]] = res

    return(invisible(env[[var]]))
}

# prints out a counter each time f is called
with_counter = function(FUN, end = "unknown") {
	counter__ <- 1
	function(...) {
		cat("Starting Run", counter__, "of", end, "\n")
		counter__ <<- counter__ + 1
		FUN(...)
	}
}

# wraps f with a timer
with_timer = function(FUN, CLEANUP = NULL) {
	function(...) {
		begin = proc.time()

		out = FUN(...)

		end = proc.time()
		cat("\tElapsed Time: ", as.numeric(end - begin)[3], "\n")

		if (is.function(CLEANUP)) {
			CLEANUP(out)
		}

		return(out)
	}
}
