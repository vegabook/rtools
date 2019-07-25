# We will produce a nice set of tools to use in R
library(xts)
library(lattice)
library(quantreg)
library(grid) # for pushing viewports onto ggplots
library(glmnet)

if(Sys.info()["sysname"] == "Windows") library("Rblpapi") # only load bberg if on Windows

regress <- function(d1, d2, hiDays = c(5, 21, 65), legnd = TRUE, greyalpha = 0.7, numlimit = NA, 
                    firstline = TRUE, secondline = TRUE, legendperiod = "days", 
                    NArange = NULL, orthogonal = T, ...) {
    #This function will plot a regression chart of the two securities specified
    if(orthogonal) lm <- orthlm # override linear model to orthogonal regression if required
    if(!is.null(NArange)) {
        ex1 <- as.numeric(d1[NArange])
        ex2 <- as.numeric(d2[NArange])
        d1[NArange] <- NA
        d2[NArange] <- NA
    }
    n1 <- as.numeric(d1)
    n2 <- as.numeric(d2)
    minlen <- min(length(n1), length(n2))
    if(!is.na(numlimit)) minlen <- min(minlen, numlimit)
    n1 <- last(n1, minlen)
    n2 <- last(n2, minlen)
    # adjust hiDays maxima
    hiDays[1] <- min(hiDays[1], minlen)
    hiDays[2] <- min(hiDays[2], minlen)
    hiDays[3] <- min(hiDays[3], minlen)
    drows <- length(n1)
    cl <- rep(addAlpha("grey30", greyalpha), drows)
    cl[(drows - hiDays[3]):drows] <- 4
    cl[(drows - hiDays[2]):drows] <- 3
    cl[(drows - hiDays[1]):drows] <- 2
    cl[drows] <- "orange"
    cx <- rep(1, drows)
    cx[drows] <- 1.5
    plot(n1, n2, pch = 16, col = cl, , cex = cx, ...)
    if(!is.null(NArange)) {
        points(ex1, ex2, col = "grey80")
    }
    lmlong <- lm(n2 ~ n1)
    lmshort <- lm(last(n2, hiDays[3]) ~ last(n1, hiDays[3]))
    rsqlong <- summary(lmlong)$r.squared
    rsqshort <- summary(lmshort)$r.squared
    selong <- last(lmlong$residuals) / sd(lmlong$residuals)
    seshort <- last(lmshort$residuals) / sd(lmshort$residuals)
    coeffshort <- coef(lmshort)
    coefflong <- coef(lmlong)
    coeffshortstring <- paste("y=", round(coeffshort[2], 3), "x", ifelse(coeffshort[1] >= 0, "+", "-"), 
                              round(abs(coeffshort[1]), 3), sep = "")
    coefflongstring <- paste("y=", round(coefflong[2], 3), "x", ifelse(coefflong[1] >= 0, "+", "-"), 
                              round(abs(coefflong[1]), 3), sep = "")
    if(firstline) abline(lmlong, lty = "solid", col = "grey", lwd = 1.5)
    if(((length(n1) / hiDays[3]) > 2) && secondline) {   # if there is enough data plot a second regression line
        abline(lmshort, col = "blue", lwd = 1.5, lty  = "dotted")
    }
    if(legnd) {
        regdirection <- mean(coef(lmlong)[2], coef(lmshort)[2])
        statstext <- c(" ",         # stats legend text, first line blank to move down from border
                     paste("rsq(", length(n1), "): ", round(rsqlong, 3)),
                     c(paste("se(", length(n1), "): ", round(selong, 3), ", val: ", round(last(lmlong$residuals), 2), sep = "")),
                     coefflongstring, 
                     paste("rsq(", hiDays[3], "): ", round(rsqshort, 3), sep = ""),
                     c(paste("se(", hiDays[3], "): ", round(seshort, 3),", val: ", round(last(lmshort$residuals), 2), sep = "")), 
                     coeffshortstring)
        legend(ifelse(regdirection > 0, "topleft", "topright"), statstext, bty = "n", bg = addAlpha("white", 0.5), 
               cex = 0.8, text.col = addAlpha("grey30", greyalpha))
        ltext = c(paste("last", length(n1), legendperiod),
                  paste("last", hiDays[3], legendperiod),    
                  paste("last", hiDays[2], legendperiod),    
                  paste("last", hiDays[1], legendperiod),    
                  "Now")
        lcol = c(cl[1], 4, 3, 2, "orange")
        if(!is.null(NArange)) {
            lcol = c(lcol, "grey80")
            ltext = c(ltext, "ommitted")
        }
        legend(ifelse(regdirection > 0, "bottomright", "bottomleft"), 
               legend = ltext, col = lcol, pch = c(rep(16, length(lcol) - 1), 1), bg = addAlpha("white", 0.5))
    }
    return(lmlong$residuals)
}

minuteCorrel <- function(inxts) {
# this will take an xts NOT of returns but of levels, that must be minute by minute
# it will then give a correlation matrix based on pairwise complete sets that have
# the correct periods to avoid gaps that are too big. 
    # find the mode
    idx <- index(inxts)
    gaps <- as.numeric(idx[-1]) - as.numeric(idx[-length(idx)]) # gaps in seconds
    themode <- as.numeric(names(table(gaps))[1])
    print(paste("Mode is", themode, "seconds"))
    if(!("POSIXct" %in% class(index(inxts)))) {
        flushprint("Need a POSIX index on an xts")
    }
    correls <- lapply(colnames(inxts), function(x) {
        sapply(colnames(inxts), function(y) {
            ddata <- na.omit(inxts[, c(x, y)])
            if(nrow(ddata) < 10) { # some small number minimum ddata necessary
                return(NA)
            }
            idx <- index(ddata)
            gaps <- as.numeric(idx[-1]) - as.numeric(idx[-length(idx)]) # gaps in seconds
            goodgaps <- (gaps == themode)
            rets <- logret(ddata)[goodgaps, ]
            return(cor(rets[, 1], rets[, 2]))
        })
    })
    correls <- do.call(cbind, correls)
    colnames(correls) <- colnames(inxts)
    return(correls)
}


makePCs <- function(sseriess, scale.unit = FALSE, equalize.vol = FALSE, decayhl = NA, colweights = NA, reorder = F, 
                    windowlen = NA) {
# returns the pcs given an eigenvector and matrix of stuff)
# scale.unit puts all means to zeros and standard deviations to one
# equalize.vol puts all series vol to same as average
    if (!is.na(decayhl)) {
        decayer <- decay(nrow(sseriess), decayhl, sumone = T)
        sseriess <- apply(sseriess, 2, function(x) x * decayer)
    }
    if(scale.unit) {
        sseriess <- scale(sseriess)
    } else if(equalize.vol) {
        vols <- apply(sseriess, 2, sd)
        newseries <- sseriess %*% diag(mean(vols) / vols)
        if("xts" %in% class(sseriess)) {
            newseries <- xts(newseries, order.by = index(sseriess))
        }
        sseriess <- newseries
    }
    if(!is.na(first(colweights))) {
        sseriess <- sseriess %*% diag(colweights)
    }
    pp <- eigen(cov(sseriess))
    eigs <- pp$values
    loadings <- pp$vectors
    # change order of factors for yield curve PCA so that level, slope, curve are in that order
    if(reorder) {
        loadrank <- apply(loadings, 2, function(x) length(rle(x > 0)$lengths))
        loadings <- loadings[, order(loadrank)]
        eigs <- eigs[order(loadrank)]
    }
    if(loadings[1, 2] > loadings[ncol(loadings), 2]) {
        loadings[, 2] <- -loadings[, 2]
    }
    testloadindex <- as.integer(round(length(pp$values) / 2))
    if(loadings[testloadindex, 1] < 0) {
        loadings[, 1] <- -loadings[, 1]
    }
    pcs <- apply(loadings, 2, function(x) sseriess %*% x)
    if("xts" %in% class(sseriess)) pcs <- xts(pcs, order.by = index(sseriess))
    # now check if we need to see the eigenproportion
    if(!is.na(windowlen)) {
        rollprop <- rollapply(sseriess, windowlen, function(x) {
            rotated <- x %*% loadings
            vars <- apply(rotated, 2, var)
            return(vars / sum(vars))
        }, by.column = FALSE)
        return(list(pcs = pcs, loadings = loadings, eigs = eigs, indata = sseriess, 
                    rollprop = rollprop))
    } else {
        return(list(pcs = pcs, loadings = loadings, eigs = eigs, indata = sseriess))
    }
}


hair <- function(inmat, recentamount = 65, 
                 haircol = "grey", recentcol = "darkgoldenrod3", yestcol = NA, mats = NA,
                 nowcol = "red", withchange = FALSE, title = "yield curve evolution", withlegend = FALSE, 
                 miny = NA, maxy = NA, xlab = "maturity", ylab = "yield") {
    if(is.na(mats)) mats = as.numeric(colnames(inmat))
    if(is.na(miny)) miny = min(inmat)
    if(is.na(maxy)) maxy = max(inmat)
    recent = last(inmat, recentamount)
    laster = as.numeric(last(inmat))
    last2 = as.numeric(first(last(inmat, 2)))
    plot(mats, laster, ylim = c(miny, maxy), axes = FALSE, col = "white", 
         type = "l", xlab = xlab, ylab = ylab, main = title)
    haircola = addAlpha(haircol, 0.1)
    recentcola = addAlpha(recentcol, 0.1)
    apply(inmat, 1, function(x) lines(mats, x, col = haircola))
    apply(recent, 1, function(x) lines(mats, x, col = recentcola))
    axis(1)
    axis(2)
    lines(mats, laster, lwd = 2, col = nowcol)
    if(!is.na(yestcol)) {
        lines(mats, last2, lwd = 2, col = yestcol)
    }
    if(withlegend) {
        if(!is.na(yestcol)) {
            legend("bottomright", fill = c(haircol, recentcol, yestcol, yestcol, nowcol), 
                   legend = c(paste("last", round(nrow(inmat)/262), "years"), paste("last", recentamount, "days"), 
                              "last bday", "now"))
        } else {
            legend("bottomright", fill = c(haircol, recentcol, nowcol), 
                   legend = c(paste("last", round(nrow(inmat)/262), "years"), paste("last", recentamount, "days"), "now"))
        }
    }
    if(withchange) {
        barplot(as.numeric(laster - last2), main = "bps change since yesterday", beside = TRUE, 
                names.arg = mats)
    }

}


colSplom <- function(df, hiDays = c(21, 65), quantiles = FALSE, greyalpha = 0.5, tranpose = F, ...) {
# colour scatter plot matrix
    if(!("data.frame" %in% class(df))) df <- as.data.frame(df)
    n <- nrow(df)
    mygrey <- addAlpha("grey40", greyalpha)
    cl <- rep(mygrey, n)
    cl[(n - hiDays[2]):n] <- "#FF555555"
    cl[(n - hiDays[1]):n] <- "#5555FF55"
    cl[n] <- "red"
    cx <- rep(0.5, n)
    cx[n] <- 1
    splom(df, axis.text.cex = 0.4, panel = function(x, y, ...) {
        panel.xyplot(x, y, ..., pch = 16, col = cl, cex = cx)
        reg <- lm(y ~ x)
        reg21 <- lm(y[(n - hiDays[2]):n] ~ x[(n - hiDays[2]):n]) #21 day regression
        regsd <- round(last(reg$residuals) / sd(reg$residuals), 1)
        regsd21 <- round(last(reg21$residuals) / sd(reg21$residuals), 1)
        panel.abline(reg)
        panel.abline(reg21, col =  "blue", lty = "dashed")
        if(quantiles) {
            panel.abline(rq(y ~ x, tau = 0.5), col = addAlpha("red", 0.5))
            panel.abline(rq(y ~ x, tau = 0.25), col = addAlpha("red", 0.5))
            panel.abline(rq(y ~ x, tau = 0.75), col = addAlpha("red", 0.5))
            panel.abline(rq(y ~ x, tau = 0.05), col = addAlpha("red", 0.5))
            panel.abline(rq(y ~ x, tau = 0.95), col = addAlpha("red", 0.5))
        }
        panel.text(max(x), min(y), paste("se l/s", regsd, ", ", regsd21, sep = ""), 
                   cex = 0.6, adj = c(1, -2.5), col = "purple")
        panel.text(max(x), min(y), paste("rsq ", format(summary(reg)$r.squared, digits = 1, nsmall = 3)), 
                   cex = 0.6, adj = c(1, 0.5))
        panel.text(max(x), min(y), paste("rsqshort ", format(summary(reg21)$r.squared, digits = 1, nsmall = 3)), 
                   cex = 0.6, adj = c(1, -1), col = "blue")
    }, ...)
}

gheatmap <- function(inmatrix, titlemain = NA) {
# ggplot-based heatmap

}


logret <- function(fxdata = fx, periods = 1, sparse = FALSE) {
#generates a period log return. If the period > 1 then spare specifies rolling or discrete series
    rets <- diff(log(fxdata), periods)[-1:-periods]
    if (sparse) {
        if(is.null(dim(fxdata))) {
            rets <- rets[rev(rev(index(rets))[rep(1:periods, length(fxdata) %/% periods) == 1])]
        } else {
            rets <- rets[rev(rev(index(rets))[rep(1:periods, nrow(fxdata) %/% periods) == 1]), ]
        }
    }
    return(rets)
}


diffret <- function(fxdata = fx, periods = 1, sparse = FALSE) {
#generates a period differnce return
    rets <- diff(fxdata, periods)[-1:-periods]
    if (sparse) {
        if(is.null(dim(fxdata))) {
            rets <- rets[rev(rev(index(rets))[rep(1:periods, length(rets) %/% periods) == 1])]
        } else {
            rets <- rets[rev(rev(index(rets))[rep(1:periods, nrow(rets) %/% periods) == 1]), ]
        }
    }
    return(rets)
}

difret <- function(fxdata = fx, periods = 1, sparse = FALSE) {
#generates a period differnce return; this one does it !!! CORRECTLY !!! as does not remove periods
    rets <- diff(fxdata, periods)
    if (sparse) rets <- rets[rev(rev(index(rets))[rep(1:periods, nrow(rets) %/% periods) == 1]), ]
    return(rets)
}


genseries <- function(returns, logrets = TRUE, starti = 1, addFirst = FALSE) {
#generates a total return series based on i, addfirst says do I add back a first date    
    i <- as.numeric(starti)[1]
    if(logrets) (for (x in returns) i <- c(i, i[length(i)] * as.numeric(exp(x)))) 
        else i <- cumsum(c(starti, returns))
    if ("xts" %in% class(returns)) { # then we must add back a first date
        if(addFirst) {
            firstDate <- as.Date(index(returns))[1] - 1
            while(weekdays(firstDate) %in% c("Saturday", "Sunday")) firstDate <- firstDate - 1
            dates <- c(firstDate, as.Date(index(returns)))
        } else {
            dates <- as.Date(index(returns))
            i <- i[-1]
        }
        return(xts(i, order.by = dates))
    } else return(i)
}


addAlpha <- function(colour, alpha) {
    cl <- col2rgb(colour)
    return(rgb(cl[[1]]/255, cl[[2]]/255, cl[[3]]/255, alpha))
}


wt.sd <- function (x, wt) {
    return(sqrt(wt.var(x, wt)))
}

wt.var <- function (x, wt) {
    s = which(is.finite(x + wt))
    wt = wt[s]
    x = x[s]
    xbar = wt.mean(x, wt)
    return(sum(wt * (x - xbar)^2) * (sum(wt)/(sum(wt)^2 - sum(wt^2))))
}

wt.mean <- function(x, wt) {
    s = which(is.finite(x * wt))
    wt = wt[s]
    x = x[s]
    return(sum(wt * x)/sum(wt))
}


decay <- function(len, halflife, sumone = TRUE, meanone = FALSE) {
#function generates an exponentially decaying series
    t <- len:1 # generate a series of numbers reverse order so biggest weights last
    lambda <- log(2) / halflife #figure out the lambda for the halflife
    w <- exp(-lambda * t) #create the weights series  
    if(sumone) w <- w / sum(w) #normalise sum to 1 if necessary
    if(meanone) w <- w / mean(w)
    return(w) 
}

howtrading <- function(avector, zwindow, zhalflife) {
# returns the genseries of a decay weighted vector of standard deviations
    sdvector <- avector/sd(avector) # get the sds
    decayer <- decay(zwindow, zhalflife)
    decayer <- decayer / last(decayer) # last one = 1
    sddecayed <- last(sdvector, zwindow) * decayer
    cumulated <- cumsum(sddecayed)
    return(cumulated)
}


cls <- function(sleep = 2) {
# clear all graphics
    graphics.off()
}

splitVector <- function(thevector, maxelements) {
# splits the vector into a list of vectors each with maxelements number of elements in it
    return(split(thevector, ceiling(seq_along(thevector) / maxelements)))
}

utctime <- function() {
    ft <- format(Sys.time(), "%Y-%m-%d %Hh%M %Z", tz = "UTC", usetz = TRUE)
    return(ft)
}

cettime <- function() {
    ft <- format(Sys.time(), "%Y-%m-%d %Hh%M", tz = "CET", usetz = TRUE)
    return(ft)
}


systime <- function() {
    ft <- format(Sys.time(), "%Y-%m-%d %Hh%M", usetx = FALSE)
    return(ft)
}


systime <- function() {
    ft <- format(Sys.time(), "%Y-%m-%d %Hh%M")
    return(ft)
}

ggstamp <- function(whatprint = cettime(), x = 0.01, y = 0.02, h = 0, 
                    font = 3, cex = 0.8, col = "grey") {
    xx <- viewport()
    pushViewport(xx)
    grid.text(whatprint, x = unit(x, "npc"), y = unit(y, "npc"), 
              hjust = h, gp = gpar(col = col, font = font, cex = cex))
    popViewport()
}


bestglmhedge <- function(cvglmobj, hedges, target, rs2tgt = 0.95) {
    for(l in cvglmobj$lambda) {
        cf <- coef(cvglmobj, s = l)[-1]
        if(!(all(cf == 0))) {
            sim <- hedges %*% cf
            linmod <- lm(sim ~ target)
            rs2 <- summary(linmod)$r.squared 
            if(rs2 > rs2tgt) {
                return(list(lambda = l, rs2 = rs2, coef = coef(cvglmobj, s = l)[-1]))
            }
        }
    }
    return(list(lambda = last(cvglmobj$lambda), rs2 = rs2, 
                coef = coef(cvglmobj, s = last(cvglmobj$lambda))[-1]))
}


hedger <- function(hedges, target, numfit = NA, rsqtarget = 0.925) {
# glmnet best hedging
    cnames <- colnames(hedges)
    minlength <- min(nrow(hedges), length(target))
    target <- last(target, minlength)
    hedges <- tail(hedges, minlength)
    target <- as.numeric(target)
    hedges <- as.matrix(hedges)
    cvfit <- cv.glmnet(hedges, target)
    se1 <- cvfit$lambda.1se
    bestfit <- coef(cvfit, s = se1)[-1]
    sparsefit <- bestglmhedge(cvfit, hedges, target, rs2tgt = rsqtarget)$coef
    if(all(sparsefit == 0)) {
        tryer <- length(cvfit$lambda)
        sparsefit = coef(cvfit, s = cvfit$lambda[tryer])[-1]
        while(sum(sparsefit != 0) > 4) {
            tryer <- tryer - 1
            sparsefit = coef(cvfit, s = cvfit$lambda[tryer])[-1]
        }
    }
    if(!(is.na(numfit))) {
        lambda <- cvfit$lambda[max(which(cvfit$nzero == numfit))]
        while(is.na(lambda)) { 
            numfit = numfit + 1
            lambda <- cvfit$lambda[max(which(cvfit$nzero == numfit))]
        }
        sparsefit <- coef(cvfit, s = lambda)[-1]
    }
    return(list(bestfit = bestfit, sparsefit = sparsefit))
}


hothot <- function(m, title = "Credit alpha matrix", xlab = NA, ylab = NA, withvalues = T) {
# do a ggplot heatmap on matrix m
    gradientscale <- scales::rescale(c(-5.5, -4, -1.9, -1, 0, 1, 2.1, 4, 5.5))
    gradientcolours <- c("mediumvioletred", "darkred", "red", "orange", "yellow", 
                         "#99BBFF", "blue", "darkblue", "mediumspringgreen")
    d <- melt(m)
    for(x in 1:nrow(d)) if(d[x, 1] == d[x, 2]) d[x, 3] <- 0
    d <- as.data.frame(d)
    g <- ggplot(d, aes(Var1, Var2))
    g <- g + geom_tile(aes(fill = value))
    if(withvalues) 
        g <- g + geom_text(aes(label = round(value, 1)), size = 3)
    g <- g + scale_fill_gradientn(colours = gradientcolours, values = gradientscale, limits = c(-5.5, 5.5))
    g <- g + theme(axis.text.x = element_text(angle = -90, vjust = 0.30))
    if(is.na(xlab)) {
        g <- g + theme(axis.title.x = element_blank())
    } else {
        g <- g + xlab(xlab)
        g <- g + theme(axis.title.x = element_text(xlab))
    }
    if(is.na(ylab)) {
        g <- g + theme(axis.title.y = element_blank())
    } else {
        g <- g + ylab(ylab) 
    }
    g <- g + ggtitle(title)
    plot(g)
}
                

matrixzs <- function(incd, days = 21, windowl = 260, mindata = 65, decayhl = NA, meanone = TRUE) {
# takes a matrix and gives final row z scores (zscores - added for searchability) for it
    if(is.na(decayhl)) decayer <- (rep(1, nrow(incd)) / ifelse(meanone, nrow(incd), 1))
        else decayer <- decay(nrow(incd), decayhl, meanone = meanone)
    usecd <- tail(incd, windowl + days)
    sds <- rollapply(usecd, windowl, function(x) 
        if(length(na.omit(x)) > mindata) wt.sd(na.omit(x), last(decayer, length(na.omit(x)))) else NA)
    means <- rollapply(usecd, windowl, function(x) 
        if(length(na.omit(x)) > mindata) wt.mean(na.omit(x), last(decayer, length(na.omit(x)))) else NA)
    lasts <- tail(usecd, days + 1)
    zs <- (lasts - means) / sds
    return(list(zs = zs, sds = sds, means = means, lasts = lasts))
}


matrixlm <- function(incd, decayhl = NA, orthogonal = F) {
# takes a matrix and gives the z scores of each column regressed against each other column
    if(is.na(decayhl)) {
        decayer <- rep(1, nrow(incd))
    } else {
        decayer <- decayer(nrow(incd), decayhl)
    }
    if(orthogonal & !is.na(decayhl)) {
        flushprint("Cannot use decay weights on orthogonal regression")
        return(-1)
    }
    apply(incd, 2, function(x) {
        apply(incd, 2, function(y) {
            if(orthogonal) {
                res = orthlm(y ~ x)$residuals
            } else {
                res = lm(y ~ x, weights = decayer)$residuals
            }
            (last(res) - wt.mean(res, decayer)) / wt.sd(res, decayer)
        })
    })
}

orthlm <- function(formula, data = NULL) {
    dta <- model.frame(formula, data = data)
    cv <- cov(dta)
    eigv <- eigen(cv)$vectors
    beta <- eigv[1, 1] / eigv[2, 1]
    intercept <- mean(dta[, 1]) - beta * mean(dta[, 2])
    residuals <- dta[, 1] - (beta * dta[, 2] + intercept)
    r.squared = cor(dta)[2, 1] ^ 2
    returner <- list(coefficients = c(intercept, beta), residuals = residuals, r.squared = r.squared)
    class(returner) <- "orthlm"
    return(returner)
}

summary.orthlm <- function(x) {
    return(x)
}

coef.orthlm <- function(x) {
    return(x$coefficients)
}
    

    








