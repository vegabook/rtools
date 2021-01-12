#THIS IS SET OF ROUTINES COMMONLY USED BY THOMAS BROWNE ANALYSES USING BLOOMBERG

library(reshape2)
library(xts)


gobb <- function() {
# connect to the terminal
    conn <<- blpConnect()
}

bdp <- function(securities, fields) {
# wrapper for bdp function that makes sure fields are valid
    fields <- toupper(gsub(" ", "_", fields))
    Rblpapi::bdp(securities, fields)
}

nafill <- function(x) {
    if("xts" %in% class(x)) {
        x0 = rev(na.locf(rev(x)))
    } else {
        x0 <- na.locf(x[nrow(x):1, ], na.rm = FALSE)
        x0 <- na.locf(x0[nrow(x):1, ], na.rm = FALSE) # twice
    }
    return(x0)
}


bbdh <- function(secs, years = 1, flds = "last_price", startDate = NULL, endDate = NULL, asDateNotPosix = FALSE, checkgoodbdp = FALSE, 
                 datasource = NULL, simplify = T) {
#this function gets secs over years from bloomberg daily data
    if(!is.null(datasource)) {
        getsecs <- sapply(secs, function(x) {
            splitsecs <- strsplit(x, " ")[[1]]
            splitsecs[1] <- paste(splitsecs[1], "@", datasource, sep = "")
            do.call(paste, as.list(splitsecs))
        })
    } else {
        getsecs <- secs
    }
    if(is.null(startDate)) startDate <- Sys.Date() - years * 365.25
    startDate <- as.Date(startDate)
    if(is.null(endDate)) endDate <- Sys.Date()
    flds <- toupper(gsub(" ", "_", flds))
    if(checkgoodbdp) {
        isgood <- sapply(getsecs, function(x) {xx = try(bdp(x, "LAST_PRICE")); !("try-error" %in% class(xx))})
        getsecs <- getsecs[isgood]
        secs <- secs[isgood]
    }
    options= c("nonTradingDayFillOptio" = "NON_TRADING_WEEKDAYS", "nonTradingDayFillMethod" = "PREVIOUS_VALUE")
    rawd <- bdh(getsecs, flds, startDate, endDate, include.non.trading.days = TRUE, verbose = FALSE)
                #options = options)
    if(class(rawd) == "data.frame") {
        rawd <- list(rawd)
        names(rawd) <- secs
    }
    rawd[sapply(rawd, nrow) == 0] <- NULL # erase zero responses
    if(length(rawd) == 0) {
        flushprint("no data")
        return(-2)
    }
    lengths <- sapply(rawd, nrow)
    if(max(lengths) != min(lengths)) {
        flushprint("different lenghts for securities")
        return(-1)
    }
    combined <- do.call(cbind, lapply(rawd, function(x) x[, 2]))
    # now re-order because new bdh can return in different order
    dates <- as.Date(as.character(rawd[[1]][, 1]))
    if(asDateNotPosix) {
        combined <- xts(combined, order.by = as.Date(dates))
    } else {
        combined <- xts(combined, order.by = as.POSIXct(as.character(dates)))
    }
    combined <- combined[, secs[secs %in% colnames(combined)]] # reorder
    nonWeekends <- !(weekdays(dates) %in% c("Saturday", "Sunday"))
    dates <- dates[nonWeekends]
    combined <- combined[nonWeekends, ]
    if(simplify) {
        colnames(combined) <- sub(" .*", "", colnames(combined)) #remove the govt, currncy bits from bb tickers
    }
    return(combined)
}


bbdt <- function(secs, days = 10, flds = "TRADE", starttime = NULL,
                 minutefreq = 1, ohlc = FALSE, fillminutes = TRUE) {
#this function gets tick data from Bloomberg
    if(is.null(starttime)) {
         starttime <- Sys.time() - days * 24 * 60 * 60
    } else {
        starttime <- strftime(starttime, format = "%Y-%m-%d %H:%M:%S.0") #needs .0
        starttime <- as.POSIXct(startime)
    }
    # create a rounded minute endtime
    endtime <- Sys.time()
    endtime <- as.POSIXlt(endtime)
    endtime$sec <- 0
    endtime <- as.POSIXct(endtime)
    # create a list of minutes
    timediff <- round(as.numeric(difftime(as.POSIXct(endtime), as.POSIXct(starttime), unit = "mins")))
    timegaps <- seq(timediff, 0, by = -minutefreq) * 60
    timelist <- endtime - timegaps
    natimes <- xts(rep(NA, length(timelist)), order.by = timelist)
    indexTZ(natimes) <- "UTC"
    alldata <- lapply(secs, function(x) {
        #get the raw data
        if("try-error" %in% class(try(rawt <- getBars(x, 
                                                      eventType = flds, 
                                                      startTime = starttime, 
                                                      endTime = endtime, 
                                                      barInterval = minutefreq)))) {
            rawt <- natimes
        } else {
            if(nrow(rawt) == 0) {
                flushprint(paste("no bar data for", secs, flds))
                rawt <- natimes
            }
            rawt[, 1] <- as.POSIXct(sub("T", " ", rawt[, 1])) #convert the first column from string to POSIXcta
            rawt <- xts(rawt[, -1], order.by = rawt[, 1]) #remove the time column and convert to xts
            indexTZ(rawt) <- "UTC"
            if(!ohlc) {
                rawt <- rawt[, 4] # can the ohlc
            }
        }
        return(rawt)
    })
    columns <- sub(" .*", "", secs)
    if(!ohlc) {
        if(fillminutes) {
            # put in a dummy variable with all the minutes so that cbind will create
            # a matrix with all possible minutes in it
            alldata[["dummy"]] <- natimes
        }
        alldata <- do.call(cbind, alldata)
        if(fillminutes) {
            alldata <- alldata[, -match("dummy", colnames(alldata))]
        }
        colnames(alldata) <- columns
    } else {
        names(alldata) <- columns
    }
    return(alldata)
}


flushcat <- function(...) {
    cat(...)
    flush.console()
}

flushprint <- function(...) {
    print(...)
    flush.console()
}

#function to add transparency to an R colour - takes param colour as a name and adds alpha [0, 1] as transparency
addAlpha <- function(colour, alpha) {
    cl <- col2rgb(colour)
    return(rgb(cl[[1]]/255, cl[[2]]/255, cl[[3]]/255, alpha))
}

consoleColour <- function(ccolour = "sienna1", points = 10) {
    temp <- tempfile()
    cat(paste("points = ", round(points), "\n", sep = ""), file = temp)
    cat(paste("normaltext = ", ccolour, "\n", sep = ""), file = temp, append = TRUE)

    loadRconsole(file = temp)
    
    #temp <- tempfile()
    #cat(paste("normaltext =", ccolour), file = temp, append = TRUE)
    #loadRconsole(file = temp)
}

wideScreen <- function(howWide=as.numeric(strsplit(system('stty size', intern=T), ' ')[[1]])[2]) {
       options(width=as.integer(howWide))
}

xts2df <- function(xtsobj) {
    if("data.frame" %in% class(xtsobj)) {
        return(xtsobj)
    }
    xtsdf <- as.data.frame(xtsobj)
	if(("POSIXct" %in% class(index(xtsobj))) | ("POSIXlt" %in% class(index(xtsobj)))) {
    	df <- data.frame("date" = as.POSIXct(rownames(xtsdf)))
	} else {
    	df <- data.frame("date" = as.Date(rownames(xtsdf)))
	}
    return(cbind(df, xtsdf))
}

rxts2df <- function(data) {
# convert any xts to df in any nested list
	if("list" %in% class(data)) {
		print(names(data))
		flush.console()
		lapply(data, rxts2df)
	} else {
		if("xts" %in% class(data)) {
			xts2df(data)
		} else {
			data
		}
	}
} 

