
logit <- function(x){ log(x/(1-x))}
expit <- function(x){1/(1+exp(-x))}

# expand the nest data for proportional hazard model
# We assume that nest data collected in standard fashion
# We generate a separate record for each day of the study plus the final 
# interval. FInally, we generate the appropriate field for survival analysis
expand.nest.data.ph <- function( nestdata, 
                              FirstFound="FirstFound", 
                              LastChecked="LastChecked",
                              LastPresent="LastPresent",
                              Fate="Fate",
                              AgeDay1="AgeDay1"){
  require(plyr)
  require(survival)
  nestdata2 <- plyr::adply(nestdata, 1, function(x){
    # expand for days when nest is known to be alive
    # create data frame because tibbles can cause problems
    x <- as.data.frame(x)
    #browser()
    orig.x <- x
    x <- x[rep(1,x[,LastPresent]-x[,FirstFound]),]  
    #browser()
    if(nrow(x)>0){
       x$Start <- x[,FirstFound][1]:(x[,LastPresent][1]-1)
       x$End   <- x$Start + 1
       x$Fail  <- 0
       #browser()
       if(AgeDay1 %in% names(orig.x)) x$NestAge <- orig.x[,AgeDay1] + x$Start -1 
    }
    
    # now for the final record where the nest fails in interval
    #browser()
    if(orig.x[,LastChecked] > orig.x[,LastPresent]){
       x2 <- orig.x
       x2$Fail  <- orig.x[,Fate]
       x2$Start <- x2[,LastPresent]
       x2$End   <- x2[,LastChecked]
       if(AgeDay1 %in% names(orig.x)) x2$NestAge  <- orig.x[,AgeDay1] + x2$Start -1
       x <- rbind(x, x2)
       #browser()
    }
    # remove any records where Start=End
    x <- x[ x$Start!=x$End,]
    # return the expanded data set
    # add the Surv
    x$Surv <- Surv(time=x$Start, time2=x$End, event=x$Fail, type="counting")
    x
  })
  nestdata2}

  
