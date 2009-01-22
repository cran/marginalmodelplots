#mmps.R part of marginalmodelplots
#by Andrew Redd
#Lisenced under GPL version 2 or newer.
library(locfit)
last<-function(x)x[length(x)]
first<-function(x)x[1]
mmplot1<-function(object,...)UseMethod("mmplot1")
mmplot1.lm<-function(object,u=predict(object), 
		 label=deparse(substitute(u)),
         locfit.control=list(nn=2/3),
         colors=c("blue","red"),...)mmplot1.glm(object,u=u,family='qgauss',link='identity',
		 label=label, locfit.control=locfit.control,colors=colors,...)
mmplot1.glm<-function(object,u=predict(object) , mean=TRUE , #sd=FALSE ,
         label=deparse(substitute(u)),
         locfit.control=list(nn=2/3),
		 subset=NULL,
         colors=c("blue","red"),
         family=NULL,link=NULL ,
		 yhat.autoquasi=TRUE ,
		 inc.legend=TRUE,...){
#adapted from mmp.glm from the alr3 package
#locfit.control is either a list of two named components, y and yhat, the components are passed onto the lp function,
# or a list of one component wich will be passed to lp for both y and yhat. A specification of h='aic' or h='gcv' will
# optimize the aic, or gcv criteria for for the fixed bandwidth, via optim.  If optimization fails nearest neighbor 
# will be used.
	if(!require(locfit))stop("locfit package is required")
	if(label=="predict(object)"){
		label <- "Linear Predictor"
	}
	data<-if(is.null(subset)) object$model else subset(object$model,subset=subset)
	y <- response <- data[,1]
	parent<-parent.frame()
	u<-tryCatch(eval(u,parent.frame()),error=function(e)eval(substitute(u,parent),data,parent))
	if(is.factor(response)){
		if(nlevels(response)!=2)stop("The number of levels in a factor must be 2.")
		y<-as.numeric(response)-1
	}
	if (is.matrix(response)) y <- response[,1]/apply(response,1,sum)
	fam <- if(is.null(family)) object$family$family else family
	link <- if(is.null(link))  object$family$link   else link
	yplot <- if((fam=='binomial') && (length(unique(y))==2)) jitter(as.numeric(y),.1) else y
	if(is.factor(response)){
		plot(u,y,xlab=label,ylab=colnames(data[1]),type='n',yaxt='n',...)
		axis(2,0:1,labels=levels(response))
	} else plot(u,y,xlab=label,ylab=colnames(data[1]),type='n',...)
	points(u,yplot)
	if(!all(c('y','yhat')%in%names(locfit.control))){  #handling for single specification for both y and yhat
		lc<-locfit.control
		locfit.control<-list()
		locfit.control$y<-lc
		locfit.control$yhat<-lc
		MatchBWToY = TRUE
	} else {
		MatchBWToY = FALSE
	}
	#y
	hopt<-numeric(0)
	if(is.numeric(locfit.control$y$h) || is.null(locfit.control$y$h)) {
		hopt<-locfit.control$y$h
	} else {
		optim.criteria<-locfit.control$y$h
		lp.par<-locfit.control$y
		gh<-function(h){
			lp.par$h=h
			lp.par$x=u
			lp.u<-do.call(lp,lp.par)
			ll.m<-locfit(y~lp.u,family=fam,link=link)
			do.call(optim.criteria,list(ll.m))[4]
		}
		optimres<-suppressWarnings(optim(diff(range(u)),gh,method='L-BFGS-B',lower=max(diff(sort(unique(u))))))
		if(optimres$convergence!=0) {
		warning(paste("optim did not converge to find the optimal bandwidth for y vs.",label,".\nThe final 		bandwidth selected was h=",optimres$par,"\n"))
		hopt<-optimres$par
		} else {
			hopt<-optimres$par
		}
	}
	nnopt<-numeric(0)
	
	if(is.numeric(locfit.control$y$nn) || is.null(locfit.control$y$nn)) {
		nnopt<-locfit.control$y$nn
	} else {
		optim.criteria<-locfit.control$y$nn
		lp.par<-locfit.control$y
		gn<-function(nn){
			lp.par$nn=nn
			lp.par$x=u
			lp.u<-do.call(lp,lp.par)
			ll.m<-locfit(y~lp.u,family=fam,link=link)
			do.call(optim.criteria,list(ll.m))[4]
		}
		optimres<-suppressWarnings(optim(.7 ,gn, method='L-BFGS-B',lower=0,upper=1))
		if(optimres$convergence!=0) {
			warning(paste("optim did not converge to find the optimal nearest neighbor for y vs.",label,".\nThe final nearest neigbor selected was nn=",optimres$par,"\n"))
			nnopt<-optimres$par
		} else {
			if(optimres$par==.7) warning(paste("optim did not move from initial position with nearest neighbor for y vs. ", label, ", using default of .7"))
			nnopt<-optimres$par
		}
	} 
	lp.par<-locfit.control$y
	if(!is.null(lp.par$h))lp.par$h<-hopt
	if(!is.null(lp.par$nn))lp.par$nn<-nnopt
	lp.par$x<-u
	lp.u<-do.call(lp,lp.par)
	locfit.y <- suppressWarnings(locfit(y ~ lp.u,family=fam,link=link))
	rtny<-list(h=lp.par$h,nn=lp.par$nn)#portion of return value
	#yhat
	yhat<-predict(object,type="response")
	yhat.fam<-if(yhat.autoquasi) paste(if(substr(fam,1,1)=="q") {""} else "q",fam,sep="") else fam
	hopt<-numeric(0)
	if(!MatchBWToY) {
		if(is.numeric(locfit.control$yhat$h) || is.null(locfit.control$yhat$h)) {
			hopt<-locfit.control$yhat$h
		} else {
			optim.criteria<-locfit.control$yhat$h
			lp.par<-locfit.control$yhat
			gh2<-function(h){
				lp.par$h=h
				lp.par$x=u
				lp.u<-do.call(lp,lp.par)
				ll.m<-locfit(yhat~lp.u,family=yhat.fam,link=link)
				do.call(optim.criteria,list(ll.m))[4]
			}
			optimres<-suppressWarnings(optim(diff(range(u)),gh2,method='L-BFGS-B',lower=max(diff(sort(unique(u))))))
			if(optimres$convergence!=0) {
				warning(paste("optim did not converge to find the optimal bandwidth for yhat vs.",label,".\nThe final bandwidth selected was h=",optimres$par,"\n"))
				hopt<-optimres$par
			} else {
				hopt<-optimres$par
			}
		} 
		nnopt<-numeric(0)
		if(is.numeric(locfit.control$yhat$nn) || is.null(locfit.control$yhat$nn)) {
			nnopt<-locfit.control$yhat$nn
		} else {
			optim.criteria<-locfit.control$yhat$nn
			lp.par<-locfit.control$yhat
			gn2<-function(nn){
				lp.par$nn=nn
				lp.par$x=u
				lp.u<-do.call(lp,lp.par)
				ll.m<-locfit(yhat~lp.u,family=fam,link=link)
				do.call(optim.criteria,list(ll.m))[4]
			}
			optimres<-suppressWarnings(optim(.7,gn2,method='L-BFGS-B',lower=0,upper=1))
			if(optimres$convergence!=0) {
				warning(paste("optim did not converge to find the optimal nearest neighbor for yhat vs.",label,".\nThe final nearest neigbor selected was nn=",optimres$par,"\n"))
				nnopt<-optimres$par
			} else {
				if(optimres$par==.7) warning(paste("optim did not move from initial position with nearest neighbor for yhat vs. ", label, ", using default of .7"))
				nnopt<-optimres$par
			}
		}
	}
	lp.par<-locfit.control$yhat
	if(!is.null(lp.par$h))lp.par$h<-hopt
	if(!is.null(lp.par$nn))lp.par$nn<-nnopt
	lp.par$x<-u
	lp.u<-do.call(lp,lp.par)
	locfit.yhat <-suppressWarnings(locfit(yhat~lp.u,family=fam,link=link))
	if (fam=='binomial' && length(unique(y))==2){
		sel<-round(yhat)!=y
		points(u[sel],yplot[sel],col='red')
	}
	new <- seq(min(u),max(u),length=200)
	pred.locfit.y<-predict(locfit.y,data.frame(lp.u=new))
	pred.locfit.yhat<-predict(locfit.yhat,data.frame(lp.u=new))
	if(mean==TRUE) {
		lines(new,pred.locfit.y,   lty=1,col=colors[1])
		lines(new,pred.locfit.yhat,lty=2,col=colors[2])}
	rtnyhat=list(h=lp.par$h,nn=lp.par$nn)
	if(inc.legend){
		legend.loc<-if(hasArg(legend.loc)) list(...)$legend.loc else xy.coords(x=first(axTicks(1)),y=last(axTicks(2)))
		legend(legend.loc, c("Raw","Fitted"),col=colors, lty=c(1,2))
	}
	invisible(list(y=rtny,yhat=rtnyhat))
}
mmplot<-function(object,exclude=NULL,layout=NULL,ask,...){
	predvars <- attr(object$terms,"predvars") #list of all variables in the model including the response
	nt <- length(predvars)-1 #number of terms excluding response
	factors <- which(sapply(2:dim(object$model)[2],function(j) 
	is.factor(object$model[[j]]))) # gives indices of the factors in the model
	if (length(factors)>0)warning("Factors were skipped") # warning if there are factors
	exclude <- c(exclude,factors) #Exclude is passed in and leave out factors
	if(is.null(layout)){
		layout <- switch(min(nt+1-length(exclude),9),
			c(1,1),c(1,2),c(2,2),c(2,2),c(3,2),c(3,2),
			c(3,3),c(3,3),c(3,3))
	}   # defines the layout if it was not defined
	op<-par(no.readonly=TRUE)
	ask <- if(missing(ask) || is.null(ask)) prod(layout)<nt-length(exclude) else ask  #sets up a default for the ask parameter
	on.exit(par(op))
	if(prod(layout) > 1)
		par(mfrow=layout,mai=c(.6,.6,.1,.1),mgp=c(2,1,0),cex.lab=1.0,cex=0.7,ask=ask) 
	else par(mfrow=layout,ask=ask) #sets the parameters
	rtn=list()
	for (j in 3:(nt+1)){
		if(is.na(match(j-2,exclude))){ 
			rtn=append(rtn, list(mmplot1(object,object$model[,j-1],label=deparse(predvars[[j]]),...)))
		}
	} #Cycle through the terms in the model skipping over the excluded ones
	rtn=append(rtn, list(mmplot1(object,,...)))
	invisible(rtn)
}
