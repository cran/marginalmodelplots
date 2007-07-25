#aicc.R
#Part of the mmps package for R
#Lisenced under the GPL version 2 or newer
mmpslider<-function(m,...)UseMethod("mmpslider")
mmpslider.lm<-function(m,pred=predict(m), bw=NULL , label=deparse(substitute(u)), colors=c('blue','red'),...)mmpslider.glm(m,pred=pred,bw=bw,label=label,colors=colors,family='qgauss',link='identity',...)
mmpslider.glm<-function(m,pred=predict(m), bw=NULL , label=deparse(substitute(u)), colors=c('blue','red'), 
family=NULL , link=NULL , ...){
#lp.par should only include extra parameters such as deg
    #validity checks
    #if(any(class(m)=="glm")) NULL else stop("m must be an object of class glm.")
    if(is.null(pred)) stop("Please provide a predictor.")
    if(!require(grid)) stop("grid package is required.  Why don't you have grid???")
    if(!require(locfit)) stop("locfit package is required")
	if((label=="predict(m)")||(label=="u")){ label <- "Linear Predictor"}

	
    #extract response and predicted values
    y <- response <- m$model[,1]
	if(is.factor(response)){
		if(nlevels(response)!=2)stop("The number of levels in a factor must be 2.")
		y<-as.numeric(response)-1
	}
	if(NCOL(y)==2){ y<-y[,1]/(y[,1]+y[,2])}
    yhat<- predict(m,type="response")
	fam <- if(is.null(family)) m$family$family else family
	link <- if(is.null(link))  m$family$link   else link
    yname <- (names(m$model)[1])

    lp.par=list(x=pred)
    #find paramteters for bw
    if(is.null(bw)){
      g<-function(h){
        lp.par$h=h
        lp.u<-do.call(lp,lp.par)
	  ll.m<-locfit(yhat~lp.u,family=fam,link=link)
        aic(ll.m)[4]
      }
      bwmin<-max(diff(sort(unique(pred))))
      optimres<-suppressWarnings(optim(diff(range(pred)),g,method='L-BFGS-B',lower=bwmin))
      if(optimres$convergence!=0) warning(paste("optim did not converge to find the optimal bandwidth for y vs.",label,".\nThe final bandwidth selected was bw=",optimres$par,"\n"))
      bwopt<-optimres$par
	bwmax<-2*bwopt
      bw<-c(bwmin,bwopt,bwmax)
    } else if((length(bw)!=3)&(length(bw)!=2))stop("if specified bw must have a length of 2 or 3 and be of either forms (min, max) or (min, starting, max)")
    if(length(bw)==2){bw<-c(bw[1],mean(bw),bw[2])}
    bwcur<-bw[2]

    # setup plotting area and grid viewports
    grid.newpage()
    lay<-grid.layout(nrow=2,ncol=1,heights=unit(c(1,1.25),c('null','cm')))
    vpmain<-viewport(layout=lay,name='main')
    vpplotarea<-viewport(layout.pos.row=1,name='plotarea')
    vpbw<-viewport(layout.pos.row=2,name='bwarea')
    pushViewport(vpTree(vpmain,vpList(vpplotarea,vpbw)))
	yplot <- if((fam=='binomial') && (length(unique(y))==2)) jitter(as.numeric(y),.1) else y
    seekViewport('plotarea')
    xscale<-1.1*(range(pred)-mean(range(pred)))+mean(range(pred))
    yscale<-1.1*(range(y)-mean(range(y)))+mean(range(y))
    mar<-c(5.1,4.1,4.1,2.1)
    vp.x<-unit(mar[2],'lines')
    vp.y<-unit(mar[1],'lines')
    vp.width<-unit(1,'npc')-unit(mar[2]+mar[4],'lines')
    vp.height<-unit(1,'npc')-unit(mar[1]+mar[3],'lines')
    pushViewport(viewport(x=vp.x,y=vp.y,width=vp.width,height=vp.height,
      default.units='native',just=c('left','bottom'),xscale=xscale,yscale=yscale,name='plot'))
    #add axes
    grid.rect()
    grid.xaxis(at=seq(min(pred),max(pred),length=5))
    if(fam=='binomial' && length(unique(y))==2) { if(is.factor(response)) grid.yaxis(at=0:1,label=levels(response)) else grid.yaxis(at=unique(y))} else grid.yaxis()
    #add labels
    grid.text(label,y=unit(-2.5,'lines'))
    grid.text(yname,x=unit(-2.5,'lines'),rot=90)
    #add points
    points<-grid.points(name='points',pred,yplot,pch=1,size=unit(3,'mm'))
    if(fam=='binomial' && length(unique(y))==2){
      sel<-y!=round(yhat)
      mispoints<-grid.points(name='misclassified points',pred[sel],yplot[sel],gp=gpar(col='red'),pch=1,size=unit(3,'mm'))
    }
    #shrink bw sel area
    seekViewport('bwarea')
    pushViewport(vpbwsel<-viewport(width=unit(.8,'npc'),height=unit(.8,'npc'),xscale=bw[-2],name='bwselarea'))
    grid.xaxis(main=F)
    grid.lines(x=unit(bw[-2],'native'),y=unit(c(.5,.5),'npc'))


    # loop through changing bw and curves
    xeval<-seq(min(pred),max(pred),length=200)
    res<-TRUE
	yhat.fam<-paste(if(substr(fam,1,1)=="q"){""}else{"q"},fam,sep="")
    while(!is.null(res)){
        h<-unclass(bwcur)
        lp.par$h<-h
        seekViewport('bwselarea')
        bwboxwidth<-grobWidth(textGrob(expression(h[y])))
        bwtext<-textGrob(bquote(h[y] == .(bwcur)),name='bwindtext',x=unit(bwcur,'native')-.5*bwboxwidth,just="left")
        bwrect<-rectGrob(name='bwindrect',x=unit(bwcur,'native'),width=bwboxwidth,gp=gpar(fill='thistle'))#,vp=vpbwsel)
        bw.indicator<-gTree(x=unit(bwcur,'native'),children=gList(bwrect,bwtext),name='bwindicator')
        grid.draw(bw.indicator)
        seekViewport('plot')
        grid.remove('ycurve',global=T,warn=F)
        grid.remove('yhatcurve',global=T,warn=F)
        #result<-sm.binomial(pred,y,h=h,display='none')
        lfx<-do.call('lp',lp.par)
        lfy<-locfit(y~lfx,family=fam,link=link)
        grid.lines(name='ycurve',x=unit(xeval,'native'), y=unit(predict(lfy,data.frame(lfx=xeval)),'native'), gp=gpar(col=colors[1],lty=1))
        lfyhat<-locfit(yhat~lfx,family=yhat.fam,link=link)
        grid.lines(name='yhatcurve',x=unit(xeval,'native'), y=unit(predict(lfyhat,data.frame(lfx=xeval)),'native'),gp=gpar(col=colors[2],lty=2))
        seekViewport('bwselarea')
        res<-grid.locator() 
        #if(unclass(res$y)<1)res<-NULL
        if(!is.null(res)){
            grid.remove('bwind*',grep=T,global=T)    
            bwcur<-res$x
            #grid.edit('bwind*',x=bwcur,grep=T,global=T,redraw=T)
            #grid.draw(bw.indicator)
        }
    }
}
