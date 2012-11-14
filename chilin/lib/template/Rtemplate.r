%# value, name, pdfName, histroyData, main, xlab, ylab


col=c('#FFB5C5','#5CACEE','#7CFC00','#FFD700','#8B475D','#8E388E','#FF6347','#FF83FA','#EEB422','#CD7054')
pch=c(21,22,24,25,21,22,24,25,21,22,24,25,21,22,24,25)

value = c(\VAR{value})
name = c(\VAR{name})

pdf('\VAR{pdfName}',height=8.5,width=8.5)
rawdata<-c(\VAR{histroyData})
fn<-ecdf(rawdata)

density <- fn(rawdata)
fndd <- cbind(rawdata,density)
fndd2<-fndd[order(fndd[,1]),]
  
%#plot(fndd2[,1],smooth(fndd2[,2]),type='l',pch=18,col='blue',main='\VAR{main}',xlab='\VAR{xlab}',ylab='\VAR{ylab}')     


%#plot(ecdf(rawdata), verticals=TRUE,col.hor='blue', col.vert='black',main='\VAR{main}',xlab='\VAR{xlab}',ylab='\VAR{ylab}')

%- if fastqcCheck
plot(fndd2[,1],smooth(fndd2[,2]),type='l',pch=18,col='blue',main='\VAR{main}',xlab='\VAR{xlab}',ylab='\VAR{ylab}')
xx=seq(0,\VAR{cutoff})    
yy = c(0,smooth(fn(xx)),0)
xx =c(0,xx,\VAR{cutoff})
polygon(xx,yy, col = 'lightpink')
        
mid = quantile(rawdata,0.5)
        
xx=seq(\VAR{cutoff},mid)   
yy = c(0,smooth(fn(xx)),0)
xx =c(\VAR{cutoff},xx,mid)
polygon(xx,yy, col = 'lightgoldenrod1')
         
xx=seq(mid,max(rawdata))
yy = c(0,smooth(fn(xx)),0)
xx =c(mid,xx,max(rawdata))  
polygon(xx,yy, col = 'palegreen')

%- endif

%- if other

plot(ecdf(rawdata),verticals=TRUE,col.hor='blue',pch='.',col.vert='black',main='\VAR{main}',xlab='\VAR{xlab}',ylab='\VAR{ylab}')
%#plot(fndd2[,1],smooth(fndd2[,2]),type='l',pch=18,col='blue',main='\VAR{main}',xlab='\VAR{xlab}',ylab='\VAR{ylab}')
xx=seq(0,\VAR{cutoff},length = 1000)    
yy = c(0,fn(xx),0)
xx =c(0,xx,\VAR{cutoff})
polygon(xx,yy, col = 'lightpink')
        
mid = quantile(rawdata,0.5)
        
xx=seq(\VAR{cutoff},mid,length = 1000)   
yy = c(0,fn(xx),0)
xx =c(\VAR{cutoff},xx,mid)
polygon(xx,yy, col = 'lightgoldenrod1')
         
xx=seq(mid,max(rawdata),length = 1000)
yy = c(0,fn(xx),0)
xx =c(mid,xx,max(rawdata))  
polygon(xx,yy, col = 'palegreen')
%- endif
len = length(value)

for (i in seq(len)){
	points(value[i],jitter(fn(value[i])),pch=pch[i],bg=col[i])
	}
legend('topleft',name,pch=pch[1:len],pt.bg=pch[1:len])

dev.off()




