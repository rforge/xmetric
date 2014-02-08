## Kettelle.R
## This function performs the Generalized Kettelle Algorithm as presented in 
## Chapter 7 of "Statistical Theory of Reliability and Life Testing, Probability Models"
##  by Richard E. Barlow and Frank Proschan.
## 
## (C) David J. Silkworth 2014
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation; either version 2, or (at your option) any
## later version.
##
## These functions are distributed in the hope that they will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the
## GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, a copy is available at
##  http://www.r-project.org/Licenses/
##

Kettelle<-function(x, limit=1e-4, data.name="", performance="FRN", show=FALSE)   {								
	numParts<-length(x[,1])							
								
	UdomAll<-matrix(data=rep(0,numParts+2),							
		nrow=1, ncol=numParts+2)						
								
		UdomAll<-data.frame(UdomAll)						
		PartNames<-NULL						
	## build lists of individual part Perf and Cost values per Qty "s"							
	## this can be referenced rather than calculated each time							
		PerfList<-list()						
		CostList<-list()						
	for (n in 1:numParts)  {							
			s=0					
			unconverged=TRUE					
			thisPerf<-NULL					
			thisCost<-NULL					
		while(unconverged)  {						
			newPerf<-FRN(s,x$Lam[n],x$TAT[n])					
			thisPerf<-c(thisPerf,newPerf)					
			thisCost<-c(thisCost,s*x$C[n])					
								
			if(s>2)  {					
				unconverged<-abs(newPerf-lastPerf)>limit				
			}else{					
				unconverged=TRUE				
			}					
			lastPerf<-newPerf					
			s=s+1					
		}						
		PerfList[[length(PerfList)+1]]<-thisPerf						
		CostList[[length(CostList)+1]]<-thisCost						
		PartNames<-c(PartNames,as.character(x[n,1]))						
	}							
	names(UdomAll)<-c(PartNames,"Perf","Cost")							
	ExAll<-UdomAll							
	thisRow<-UdomAll							
								
	j=1							
	## Fill the UdomAll for the first part with no others							
	for(k in 2:length(CostList[[j]])-1)   {							
		thisRow[j,1]=k						
		for(f in 2:numParts)  {						
			thisRow[j,f]=0					
		}						
		thisRow$Perf<-PerfList[[j]][k+1]						
		thisRow$Cost<-CostList[[j]][k+1]						
		UdomAll<-rbind(UdomAll,thisRow)						
	}							
	UdomLast<-UdomAll							
								
	for(f in 1:(numParts-1))   {							
								
		for(j in 2:length(UdomLast[,1]-1))  {						
								
			for(k in 1:length(CostList[[f+1]]))   {					
				thisRow<-UdomLast[j,]				
				thisRow[1,f+1]<-k-1				
				thisRow$Cost[1]<-thisRow$Cost[1]+CostList[[f+1]][k]				
				thisRow$Perf[1]<-thisRow$Perf[1]+PerfList[[f+1]][k]				
								
				if(thisRow$Cost[1] %in% UdomAll$Cost)   {				
					pos<-match(thisRow$Cost[1],UdomAll$Cost)			
					if(UdomAll$Perf[pos]<thisRow$Perf[1])  {			
						ex<-UdomAll[pos,]		
						UdomAll[pos,]<-thisRow		
						ExAll<-rbind(ExAll,ex)		
	##  Does this newly identified dom sequence dominate anything elseof higher cost?							
						if(max(UdomAll$Cost)>thisRow$Cost)  {		
						majorDom<-UdomAll[(pos+1):length(UdomAll[,1]),]		
							if(min(majorDom$Perf)<thisRow$Perf)  {	
								ex<-majorDom[sapply(majorDom$Perf, function(x) x<thisRow$Perf),]
								ExAll<-rbind(ExAll,ex)
								majorDom<-majorDom[sapply(majorDom$Perf, function(x) x>thisRow$Perf),]
								UdomAll<-rbind(UdomAll[1:pos,],majorDom)
							}	
						}		
					}else{			
	## this sequence is dominated by a previous sequence at this Cost							
						ExAll<-rbind(ExAll,thisRow)		
					}			
				}else{				
		## this is where additional activity must take place						
		## need to check if this sequence is dominated by a lower cost						
		## need to error trap here for case where this is a new high Cost						
					if(max(UdomAll$Cost)>thisRow$Cost)  {			
						pos<-min(which(UdomAll$Cost>thisRow$Cost))		
						minorDom<-UdomAll[1:(pos-1),]		
						if(max(minorDom$Perf)>thisRow$Perf)  {		
		## this is a new Cost sequence that is dominated by a lesser cost sequence						
							ExAll<-rbind(ExAll,thisRow)	
			##  this must terminate all processing of thisRow					
						}else{		
		## if not, it is dominant, but did it dominate something elseof higher cost?						
		##  Does this newly identified Udom alloccation dominate anything elseof higher cost?						
							majorDom<-UdomAll[pos:length(UdomAll[,1]),]	
							if(min(majorDom$Perf)<thisRow$Perf)  {	
								ex<-majorDom[sapply(majorDom$Perf, function(x) x<thisRow$Perf),]
								majorDom<-majorDom[sapply(majorDom$Perf, function(x) x>thisRow$Perf),]
								ExAll<-rbind(ExAll,ex)
							}	
							UdomAll<-rbind(UdomAll[1:(pos-1),],thisRow,majorDom)	
							## it is impossible to state a domination over anything at lower cost	
						}		
					}else{			
				## at the end of all logic this must be an Undominated Allocation				
						UdomAll<-rbind(UdomAll,thisRow)		
					}			
				}				
			## close incremental new part					
			}					
		## close UdomLast						
		}						
		UdomLast<-UdomAll						
								
	## close all parts							
	}							
								
	title<-"Kettelle Part Allocation Optimization"							
								
	if(data.name!="") {							
		title<-c(title, data.name)						
	}							
	if(performance=="Fill Rate")   {							
		ExAll$Perf<-ExAll$Perf/sum(x[,2])						
		UdomAll$Perf<-UdomAll$Perf/sum(x[,2])						
	  }else{							
	performance="FRN"							
	}

	outnames<-names(UdomAll)
	outcols<-length(outnames)
	outnames[outcols-1]<-performance
	names(UdomAll)<-outnames

								
	if(show==TRUE)  {							
		plot(ExAll$Cost,ExAll$Perf,						
		xlab="Cost",ylab=performance,						
		pch=19, col="blue", cex=0.5,						
		main=title						
								
		)						
		points(UdomAll$Cost,UdomAll[,outcols-1],						
		pch=21, bg="red")						
	}							
	return(UdomAll)							
}								
