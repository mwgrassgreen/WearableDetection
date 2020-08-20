#### source functions to implement RHR-Diff offline detection and CuSum online detection for short term data
#### Author: Meng Wang
#### Email: mengw1@stanford.edu
#### Date: 2020 August


########################################################
####    to obtain residuals and CuSum statistics    ####
####            based on sliding window             ####
########################################################

get.stats.fn = function(dir.hr, dir.step, start.day="2020-02-01",  smth.k.par = 10, rest.min.par = 10, base.num.par=28, resol.tm.par=60, test.r.fixed.par=FALSE, test.r.par=NA, res.quan.par=0.9, pval.thres.par=0.01) {
	
	formatted.dat = dat.format.fn(dir.hr, dir.step)
	dat.hr.1 = formatted.dat$dat.hr
	dat.step.1 = formatted.dat$dat.step
	# head(dat.hr.1)
	# head(dat.step.1)
	
	hr.dates = unique(dat.hr.1[, 2])
	gap.dates = as.numeric(diff(as.Date(hr.dates)))
	hr.date.start = hr.dates[-1][gap.dates >= 60]
	hr.date.start = hr.date.start[length(hr.date.start)]
	if (length(hr.date.start) == 1) {
	   dat.hr.1 = dat.hr.1[as.numeric(as.Date(dat.hr.1[,2]) - as.Date(hr.date.start)) >= 0, ]
	}
	
	if (!is.null(start.day)) {
		hr.dates = unique(dat.hr.1[, 2])
		if (as.numeric(as.Date(hr.dates[1]) - as.Date(start.day)) < 0) {
			dat.hr.1 = dat.hr.1[as.numeric(as.Date(dat.hr.1[,2]) - as.Date(start.day)) >= 0, ]
		}
	}
		
	#####	
	day.nm = unique(dat.hr.1[,"date"])
	length(day.nm) 
	timerange1 = "20160106 0000/20160106 2359"
	min.nm = format(timeBasedSeq(timerange1), "%H:%M:%S")
		
	hr.mx = matrix(NA, length(min.nm), length(day.nm))
	rownames(hr.mx) = min.nm
	colnames(hr.mx) = day.nm
	for (i in 1:length(day.nm)) {
		day.1 = day.nm[i]
		day.hr = dat.hr.1[dat.hr.1[,"date"] == day.1, "hr"]
		names(day.hr) = dat.hr.1[dat.hr.1[,"date"] == day.1, "time"]
		hr.mx[names(day.hr), day.1] = day.hr
	}
	head(hr.mx)
	
	hr.mx = apply(hr.mx, 2, as.numeric)
	rownames(hr.mx) = min.nm
	colnames(hr.mx) = day.nm
		
	if (ncol(hr.mx) > base.num.par) {
	
		# to initialize baseline 
	    day.train.0 = colnames(hr.mx)[1:base.num.par]
		dat.train.0 = hr.mx[, day.train.0]
		X = dat.train.0 
		rownames(X) = rownames(dat.train.0)
		x = c(X)
		K = smth.k.par
		x.conv = c(rep(NA, K-1), rollmean(x, k=K, align="right"))
		names(x.conv) = paste(rep(colnames(X), each=nrow(X)), rep(rownames(X),ncol(X)), sep=" ")
		X.conv = matrix(x.conv, nrow=nrow(X))
		rownames(X.conv) = rownames(X)
		colnames(X.conv) = colnames(X)
	
	    X.conv.f = X.conv
		X.conv.f[is.na(X)] = NA
		for (j in 1:length(day.train.0)) {
			day.1 = day.train.0[j]
			day.step.1 = dat.step.1[dat.step.1[,"date"] == day.1, c("time", "steps")]
			ind = day.step.1[,2]
			names(ind) = day.step.1[,1]
			ind.1 = names(ind)[ind > 0]
			ind.ext = unique(c(unlist(sapply(ind.1, function(tm) rownames(X.conv)[grep(tm, rownames(X.conv))+0:rest.min.par] ))))
			ind.ext = ind.ext[!is.na(ind.ext)]
			X.conv.f[ind.ext ,j] = NA
		}
		
		result = stats.fn(hr.mx, dat.step.1, X.conv.f, base.num=base.num.par, resol.tm.1=resol.tm.par, test.r.fixed.1=test.r.fixed.par, test.r.1=test.r.par, res.quan.1=res.quan.par, pval.thres=pval.thres.par) 
			         
	} else {
		result = NULL
		print("number of total recorded days is less than baseline days")
	}
    
    return(result)

}


#### preformatting the dataset
dat.format.fn = function(dir.hr, dir.step) {
	dat.hr = read.csv(dir.hr)
	dat.hr = as.matrix(dat.hr)
	#head(dat.hr)
	
	# to summarize hear rate in one minute resolution
	hr.med = tapply(as.numeric(dat.hr[, "heartrate"]), substring(dat.hr[,"datetime"], 1, 16), median, na.rm=TRUE)
	dat.hr.1 = cbind(rep(dat.hr[1, "user"], length(hr.med)), paste(names(hr.med), ":00", sep=""), hr.med )
	dat.hr.1 = cbind(dat.hr.1[,1], substring(dat.hr.1[,2], 1, 10), substring(dat.hr.1[,2], 12), dat.hr.1[,3])
	colnames(dat.hr.1) = c("user", "date", "time", 'hr')
	dat.hr.1[,"hr"] = as.numeric(dat.hr.1[,"hr"])
	#head(dat.hr.1)
	
	dat.step = read.csv(dir.step)
	head(dat.step)
	dat.step.1 = dat.step[, c("user", "datetime", "steps")]
	dat.step.1 = cbind(dat.step.1[,1], substring(dat.step.1[,2], 1, 10), substring(dat.step.1[,2], 12), dat.step.1[,2:3])
	colnames(dat.step.1) = c("user", "date", "time", 'datetime', 'steps')
	#head(dat.step.1)

    return(list(dat.hr=dat.hr.1, dat.step=dat.step.1))
}



#### to obtain daily resting heart rate
daily.conv.filter.fn = function(x, day.1, dat.step.1, K.conv=10, K.rest=10) {
	x.conv = c(rep(NA, K.conv-1), rollmean(x, k=K.conv, align="right"))
	day.step.1 = dat.step.1[dat.step.1[,"date"] == day.1, c("time", "steps")]
	ind = day.step.1[,2]
	names(ind) = day.step.1[,1]
	ind.1 = names(ind)[ind > 0]
	ind.ext = unique(c(unlist(sapply(ind.1, function(tm) names(x)[grep(tm, names(x))+0:K.rest] ))))
	ind.ext = ind.ext[!is.na(ind.ext)]
	x.conv.f = x
	x.conv.f[ind.ext] = NA
	names(x.conv.f) = paste(rep(day.1, length(x)), names(x), sep=" ")
    return(x.conv.f) 
}

#### to obtain baseline residuals 
base.stand.res.fn = function(X.conv.f.base, resol.tm=60, test.r.fixed=FALSE, test.r=NULL, res.quan=0.9){
	cnt.tm = nrow(X.conv.f.base)/resol.tm	
	X.bin = sapply(1:cnt.tm, function(k) {
		X.sub = X.conv.f.base[(1+(k-1)*resol.tm):(k*resol.tm), ]
		c(X.sub)
	} )
	dim(X.bin)
	
	x.cen = apply(X.bin, 2, mean, na.rm=TRUE)
	x.sd = apply(X.bin, 2, sd, na.rm=TRUE)
	
	z.mx = matrix(NA, cnt.tm, ncol(X.conv.f.base))
	for (j in 1:ncol(X.conv.f.base)) {
		x = X.conv.f.base[,j]
		x.bin = matrix(x, resol.tm, cnt.tm)
		x.avg = colMeans(x.bin, na.rm=TRUE)
	    z.mx[,j] = (x.avg - x.cen)/x.sd		
	}
	
	z.abs.thres=4
	z.v = c(z.mx)
	z.v = pmin(z.v, z.abs.thres)
	z.v = pmax(z.v, - z.abs.thres)
	names(z.v) = paste(rep(colnames(X.conv.f.base), each=cnt.tm), rep(rownames(X.conv.f.base)[seq(1, nrow(X.conv.f.base), by=resol.tm)], ncol(X.conv.f.base)), sep=" ")
	
	if (test.r.fixed ==FALSE) {
		test.r = quantile(z.v, res.quan, na.rm=TRUE)/2
	}	
	cusum.v = cusum.fn(z.v, test.r)
	return(list(z.v=z.v, test.r=test.r, cusum.v=cusum.v))
}


#### to update residuals along with time
update.stand.res.fn = function(X.conv.f.base, x.conv.f, resol.tm=60, z.abs.thres=4){
	cnt.tm = nrow(X.conv.f.base)/resol.tm
	X.bin = sapply(1:cnt.tm, function(k) {
		X.sub = X.conv.f.base[(1+(k-1)*resol.tm):(k*resol.tm), ]
		c(X.sub)
	} )
	dim(X.bin)
	
	x.cen = apply(X.bin, 2, mean, na.rm=TRUE)
	x.sd = apply(X.bin, 2, sd, na.rm=TRUE)
	
	x.bin = matrix(x.conv.f, resol.tm, cnt.tm)
	x.avg = colMeans(x.bin, na.rm=TRUE)
	z.v = (x.avg - x.cen)/x.sd
	z.v = pmin(z.v, z.abs.thres)
	z.v = pmax(z.v, - z.abs.thres)		
	names(z.v) = paste(rep(substring(names(x.conv.f)[1], 1, 10), cnt.tm), rownames(X.conv.f.base)[seq(1, nrow(X.conv.f.base), by=resol.tm)], sep=" ")
	
	return(z.v)
}

### to obtain initial CuSum stats
cusum.fn = function(z.v, test.r=1) {
	z.t = z.v[!is.na(z.v)]
	w.t = as.numeric(length(z.t))
	w.t[1] = 0
	for (i in 1:length(z.t)) {
		w.t[i+1] = max(0, w.t[i]+z.t[i]-test.r) 
	}
	w.t = w.t[-1]
	names(w.t) = names(z.t)
	return(cusum.stats = w.t)
}

#### to update CuSum stats along with time
udpate.cusum.fn = function (x.res, cusum.t, test.r=1) {
	z.t = x.res[!is.na(x.res)]
	w.t = as.numeric(length(z.t))
	w.t[1] = cusum.t[length(cusum.t)]
	for (i in 1:length(z.t)) {
		w.t[i+1] = max(0, w.t[i]+z.t[i]-test.r) 
	}
	w.t = w.t[-1]
	names(w.t) = names(z.t)
	return(cusum.stats = w.t)
}


#### to obtain emprical null distribution for test stats	
null.ecdf.fn = function(base.test) {
	 ind.0 = (1:(length(base.test)-1))[base.test[1:(length(base.test)-1)] == 0]
	 ind.0 = ind.0[base.test[ind.0 + 1] > 0]
	 if (length(ind.0) == 0) {
	 	null.ecdf= NULL
	 } else{
		 max.step = ifelse(length(ind.0) == 1, 1, max(diff(ind.0)))
		 null.ecdf = list(1:max.step)			
		 for (step.k in 1:max.step) {
			 max.test = sapply(1:length(ind.0), function(j) max(base.test[ind.0[j]:(ind.0[j] + step.k)]))
			 max.test = max.test[!is.na(max.test)]
		     null.ecdf[[step.k]] = ecdf(max.test)
		 }	 	
	 }
	 
     return(null.ecdf)
 } 
 
#### run index for CuSum stats
base.test.ind.fn = function(base.test) {
	 base.test.ind = rep(0, length(base.test))
	 base.test.ind[1] = ifelse(base.test[1] == 0, 0, 1)
	 for (i in 2:length(base.test.ind)) {
	 	 if (base.test[i] > 0) {
	 	 	 base.test.ind[i] = base.test.ind[i-1] + 1
	 	 }
	 }
	 names(base.test.ind) = names(base.test)
	 return(base.test.ind)			 	
}			

#### to obtain residuals and CuSum stats compared to sliding window baseline
stats.fn = function (hr.mx, dat.step.1, X.conv.f, base.num=28, resol.tm.1=60, test.r.fixed.1=FALSE, test.r.1=NA, res.quan.1=0.9, pval.thres=0.01) {
	
	# to obtain initial baseline residuals, CuSum stats, null distribution and p-value
	X.conv.f.base = X.conv.f
	base.info = base.stand.res.fn(X.conv.f.base, resol.tm.1, test.r.fixed.1, test.r=test.r.1, res.quan.1)
    test.r.base = base.info$test.r
	base.res = base.info$z.v		
    test.r.seq = rep(test.r.base, length(base.res))
    names(test.r.seq) = names(base.res)   	
	base.test = base.info$cusum.v
	null.ecdf.test = null.ecdf.fn(base.test)			
    base.test.ind = base.test.ind.fn(base.test)
   
    base.days = colnames(X.conv.f.base)
	day.test = colnames(hr.mx)[(base.num+1):ncol(hr.mx)]
	res.t = base.res
	
	score.thres = 1 - pval.thres
	test.t = base.test
	test.t.ind = base.test.ind
	risk.score = rep(NA, length(base.test))
	names(risk.score) = base.test
	for (j in 1:length(base.test.ind)) {
		if (base.test.ind[j] >= 2 & base.test.ind[j] <= length(null.ecdf.test)) {
			risk.score[j] =  null.ecdf.test[[base.test.ind[j]]](base.test[j] - 0.001)
		}
		if (base.test.ind[j] > length(null.ecdf.test)) {
			risk.score[j] = 1
		}				
	}
	names(risk.score) = names(base.test)	
		
	# to obtain residuals, CuSum stats, and p-value along with time
	for (k in 1:length(day.test)) {
		#print(k)
		day.1 = day.test[k]
		x = hr.mx[, day.1]
	
		# moving average from past 10mins
		x.conv.f = daily.conv.filter.fn(x, day.1, dat.step.1, K.conv=10, K.rest=10) 
			
		# to obtain standardized residuals based on baseline
	    x.res = update.stand.res.fn(X.conv.f.base, x.conv.f, resol.tm=resol.tm.1, z.abs.thres=4)

        if (sum(!is.na(x.res)) >= 2 ) { 
            res.t = c(res.t, x.res)
	        x.test.r = rep(test.r.base, length(x.res))
	        names(x.test.r) = names(x.res)
            test.r.seq = c(test.r.seq, x.test.r)
            
			# to obtain cusum stats
 			x.test = udpate.cusum.fn(x.res, test.t, test.r=test.r.base)			
			test.t = c(test.t, x.test)
			
			x.test.ind = rep(0, length(x.test))
			names(x.test.ind) = names(x.test)
			x.test.ind[1] = ifelse(x.test[1] == 0, 0, 1)
			if (test.t.ind[length(test.t.ind)] > 0 & x.test.ind[1] == 1) {
				x.test.ind[1] = 1 + test.t.ind[length(test.t.ind)]
			}
			for (i in 2:length(x.test)) {
				if (x.test[i] > 0) {
				  x.test.ind[i]	= x.test.ind[i-1] + 1					
				}				
			}
			test.t.ind = c(test.t.ind, x.test.ind)

			x.risk.score = rep(NA, length(x.test))
			names(x.risk.score) = names(x.test)
			for (j in 1:length(x.test.ind)) {
				if (x.test.ind[j] >= 2 & x.test.ind[j] <= length(null.ecdf.test)) {						x.risk.score[j] = null.ecdf.test[[x.test.ind[j]]]( test.t[(1:length(test.t))[names(test.t.ind) == names(x.test[j])] ]  - 0.001)
				}
				if (x.test.ind[j] > length(null.ecdf.test)) {
					x.risk.score[j] = 1
				}				
			}
			risk.score = c(risk.score, x.risk.score)	
						
			if (sum(x.risk.score > score.thres, na.rm=TRUE) == 0) {

				base.days = c(base.days, day.1)
			    X.conv.f.base = cbind(X.conv.f.base, x.conv.f)
			    X.conv.f.base = X.conv.f.base[, max(1, ncol(X.conv.f.base) - base.num+1):ncol(X.conv.f.base)]
				base.info = base.stand.res.fn(X.conv.f.base, resol.tm.1, test.r.fixed.1, test.r=test.r.1, res.quan.1)

				base.res = base.info$z.v
				base.test = base.info$cusum.v
				test.r.base = base.info$test.r				
				null.ecdf.test = null.ecdf.fn(base.test)						
			    base.test.ind = base.test.ind.fn(base.test)
			    			    			       
			}
			
		} 
	
	}
		
	return(list(res.t = res.t, test.t=test.t, test.t.ind=test.t.ind, test.pval=1-risk.score, test.r.seq= test.r.seq))
}


#########################################
####    CuSum online detection       ####
#########################################
cusum.detection.fn = function (id, test.t, test.t.ind, test.pval, pval.thres=0.01, max.hour=24, dur.hour=48) {
    risk.score = 1 - test.pval 
    test.sum = cbind(test.t, test.t.ind, risk.score)
    
	ind = rep(0, length(test.t.ind))
	names(ind) = names(test.t.ind)
	ind[test.t.ind > 0] = 1
	ind.mx = with(rle(ind), data.frame(number = values, start = cumsum(lengths) - lengths + 1, end = cumsum(lengths))[order(values),])
    ind.mx = ind.mx[ind.mx$number == 1, ]
    
    score.thres = 1 - pval.thres
    risk.score[is.na(risk.score)] = 0
    test.alarming.tm = c()
    run.start = c()
    run.max = c()
    run.num = c()
    run.max.val = c()
    run.len = c()
    for (i in 1:nrow(ind.mx)) {
    	x.risk = risk.score[ind.mx[i, 2]:ind.mx[i,3]]
    	x.test = test.t[ind.mx[i, 2]:ind.mx[i,3]]
    	if (max(x.risk, na.rm=TRUE) > score.thres) {
    		tm = names(x.risk)[min((1:length(x.risk))[x.risk > score.thres])]   
    		test.alarming.tm = c(test.alarming.tm, tm)
    		run.start = c(run.start, names(x.risk)[1])
    		run.max = c(run.max, names(which.max(x.test)))
    		run.num = c(run.num, which.max(x.test))
    		run.max.val = c(run.max.val, max(x.test, na.rm=TRUE))
    		run.len = c(run.len, length(x.risk))		
    	}
    }
    if (!is.null(test.alarming.tm)) {
	    alarm.tab = cbind(test.alarming.tm, run.start, run.max, run.num, round(run.max.val,2), run.len)
	    alarm.tab = as.matrix(alarm.tab, ncol=6)
	    if (sum(as.numeric(alarm.tab[,6]) > dur.hour & as.numeric(alarm.tab[,4]) > max.hour) == 0 ) {
	    	alarm.sum = matrix(c(id, rep("N/A", 7 )), ncol=8)
            colnames(alarm.sum) = c("ParticipantID", "alarming time", "alarming hour", "duration hours to max", "total hours", "max CuSum statistics", "CuSum statistics at alarm", "p-value")

	    } else {
		    alarm.tab = alarm.tab[as.numeric(alarm.tab[,6]) > dur.hour & as.numeric(alarm.tab[,4]) > max.hour, ]
		    alarm.tab = matrix(alarm.tab, ncol=6)
		    colnames(alarm.tab) = c("alarming.time", "anormaly_starting_time", "anormaly_max_time", "duration_hours_to_max", "max_test_stats", "total_hours")	    
		    alarm.sum = cbind(rep(id, nrow(alarm.tab)), alarm.tab[, "alarming.time"], test.sum[alarm.tab[, "alarming.time"], 2],          test.sum[alarm.tab[, "anormaly_max_time"], 2], alarm.tab[,"total_hours"], alarm.tab[,"max_test_stats"],  round(as.numeric(test.sum[alarm.tab[, "alarming.time"], 1]), 2), round(1 - as.numeric(test.sum[alarm.tab[, "alarming.time"], 3]), 2) )
		    colnames(alarm.sum) = c("ParticipantID", "alarming time", "alarming hour", "duration hours to max", "total hours", "max CuSum statistics", "CuSum statistics at alarm", "p-value")
		    rownames(alarm.sum) = NULL
	    }
    } else {
    	alarm.sum = matrix(c(id, rep("N/A", 7)), ncol=8)
        colnames(alarm.sum) = c("ParticipantID", "alarming time", "alarming hour", "duration hours to max", "total hours", "max CuSum statistics", "CuSum statistics at alarm", "p-value")
    }
    
    return(alarm.sum)
}



#########################################
####    RHR-Diff offline detection   ####
#########################################
       
# reference: https://github.com/mwgrassgreen/RankScan
# paper: Arias-Castro, Ery, Rui M Castro, Ervin Tánczos, and Meng Wang. 2018. "Distribution-free detection of structured anomalies: Permutation and rank-based scans." Journal of the American Statistical Association 113 (522): 789-801.


#### rank scan detection
rhr.diff.detection.fn = function (id, res.t, alpha=0.05, B=1000) {
	timerange1 = "20160106 0000/20160106 2359"
	min.nm = format(timeBasedSeq(timerange1), "%H:%M:%S")
	hour.nm = min.nm[(1:length(min.nm)) %% 60 == 1]
	
	nm.dum = seq(as.Date(names(res.t)[1]), as.Date(names(res.t)[length(res.t)]), by=1)
    nm.dum = paste(rep(nm.dum, each=24), rep(hour.nm, length(nm.dum)), sep=" ")
    res.t.1 = rep(0, length(nm.dum))
    names(res.t.1) = nm.dum
    res.t.1[names(res.t)] = res.t
    res.t.1[is.na(res.t.1)] = 0

	z.v = res.t.1
	n = length(z.v)
	Q = floor(log2(n))
	null.result = crt.val.rank(n, Q, alpha, B) 
	d = null.result$crt.val
	null.B = null.result$test.null
	#d = 43120.513817 # n=2^15
	a = rank(z.v)
	names(a) = NULL
	q = floor(log2(length(z.v)))
	result = rank_scan_est(a, q, d)
	pval = sapply(abs(result[,3]), function(x) sum(null.B >= x)/(1+B))	
	detect = cbind(names(z.v)[result[,1]], names(z.v)[result[,2]], round(result[,3],3), round(pval, 3))
	if (ncol(detect) == 1) {
		detect.pos = matrix(c(id, rep("N/A", 4)), ncol=5)
	} else {
		detect = detect[!duplicated(detect[,3]), ]
		detect = matrix(detect, ncol=4)	
        hr.diff = difftime(detect[,2], detect[,1], unit="hours")
		detect = detect[as.numeric(hr.diff) >= 24, ] # to remove the detected interval < 24 hours
		detect = matrix(detect, ncol=4)
		if (sum(detect[,3] > 0) > 0) {
			detect.pos = matrix(detect[detect[,3] > 0, ],ncol=4)
			detect.pos = matrix(cbind(rep(id, nrow(detect.pos)), detect.pos ), ncol=5)
		} else {
			detect.pos = matrix(c(id, rep("N/A", 4)), ncol=5)
		}
	}
	colnames(detect.pos) = c("ParticipantID", "staring time of detected interval", "ending time of detected interval", "rank scan test statistics", "p-value")
	return(detect.pos)
}



# to get critival value for rank scan test
# only scan intervals in length 2^(1:(log2(N)-1))

GetMaxAtk <- function(X, mu.X, k){
            # k--the length of signal interval candidate	
	        N <- length(X)
	        Intv <- c(rep(1, k), rep(0, N-k))
	        return(max(convolve(X - mu.X, Intv)[1:(N-k+1)])/sqrt(k))
}

GetRankScan <- function(Rank.X, Q){
	        N <- length(Rank.X)
	        # scan intervals in dyalic lengths from 2 to [N/2]
	        #Q <- floor(log2(N))
	        max.at.k <- Q
	        for (q in 1:(Q-1)){
	        	max.at.k[q] <- GetMaxAtk(Rank.X, (N+1)/2, 2^q)
	        }
	        stats <- max(max.at.k)
            return(stats)
}


crt.val.rank = function (N, Q, alpha, B) {
	rank.scan.null <- rep(NA, B)
	for (b in 1:B){
		Rank.X.null <- sample(1:N)			
		rank.scan.null[b] <- GetRankScan(Rank.X.null, Q)
	}
	
	crt.val <- quantile(rank.scan.null, 1 - alpha)
    return(list(crt.val=crt.val, test.null=rank.scan.null))
}



#### to get identified anomalous intervals from rank scans

rank_scan_est = function(a, q, d){
	
	# a is the rank sequence of the orignial data X
	# q is the log_2 of maximum length of the candidate intervals
	# d is the threshold for normalized sum of ranks in one interval i.e., 1/sqrt(|S|) sum_{v in S} (R_v - (N+1)/2) in the notation of (Arias-Castro et al 2018)
	# reference for reporting the identified intervals: Jeng, X. J., Cai, T. T., and Li, H. (2010), “Optimal Sparse Segment Identification With Application in Copy Number Variation Analysis” 

	
	b = length(a)
	x = rep(0, b*(2^q+1)) # store all the intervals in length 2^(1:q)
	for (i in 1:(b-1)) {
		#print(i)
		# i is the starting index of a candidate interval
		for (j in pmin(b, i+2^(1:q)-1)) {
		    # j is its ending index
		    #print(j)
			x[(i-1)*2^q + j] = sum(a[i:j] - (b+1)/2)/sqrt(j-i+1)
			#print(x[(i-1)*2^q + j])
		}
	}
	
	
	k= which(abs(x)>d);
	i = ceiling(k/(2^q+1)); 
	j = k - (i-1)*2^q;
	list = cbind(i,j, x[k]);
	
	start = rep(0,1);
	end = rep(0,1);
	Rank_scan = rep(0,1);
	t=1;
	
	while (length(list)> 3) {
	ind = which(abs(list[,3]) == max(abs(list[,3])));
	len.ind = length(ind)
	start[t:(t+len.ind-1)] = list[ind,1];
	end[t:(t+len.ind-1)] = list[ind,2];
	Rank_scan[t:(t+len.ind-1)] = list[ind,3];
	II = c()
	for (l in 1:len.ind) {
		s = t+l-1
		II = c(II, which(list[,1]<=end[s] & list[,2]>=start[s]))
	}
	#II = which(list[,1]<=end[t] & list[,2]>=start[t]);
	list = list[-II,];
	t = t+len.ind; 
	} 
	
	if(length(list)==3) {
	start[t] = list[1];
	end[t] = list[2];
	Rank_scan[t] = list[3];
	}
	
	peaks = cbind(start, end, Rank_scan); 
	return(peaks)
}

#########################################
####    plot of detection result     ####
#########################################

result.plot.fn = function (id, sym.date=NA, diag.date=NA, res.t, cusum.t, cusum.test.r, offline.result, cusum.alarm.result) {

    pdf(file="./figure/detectoin_plot.pdf", width=17, height=8)
    par(mfrow=c(2,1))
    par(mai=c(0.5, 1, 0.5, 0.1))		
    
    timerange1 = "20160106 0000/20160106 2359"
	min.nm = format(timeBasedSeq(timerange1), "%H:%M:%S")
	hour.nm = min.nm[(1:length(min.nm)) %% 60 == 1]
    nm.dum = seq(as.Date(names(res.t)[1]), as.Date(names(res.t)[length(res.t)]), by=1)
    nm.dum = paste(rep(nm.dum, each=24), rep(hour.nm, length(nm.dum)), sep=" ")
    res.t.1 = rep(0, length(nm.dum))
    names(res.t.1) = nm.dum
    res.t.1[names(res.t)] = res.t
    res.t.1[is.na(res.t.1)] = 0
	plot(res.t.1, type="l", xaxt="n", xlab="", ylab="residuals", cex.lab=2, cex.main=2, main=paste("offline detection for", id), ylim=c(-4, 4)) #  ylim=c(min(-2.5, min(res.t.1), max(2.5, max(res.t.1))))
	pos = (1:length(res.t.1))[!duplicated(substring(names(res.t.1), 1, 10))]
	names(pos) = unique(substring(names(res.t.1), 1, 10))
	#axis(1, pos,  names(pos), las=2, cex.axis=1)
	axis(1, pos,  labels=NA, las=2, cex.axis=1)	
	test.r.thres =  cusum.test.r
    test.r.thres.1 =  test.r.thres[names(res.t.1)]
    plot(stepfun(1:(length(res.t.1)), c(test.r.thres.1,test.r.thres.1[length(test.r.thres.1)] )), pch=".", col='darkgreen', add=TRUE,lwd=2)
    abline(h=c(0),lty=c(2), lwd=2, col='darkgreen')  
	abline(v=(1:length(res.t.1))[!duplicated(substring(names(res.t.1), 1, 10))], lty=2, col="grey")
	if (sum(diag.date %in% names(pos)) > 0) {
		#axis(1, pos[diag.date],  diag.date, col.axis="purple", las=2, cex.axis=1)
	    abline(v=pos[diag.date], col="purple", lwd=2, lty=2)			
	}
	if (sum(sym.date %in% names(pos)) > 0 ) {
	   #axis(1, pos[sym.date],  sym.date, col.axis="red", las=2, cex.axis=1)
	   abline(v=pos[sym.date], col="red", lwd=2, lty=2)
	}
	line.pos=-2
	z.v = res.t.1
	detect = offline.result[, c("staring time of detected interval", "ending time of detected interval")]
	detect = matrix(detect, ncol=2)
	if (detect[1] != "N/A"){
		for (i in 1:nrow(detect)) {
				arrows((1:length(z.v))[names(z.v) == detect[i,1] ], line.pos, (1:length(z.v))[names(z.v) == detect[i,2] ], line.pos, col="red", lwd=2, length=0.1)
		        arrows((1:length(z.v))[names(z.v) == detect[i,2] ], line.pos, (1:length(z.v))[names(z.v) == detect[i,1] ], line.pos, col="red", lwd=2, length=0.1)		

	    }	
	}


    ##### 
	test.t = cusum.t
	test.alarming.tm = cusum.alarm.result[,"alarming time"]
	plot(test.t, type="l", xaxt="n", xlab="", ylab="CuSum statistics", cex.lab=2, cex.main=2, main=paste("online detection for", id))
	points(test.t, pch=".", cex=3)
	abline(v=(1:length(test.t))[!duplicated(substring(names(test.t), 1, 10))], lty=2, col="grey")
	pos = (1:length(test.t))[!duplicated(substring(names(test.t), 1, 10))]
	names(pos) = unique(substring(names(test.t), 1, 10))
	#axis(1, pos, names(pos), las=2, cex.axis=1)
	axis(1, pos,  labels=NA, las=2, cex.axis=1)
	if (sum(diag.date %in% names(pos)) > 0 ) {
		#axis(1, pos[diag.date],  diag.date, col.axis="purple", las=2, cex.axis=1)
	    abline(v=pos[diag.date], col="purple", lwd=2, lty=2)			
	}
	if (sum(sym.date %in% names(pos)) > 0 ) {
	   #axis(1, pos[sym.date],  sym.date, col.axis="red", las=2, cex.axis=1)
	   abline(v=pos[sym.date], col="red", lwd=2, lty=2)
	}
	
	pos = 1:length(test.t)
	names(pos) = names(test.t)
	abline(v=pos[test.alarming.tm], col="blue", lwd=2, lty=2 )
    dev.off()
}

