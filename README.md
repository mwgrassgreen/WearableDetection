# WearableDetection
WearableDetection is a R package for anomaly detection in heart rates from wearble device. Current version provides algorithms of RHR-Diff (resting heart rate difference offline detection) method and CuSum online detection method for wearable data from FitBit smartwatch.


## Dependence
* [R](https://www.r-project.org/) (version >= 3.3.0)

* R package: "xts"

* R function file: "detection_fn.R"

## Usage
`source("./R/detection_fn.R")`

`library("xts")`

To obtain residuals and test statistics based on moving baseline sliding window:

`stats.result = get.stats.fn(dir.hr, dir.step, start.day=NULL, smth.k.par = 10, rest.min.par = 10, base.num.par=28, resol.tm.par=60, test.r.fixed.par=FALSE, test.r.par=NA, res.quan.par=0.9, pval.thres.par=0.01)`

--dir.hr  The working directory for heart rate data.

--dir.step The working directory for step data.

--start.day The starting day for detection.

--smth.k.par The moving average time windonw (default: 10min).

--rest.min.par The resting time after nonzero steps (default: 10min).

--base.num.par The baseline sliding window (default: 28days).

--resol.tm.par The resoluation time for summarized statistics (default: 60min).

--test.r.fixed.par, test.r.par, res.quan.par The threshold in the CuSum statistics. If test.r.fixed.par is FALSE, the threshold is based on the parameter res.quan.par; otherwise, the threshold is based on the parameter test.r.par (default: test.r.fixed.par=FALSE, test.r.par=NA, res.quan.par=0.9).

--pval.thres.par The threshold for p-value (default: 0.01).

The output from get.stats.fn:

`res.t = stats.result$res.t` sequence of the residuals

`cusum.t = stats.result$test.t` sequence of the CuSum statistics

`cusum.t.ind = stats.result$test.t.ind` sequence of run index of the CuSum statistics

`cusum.test.r = stats.result$test.r.seq`sequence of thresholds of the CuSum statistics 

`cusum.pval = stats.result$test.pval` sequence of p-value for CuSum statistics

To implement offline detection:

`offline.result = rhr.diff.detection.fn(res.t, alpha=0.05)`

`write.csv(offline.result, file="./result/RHRDiff_offline_detection.csv" )`

Reference: https://github.com/mwgrassgreen/RankScan

--alpha The significance level to control FWER.

To implement online detection:

`cusum.alarm.result = cusum.detection.fn(cusum.t, cusum.t.ind, cusum.pval, pval.thres=0.01, max.hour=24, dur.hour=48)`

`write.csv(cusum.alarm.result, file="./result/CuSum_online_detection.csv" )`

--pval.thres The threshold for p-value (default: 0.01).

--max.hour The threshold for the maximum hour (default: 24hours).

--dur.hour The threshold for duration hours (default: 48hours).


Visualization of the detection results:

`result.plot.fn(id, sym.date=NA, diag.date=NA,  res.t, cusum.t, cusum.test.r, offline.result, cusum.alarm.result)`

--id The participant id.


## Example 

Input data from one participant:

`dir.hr = "./data/AHYIJDV_hr.csv" ` raw heart rate data

`dir.step = "./data/AHYIJDV_step.csv" ` step data

`id = "AHYIJDV" `

Output of detection results:

`./result/RHRDiff_offline_detection.csv` offline detection result

`./result/CuSum_online_detection.csv` online detection result

`./figure/detetion_plot.pdf` detection plot

The time to run the whole procedure in this example is approximate 16s.
