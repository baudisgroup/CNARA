#### scripts and sample code

#install.packages("CNAQA", dependencies=c("robfilter", "DNAcopy", "kernlab", "e1071"))
#install the 4 dependencies by yourself at the moment. Will be installed automatically when the package "CNAQA" is built.
#Quality assessment for a series. Output QA metrics, classification and decision value from SVM and case diagnosis.
#Specify the corresponding files/directories/labels.
#setwd("/Users/XXX(input yours)/CNAQA")

workPath <- getwd()
samplePath <- file.path(workPath, "sample")
series <- "sample_series"
probeFileName <- "CNprobes.tab"
segFileName <- "segments.tab"

trainingFile <- file.path(workPath, "trainingSet.txt")
outputFile <- file.path(workPath, "output.txt")

arrayDir <- file.path(samplePath, series)
arrayList <- list.files(arrayDir)
noArray <- length(arrayList)

tableHeader <- paste("Sample", "Series", "Speak", "Breakpoint_Step", "Breakpoint_CBS", "Spread", "Classification_Label", "Decision_Value", "Good_Poor", "Case_Diagnosis", sep = "\t")
write.table(tableHeader, file = outputFile, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

classifier <- trainSVM(trainingFile)

for (arr in 1:noArray) {
	print(paste0("processing ", arrayList[arr], " (", arr, " of ", noArray, ")"))
	probeFile <- file.path(arrayDir, arrayList[arr], probeFileName)
	segFile <- file.path(arrayDir, arrayList[arr], segFileName)

	newCNProbe <- readProbe(probeFile=probeFile, sampleID=arrayList[arr])
	newSpeakCNAno <- calSpeakCNAno(newCNProbe)

	#plot S graph
	#plot(quality(newSpeakCNAno), xlab = "number of iterations", ylab = "S", main=arrayList[arr])

	segNumberCBS <- calCBSBreakpoints(newCNProbe)
	segSpread <- calSpread(newCNProbe, segFile=segFile)
	#segSpread <- calSpread(newCNProbe) # if you don't have segmentation file

	CNAno <- numberOfCNA(newSpeakCNAno)
	Speak <- speak(newSpeakCNAno)

	#create a new object "Metrics" for the copy number profile for quality assessment
	CNProfileMetrics <- createMetrics(sampleID=arrayList[arr], speak=Speak, numberOfCNA=CNAno, cbsBreakpoints=segNumberCBS, spread=segSpread)

	#quality assessment for the copy number profile
	assessment <- assessQuality(CNProfileMetrics, svmClassifier=classifier)

	tmp <- paste(arrayList[arr], series, Speak, CNAno, segNumberCBS, segSpread, assessment$label, assessment$decision.values, assessment$flag, assessment$caseDiag, sep = "\t")
	write.table(tmp, file = outputFile, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

	### write QA log file for each sample
    #sampleLog <- file.path(arrayDir, arrayList[arr], "QA.tab")
	#for (entry in 1:11)
	#	write.table(paste(strsplit(tableHeader, "\t")[[1]][entry], strsplit(tmp, "\t")[[1]][entry], sep="\t"), file = sampleLog, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

	#plot a best fit and counter-fit for a copy number profile. Call calSpeakCNAno again, can be omitted if you just want the assessment without the plot.
 	#plotFCF(calSpeakCNAno(newCNProbe, iterations=CNAno))
}