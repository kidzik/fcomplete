## THIS IS AN INTERNAL SCRIPT FOR DATA PREPARATION
if (!("all.data" %in% ls())){
  all.data = read.csv("/home/kidzik/Dropbox/DATA/CP/alldata.csv")
  gait.cycles = t(read.csv("/home/kidzik/Dropbox/DATA/CP/G_avg_CP.csv"))
  gdi = read.csv("/home/kidzik/Dropbox/DATA/CP/gdi.csv")
  pcas = prcomp(gait.cycles)
  all.data = cbind(all.data, pcas$x[,1:10])
  all.data = merge(all.data, gdi,by = c("Patient_ID","examid","side"))
  trialInfo = read.csv("/home/kidzik/Dropbox/DATA/CP/trialInfo_CP.csv")
}

all.data.subset = all.data[all.data$side == "L",]
all.data.subset = all.data.subset[all.data.subset$age < 25,] # remove outliers
all.data.subset = all.data.subset[all.data.subset$age > 3,] # remove outliers
all.data.subset = all.data.subset[all.data.subset$bmi > 5,] # remove outliers
all.data.subset = all.data.subset[all.data.subset$bmi < 25,] # remove outliers

cmeans = colMeans(all.data.subset[,300:309])
all.data.subset$score = sqrt(rowSums(t(t(all.data.subset[,300:309]) - cmeans)**2))

all.data.filtered = all.data.subset[,c("Patient_ID", "speed","cadence","age","bmi","height","KneeFlex_maxExtension","O2cost","GDI","PC1","PC2","PC3")]
all.data.filtered = all.data.filtered[order(all.data.filtered$Patient_ID, all.data.filtered$age),]
TBL = table(all.data.filtered$Patient_ID)
all.data.filtered = all.data.filtered[all.data.filtered$Patient_ID %in% names(TBL[TBL > 2]),]
all.data.filtered = na.omit(all.data.filtered)
all.data.filtered = all.data.filtered[!is.nan(all.data.filtered$GDI),]
pats = table(all.data.filtered$Patient_ID)
pats = names(pats[pats>=3])

impairement.info = trialInfo[trialInfo["side"]=="L",c("Patient_ID","dxmod")]
impairement.info = impairement.info[!duplicated(impairement.info["Patient_ID"]),]
save(all.data.filtered, impairement.info, file="data/gdi.Rda")
