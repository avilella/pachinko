myarg <- commandArgs()
basedir = myarg[length(myarg) - 1]
outfile = myarg[length(myarg)]
print(basedir)
print(outfile)
filelist = list.files(basedir)
for (i in 1:length(filelist)) {
  if (i==1) {
    print(paste("first:",filelist[i]))
#    if (file.info(filelist[i])$size != 0) {
#      print "non-empty"
#    }
    df = read.table(filelist[i], colClasses = c("character","NULL","character","NULL"))
    colnames(df)[1] = paste("C",i,sep='')
  } else {
    print(paste("other:",filelist[i]))
    df = merge(df,read.table(filelist[i], colClasses = c("character","NULL","character","NULL")),by.x="V3",by.y="V3",all=TRUE)
    colnames(df)[i+1] = paste("C",i,sep='')
  }
}
write.table(df[order(df$C1,df$C2),], quote=FALSE, row.names=FALSE, col.names=FALSE, file = outfile,sep=',')
q()
