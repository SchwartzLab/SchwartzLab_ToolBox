# Title: Aligning efficiency summary
# Author: Miguel Angel Garcia-Campos https://github.com/AngelCampos

ls <- list.files()
files <- ls[grep("Log.final", ls)]

for (i in files){
  x <- read.delim(file = i, header = F, sep = "\t", stringsAsFactors = F, row.names = 1)
    assign(i, x)
} 
rm(ls, x, i)

nRow <- nrow(eval(parse(text= ls()[ls() != "files"][1])))
rNames <- rownames(eval(parse(text= ls()[ls() != "files"][1])))
summary <- matrix(NA, nrow = nRow, ncol = length(files))
colnames(summary) <- files; rownames(summary) <- rNames

for (i in files){
  summary[,i] <- eval(parse(text= i))[,1]
}

write.table(file = "mappingSummary.txt", x = summary, sep = "\t", quote = F, col.names = NA)
