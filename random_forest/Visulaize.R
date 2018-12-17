require(ggplot2)
require(reshape2)
DE        = read.table(file="/Users/alvajames/Documents/Sensitivity_JHO/Sensitivity_Journal_comemnts_October_2018/Transposed_Allde.tsv",sep='\t',header=T,stringsAsFactors=FALSE)
Types    = read.table(file="/Users/alvajames/Documents/Sensitivity_JHO/Sensitivity_Journal_comemnts_October_2018/Samples",sep='\t',header=T,stringsAsFactors=FALSE)
DE_sample = merge(DE, Types, on=sample_name)
dim(DE_sample) #  82 1237
#mat_scaled = t(apply(DE_sample[,2:1236], 1, function(x) (x - mean(x)) / sd(x)))
DE_sample=DE_sample[!grepl("LH", DE_sample$Subtype),]
DE_sample=DE_sample[!grepl("unassigned", DE_sample$Subtype),]
df_train =  t(apply(DE_sample[,2:1236], 1, function(x){as.numeric(x)}))
df_train = t(apply(df_train, 1, function(x) (x - mean(x)) / sd(x)))
df_train = as.data.frame(df_train)
colnames(df_train) = colnames(DE_sample)[2:1236]
df_train$Subtype = DE_sample$Subtype
df_train <- melt(df_train, id.vars=c("Subtype"))

df_valid =  t(apply(vali_t_ann[,2:1236], 1, function(x){as.numeric(x)}))
df_valid = t(apply(df_valid, 1, function(x) (x - mean(x)) / sd(x)))
df_valid = as.data.frame(df_valid)
colnames(df_valid) = colnames(vali_t_ann)[2:1236]
df_valid$Subtype = vali_t_ann$subgroup
df_valid <- melt(df_valid, id.vars=c("Subtype"))

df <- rbind(df_train, df_valid)
df$set <- c(rep("train", dim(df_train)[1]), rep("test", dim(df_valid)[1]))
ggplot(df, aes(Subtype,as.numeric(value))) + geom_boxplot(aes(colour=set)) #+ ylim(c(0,5))

varImpPlot(rf, n.var=10)

var_imp = data.frame(importance(rf,type=2))
top_genes = rownames(var_imp)[order(var_imp[,1],decreasing=TRUE)[1:50]]
ggplot(df[which(df$variable %in% top_genes),], aes(Subtype,as.numeric(value))) + geom_boxplot(aes(colour=set)) #+ ylim(c(0,5))