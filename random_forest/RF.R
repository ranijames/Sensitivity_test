# Random forest on the validation cohort

library(randomForest)
Validation_cohort    = read.table(file="/Users/alvajames/Documents/Sensitivity_JHO/Sensitivity_Journal_comemnts_October_2018/BCP_ALL_DE_on_Berlin_count_noDUO",sep='\t',header=T,stringsAsFactors=FALSE)
# Removing the un-wamted columns
Validation_cohort$Berlin2492<-NULL
Validation_cohort$Berlin2493<-NULL
Validation_cohort$Berlin2494<-NULL
Info=read.table("/Users/alvajames/Documents/Sensitivity_JHO/Sensitivity_Journal_comemnts_October_2018/Sample_group3_subtype.tsv", header=TRUE, sep="\t",check.names=FALSE)
Annotation   = read.table(file="/Users/alvajames/Documents/Sensitivity_JHO/Sensitivity_Journal_comemnts_October_2018/lncRNAs_annotation",sep='\t',header=T,stringsAsFactors=FALSE)
Validation_cohort_annon=merge(Validation_cohort,Annotation, on='Geneid')
validation_class=merge(Validation_cohort_annon,All_de[ , c("Gene_sym", "Subtype")],on='Gene_sm')
validation_class=validation_class[names(validation_class) %in% Info$sample_name]
row.names(validation_class) <- validation_class$Gene_sym
validation_class$Gene_sym <-NULL
dim(validation_class) #1235 49


table(Subtype)
# DUX4  NHHeH Phlike 
# 565    287    383 
set.seed(1000)
#Select the trianing data for 70%
1235* 0.7 
# 864.5
## Pull random 865 rows from 1235 rows
train     = sample(1:1235,865,replace = FALSE)
traindata = validation_class[train,]
Subtype   = factor(traindata$Subtype)
traindata = t(apply(traindata[,1:47], 1, function(x) (x - mean(x)) / sd(x)))

testdata  = validation_class[-train,]

#Chaniging the labels, here Subtypes into factors as randomforest excepts only factors. This is response vector
traindata$Subtype = factor(traindata$Subtype) 
fit       = randomForest(Subtype ~. , data= traindata,ntree=1400,importance=TRUE,proximity=TRUE)
fit4       = randomForest(Subtype ~. , data= traindata,ntree=1200,importance=TRUE,proximity=TRUE)

predictions <- as.data.frame(predict(fit, testdata, type = "prob"))
predictions$predict <- names(predictions)[1:3][apply(predictions[,1:3], 1, which.max)]
testdata$Subtype=factor(testdata$Subtype)
predictions$observed <- testdata$Subtype

### Sensitivity and specificity
RFPredictiction_vali <- prediction(predictions , Subtype)
roc.dux <- roc(ifelse(predictions$observed=="DUX4","DUX4","non-DUX4"), as.numeric(predictions$DUX4))
plot(roc.dux, col = "gray60")

# others
roc.Phlike <- roc(ifelse(predictions$observed=="Phlike", "Phlike", "non-Phlike"), as.numeric(predictions$Phlike))
roc.NHeH <- roc(ifelse(predictions$observed=="NHHeH", "NHHeH", "non-NHHeH"), as.numeric(predictions$NHHeH))
lines(roc.Phlike, col = "blue")
lines(roc.NHeH, col = "red")

############# Using Random forest on Discover chohort to build a model
DE_sample$Subtype=factor(DE_sample$Subtype)
vali_t_ann$subgroup=factor(vali_t_ann$subgroup)
Modelon_disCH    =randomForest(DE_sample$Subtype ~. , data= DE_sample[,2:1236],ntree=1000,importance=TRUE,proximity=TRUE)

###################### Iteration  ###############################

N = 10
sens = spec =as.data.frame(matrix(NA, nrow=1, ncol=3))
for(i in 1:N){
  set.seed(250+i)
  train     = sample(1:dim(DE_sample)[1],2*dim(DE_sample)[1]/3,replace = FALSE)
  traindata = DE_sample[train,]
  Subtype   = factor(traindata$Subtype)
  #traindata = t(apply(traindata[,1:47], 1, function(x) (x - mean(x)) / sd(x)))
  testdata  = DE_sample[-train,]
  rf = randomForest(Subtype ~. , data=traindata[,2:1236],ntree=1000,importance=TRUE,proximity=TRUE)
  
  pred_valid <- as.data.frame(predict(rf, testdata, type = "prob"))
  predicted <- as.factor(names(pred_valid)[1:3][apply(pred_valid[,1:3], 1, which.max)])
  observed <- as.factor(testdata$Subtype)
  conf_mat <- confusionMatrix(predicted,observed)
  #errors=rbind(errors,rf_valid$confusion[,'class.error'])
  sens <- rbind(sens,conf_mat$byClass[,"Sensitivity"])
  spec <- rbind(spec,conf_mat$byClass[,"Specificity"])
}
sens=sens[-1,]
spec=spec[-1,]
sens_mean <- apply(sens,2,mean)
spec_mean <- apply(spec,2,mean)

# validation
traindata = DE_sample
testdata  = vali_t_ann
rf = randomForest(traindata$Subtype ~. , data=traindata[,2:1236],ntree=1000,importance=TRUE,proximity=TRUE)
pred_valid <- as.data.frame(predict(rf, testdata[,2:1236], type = "prob"))
predicted <- as.factor(names(pred_valid)[1:3][apply(pred_valid[,1:3], 1, which.max)])
observed <- as.factor(testdata$subgroup)
conf_mat <- confusionMatrix(predicted,observed)

require(ggplot2)
require(reshape2)
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

###################### END of Iteration ###############################

Predict_validation <- as.data.frame(predict(Modelon_disCH, vali_t_ann[,2:1236], type = "prob"))

Predict_validation$predict <- names(Predict_validation)[1:3][apply(Predict_validation[,1:3], 1, which.max)]
Predict_validation$observed <- vali_t_ann$subgroup

Predict_validation$predict <- names(Predict_validation)[1:3][apply(Predict_validation[,1:3], 1, which.max)]
Predict_validation$observed <- vali_t_ann$subgroup
colnames(Predict_validation) =c("DUX4","NHHeH","PHlike","predict","observed")

### Sensitivity and specificity
RFPredictiction <- prediction(Predict_validation[,2] , vali_t_ann$subgroup)

ss <- performance(Predict_validation,)



## Plotting ROC curve
roc.dux <- roc(ifelse(Predict_validation$observed=="DUX4","DUX4","non-DUX4"), as.numeric(Predict_validation$DUX4))
plot(roc.dux, col = "gray60")

# others
roc.Phlike <- roc(ifelse(Predict_validation$observed=="Ph-like", "Ph-like", "non-Phlike"), as.numeric(Predict_validation$PHlike))
roc.NHeH <- roc(ifelse(Predict_validation$observed=="NH-HeH", "NH-HeH", "non-NHHeH"), as.numeric(Predict_validation$NHHeH))
lines(roc.Phlike, col = "blue")
lines(roc.NHeH, col = "red")


matplot(data.frame(roc.dux $sensitivities, roc.dux $specificities), x = roc.dux $thresholds, type='l', xlab = 'threshold', ylab='TPR, TNR')
legend('bottomright', legend=c('TPR', 'TNR'), lty=1:2, col=1:2)
matplot(data.frame(roc.Phlike $sensitivities, roc.Phlike $specificities), x = roc.Phlike $thresholds, type='l', xlab = 'threshold', ylab='TPR, TNR')
legend('bottomright', legend=c('TPR', 'TNR'), lty=1:2, col=1:2)

label <- ifelse(Predict_validation > 0.57, 'DUX4','Ph-like')
confusionMatrix(label, Predict_validation$observed, positive = 'versicolor')
accuracyData <- confusionMatrix(Predict_validation[,5],vali_t_ann$subgroup)
accuracyData

