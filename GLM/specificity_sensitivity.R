##Predict the model, transfose the data and then use it in next steps
library(magrittr)
library(caret)
library(data.table)
library(glmnet)
library(randomForest)
library(pROC)


All_DE_FPKM = read.table(file="/Users/alvajames/Documents/Sensitivity_JHO/Sensitivity_Journal_comemnts_October_2018/All_DE_FPKM",sep='\t',header=T,stringsAsFactors=FALSE)
All_DE_FPKM$Geneid=NULL
All_DE_FPKM=All_DE_FPKM[1:83]
row.names(All_DE_FPKM) =All_DE_FPKM$Gene_sym
All_DE_FPKM$Gene_sym=NULL
All_DE_FPKM_t=t(All_DE_FPKM)
write.table(All_DE_FPKM_t,file="/Users/alvajames/Documents/Sensitivity_JHO/Sensitivity_Journal_comemnts_October_2018/All_DE_FPKM_transposed",sep='\t')
DE_sample_FPKM=read.table(file="/Users/alvajames/Documents/Sensitivity_JHO/Sensitivity_Journal_comemnts_October_2018/All_DE_FPKM_transposed",sep='\t',header=T,stringsAsFactors=FALSE)
vali_t_ann=read.table(file="/Users/alvajames/Documents/Sensitivity_JHO/Sensitivity_Journal_comemnts_October_2018/Validation_trans.tsv",sep='\t',header=T,stringsAsFactors=FALSE)

All_de            = read.table(file="/Users/alvajames/Documents/Sensitivity_JHO/Sensitivity_Journal_comemnts_October_2018/All_de_nodup",sep='\t',header=T,stringsAsFactors=FALSE)
Types    = read.table(file="/Users/alvajames/Documents/Sensitivity_JHO/Sensitivity_Journal_comemnts_October_2018/Samples",sep='\t',header=T,stringsAsFactors=FALSE)
Types$Relapsed.days = NULL
Types$Age = NULL
Validation_cohort    = read.table(file="/Users/alvajames/Documents/Sensitivity_JHO/Sensitivity_Journal_comemnts_October_2018/BCP_ALL_DE_on_Berlin_count_noDUO",sep='\t',header=T,stringsAsFactors=FALSE)
Validation_cohort$Berlin2492<-NULL
Validation_cohort$Berlin2493<-NULL
Validation_cohort$Berlin2494<-NULL
Info=read.table("/Users/alvajames/Documents/Sensitivity_JHO/Sensitivity_Journal_comemnts_October_2018/Sample_group3_subtype.tsv", header=TRUE, sep="\t",check.names=FALSE)
Annotation   = read.table(file="/Users/alvajames/Documents/Sensitivity_JHO/Sensitivity_Journal_comemnts_October_2018/lncRNAs_annotation",sep='\t',header=T,stringsAsFactors=FALSE)
Validation_cohort_annon=merge(Validation_cohort,Annotation, on='Geneid')
Validation_cohort_annon<-Validation_cohort_annon[,c(146, 1:145)]
Validation_cohort_annon$Geneid<-NULL
validation_class=merge(Validation_cohort_annon,All_de[ , c("Gene_sym", "Subtype")],on='Gene_sm')

row.names(Validation_cohort_annon) <- Validation_cohort_annon$Gene_sym

validation_class=validation_class[names(validation_class) %in% Info$sample_name]

vali_t =t(Validation_cohort_annon)
vali_t=as.data.frame(vali_t)
library(data.table)
setDT(vali_t, keep.rownames = TRUE)
colnames(vali_t)[1] <- "sample_name"
vali_t_ann=merge(vali_t,Info,by='sample_name')
dim(vali_t_ann[,2:1236])
write.table(vali_t_ann,file="/Users/alvajames/Documents/Sensitivity_JHO/Sensitivity_Journal_comemnts_October_2018/Validation_trans",sep='\t')
DE        = read.table(file="/Users/alvajames/Documents/Sensitivity_JHO/Sensitivity_Journal_comemnts_October_2018/Transposed_Allde.tsv",sep='\t',header=T,stringsAsFactors=FALSE)
DE_sample = merge(DE, Types, on=sample_name)
dim(DE_sample) #  82 1237
#mat_scaled = t(apply(DE_sample[,2:1236], 1, function(x) (x - mean(x)) / sd(x)))
DE_sample=DE_sample[!grepl("LH", DE_sample$Subtype),]
DE_sample=DE_sample[!grepl("unassigned", DE_sample$Subtype),]
table(DE_sample$Subtype)

## Perform penalized cross validation using glm

### 10-fold lasso model
setcolorder(vali_t_ann[,2:1236], names(DE_sample[,2:1236]))

traintest=rbind(DE_sample[,2:1236],vali_t_ann[,2:1236])

X = sparse.model.matrix(as.formula(paste("y ~", paste(colnames(DE_sample[,2:1236]), sep = "", collapse =" +" ))), data = traintest)

cv.lassoModel <- cv.glmnet(x=data.matrix(DE_sample[,2:1236]),y=DE_sample$Subtype,standardize=TRUE,alpha=0,nfolds=10,family="multinomial")

##### replace the DE_sample with validation cohort ****
response<-predict(cv.lassoModel, data.matrix(t(Validation_cohort_annon)), s= "lambda.min", type="class")

## Get confusion matrix for Validation and discovry cohort
newX <- model.matrix(~.-DE_sample$Subtype,data=DE_sample[,2:1236])

discovery_conf_matrix=confusionMatrix(table(predict(cv.lassoModel,newx =newX )))



full <- glm(DE_sample$Subtype ~ colnames(DE_sample[,2:1236]), data=DE_sample[,2:1236])

small.lambda.index <- which(cv.lassoModel$lambda == cv.lassoModel$lambda.min)

#Plot variable deviances vs. shrinkage parameter, λ (lambda)
plot(cv.lassoModel)

## getting the best predictors:
idealLambda <- cv.lassoModel$lambda.min
idealLambda1se <- cv.lassoModel$lambda.1se
print(idealLambda); print(idealLambda1se)

#Derive coefficients for each gene to each subtype
co <- coef(cv.lassoModel, s=idealLambda)
co

co.se <- coef(cv.lassoModel, s=idealLambda1se, exact=TRUE)

co.se
###### finding the DUX4 variables genes with least coefficeint
DUX_preditcors <- as.data.frame(as.matrix(co$DUX4))
rownames(co$DUX4)
rownames(co$`Ph-like`)[which(co$`Ph-like`[,1]>0)]

###
finalLasso <- glm(DE_sample$Subtype~., data=DE_sample[,2:1236], family= poisson)
summary(finalLasso)

### Final comparision between model and already classified classes(subtypes)
cbind(response, DE_sample$Subtype)

require(MASS)
fit1_LM <- stepAIC(cv.lassoModel, direction = 'backward')

library(ROCR)

prediction(predict(cv.lassoModel, data.matrix(vali_t_ann[,2:1236]), s= "lambda.min", type="class"),vali_t_ann$subgroup)

pred <- prediction(predict(cv.lassoModel, type="response",newx=newX), DE_sample$Subtype)
ss <- performance(response, "sens", "spec")

############ Apply the CV glm for 1000 times for Discovery cohort #######
N = 1000
sens = spec =as.data.frame(matrix(NA, nrow=1, ncol=3))
for(i in 1:N){
  set.seed(250+i)
  train     = sample(1:dim(DE_sample)[1],2*dim(DE_sample)[1]/3,replace = FALSE)
  traindata = DE_sample[train,]
  Subtype   = factor(traindata$Subtype)
  #traindata= t(apply(traindata[,1:47], 1, function(x) (x - mean(x)) / sd(x)))
  testdata  = DE_sample[-train,]
  CV        = cv.glmnet(x=data.matrix(traindata[,2:1236]),y=Subtype,standardize=TRUE,alpha=0,nfolds=3,family="multinomial")
  response_predict=predict(CV, data.matrix(testdata[,2:1236]),type="response")
  predicted <- as.factor(colnames(response_predict)[apply(response_predict, 1, which.max)])
  observed  = as.factor(testdata$Subtype)
  conf_matD = confusionMatrix(predicted,observed)
  #errors=rbind(errors,rf_valid$confusion[,'class.error'])
  sens = rbind(sens,conf_matD$byClass[,"Sensitivity"])
  spec = rbind(spec,conf_matD$byClass[,"Specificity"])
}
sens=sens[-1,]
spec=spec[-1,]
sens_mean <- apply(sens,2,mean)
spec_mean <- apply(spec,2,mean)

# validation
traindata = DE_sample
testdata  = vali_t_ann
CVV = cv.glmnet(x=data.matrix(traindata[,2:1236]),y=traindata$Subtype,standardize=TRUE,alpha=0,nfolds=10,family="multinomial")
observed <- as.factor(testdata$subgroup)
predicted <- as.data.frame(predict(CVV, as.matrix(testdata[,2:1236]), s="lambda.min", type="class"))
predicted <- factor(predicted[,1], levels=levels(observed))
conf_mat_validation <- confusionMatrix(predicted,observed)

# validation - binomial
traindata = DE_sample
subtype2 = rep("DUX4",length(traindata$Subtype))
subtype2[-which(traindata$Subtype=="DUX4")] = "PHnh"
subtype2 = as.factor(subtype2)
testdata  = vali_t_ann
CVV = cv.glmnet(x=data.matrix(traindata[,2:1236]),y=subtype2,standardize=TRUE,alpha=0,nfolds=10,family="binomial")
observed <- testdata$subgroup
observed[-which(testdata$subgroup=="DUX4")] = "PHnh"
observed <- as.factor(observed)
predicted <- as.data.frame(predict(CVV, as.matrix(testdata[,2:1236]), s="lambda.min", type="class"))
conf_mat_validation <- confusionMatrix(predicted[,1],observed)


## Train and test on validation data
dim(vali_t_ann)
#Select the trianing data for 70%
47* 0.7 
# 864.5
## Pull random 865 rows from 1235 rows
train     = sample(1:47,33,replace = FALSE)
traindata = vali_t_ann[train,]
Subtype   = factor(traindata$subgroup)
traindata = t(apply(traindata[,2:1236], 1, function(x) (x - mean(x)) / sd(x)))

testdata  = vali_t_ann[-train,]
dim(testdata)
CVV = cv.glmnet(x=traindata,y=Subtype,standardize=TRUE,alpha=0,nfolds=10,family="multinomial")

pred_valid <- as.data.frame(predict(rf, testdata[,2:1236], type = "prob"))
predicted <- as.factor(names(pred_valid)[1:3][apply(pred_valid[,1:3], 1, which.max)])
observed <- as.factor(testdata$subgroup)
conf_mat_validationonly <- confusionMatrix(predicted,observed)
