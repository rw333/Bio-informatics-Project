setwd("/Users/richardwen/Desktop/informatics")

#read in all the tumor file names 
inf_anti_files=file("inf_anti_files.txt",open="r")
anti_file_names=readLines(inf_anti_files)
close(inf_anti_files)

#reading in a viral genome expressions and return the data table 
get_anti_exp = function(anti_files)
{
  anti_exp = list()
  i = 0
  for (f_name in anti_files)
  {
    i=i+1
    inf_anti = as.data.frame(read.csv(paste("inf_anti_data/",f_name,sep=""),header = TRUE))
    colnames(inf_anti)=unlist(inf_anti[1,])
    inf_anti = inf_anti[-1,]
    anti_exp[i] = list(inf_anti$Description)
  }
  return (anti_exp)
}

#functions starts here, need to calculate score
get_match_score = function(mutations,gene_exps)
{
  #keep track of number of matches
  matches = rep(0,length(gene_exps))
  i = 1
  #go through each antigen in the file 
  for (anti in gene_exps)
  {
    total = 0
    for (mut in mutations)
    {
      #compare mut to this specific specific antigen 
      score = sum(grepl(mut,anti))
      if (score>0)
        total = 1
    }
    matches[i]=total
    i=i+1
  }
  return(matches)
}
#read in tumor data
tumor_samples = read.csv("tumor_data.csv",header = T)

#read in all the antigen expressions
anti_gene_exps = get_anti_exp(anti_file_names)

tetras=colnames(tumor_samples)[3:length(colnames(tumor_samples))]

count_score=function(r)
{
  muts = tetras[r==1]
  return(get_match_score(muts,anti_gene_exps))
}

#results = apply(tumor_samples[,-(1:2)], 1, count_score)

results2 = apply(tumor_samples[,-(1:2)], 1, count_score)


lb_scores = as.data.frame(results[tumor_samples[,2]==1])
nb_scores = as.data.frame(results[tumor_samples[,2]==0])
lb_scores[,2]=1
nb_scores[,2]=0
colnames(lb_scores)=c("Matches","Benefit")
colnames(nb_scores)=c("Matches","Benefit")
allscores = rbind(lb_scores,nb_scores)

model = wilcox.test(Matches~Benefit,data = allscores,exact=FALSE)
print(model$p.value)

#create the bigger table
lb_scores = as.data.frame(results2[tumor_samples[,2]==1,])
nb_scores = as.data.frame(results2[tumor_samples[,2]==0,])
lb_scores$Benefit=1
nb_scores$Benefit=0
final_dat = rbind(lb_scores,nb_scores)
antigens = c()
for (v in anti_file_names)
{
  antigens = c(antigens,substr(v,1,nchar(v)-4))
}
colnames(final_dat)=c(antigens,"Benefit")

#let the logistic regression begin
library("glmnet")
X=data.matrix(final_dat[,1:59])
Y=final_dat[,60]
model = cv.glmnet(X,Y,family = "binomial",nfolds=20)
plot(model)
plot(model$glmnet.fit,"lambda")
lam = exp(-2.7) 
Y = as.factor(Y)
mod_l1 = glmnet(X,Y,family ="binomial",lambda = lam)
pred_l1 = predict(mod_l1,X,type="response")
roc_curve_l1 = roc.curve(pred_l1[Y==0],pred_l1[Y==1],curve = T)
roc_curve_l1_2 = roc(Y,pred_l1) 
plot(roc_curve_l1)
plot(roc_curve_l1_2)

#do the final train/test modeling
set.seed(123)
smp_size = floor(.80*nrow(final_dat))
roc_area = c()
pr_area = c()
ba = c()
for (i in 1:40)
{
  train_ind = sample(seq_len(nrow(final_dat)), size = smp_size) 
  X_train = X[train_ind,]
  Y_train = Y[train_ind]
  X_test = X[-train_ind,]
  Y_test = Y[-train_ind]
  mod_l1_2 = glmnet(X_train,Y_train,family ="binomial",lambda = lam)
  pred_l1_2 = predict(mod_l1_2,X_train,type="response")
  pred_test = predict(mod_l1_2,X_test,type="response")
  
  roc_curve = roc.curve(pred_test[Y_test==0],pred_test[Y_test==1],curve = T)
  pr_curve = pr.curve(pred_test[Y_test==0],pred_test[Y_test==1],curve = T)
  
  roc_curve_2 = roc(Y_test,pred_test)
  pr_area = c(pr_area,pr_curve$auc.integral)
  roc_area = c(roc_area,roc_curve_2$auc)
  ba = c(ba, max(roc_curve_2$sensitivities+roc_curve_2$specificities)/2)
}
print(mean(roc_area))
print(mean(pr_area))
print(mean(ba))

#decision tree
library(rpart)
#new_dat = final_dat[,c(vip_v,60)]
mod_tree = rpart(Benefit~.,data=final_dat,method="class")
plot(mod_tree, uniform=TRUE, main="Classification Tree for Paitents")
text(mod_tree, use.n=TRUE, all=TRUE, cex=.8)
pred_tree = predict(mod_tree,as.data.frame(X))
roc_tree = roc(as.numeric(Y),pred_tree[,2])
plot(roc_tree)
table(Y,pred_tree[,2]>=.476) 

#svm
library(e1071)
mod_svm = svm(X,as.factor(as.numeric(Y)),cost = 10000)
pred_svm = predict(mod_svm,X)
roc_svm = roc(as.factor(as.numeric(Y)),pred_svm)
plot(roc_svm)
table(as.numeric(Y),pred_svm)
