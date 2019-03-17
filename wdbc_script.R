
# check to see if the required packages are installed and if not, install them
if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(caretEnsemble)) install.packages("caretEnsemble", dependencies = TRUE,
                                             repos = "http://cran.us.r-project.org")
if(!require(captioner)) install.packages("captioner", repos = "http://cran.us.r-project.org")
if(!require(corrplot)) install.packages("corrplot", repos = "http://cran.us.r-project.org")

# # read the data file from Google Drive
# wdbc_data <- read.csv(file = "https://drive.google.com/uc?export=download&id=1suhP_ATOyeiQVMGmWyEudF5JZNFnmVTC",
#                       header = FALSE, sep = ",")

# read the data file from the GitHub repository
wdbc_data <- read.csv(file = "https://raw.githubusercontent.com/lschoch/wdbc/master/wdbc.data.csv",
                      header = FALSE, sep = ",")

# add column names
names(wdbc_data) <- c("ID", "diagnosis", "radius", "texture", "perimeter", "area", "smoothness",
                      "compactness", "concavity", "concave_points", "symmetry", "fractal_dimension",
                      "radius_se", "texture_se", "perimeter_se", "area_se", "smoothness_se", "compactness_se",
                      "concavity_se", "concave_points_se", "symmetry_se", "fractal_dimension_se",
                      "radius_worst", "texture_worst", "perimeter_worst", "area_worst", "smoothness_worst",
                      "compactness_worst", "concavity_worst", "concave_points_worst", "symmetry_worst",
                      "fractal_dimension_worst")
# reorder diagnosis levels so that "M" will be the positive value
levels(wdbc_data$diagnosis) = c("M", "B")

# barplot of diagnosis
wdbc_data %>% ggplot(aes(diagnosis)) + geom_bar(fill="light green") + coord_cartesian(ylim = c(0, 400)) + geom_text(stat='count',aes(label=..count..),vjust=-1)

# correlation plot
corrplot(cor(wdbc_data[ ,3:32]), method = "square")

# create the principal components
prin_comp <- prcomp(wdbc_data[ ,3:32], scale=TRUE)
# compute the standard deviation of each principal component
std_dev <- prin_comp$sdev
#compute the variance of each principal component
var <- std_dev^2
#compute the proportion of variance 
prop_var <- var/sum(var)

#cumulative scree plot
plot(cumsum(prop_var), xlab = "Principal Component",
ylab = "Cumulative Variance",
type = "b")

# partition the data, 80% as train_set and 20% as test_set
set.seed(2019)
partition_index <- createDataPartition(y = wdbc_data$diagnosis, times = 1, p = 0.2, list = FALSE)
train_set <- wdbc_data[-partition_index, ]
test_set <- wdbc_data[partition_index, ]

# create the target and predictors datasets
y <- train_set$diagnosis # target set
x <- train_set[ ,3:32] # predictors set

# create a list of models
model_list <- c('lda', 'rpart', 'knn', 'svmRadial', 'Rborist','nnet')

# train the list of models using the caretList function
control <- trainControl(method="repeatedcv", 
number=30,
index=createResample(y, 30),
repeats=3, 
savePredictions='final', 
classProbs=TRUE)
set.seed(2019)
models <- caretList(x,y,
trControl=control,
methodList=model_list,
preProcess = c('center','scale'))
# get results by resampling the models
results <- resamples(models)

# train the list of models using pca
set.seed(2019)
models_pca <- caretList(x,y,
trControl=control,
methodList=model_list,
preProcess = c('center','scale','pca'))
# get results by resampling the models
results_pca <- resamples(models_pca)

# create an object to contain the accuracy columns of the results object and round to 4 digits
res <- round(results$values[ ,seq(2,12,2)], digits = 4)
# rename the columns
names(res) <- model_list
# create an object to contain the accuracy columns of the results_pca object and round to 4 digits
res_pca <- round(results_pca$values[ ,seq(2,12,2)], digits = 4)
# rename the columns
names(res_pca) <- model_list

# obtain the average accuracy and sd for each model that was computed without pca
res_new <- data.frame(round(colMeans(res), digits = 4)) %>%  # calculate avg accuracy for each model
rownames_to_column("model") %>% # create a column for model
mutate('sd_nopca' = round(apply(res, 2, sd), digits = 4)) # add a column for standard deviation
names(res_new) <- c("model", "avg_nopca", "sd_nopca") # rename the columns

# obtain the average accuracy and sd for each model that was computed with pca
res_pca_new <- data.frame(round(colMeans(res_pca), digits = 4)) %>%  # calculate avg accuracy for each model
rownames_to_column("model") %>% # create a column for model
mutate('sd_pca' = round(apply(res_pca, 2, sd), digits = 4)) # add a column for standard deviation
names(res_pca_new) <- c("model", "avg_pca", "sd_pca") # rename the columns

# combine res_new and res_pca_new (building a table to be printed by the kable function)
res_combined <- res_new %>% left_join(res_pca_new, by = 'model')

# rename columns of res_pca so they won't match column names of res
names(res_pca) <- c("lda_pca", "rpart_pca", "knn_pca", "svmRadial_pca", "Rborist_pca", "nnet_pca")
# combine res and res_pca to calculate p values
res_p <- cbind(res, res_pca)
#calculate the vector of p values
dat <- seq(1, 6) # to select colums in res_p for t.test
p_vals <- sapply(dat, function(d) {
  pv <- t.test(res_p[ , d], res_p[ , d+6], paired = TRUE)$p.value 
  return(round(pv, digits = 3))
})
# add p values as a column in res_combined
res_combined <- res_combined %>% mutate("p_val" = p_vals)
# select the columns for the table to be printed by the kable function
res_combined_kable <- res_combined %>% select(-sd_nopca, -sd_pca)
# create character vector to name columns in printed table
col_names <- c("model", "accuracy-no pca", "accuracy-pca", "p-value")

# evaluate effect size (Cohen's d) for svmRadial and Rborist (p-values were < 0.05)
# create function to calculate Cohen's d
cohenD <- function(x,y) {
  lx <- length(x) - 1
  ly <- length(y) - 1
  md <- abs(mean(x) - mean(y))
  csd <- lx*var(x) + ly*var(y)
  csd <- csd/((lx + ly))
  csd <- sqrt(csd)
  return(md/csd)
}
# calculate effect size
svm_cohend <- round(cohenD(res$svmRadial, res_pca$svmRadial_pca), digits = 3)
Rbor_cohend <- round(cohenD(res$Rborist, res_pca$Rborist_pca), digits = 3)

# create an ensemble of the models using the caretEnsemble function 
stackControl <- trainControl(method="repeatedcv", number=30, repeats=3, savePredictions=TRUE, classProbs=TRUE)
set.seed(2019)
ensemble <- caretEnsemble(models, metric="Accuracy", preProcess = c('center','scale'), 
                          trControl=stackControl)
# use the ensemble to predict diagnoses for the test_set
pred <- predict(ensemble, newdata = test_set[ ,3:32])
# check accuracy of test_set predictions
cm <- confusionMatrix(pred, test_set$diagnosis)

# create an ensemble of the models that were trained using pca 
ensemble_pca <- caretEnsemble(models_pca, metric="Accuracy", preProcess = c('center','scale'), 
                              trControl=stackControl)
# use the ensemble to predict diagnoses for the test_set
pred_pca <- predict(ensemble_pca, newdata = test_set[ ,3:32])
# check accuracy of test_set predictions
cm_pca <- confusionMatrix(pred_pca, test_set$diagnosis)

# create and print table of individual model and ensemble accuracies
res_avg <- data.frame(round(colMeans(res), digits = 4)) # accuracy for each model = average of resamples
# add ensemble accuracy
res_avg <- rbind(res_avg, "ensemble" = round(ensemble$error$Accuracy, digits = 4)) 
res_avg_pca <- data.frame(round(colMeans(res_pca), digits = 4)) # avg accuracy for models with pca
# add ensemble accuracy
res_avg_pca <- rbind(res_avg_pca, "ensemble" = round(ensemble_pca$error$Accuracy, digits = 4))
# combine res_avg and res_avg_pca
res_avg_combined <- cbind(res_avg, res_avg_pca)

# create table of ensemble accuracy results, with and without pca
ens_acc_train <- round(ensemble$error$Accuracy, digits = 4) # accuracy for train set, no pca
ens_acc_test <- round(cm$overall[[1]], digits = 4) # accuracy for test set, no pca
ens_pca_acc_train <- round(ensemble_pca$error$Accuracy, digits = 4) # accuracy for train set, with pca
ens_pca_acc_test <- round(cm_pca$overall[[1]], digits = 4) # accuracy for test set, with pca
# create matrix of accuracy results
ens_res <- matrix(c(ens_acc_train,ens_acc_test,ens_pca_acc_train,ens_pca_acc_test), ncol = 2)
colnames(ens_res) <- c("no pca", "pca")
rownames(ens_res) <- c("accuracy-train","accuracy-test")

######################################################################################################
#   accuracies: 
#    training data, 
#    with and without pca
knitr::kable(res_avg_combined, col.names = c("accuracy-no pca", "accuracy-pca"), align = c("c","c"))
#
######################################################################################################
#   accuracies: 
#    ensembles 
knitr::kable(as.table(ens_res), align = c("c", "c"))
#
######################################################################################################