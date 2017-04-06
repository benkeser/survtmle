
# library(devtools); library(roxygen2)
# setwd("~/Dropbox/R/")
# document("survtmle")
# build("survtmle")
# install("survtmle")
# library(survtmle)


# library(foreach); library(doParallel); library(parallel); library(caret)
# library(rpart); library(randomForest)

# # make toy data
# library(dplyr)
# set.seed(6401^2)
# data(rv144)
# rv144 <- rv144[1:5000,]
# "%ni%" <- Negate("%in%")

# rv144_v2 <- rv144 %>%
#   dplyr::mutate(
#     ftime = I(ftime),
#     ftype0 = I(ftype),
#     ftype1 = sample(ftype),
#     ftype2 = sample(ftype),
#     ftype3 = sample(ftype),
#     vax = I(vax),
#     male = I(male),
#     year04 = I(year04),
#     year05 = I(year05),
#     medRisk = I(medRisk),
#     highRisk = I(highRisk),
#     medAge = I(medAge),
#     highAge = I(highAge)
#   ) %>%
#  dplyr::select(which(colnames(.) %ni% "ftype"))

# rm(rv144)

# # function inputs
# folds = 2
# listPkgs <- c("caret", "survtmle", "rpart", "randomForest")
# libSL <- c("SL.glm", "SL.mean")
# sites <- list("ftype0", "ftype1")
# baseline <- c("highRisk")

# survtmle_adaptive(data = rv144_v2,
#                               ftime = "ftime",
#                               vax = "vax",
#                               covars = baseline,
#                               method = "hazard",
#                               t0 = 4, 
#                               glm.ctime = "t + highRisk",
#                               glm.ftime = "t + highRisk",
#                               glm.trt = "1",
#                               sites = sites,
#                               nfolds = 2,
#                               parallel = TRUE,
#                               cv_tmle = TRUE)
