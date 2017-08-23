
# preliminaries
library(microbenchmark)

# set seed and constants
set.seed(341796)
t_0 <- 15

## get correct version of `survtmle`
if ("survtmle" %in% installed.packages()) {
    remove.packages("survtmle")
}
suppressMessages(devtools::install_github("benkeser/survtmle", ref = "speedglm"))
library(survtmle)

# functions for this simulation
get.ftimeForm <- function(trt, site){
	form <- "-1"
	for(i in unique(trt)){
		for(s in unique(site)){
			form <- c(form, 
			  paste0("I(trt==",i,"& site == ",s," & t==",
			         unique(ftime[ftype>0 & trt==i & site == s]),")",
			         collapse="+"))
		}
	}
	return(paste(form,collapse="+"))
}

get.ctimeForm <- function(trt, site){
	form <- "-1"
	for(i in unique(trt)){
		for(s in unique(site)){
			form <- c(form, 
			  paste0("I(trt==",i,"& site == ",s," & t==",
			         unique(ftime[ftype==0 & trt==i & site == s]),")",
			         collapse="+"))
		}
	}
	return(paste(form,collapse="+"))
}

# simulate data
n <- 100
trt <- rbinom(n, 1, 0.5)

# e.g., study site
adjustVars <- data.frame(site = (rbinom(n,1,0.5) + 1))
ftime <- round(1 + runif(n, 1, 350) - trt + adjustVars$site)
ftype <- round(runif(n, 0, 1))

glm.ftime <- get.ftimeForm(trt = trt, site = adjustVars$site)
glm.ctime <- get.ctimeForm(trt = trt, site = adjustVars$site)

system.time(
    fit <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
                    glm.trt = "1", glm.ftime = glm.ftime, glm.ctime = glm.ctime,
                    method = "hazard", t0 = t_0)
)

suppressMessages(
    m4 <- microbenchmark(unit = "s",
        fit <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
                        glm.trt = "1", glm.ftime = glm.ftime, glm.ctime = glm.ctime,
                        method = "hazard", t0 = t_0)
    )
)
summary(m4)

# simulate data
n <- 1000
trt <- rbinom(n, 1, 0.5)

# e.g., study site
adjustVars <- data.frame(site = (rbinom(n,1,0.5) + 1))
ftime <- round(1 + runif(n, 1, 350) - trt + adjustVars$site)
ftype <- round(runif(n, 0, 1))

glm.ftime <- get.ftimeForm(trt = trt, site = adjustVars$site)
glm.ctime <- get.ctimeForm(trt = trt, site = adjustVars$site)

system.time(
    fit <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
                    glm.trt = "1", glm.ftime = glm.ftime, glm.ctime = glm.ctime,
                    method = "hazard", t0 = t_0)
)

suppressMessages(
    m5 <- microbenchmark(unit = "s",
        fit <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
                        glm.trt = "1", glm.ftime = glm.ftime, glm.ctime = glm.ctime,
                        method = "hazard", t0 = t_0)
    )
)
summary(m5)

# simulate data
#n <- 5000
#trt <- rbinom(n, 1, 0.5)

# e.g., study site
#adjustVars <- data.frame(site = (rbinom(n,1,0.5) + 1))
#ftime <- round(1 + runif(n, 1, 350) - trt + adjustVars$site)
#ftype <- round(runif(n, 0, 1))

#glm.ftime <- get.ftimeForm(trt = trt, site = adjustVars$site)
#glm.ctime <- get.ctimeForm(trt = trt, site = adjustVars$site)

#system.time(
#    fit <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#                    glm.trt = "1", glm.ftime = glm.ftime, glm.ctime = glm.ctime,
#                    method = "hazard", t0 = t_0)
#)

#m6 <- microbenchmark(unit = "s",
#    fit <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
#                    glm.trt = "1", glm.ftime = glm.ftime, glm.ctime = glm.ctime,
#                    method = "hazard", t0 = t_0)
#)
#summary(m6)
