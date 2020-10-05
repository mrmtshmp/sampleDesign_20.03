#' Comparison of changes in RA activity between treatment groups, with adjusting sampling (assignment) unbalance.
#' PI: Dr Endoh
#' date created: 2020/09/28
#' ---

setwd('./src/R')

# models ------------

vars.1 <- c(
  "Age","Duration_month","ACPA",
  "MTX.init","PSL.init",
  "Tender.init",
  "Swollen.init",
  "PR_VAS.init",
  "MR_VAS.init"
  )

vars <- "vars.1"

fml.ps_model <- sprintf(
  "treatment ~ %s", 
  paste(eval(parse(text = vars)),collapse = "+")
  ) 

# subroutines-------

dir.sub <- './sub'
Bibtex <- FALSE

list.files.dir.sub <- list.files(path = dir.sub)

for(i in 1:length(list.files.dir.sub))
  source(sprintf("%s/%s", dir.sub, list.files.dir.sub[i]))

sink('sessionInfo.')
sessionInfo()
sink()

# Load data ---------------------------------------------------------------

load(file = sprintf("%s/%s", dir.ADS, fn.ADS))

data$tratment <- factor(data$treatment)

## Construct a table
tabUnmatched <-
  CreateTableOne(
    vars = eval(parse(text=vars)), 
    strata = "treatment", 
    data = data, 
    test = FALSE)
## Show table with SMD
sink(
  sprintf("%s/TableOne.txt", dir.output)
  )
print(tabUnmatched, smd = TRUE)
sink()

quartz(family = 'Arial',type = 'pdf',file = sprintf("%s/cov_rel.pairwise.pdf", dir.output))
GGally::ggpairs(data[,eval(parse(text=vars))])
dev.off()

# Propensity score model ---------------

propensityScoreModel <-
  glm(
    gsub("treatment","factor\\(treatment\\)",fml.ps_model),
    family  = binomial(link = "logit"),
    data = data, na.action = na.exclude
    )

# propensityScoreModel <-
#   brglm::brglm(
#     gsub("treatment","factor\\(treatment\\)",fml.ps_model),
#     family = binomial(link="logit"),
#     data = data %>% data.frame(),
#     na.action = na.exclude
#   )

# propensityScoreModel <- 
#   rpart::rpart(
#     fml.ps_model,
#     data = data
#     )

data.propensityScores <-
  data %>% data.frame()
  
data.propensityScores$propensity_score <- 
  predict(
    propensityScoreModel,
    # data=data, 
    type= "response",
    # type= "prob",
    na.action = na.exclude()
    ) %>% unlist()

data.propensityScores_IPW <-
  IPW_weights(
    treatment = data.propensityScores$treatment,
    propensity_score = data.propensityScores$propensity_score,
    dat = data.propensityScores
    )

res.roc.propensity_score <- roc(
  response = 
    as.factor(data.propensityScores_IPW$treatment),
  predictor = 
    data.propensityScores_IPW$propensity_score
  )

# Plot distribution of propensity score and weighted counts. ----------------
#'
ggdata.propensityScores <-
  ggplot(
    data =
      data.propensityScores_IPW,
    aes(
      x = propensity_score
      )
    )

ggdata.propensityScores.weighted_count <-
  ggplot(
    data =
      data.propensityScores_IPW %>%
      pivot_longer(
        cols =
          c(starts_with("w_at")),
        values_to = "weight",
        names_to = "target_pop"
      ) %>%
      dplyr::filter(target_pop=="w_ato"),
    aes(
      x = propensity_score,
      weight = weight
    )
  )

pdf(
  file = 
    sprintf("%s/IPWcount.pdf", dir.output),
  width = 21
  )
plot(
  ggdata.propensityScores + 
    geom_density(
      aes(
        fill = 
          as.factor(treatment)
        ),
      bw="SJ",
#      binwidth = FD,
      alpha=0.5,
      position="identity"
      ) +
  geom_point(
    aes(
      y=as.numeric(treatment),
      x=propensity_score,
      color=as.factor(treatment),
      size=2
      )
    ) +
    theme_bw()
  )
plot(
  ggdata.propensityScores.weighted_count + 
    geom_density(
      aes(fill=as.factor(treatment)),
      bw="SJ",
#      binwidth = FD,
      alpha=0.5#,
#      position="dodge"
      ) + 
    facet_grid(~target_pop) + 
    theme_bw()
    )
plot(
  res.roc.propensity_score,
  print.thres=TRUE 
  )
legend(
  x = 0.6, y=0.5,cex = 0.7, 
  # lwd = c(2,2,0), lty = 1:2,
  legend = c(
    sprintf(
      "AUC = %s (0.95CI: %s, %s)",
      round(auc(res.roc.propensity_score),3),
      round(ci(auc(res.roc.propensity_score))[1],3),
      round(ci(auc(res.roc.propensity_score))[3],3)
    )
  ),
  bty = "n"
)
dev.off()


# Balancing assessment ----------------------------------------------------

#' The standardized mean differences between the two patients' groups of *treatment* 
#' Reference:
#' https://cran.r-project.org/web/packages/tableone/vignettes/smd.html
#' [accessed:2020/08/31])

res.svydesign.weighted <- 
  survey::svydesign(
    ids = ~ 1, 
    data = data.propensityScores_IPW %>% dplyr::filter(!is.na(propensity_score)),
    weights = ~ w_ato
    )
  
## Construct a table
tabWeighted.weighted <- 
  svyCreateTableOne(
    vars = eval(parse(text = vars)),
    strata = "treatment", 
    data = res.svydesign.weighted, 
    test = FALSE
    )
## Show table with SMD

sink(
  sprintf("%s/TableOne._weighted.txt",dir.output)
  )
print(tabWeighted.weighted, smd = TRUE)
sink()


## Construct a data frame containing variable name and SMD from all methods
dataPlot <- data.frame(
  variable  = rownames(ExtractSmd(tabUnmatched)),
  rawdata = as.numeric(ExtractSmd(tabUnmatched)),
  weighted_data = as.numeric(ExtractSmd(tabWeighted.weighted))
  )

## Create long-format data for ggplot2
dataPlotMelt <-
  melt(
    data          = dataPlot,
    id.vars       = c("variable"),
    variable.name = "Method",
    value.name    = "SMD"
    )

## Order variable names by magnitude of SMD
varNames <- unique(
  as.character(dataPlot$variable)[
    order(dataPlot$rawdata)
    ]
  )


## Order factor levels in the same order
dataPlotMelt$variable <- 
  factor(
    dataPlotMelt$variable,
    levels = varNames
    )

## Plot using ggplot2

quartz(
  family = "Arial",type = "pdf",
  file =   sprintf("%s/smd.pdf",dir.output)
  )
ggplot(
  data = dataPlotMelt,
  mapping = aes(x = variable, y = SMD, group = Method, color = Method)
  ) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.1, color = "black", size = 0.1) +
  coord_flip() +
  theme_bw() + theme(legend.key = element_blank())
dev.off()

# Make long data ----------------------------------------------------------

ADS <- data %>%
  pivot_longer(
    cols =ends_with("oc_init")
  ) %>%
  mutate(
    name=gsub("(.+)\\.oc_init","\\1",name)
  ) %>%
  dplyr::rename(
    value.pre=value
  ) %>%
  left_join(
    data %>%
      pivot_longer(
        cols =ends_with("oc_post")
      ) %>%
      mutate(
        name=gsub("(.+)\\.oc_post","\\1",name)
      ) %>%
      dplyr::select(ID, name, value),
    by = c("ID","name") 
  ) %>%
  mutate(
    value.change = value.pre - value
  ) %>%
  left_join(
    data.propensityScores_IPW[,grep("^ID$|^w_.+",colnames(data.propensityScores_IPW))]
  )


## Weighted model (use svyglm function)

#' From help(svyglm):
#' Note
#' svyglm always returns 'model-robust' standard errors;
#'  the Horvitz-Thompson-type standard errors used 
#'  everywhere in the survey package are a generalisation
#'   of the model-robust 'sandwich' estimators. 
#'  In particular, a quasi-Poisson svyglm will return
#'   correct standard errors for relative risk regression models.

list.res.glm <-
  dlply(
    ADS,
    .(name),
    function(D){
      res.svydesign.weighted <- 
        survey::svydesign(
          ids = ~ 1, 
          data = D %>% dplyr::filter(!is.na(w_ato)),
          weights = ~ w_ato
        )
      res.glmWeighted <-
        svyglm(
          formula = value.change ~ treatment,
          # family  = quasibinomial(link = "logit"),                  
          design  = res.svydesign.weighted
        )
      }
    )


res.svydesign.weighted.binomial_outcome <- 
  survey::svydesign(
    ids = ~ 1, 
    data = data.propensityScores_IPW %>% dplyr::filter(!is.na(w_ato)),
    weights = ~ w_ato
  )
res.glmWeighted.binomial_outcome <-
  svyglm(
    formula = as.factor(outcome) ~ treatment,
    family  = quasibinomial(link = "logit"),                  
    design  = res.svydesign.weighted.binomial_outcome
    )

sink(sprintf("%s/res.glm.txt", dir.output))
print(summary(res.glmWeighted.binomial_outcome))
print(exp(res.glmWeighted.binomial_outcome$coefficients))
print(exp(confint(res.glmWeighted.binomial_outcome)))

lapply(list.res.glm, function(L){list(summary(L),confint(L))})
sink()

#```