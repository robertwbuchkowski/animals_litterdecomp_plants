# Running the stable isotope mixing models for my plant data

#Load library
library(MixSIAR)

# Get custom plotting function:
source("plot_data_iso_one_custom.R")


#*************************************************#
#           Various Mixing Models below           #
#                                                 #
#*************************************************#

# NOTE: Mixing models require significant time to run

# Forbs Run only fractional 15N -------------------------------------------
mix <- load_mix_data(filename="forb_consumers_fraction.csv", 
                     iso_names=c("d15N"), 
                     factors="LitterType", 
                     fac_random=c(FALSE), 
                     fac_nested=c(FALSE), 
                     cont_effects="Isopods")

# Load the source data
source <- load_source_data(filename="forb_source.csv",
                           source_factors="LitterType", 
                           conc_dep=F, 
                           data_type="means", 
                           mix)

# Load the discrimination/TDF data
discr <- load_discr_data(filename="forb_discrimination.csv", mix)

colnames(source$S_MU) <-c("d15N", "LitterType")
colnames(source$S_SIG) <-c("d15N", "LitterType")

plot_data_one_iso(filename="isospace_plot_ch2", plot_save_pdf=F, plot_save_png=FALSE, mix,source,discr)

calc_area(source=source,mix=mix,discr=discr)

plot_prior(alpha.prior=1,source, filename="prior_plot_forbs")

model_filename <- "MixSIAR_forbs_fractional.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

#run <- list(chainLength=200000, burn=150000, thin=50, chains=3, calcDIC=TRUE)

jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics_forbs_fractional",
                       sup_post = T,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density_forbs_fractional",
                       sup_pairs = T,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot_forbs_fractional",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot_forbs_fractional",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics_forbs_fractional",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

jags_forbs_fractional <- run_model(run="long", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

output_JAGS(jags_forbs_fractional, mix, source, output_options)

# Forbs discrimination soil Run only fractional 15N (FINAL MODEL USED) -------------------------------------------
mix <- load_mix_data(filename="forb_consumers_fraction.csv", 
                     iso_names=c("d15N"), 
                     factors="LitterType", 
                     fac_random=c(FALSE), 
                     fac_nested=c(FALSE), 
                     cont_effects="Isopods")

# Load the source data
source <- load_source_data(filename="forb_source_discrim.csv",
                           source_factors="LitterType", 
                           conc_dep=F, 
                           data_type="means", 
                           mix)

# Load the discrimination/TDF data
discr <- load_discr_data(filename="forb_discrimination.csv", mix)

colnames(source$S_MU) <-c("d15N", "LitterType")
colnames(source$S_SIG) <-c("d15N", "LitterType")

plot_data_one_iso_custom(filename="forbsdiscr_plot_ch2", plot_save_pdf=T, plot_save_png=FALSE, mix,source,discr)

#calc_area(source=source,mix=mix,discr=discr)

plot_prior(alpha.prior=1,source, filename="prior_plot_forbs_discr")

model_filename <- "MixSIAR_forbsdiscr_fractional.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

#run <- list(chainLength=200000, burn=150000, thin=50, chains=3, calcDIC=TRUE)

jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics_forbsdiscr_fractional",
                       sup_post = T,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density_forbsdiscr_fractional",
                       sup_pairs = T,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot_forbsdiscr_fractional",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot_forbsdiscr_fractional",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics_forbsdiscr_fractional",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

jags_forbsdiscr_fractional <- run_model(run="long", mix, source, discr, model_filename, 
                                   alpha.prior = 1, resid_err, process_err)

output_JAGS(jags_forbsdiscr_fractional, mix, source, output_options)

# Grass Run only fractional 15N -------------------------------------------
mix <- load_mix_data(filename="grass_consumers_fraction.csv", 
                     iso_names=c("d15N"), 
                     factors="LitterType", 
                     fac_random=c(FALSE), 
                     fac_nested=c(FALSE), 
                     cont_effects="Isopods")

# Load the source data
source <- load_source_data(filename="grass_source.csv",
                           source_factors="LitterType", 
                           conc_dep=F, 
                           data_type="means", 
                           mix)

# Load the discrimination/TDF data
discr <- load_discr_data(filename="grass_discrimination.csv", mix)

colnames(source$S_MU) <-c("d15N", "LitterType")
colnames(source$S_SIG) <-c("d15N", "LitterType")

plot_data_one_iso_custom(filename="grass_isospace_plot", plot_save_pdf=T, plot_save_png=FALSE, mix,source,discr)

plot_prior(alpha.prior=1,source, filename="prior_plot_grass")

model_filename <- "MixSIAR_grass_fractional.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

#run <- list(chainLength=200000, burn=150000, thin=50, chains=3, calcDIC=TRUE)

jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics_grass_fractional",
                       sup_post = T,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density_grass_fractional",
                       sup_pairs = T,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot_grass_fractional",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot_grass_fractional",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics_grass_fractional",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

jags_grass_fractional <- run_model(run="long", mix, source, discr, model_filename, 
                                   alpha.prior = 1, resid_err, process_err)

output_JAGS(jags_grass_fractional, mix, source, output_options)


# Grass Discrim Run only fractional 15N -------------------------------------------
mix <- load_mix_data(filename="grass_consumers_fraction.csv", 
                     iso_names=c("d15N"), 
                     factors="LitterType", 
                     fac_random=c(FALSE), 
                     fac_nested=c(FALSE), 
                     cont_effects="Isopods")

# Load the source data
source <- load_source_data(filename="grass_source_discrim.csv",
                           source_factors="LitterType", 
                           conc_dep=F, 
                           data_type="means", 
                           mix)

# Load the discrimination/TDF data
discr <- load_discr_data(filename="grass_discrimination.csv", mix)

colnames(source$S_MU) <-c("d15N", "LitterType")
colnames(source$S_SIG) <-c("d15N", "LitterType")

plot_data_one_iso_custom(filename="grassdiscrim_isospace_plot", plot_save_pdf=T, plot_save_png=FALSE, mix,source,discr)

plot_prior(alpha.prior=1,source, filename="prior_plot_grassdiscrim")

model_filename <- "MixSIAR_grassdiscrim_fractional.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

#run <- list(chainLength=200000, burn=150000, thin=50, chains=3, calcDIC=TRUE)

jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics_grassdiscrim_fractional",
                       sup_post = T,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density_grassdiscrim_fractional",
                       sup_pairs = T,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot_grassdiscrim_fractional",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot_grassdiscrim_fractional",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics_grassdiscrim_fractional",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

jags_grassdiscrim_fractional <- run_model(run="long", mix, source, discr, model_filename, 
                                   alpha.prior = 1, resid_err, process_err)

output_JAGS(jags_grassdiscrim_fractional, mix, source, output_options)


# Grass Discrim Run only fractional 15N, no Isopods (FINAL MODEL USED) -------------------------------------------
mix <- load_mix_data(filename="grass_consumers_fraction.csv", 
                     iso_names=c("d15N"), 
                     factors="LitterType", 
                     fac_random=c(FALSE), 
                     fac_nested=c(FALSE), 
                     cont_effects=NULL)

# Load the source data
source <- load_source_data(filename="grass_source_discrim.csv",
                           source_factors="LitterType", 
                           conc_dep=F, 
                           data_type="means", 
                           mix)

# Load the discrimination/TDF data
discr <- load_discr_data(filename="grass_discrimination.csv", mix)

colnames(source$S_MU) <-c("d15N", "LitterType")
colnames(source$S_SIG) <-c("d15N", "LitterType")

plot_data_one_iso_custom(filename="grassdiscrimn_oiso_isospace_plot", plot_save_pdf=T, plot_save_png=FALSE, mix,source,discr)

plot_prior(alpha.prior=1,source, filename="prior_plot_grassdiscrim_noiso")

model_filename <- "MixSIAR_grassdiscrim_nosio_fractional.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

#run <- list(chainLength=200000, burn=150000, thin=50, chains=3, calcDIC=TRUE)

jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics_grassdiscrim_noiso_fractional",
                       sup_post = T,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density_grassdiscrim_noiso_fractional",
                       sup_pairs = T,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot_grassdiscrim_noiso_fractional",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot_grassdiscrim_noiso_fractional",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics_grassdiscrim_noiso_fractional",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

jags_grassdiscrim_noiso_fractional <- run_model(run="long", mix, source, discr, model_filename, 
                                          alpha.prior = 1, resid_err, process_err)

output_JAGS(jags_grassdiscrim_noiso_fractional, mix, source, output_options)

#Plot plant outputs
par(mfrow=c(1,2))

library(R2jags)

attach.jags(jags_forbsdiscr_fractional)
postlitter = data.frame(p.fac1[,,1])

postlitter2 = data.frame(cbind(rep(c(2,3,4,1),each=3000),c(postlitter[,1],postlitter[,2],postlitter[,3],postlitter[,4])))
names(postlitter2) = c("Litter", "Value")

boxplot(postlitter2[,2]~postlitter2[,1], ylab="Proportion Litter Nitrogen",main="Forbs",
        names=c("Control", "Herbivory", "Fertilization",
                "Both"))
text(x=0.5, y=0.15, label="a")

attach.jags(jags_grassdiscrim_fractional)
postlitter = data.frame(p.fac1[,,1])

gpostlitter2 = data.frame(cbind(rep(c(2,3,4,1),each=3000),c(postlitter[,1],postlitter[,2],postlitter[,3],postlitter[,4])))
names(gpostlitter2) = c("Litter", "Value")


boxplot(gpostlitter2[,2]~gpostlitter2[,1], ylab="Proportion Litter Nitrogen",main="Grass",
        names=c("Control", "Herbivory", "Fertilization",
                "Both"))
text(x=0.5, y=0.95, label="b")
text(x=1, y=0.6, label="#")

# Isopod Run only fractional 15N -------------------------------------------
mix <- load_mix_data(filename="isopod_consumers_fraction.csv", 
                     iso_names=c("d15N"), 
                     factors="LitterType", 
                     fac_random=c(FALSE), 
                     fac_nested=c(FALSE),
                     cont_effects=NULL)

# Load the source data
source <- load_source_data(filename="isopod_source.csv",
                           source_factors="LitterType", 
                           conc_dep=F, 
                           data_type="means", 
                           mix)

# Load the discrimination/TDF data
discr <- load_discr_data(filename="isopod_discrimination.csv", mix)

colnames(source$S_MU) <-c("d15N", "LitterType")
colnames(source$S_SIG) <-c("d15N", "LitterType")

plot_data_one_iso_custom(filename="isopod_isospace_plot", plot_save_pdf=T, plot_save_png=FALSE, mix,source,discr)

plot_prior(alpha.prior=1,source, filename="prior_plot_isopod")

model_filename <- "MixSIAR_isopod_fractional.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

#run <- list(chainLength=200000, burn=150000, thin=50, chains=3, calcDIC=TRUE)

jags.1 <- run_model(run="test", mix, source, discr, model_filename, 
                    alpha.prior = 1, resid_err, process_err)

output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics_isopod_fractional",
                       sup_post = T,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density_isopod_fractional",
                       sup_pairs = T,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot_isopod_fractional",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = F,
                       plot_xy_name = "xy_plot_isopod_fractional",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics_isopod_fractional",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

jags_isopod_fractional <- run_model(run="long", mix, source, discr, model_filename, 
                                   alpha.prior = 1, resid_err, process_err)

output_JAGS(jags_isopod_fractional, mix, source, output_options)


