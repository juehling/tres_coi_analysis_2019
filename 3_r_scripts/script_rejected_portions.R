


bet1 <- betareg(percent_aquatic_ra ~ age*exp_treat + site_box_year, data = aquatic, link = "logit")


bet2 <- betareg(fCover.m$covert ~ fCover.m[, 1] | fCover.m[, 1], link = "logit")
bet3 <- betareg(1 - fCover.m$covert ~ fCover.m[, 1] - 1 | fCover.m[, 1], link = "log")


m1 <- gamlss(percent_aquatic_ra ~ age + exp_treat + age:exp_treat + re(random = ~1|site_box_year), family= BEINF, data = aquatic, trace = FALSE)



#############
# PLOT: Aquatic across age (nestling vs. adult) and sites
#############

p <- ggplot(aquatic_prov) +
  geom_boxplot((aes(x=age, y=percent_aquatic_ra, fill=age))) + 
  facet_wrap(~ site, ncol = 2) +
  xlab("Age") +
  ylab("Percent Aquatic Insects in Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))

ggsave(here("3_r_scripts/lh_age_site.png"), p, width = 7, height = 7, device = "png")

#############
# PLOT: Aquatic across age (nestling days 6, 12, and 15, adult) and sites
#############

p <- ggplot(aquatic_prov) +
  geom_boxplot((aes(x=age_nestad, y=percent_aquatic_ra, fill=age_nestad))) + 
  facet_wrap(~ site, ncol = 2) +
  xlab("Age") +
  ylab("Percent Aquatic Insects in Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))

# set order for ages for neatness in plots
age_order = c("6", "12", "15", "Adult")
p$data$age_nestad <- as.character(p$data$age_nestad)
p$data$age_nestad <- factor(p$data$age_nestad, levels=age_order)

ggsave(here("3_r_scripts/lh_age_days_site.png"), p, width = 10, height = 8, device = "png")

#############
# PLOT: Aquatic vs. terrestrial across days
#############

p <- ggplot(aquatic, aes(x=as.numeric(cap_doy), y=percent_aquatic_ra)) +
  geom_point() + geom_smooth(method = "loess") +
  facet_wrap(~ age + site, ncol = 2) +
  xlab("Day of Year") +
  ylab("Percent aquatic insects in diet") +
  theme_classic() +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))

ggsave(here("3_r_scripts/lh_doy.png"), p, width = 10, height = 8, device = "png")

#############
# PLOT: Aquatic vs. terrestrial across days, separated by adults vs. nestlings
#############

p <- ggplot(aquatic, aes(x=as.numeric(cap_doy), y=percent_aquatic_ra)) +
  geom_point() + geom_smooth(method = "loess") +
  facet_wrap(~ age, ncol = 1) +
  xlab("Day of Year") +
  ylab("Percent aquatic insects in diet") +
  theme_classic() +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))

ggsave(here("3_r_scripts/lh_doy_age.png"), p3, width = 10, height = 8, device = "png")

################################################################################
# Calculate percent aquatic in each nest ----



# Filter data just to adult and nestling captures during provisioning
# (i.e. take out first adult captures during provisioning)
aquatic_ad <- aquatic_wild[aquatic_wild$cap_num != "1" & aquatic_wild$age == "Adult" ,]
aquatic_wild_nest <- aquatic_wild[aquatic_wild$age == "Nestling" ,]
aquatic_wild_prov <- rbind(aquatic_wild_ad, aquatic_wild_nest)

# Pull out nestlings to average across nest
aquatic_wild_prov_n <- aquatic_wild_prov[aquatic_wild_prov$age == "Nestling" ,]

# Pull out adults
aquatic_wild_prov_a <- aquatic_wild_prov[aquatic_wild_prov$age == "Adult" ,]

######## Relative abundance calculations

# Create an empty dataframe to store relative abundance information for each nest
nests <- unique(aquatic_wild_prov_n$site_box_year)
aquatic_wild_prov_npooled <- data.frame(matrix(NA, ncol = 2, nrow = length(nests)))
aquatic_wild_prov_npooled <- data.frame()

# Calculate average aquatic percentage with relative abundance for the nestlings in each nest
for (i in 1:length(nests)){
  nest <- nests[i] # identify nest
  list <- aquatic_wild_prov_n[aquatic_wild_prov_n$site_box_year == nest ,] # pull out all nestlings from that nest
  mean_aq <- mean(list$percent_aquatic_ra) # average % aquatic insects
  aquatic_wild_prov_npooled[i,1] <- nest # save the nest name
  aquatic_wild_prov_npooled[i,2] <- mean_aq # save the total aquatic relative abundance
}

# Name columns
names(aquatic_wild_prov_npooled)[names(aquatic_wild_prov_npooled) == "V1"] <- "site_box_year"
names(aquatic_wild_prov_npooled)[names(aquatic_wild_prov_npooled) == "V2"] <- "percent_aquatic_ra_nestlings_mean"

# Combine adults with pooled nestling percentages
aquatic_wild_prov_a <- merge(aquatic_wild_prov_a, aquatic_wild_prov_npooled, by = "site_box_year")

######## Presence/absence calculations

# Create an empty dataframe to store presence/absence information for each nest
nests <- unique(aquatic_wild_prov_n$site_box_year)
aquatic_wild_prov_npooled <- data.frame(matrix(NA, ncol = 2, nrow = length(nests)))
aquatic_wild_prov_npooled <- data.frame()

# Calculate average aquatic percentage with presence/absence for the nestlings in each nest
for (i in 1:length(nests)){
  nest <- nests[i] # identify nest
  list <- aquatic_wild_prov_n[aquatic_wild_prov_n$site_box_year == nest ,] # pull out all nestlings from that nest
  mean_aq <- mean(list$percent_aquatic_pa) # average % aquatic insects using presence/absence
  aquatic_wild_prov_npooled[i,1] <- nest # save the nest name
  aquatic_wild_prov_npooled[i,2] <- mean_aq # save the total aquatic percent
}

# Name columns
names(aquatic_wild_prov_npooled)[names(aquatic_wild_prov_npooled) == "V1"] <- "site_box_year"
names(aquatic_wild_prov_npooled)[names(aquatic_wild_prov_npooled) == "V2"] <- "percent_aquatic_pa_nestlings_mean"

# Combine adults with pooled nestling percentages
aquatic_wild_prov_a <- merge(aquatic_wild_prov_a, aquatic_wild_prov_npooled, by = "site_box_year")


################################################################################
# Nestling condition vs. percent aquatic

aquatic_day12nest <- aquatic[aquatic$age_nestad == "12" ,]

p <- ggplot(aquatic_day12nest, aes(x=percent_aquatic_ra, y=mass)) + geom_point()
p

mod_nestcon <- lm(mass ~ percent_aquatic_ra*site, data = aquatic_day12nest)
summary(mod_nestcon)

p <- ggplot(aquatic_day12nest, aes(x=percent_aquatic_pa, y=mass)) + geom_point()

p

################################################################################
# Plot Patterns: Ordinations/Beta Diversity ----

# Plot some ordinations using relative abundance and presence/absence
ord_data_ra <- ordinate(coi_ra2, method = "PCoA", distance = "bray")
ord_data_pa <- ordinate(coi_pa, method = "PCoA", distance = "bray")

# Site -- relative abundance
p1_ra <- plot_ordination(coi_ra2, ord_data_ra, color = "site", title = "Bray-Curtis PCoA") + 
  geom_point(size = 2, alpha = .4) + theme_classic() +
  stat_ellipse(aes(group = site), level = 0.9) +
  theme(axis.title = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14))
ggsave(here("3_r_scripts/site_ordinate_ra.png"), p1_ra, width = 7, height = 5.4, device = "png")

# Site -- presence/absence
p1_pa <- plot_ordination(coi_pa, ord_data_pa, color = "site", title = "Bray-Curtis PCoA") + 
  geom_point(size = 2, alpha = .4) + theme_classic() +
  stat_ellipse(aes(group = site), level = 0.9) +
  theme(axis.title = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14))
ggsave(here("3_r_scripts/site_ordinate_pa.png"), p1_pa, width = 7, height = 5.4, device = "png")

# Age -- relative abundance
p2_ra <- plot_ordination(coi_ra2, ord_data_ra, color = "age", title = "Bray-Curtis PCoA") + 
  geom_point(size = 2, alpha = .4) + theme_classic() +
  stat_ellipse(aes(group = age), level = 0.9)
ggsave(here("3_r_scripts/age_ordinate_ra.png"), p2_ra, width = 7, height = 5.4, device = "png")

# Age -- presence absence
p2_pa <- plot_ordination(coi_pa, ord_data_pa, color = "age", title = "Bray-Curtis PCoA") + 
  geom_point(size = 2, alpha = .4) + theme_classic() +
  stat_ellipse(aes(group = age), level = 0.9)
ggsave(here("3_r_scripts/age_ordinate_pa.png"), p2_pa, width = 7, height = 5.4, device = "png")

################################################################################
# Captive nestling gut passage samples ----

coi_pa_cap <- subset_samples(coi_pa, site == "Gut_passage")

# All genera
p <- plot_bar(coi_pa_cap, "family") + theme_classic() + theme(axis.text.x = element_text(angle = 90, size = 12))
p <- p + facet_grid(~ age.1)
age_order = c("6", "7", "8", "9", "10", "11", "12")
p$data$age.1 <- as.character(p$data$age.1)
p$data$age.1 <- factor(p$data$age.1, levels=age_order)

ggsave(here("3_r_scripts/family_bar.png"), width = 10, height = 4.5, device = "png")

# Try just certain days for readability

coi_pa_cap_12 <- subset_samples(coi_pa_cap, age.1 == "12")
p <- plot_bar(coi_pa_cap_12, "family") + theme_classic() + theme(axis.text.x = element_text(angle = 90, size = 12))

coi_pa_cap_11 <- subset_samples(coi_pa_cap, age.1 == "11")
p <- plot_bar(coi_pa_cap_11, "family") + theme_classic() + theme(axis.text.x = element_text(angle = 90, size = 12))


################################################################################
# Plot Patterns: Select objects to plot ----

# Note that we will use coi_ra2 for plots for relative abundance, and 
# coi_pa for plots for presence/absence.

################################################################################
# Plot Patterns: Richness ----

# Richness by sample type. See help for many more options of different alpha metrics
# Cannot use coi_ra2 here because plot_richness only accepts integers

# First, look at adult r

coi_adults_comp_F <- subset_samples(coi_pa, cap_num == "3" & sex == "F")
coi_adults_comp_M <- subset_samples(coi_pa, sex == "M")
coi_adults_comp <- merge_phyloseq(coi_adults_comp_F, coi_adults_comp_M)

p <- plot_richness(coi_pa, x = "age_nestad", measures = c("Observed", "InvSimpson", "Shannon")) + 
  geom_boxplot(alpha = 0.3, aes(fill = age_nestad))
# set order for ages for neatness in plots 
age_order = c("6", "12", "15", "Adult")
p$data$age_nestad <- as.character(p$data$age_nestad)
p$data$age_nestad <- factor(p$data$age_nestad, levels=age_order)
ggsave(here("3_r_scripts/age_alpha.png"), p, width = 10, height = 4, device = "png")

# Richness by age (nestling days and adult years) and site
p <- plot_richness(coi_pa, x = "age_nestad", measures = c("Observed")) +
  geom_boxplot(alpha = 0.3, aes(fill = age_nestad)) + facet_wrap(~ site) +
  xlab("Age") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, size = 12)) + 
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(axis.title = element_text(size = 16)) +
  theme(strip.text = element_text(size = 16)) +
  guides(fill=guide_legend(title="Age"))
# set order for sites for neatness in plots 
site_order = c("Turkey_Hill", "Unit_4", "Unit_1", "Unit_2")
p$data$site <- as.character(p$data$site)
p$data$site <- factor(p$data$site, levels=site_order)
# set order for ages for neatness in plots 
age_order = c("6", "12", "15", "Adult")
p$data$age_nestad <- as.character(p$data$age_nestad)
p$data$age_nestad <- factor(p$data$age_nestad, levels=age_order)
ggsave(here("3_r_scripts/age_alpha_site.png"), p, width = 9, height = 6.5, device = "png")

# Richness by age (just nestling vs. adult) and site

# For a more direct comparison across sites, only use nestlings day 12
coi_pa_comp <- subset_samples(coi_pa, age_nestad == "Adult" | age_nestad == "12")

p <- plot_richness(coi_pa_comp, x = "age", measures = c("Observed")) +
  geom_boxplot(alpha = 0.3, aes(fill = age)) + facet_wrap(~ site) +
  xlab("Age") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, size = 12)) + 
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(axis.title = element_text(size = 16)) +
  theme(strip.text = element_text(size = 16)) +
  guides(fill=guide_legend(title="Age"))
# set order for sites for neatness in plots 
site_order = c("Turkey_Hill", "Unit_4", "Unit_1", "Unit_2")
p$data$site <- as.character(p$data$site)
p$data$site <- factor(p$data$site, levels=site_order)
ggsave(here("3_r_scripts/age_nva_alpha_site.png"), p, width = 7, height = 7, device = "png")

# Richness by adult capture number and site
coi_adults <- subset_samples(coi_pa, cap_num != "")
p <- plot_richness(coi_adults, x = "cap_num", measures = c("Observed")) + 
  geom_boxplot(alpha = 0.3, aes(fill = cap_num)) + facet_wrap(~ site)
ggsave(here("3_r_scripts/capnum_alpha_site.png"), p, width = 9, height = 6.5, device = "png")

# Richness between males and females (use only third capture for females because 
# this was approximate time males were caught)
coi_adults_comp_F <- subset_samples(coi_pa, cap_num == "3" & sex == "F")
coi_adults_comp_M <- subset_samples(coi_pa, sex == "M")
coi_adults_comp <- merge_phyloseq(coi_adults_comp_F, coi_adults_comp_M)
p <- plot_richness(coi_adults_comp, x = "sex", measures = c("Observed")) + 
  geom_boxplot(alpha = 0.3, aes(fill = sex)) + facet_wrap(~ site)
ggsave(here("3_r_scripts/alpha_site_sex.png"), p, width = 9, height = 6.5, device = "png")

# Richness by day of year
# First look at richness with a series of box plots
p <- plot_richness(coi_pa, x = "cap_doy", measures = c("Observed")) + 
  geom_boxplot(alpha = 0.3, aes(fill = age.1)) + facet_wrap(~ site)
# Now look at richness as a loess curve
richness <- data.table(p$data)
richness$cap_doy <- as.numeric(richness$cap_doy)
p <- ggplot(richness, mapping = aes(x = cap_doy, y = value)) + geom_point() + geom_smooth(method = "loess") +
  xlab("Capture day of year") + ylab("Alpha Diversity Measure") + facet_wrap(~ age, ncol = 1) +
  theme_classic() +
  theme(axis.text = element_text(size = 14)) + 
  theme(axis.title = element_text(size = 16)) +
  theme(strip.text = element_text(size = 16))
ggsave(here("3_r_scripts/age_site_richness.png"), p, width = 8.7, height = 6.8, device = "png")

################################################################################
# Model Patterns: Richness ----

# Based on our dataset, it seems like the most reasonable comparison to make is
# between nestling day 12 samples and adults (3rd capture) across sites

coi_pa_comp_div <- data.frame(
  "Observed" = phyloseq::estimate_richness(coi_pa_comp, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(coi_pa_comp, measures = "Shannon"),
  "sampleID" = phyloseq::sample_data(coi_pa_comp)$sampleID,
  "age" = phyloseq::sample_data(coi_pa_comp)$age,
  "site" = phyloseq::sample_data(coi_pa_comp)$site)

mod_adnest_richness <- lm(Observed ~ age + site, data = coi_pa_comp_div)
qqPlot(residuals(mod_adnest_richness))
plot(mod_adnest_richness)
summary(mod_adnest_richness)
# Don't think this model is appropriate, but there does not seem to be a difference here.