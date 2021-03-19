#visualization of the "demo" inputs on a fake landscape

library(tidyverse)
library(cowplot)


# plasticity on viz -------------------------------------------------------
setwd("./inputs_demo/output1616171118")

#create list of the outputs created after running CDMetaPOP
output_files = list.files(pattern = paste("ind*", sep=''), 
                          full.names = TRUE, 
                          recursive = TRUE, 
                          include.dirs = TRUE) %>% 
  map_df(function(x) read_csv(x, col_types = cols(.default = "c")) %>% mutate(filename=paste(dirname((x)),basename(x),sep="/")))


#break up the column with the file name
data_df <- separate(data = output_files, col = filename, into = c('junk', 'junk1', 'year'), sep = "/")

data_df <- separate(data = data_df, col = year, into = c('ind', 'year'), sep = "d")
data_df <- separate(data = data_df, col = year, into = c('year', 'junk2'), sep = ".cs")

#set column data type and remove unneeded columns
data_df$junk <- NULL
data_df$junk1 <- NULL
data_df$junk2 <- NULL
data_df$ID <- NULL
data_df$sex <- NULL
data_df$size <- NULL
data_df$mature <- NULL
data_df$newmature <- NULL
data_df$layeggs <- NULL
data_df$capture <- NULL
data_df$infection <- NULL
data_df$Hindex <- NULL
data_df$Species <- NULL
data_df$recapture <- NULL
data_df$SubPatchID <- NULL
data_df$CDist <- NULL
data_df$ind <- NULL
data_df$ClassFile <- NULL
data_df$L0A0 <-  as.numeric(data_df$L0A0)
data_df$L0A1 <- as.numeric(data_df$L0A1)
data_df$L1A0 <-  as.numeric(data_df$L1A0)
data_df$L1A1 <- as.numeric(data_df$L1A1)
data_df$age <- as.factor(data_df$age)
data_df$year <- as.numeric(data_df$year)
data_df$XCOORD <- as.numeric(data_df$XCOORD)
data_df$YCOORD <- as.numeric(data_df$YCOORD)

#summarise to go from idnividuals to patch statistics. Then calculate the proportion of individuals in a given year based on the total number of individuals in that year

data_df_freq <- data_df %>%
  group_by(PatchID, XCOORD, YCOORD, year, .drop = FALSE) %>%
  summarise(n = n()) %>%
  group_by(year) %>%
  mutate(year_tot = sum(n)) %>%
  mutate(year_prop = n/year_tot)

#filter to only the initial and last year
data_df_plot <- data_df_freq %>%
  filter(year == "5" | year == "45")

data_df_plot$Temperature <- c("Cold", "Cold", "Trigger", "Trigger", "Avoid", "Avoid")
data_df_plot$Temperature <- as.factor(data_df_plot$Temperature)
data_df_plot <- rename(data_df_plot, Population = n)


#create plot for pop count
space1 <- data_df_plot %>% ggplot(aes(x = XCOORD, y = abs(YCOORD))) + 
  geom_point(aes(color = Temperature, size = Population)) +
  facet_wrap(~year) +
  scale_color_manual(values = c("red", "blue", "orange"))+
  xlim(0, 1000) +
  xlab("X") +
  ylab("Y") +
  scale_size(range = c(1,20))+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  theme_bw()



#next, scripts to create the figure with allele frequencies in each location

#first need the population. can do per patch, but note frequencies very high for 2 in hot due to low sample size
data_df_patch_pop <- data_df %>%
  group_by(year, .drop = FALSE) %>%
  summarise(pop = n())

#now to work with plastic allele, counting up
data_L0A0 <- data_df %>%
  group_by(PatchID, XCOORD, YCOORD, year, .drop = FALSE) %>%
  count(L0A0, name = "L0A0_count")

data_L0A0_spread <- data_L0A0 %>%
  spread(L0A0, L0A0_count)

#renaming last couple of columns for when merging
colnames(data_L0A0_spread) <- c("PatchID", "XCOORD", "YCOORD", "year", "0.x", "1.x", "2.x")


data_L0A1 <- data_df %>%
  group_by(PatchID, XCOORD, YCOORD, year, .drop = FALSE) %>%
  count(L0A1, name = "L0A1_count")

data_L0A1_spread <- data_L0A1 %>%
  spread(L0A1, L0A1_count)

colnames(data_L0A1_spread) <- c("PatchID", "XCOORD", "YCOORD","year","0.y", "1.y", "2.y")

data_df_L0 <- left_join(x=data_L0A0_spread, y=data_L0A1_spread)


#next sum together across

data_df_L0$L0A0 <- data_df_L0$'0.x' + data_df_L0$'0.y'
#this will throw an eroor if all switches flipped
data_df_L0$L0A1 <- data_df_L0$'1.x' + data_df_L0$'1.y'
data_df_L0$L0A2 <- data_df_L0$'2.x' + data_df_L0$'2.y'

#remove old rows
data_df_L0$'0.x' <- NULL
data_df_L0$'0.y' <- NULL
data_df_L0$'1.x' <- NULL
data_df_L0$'1.y' <- NULL
data_df_L0$'2.x' <- NULL
data_df_L0$'2.y' <- NULL

#fill NA

data_df_L0[is.na(data_df_L0)] = 0

#add pop total column

data_df_L0$pop <- rep(c(data_df_patch_pop$pop))

data_df_L0_gathered <- gather(data_df_L0, key = "Allele", counts, L0A0:L0A2)

#now to do the same thing with a neutral locus for comparison

data_df_L1 <- data_df %>%
  group_by(PatchID, XCOORD, YCOORD, year, .drop = FALSE) %>%
  summarise_at(vars(L1A0:L1A1), sum, na.rm = TRUE)


#add in the per patch pop
data_df_L1$pop <- rep(data_df_patch_pop$pop)

#gathering the columns of alleles
data_df_L1_gathered <- gather(data_df_L1, key = "Allele", counts, L1A0:L1A1)

#combine, and then calculate proportion

data_df_gen <- rbind(data_df_L1_gathered, data_df_L0_gathered)

#seperate allele and locus values for figure

data_df_gen <- separate(data = data_df_gen, col = Allele, into = c('Locus', 'Allele'), sep = 3)

#change L1 and L0 to names that are descritive
data_df_gen <- data_df_gen %>% 
  mutate(Locus = replace(Locus, Locus == 'L0A', "Plastic")) %>%
  mutate(Locus = replace(Locus, Locus == 'L1A', "Neutral"))

data_df_gen$Locus <- as.factor(data_df_gen$Locus)

data_df_gen$Temperature <- rep(c("Cold", "Cold","Cold","Cold","Cold","Cold","Cold","Cold","Cold","Cold","Cold","Cold", "Trigger", "Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger", "Avoid", "Avoid", "Avoid","Avoid","Avoid","Avoid","Avoid","Avoid","Avoid","Avoid","Avoid","Avoid"), times = 5)
data_df_gen$Temperature <- as.factor(data_df_gen$Temperature)

#show  alleles through time
data_df_gen %>% 
  ggplot(aes(x = year, y = counts)) + 
  geom_line(aes(color = Temperature)) +
  facet_grid(rows = vars(Locus), cols = vars(Allele)) +
  scale_color_manual(values = c("red", "blue", "orange"))+
  xlim(0, 50) +
  xlab("X") +
  ylab("Y") +
  guides(colour = guide_legend(override.aes = list(size=10)))+
  theme_bw()

data_df_gen %>% 
  ggplot(aes(x = year, y = counts)) + 
  geom_line(aes(color = Temperature, linetype = Locus), size = 1.1) +
  facet_wrap(~Allele) +
  scale_color_manual(values = c("red", "blue", "orange"))+
  xlim(0, 50) +
  xlab("X") +
  ylab("Y") +
  guides(colour = guide_legend(override.aes = list(size=10)))+
  theme_bw()


#now to show the patches population growths

data_df_patch_pop <- data_df %>%
  group_by(PatchID, year, .drop = FALSE) %>%
  summarise(pop = n())

data_df_patch_pop$Temperature <- c("Cold", "Cold","Cold","Cold","Cold","Cold","Cold","Cold","Cold","Cold","Cold","Cold", "Trigger", "Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger", "Avoid", "Avoid", "Avoid","Avoid","Avoid","Avoid","Avoid","Avoid","Avoid","Avoid","Avoid","Avoid")
data_df_patch_pop$Temperature <- as.factor(data_df_patch_pop$Temperature)


pop_time_1 <- data_df_patch_pop %>% 
  ggplot(aes(x = year, y = as.numeric(pop))) + 
  geom_line() +
  facet_wrap(~Temperature) +
  ggtitle("Plasticity On")+
  xlab("Year") +
  ylab("Population") +
  theme_bw()



# plasticity off viz ------------------------------------------------------

setwd("../../inputs_demo/no_plast1616196029")
output_files = list.files(pattern = paste("ind*", sep=''), 
                          full.names = TRUE, 
                          recursive = TRUE, 
                          include.dirs = TRUE) %>% 
  map_df(function(x) read_csv(x, col_types = cols(.default = "c")) %>% mutate(filename=paste(dirname((x)),basename(x),sep="/")))


#break up the column with the file name
data_df <- separate(data = output_files, col = filename, into = c('junk', 'junk1', 'year'), sep = "/")

data_df <- separate(data = data_df, col = year, into = c('ind', 'year'), sep = "d")
data_df <- separate(data = data_df, col = year, into = c('year', 'junk2'), sep = ".cs")

#set column data type and remove unneeded columns
data_df$junk <- NULL
data_df$junk1 <- NULL
data_df$junk2 <- NULL
data_df$ID <- NULL
data_df$sex <- NULL
data_df$size <- NULL
data_df$mature <- NULL
data_df$newmature <- NULL
data_df$layeggs <- NULL
data_df$capture <- NULL
data_df$infection <- NULL
data_df$Hindex <- NULL
data_df$Species <- NULL
data_df$recapture <- NULL
data_df$SubPatchID <- NULL
data_df$CDist <- NULL
data_df$ind <- NULL
data_df$ClassFile <- NULL
data_df$L0A0 <-  as.numeric(data_df$L0A0)
data_df$L0A1 <- as.numeric(data_df$L0A1)
data_df$L1A0 <-  as.numeric(data_df$L1A0)
data_df$L1A1 <- as.numeric(data_df$L1A1)
data_df$age <- as.factor(data_df$age)
data_df$year <- as.numeric(data_df$year)
data_df$XCOORD <- as.numeric(data_df$XCOORD)
data_df$YCOORD <- as.numeric(data_df$YCOORD)

#summarise to go from idnividuals to patch statistics. Then calculate the proportion of individuals in a given year based on the total number of individuals in that year

data_df_freq <- data_df %>%
  group_by(PatchID, XCOORD, YCOORD, year, .drop = FALSE) %>%
  summarise(n = n()) %>%
  group_by(year) %>%
  mutate(year_tot = sum(n)) %>%
  mutate(year_prop = n/year_tot)

#filter to only the initial and last year
data_df_plot <- data_df_freq %>%
  filter(year == "5" | year == "45")

data_df_plot$Temperature <- c("Cold", "Cold", "Trigger", "Trigger", "Avoid", "Avoid")
data_df_plot$Temperature <- as.factor(data_df_plot$Temperature)
data_df_plot <- rename(data_df_plot, Population = n)


#create plot for pop count
space1 <- data_df_plot %>% ggplot(aes(x = XCOORD, y = abs(YCOORD))) + 
  geom_point(aes(color = Temperature, size = Population)) +
  facet_wrap(~year) +
  scale_color_manual(values = c("red", "blue", "orange"))+
  xlim(0, 1000) +
  xlab("X") +
  ylab("Y") +
  scale_size(range = c(1,20))+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  theme_bw()



#next, scripts to create the figure with allele frequencies in each location

#first need the population. can do per patch, but note frequencies very high for 2 in hot due to low sample size
data_df_patch_pop <- data_df %>%
  group_by(year, .drop = FALSE) %>%
  summarise(pop = n())

#now to work with plastic allele, counting up
data_L0A0 <- data_df %>%
  group_by(PatchID, XCOORD, YCOORD, year, .drop = FALSE) %>%
  count(L0A0, name = "L0A0_count")

data_L0A0_spread <- data_L0A0 %>%
  spread(L0A0, L0A0_count)

#renaming last couple of columns for when merging
colnames(data_L0A0_spread) <- c("PatchID", "XCOORD", "YCOORD", "year", "0.x", "1.x", "2.x")


data_L0A1 <- data_df %>%
  group_by(PatchID, XCOORD, YCOORD, year, .drop = FALSE) %>%
  count(L0A1, name = "L0A1_count")

data_L0A1_spread <- data_L0A1 %>%
  spread(L0A1, L0A1_count)

colnames(data_L0A1_spread) <- c("PatchID", "XCOORD", "YCOORD","year","0.y", "1.y", "2.y")

data_df_L0 <- left_join(x=data_L0A0_spread, y=data_L0A1_spread)


#next sum together across

data_df_L0$L0A0 <- data_df_L0$'0.x' + data_df_L0$'0.y'
#this will throw an eroor if all switches flipped
data_df_L0$L0A1 <- data_df_L0$'1.x' + data_df_L0$'1.y'
data_df_L0$L0A2 <- data_df_L0$'2.x' + data_df_L0$'2.y'

#remove old rows
data_df_L0$'0.x' <- NULL
data_df_L0$'0.y' <- NULL
data_df_L0$'1.x' <- NULL
data_df_L0$'1.y' <- NULL
data_df_L0$'2.x' <- NULL
data_df_L0$'2.y' <- NULL

#fill NA

data_df_L0[is.na(data_df_L0)] = 0

#add pop total column

data_df_L0$pop <- rep(c(data_df_patch_pop$pop))

data_df_L0_gathered <- gather(data_df_L0, key = "Allele", counts, L0A0:L0A2)

#now to do the same thing with a neutral locus for comparison

data_df_L1 <- data_df %>%
  group_by(PatchID, XCOORD, YCOORD, year, .drop = FALSE) %>%
  summarise_at(vars(L1A0:L1A1), sum, na.rm = TRUE)


#add in the per patch pop
data_df_L1$pop <- rep(data_df_patch_pop$pop)

#gathering the columns of alleles
data_df_L1_gathered <- gather(data_df_L1, key = "Allele", counts, L1A0:L1A1)

#combine, and then calculate proportion

data_df_gen <- rbind(data_df_L1_gathered, data_df_L0_gathered)

#seperate allele and locus values for figure

data_df_gen <- separate(data = data_df_gen, col = Allele, into = c('Locus', 'Allele'), sep = 3)

#change L1 and L0 to names that are descritive
data_df_gen <- data_df_gen %>% 
  mutate(Locus = replace(Locus, Locus == 'L0A', "Plastic")) %>%
  mutate(Locus = replace(Locus, Locus == 'L1A', "Neutral"))

data_df_gen$Locus <- as.factor(data_df_gen$Locus)

data_df_gen$Temperature <- rep(c("Cold", "Cold","Cold","Cold","Cold","Cold","Cold","Cold","Cold","Cold","Cold","Cold", "Trigger", "Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger", "Avoid", "Avoid", "Avoid","Avoid","Avoid","Avoid","Avoid","Avoid","Avoid","Avoid","Avoid","Avoid"), times = 5)
data_df_gen$Temperature <- as.factor(data_df_gen$Temperature)

#show  alleles through time
data_df_gen %>% 
  ggplot(aes(x = year, y = counts)) + 
  geom_line(aes(color = Temperature)) +
  facet_grid(rows = vars(Locus), cols = vars(Allele)) +
  scale_color_manual(values = c("red", "blue", "orange"))+
  xlim(0, 50) +
  xlab("X") +
  ylab("Y") +
  guides(colour = guide_legend(override.aes = list(size=10)))+
  theme_bw()

data_df_gen %>% 
  ggplot(aes(x = year, y = counts)) + 
  geom_line(aes(color = Temperature, linetype = Locus), size = 1.1) +
  facet_wrap(~Allele) +
  scale_color_manual(values = c("red", "blue", "orange"))+
  xlim(0, 50) +
  xlab("X") +
  ylab("Y") +
  guides(colour = guide_legend(override.aes = list(size=10)))+
  theme_bw()


#now to show the patches population growths

data_df_patch_pop <- data_df %>%
  group_by(PatchID, year, .drop = FALSE) %>%
  summarise(pop = n())

data_df_patch_pop$Temperature <- c("Cold", "Cold","Cold","Cold","Cold","Cold","Cold","Cold","Cold","Cold","Cold","Cold", "Trigger", "Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger","Trigger", "Avoid", "Avoid", "Avoid","Avoid","Avoid","Avoid","Avoid","Avoid","Avoid","Avoid","Avoid","Avoid")
data_df_patch_pop$Temperature <- as.factor(data_df_patch_pop$Temperature)

pop_time_2 <- data_df_patch_pop %>% 
  ggplot(aes(x = year, y = as.numeric(pop))) + 
  geom_line() +
  facet_wrap(~Temperature) +
  ggtitle("Plasticity Off")+
  xlab("Year") +
  ylab("Population") +
  theme_bw()

plot_grid(pop_time_1, pop_time_2, labels = c('A', 'B'), ncol = 1)
