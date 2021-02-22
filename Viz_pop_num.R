#script for figures for the presentation
#do after Pop_occupancy_3_site.R


library(tidyverse)
library(cowplot)
library(ggpubr)
library(plyr)

pop <- data_df %>%
  group_by(year, loc, run) %>%
  tally()

pop$year <- as.numeric(pop$year)

pop$loc <- as.factor(pop$loc)
pop$loc <- revalue(pop$loc, c("cool"="Keithley/Mann", "dry"="Dry Creek", "desert" = "Jack's Creeks"))

pop$run <- as.factor(pop$run)
pop$run <- revalue(pop$run, c("null"="Null", "plast"="Plasticity: Habitat Selection", "sel" = "Local Adaptation: Thermal Tolerance", "sel-plast"="Combination"))


pop$Model <- pop$run

ggplot(pop, aes(x=year, y=n, color = Model)) +
  geom_line(size=2) +
  labs(x = "Time (years)", y = "Population Size") +
  facet_grid(~loc) +
  scale_color_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2") +
  xlim(20,100) +
  scale_y_continuous(labels = scales::comma, limits = c(0, 100000)) +
  theme_classic(base_size = 16)

ggsave(filename = "pop1.png", plot = last_plot(), dpi = 400, width = 10, height = 8, units = "in")




#convert chr to factor
data_df$L0A0 <- as.numeric(data_df$L0A0)
data_df$L0A1 <- as.numeric(data_df$L0A1)
data_df$L1A0 <- as.numeric(data_df$L1A0)
data_df$L1A1 <- as.numeric(data_df$L1A1)
data_df$L2A0 <- as.numeric(data_df$L2A0)
data_df$L2A1 <- as.numeric(data_df$L2A1)

#going to first see what the total should be
Population <- data_df %>%
  group_by(year, loc, run) %>%
  tally()

Population$year <- as.numeric(Population$year)

#going to make seperate dataframes here for each locus
#sum gives the total number of alleles in the population
data_df_L0 <- data_df %>%
  group_by(year, loc, run) %>%
  summarise_at(vars(L0A0:L0A1), sum, na.rm = TRUE)


#plastic region a bit different, need the zeroes

#first count the homozygotes
data_df %>%
  group_by(year, L1A0, L0A1, loc, run) %>%
  tally()

data_L1A1 <- data_df %>%
  group_by(year, loc, run) %>%
  dplyr::count(L1A1, name = "L1A1_count")

data_L1A1_spread <- data_L1A1 %>%
  spread(L1A1, L1A1_count)

data_L1A0 <- data_df %>%
  group_by(year, loc, run) %>%
  dplyr::count(L1A0, name = "L1A0_count")

data_L1A0_spread <- data_L1A0 %>%
  spread(L1A0, L1A0_count)


data_df_L1 <- merge(x=data_L1A0_spread, y=data_L1A1_spread, by=c("year","loc","run"), all = TRUE)

head(data_df_L1)

#next sum together across

data_df_L1$L1A0 <- data_df_L1$'0.x' + data_df_L1$'0.y'
#this will throw an error if all switches flipped
data_df_L1$L1A1 <- data_df_L1$'1.x' + data_df_L1$'1.y'
data_df_L1$L1A2 <- data_df_L1$'2.x' + data_df_L1$'2.y'

#remove old rows
data_df_L1$'0.x' <- NULL
data_df_L1$'0.y' <- NULL
data_df_L1$'1.x' <- NULL
data_df_L1$'1.y' <- NULL
data_df_L1$'2.x' <- NULL
data_df_L1$'2.y' <- NULL




data_df_L2 <- data_df %>%
  group_by(year, loc, run) %>%
  summarise_at(vars(L2A0:L2A1), sum, na.rm = TRUE)

#bringing together for comparisons
df_all <- merge(x = data_df_L0, y = data_df_L1, by=c("year", "loc", "run"), all = TRUE)
df_all <- merge(x=df_all, y=data_df_L2, by=c("year", "loc", "run"), all = TRUE)

df_all$year <- as.numeric(df_all$year)



#gathering the columns of alleles
df_all_plot <- gather(df_all, key = "Allele", value, L0A0:L2A1)

#going split allele column to include locus for better graphing

df_all_plot <- separate(data = df_all_plot, col = Allele, into = c('Locus', 'Allele'), sep = 2)

#plasticity locus only plot
df_plastic_plot <- subset(df_all_plot, df_all_plot$Locus != "L0" & df_all_plot$Locus != "L2")

head(df_plastic_plot)



#graph



#rename loci then plot
df_plot_final <- df_all_plot %>% 
  mutate(Locus = replace(Locus, Locus == 'L0', "Adaptive")) %>%
  mutate(Locus = replace(Locus, Locus == 'L2', "Neutral")) %>%
  mutate(Locus = replace(Locus, Locus == 'L1', "Plastic"))

df_plot_final$Locus <- as.factor(df_plot_final$Locus)

df_plot_final$Locus <- factor(df_plot_final$Locus,levels(df_plot_final$Locus)[c(1,3,2)])



df_plot_final <- merge(x=df_plot_final, y=Population, by=c("year", "loc", "run"), all = TRUE)
df_plot_final$value <- df_plot_final$value / (2*df_plot_final$n)


df_plot_final$run <- revalue(df_plot_final$run, c("null"="Null", "plast"="Plasticity", "sel" = "Local Adaptation", "sel-plast"="Combination"))

for (i in unique(df_plot_final$loc)) {
  loc = i 
  
  data_plot <- df_plot_final %>%
    filter(loc == i)

my_plot <- ggplot(data=data_plot, aes(x = year, y = value, color=as.factor(Allele), fill=as.factor(Allele), shape=as.factor(Allele))) +
  geom_point(size=4, alpha=0.55) +
  geom_line(size=1.25)+
  facet_grid(rows = vars(run), cols = vars(Locus)) +
  labs(x = "Time (years)", y = "Allele Proportion", color="Allele", fill="Allele", shape="Allele") +
  ggtitle(paste("Location",i,sep = " ")) +
  scale_color_brewer(palette = "Dark2")+
  theme_bw(base_size=14)

ggsave(my_plot, filename=paste("frequencies",i,".png",sep=""),dpi = 400, width = 10, height = 8, units = "in")
}


#as opposed to comparing frequencies in a site, can also compare across sites by making df separated by runs

for (i in unique(df_plot_final$run)) {
  run = i 
  
  data_plot <- df_plot_final %>%
    filter(run == i)
  
  my_plot <- ggplot(data=data_plot, aes(x = year, y = value, color=as.factor(Allele), fill=as.factor(Allele), shape=as.factor(Allele))) +
    geom_point(size=4, alpha=0.55) +
    geom_line(size=1.25)+
    facet_grid(rows = vars(loc), cols = vars(Locus)) +
    labs(x = "Time (years)", y = "Allele Proportion", color="Allele", fill="Allele", shape="Allele") +
    ggtitle(paste("Run",i,sep = " ")) +
    scale_color_brewer(palette = "Dark2")+
    theme_bw(base_size=14)
  
  ggsave(my_plot, filename=paste("site_compare",i,".png",sep=""),dpi = 400, width = 10, height = 8, units = "in")
}


#as opposed to comparing frequencies in a site, can also compare across sites by making df separated by runs