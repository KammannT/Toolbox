# Donor-to-Donor expression variation (DEV)** metric adopted from: Strunz et al., Chronic hepatitis C virus infection irreversibly impacts human natural killer cell repertoire diversity, Nat Commun. 2018 Jun 11;9(1):2275. doi: 10.1038/s41467-018-04685-9.<br>
# DEV gives variation of protein expression between replicates (= study subjects or donors), e.g. answering the question:
# "How much does the expression of a molecule differ between multiple samples?"
# The average DEV for a set of parameters (= proteins) can be compared between different variables, e.g. stimulations, tissues, timepoints and answers the question:
# "Is the variability of proteins similar between different levels of the variable (e.g. treatment groups)?*

# List of packages to check and install if necessary
packages <- c('tidyverse', 'reshape2', 'see', 'ggrepel')
suppressMessages(for(pack in packages){if(!require(pack, character.only=TRUE)){
  install.packages(p)
  library(p, character.only=TRUE)}})

donors   <- paste0('D', 1:6) # donor variable for Donor-to-Donor expression variation analysis
variable <- c('A', 'B', 'C') # variable for more complex dataset, e.g. different stimulations or parameters.

DEVdata <- data.frame( #dummy data frame
  Donor = sample(donors, length(donors)*length(variable), replace = TRUE),
  Variable = sample(variable, length(donors)*length(variable), replace = TRUE),
  ProteinA = rnorm(length(donors)*length(variable), mean=80, sd=10),
  ProteinB = rnorm(length(donors)*length(variable), mean=40, sd=10),
  ProteinC = rnorm(length(donors)*length(variable), mean=20, sd=15),
  ProteinD = rnorm(length(donors)*length(variable), mean=20, sd=20),
  ProteinE = rnorm(length(donors)*length(variable), mean=60, sd=60))

DEVdata <- melt(DEVdata, variable.name='Protein', value.name='Expression')

# drop NA values and perform calculation for each parameter (=protein)
DEV_data <- DEVdata %>% drop_na() %>% group_by(Protein) %>%
  # calculate median for each protein across replicates and variable levels
  mutate(Median=median(Expression)) %>%
  # normalize each value to the global median
  mutate(Norm=Expression/Median) %>% ungroup %>%
  group_by(Variable, Protein) %>%
  # calculate standard deviation of normalized protein expression
  mutate(DEV=sd(Norm)) %>% ungroup
print(DEV_data)

DEV_data <- DEV_data %>% select(-Donor, -Expression, -Norm, -Median) %>% distinct() # remove duplicate rows from donors
print(head(DEV_data))

ggplot(data=DEV_data, aes(x=Variable, y=DEV, fill=Variable))+
  geom_violindot(dots_size=1, show.legend=FALSE)+
  ggrepel::geom_label_repel(data=subset(DEV_data, DEV>1), aes(label=Protein), label.size=0.2, fill='white')+
  labs(y='Donor-to-donor expression variation (DEV)', x='Variable (e.g. subset or treatment group)')+
  theme_classic()+
  scale_y_continuous(limits=c(0, (max(DEV_data$DEV)+max(DEV_data$DEV)*0.1)), n.breaks=3, expand=c(0,0))
