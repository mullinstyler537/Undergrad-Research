#First Beginner R Project for Ecology

# Installing the packages tidyverse and vegan
install.packages("tidyverse")
install.packages("vegan")

#Loading the packages into the current session
library(tidyverse)
library(vegan)

#Load the built in datasets
data(varespec)
data(varechem)

#Calculate the Shannon diversity for each of the 24 sites
shannon_div <- diversity(varespec, index = "shannon")

#Take that previous result and add it as a column inside the soil dataset
varechem$Biodiversity <- shannon_div

#Check the dataset to ensure the "Biodiversity" column was successfully added to the end
head(varechem)

# Build the plot by each layer
ggplot(data = varechem, aes(x = Al, y = Biodiversity)) +
  geom_point(color = "darkgreen", size = 3, alpha = 0.7) + # Adds data points
  geom_smooth(method = "lm", color = "darkred", se = TRUE) + # Adds a linear trendline with a gray error ribbon
  labs(
    title = "Impact of Soil Aluminum Levels on Understory Biodiversity",
    subtitle = "Data source: vegan package (Pine forest lichen pastures)",
    x = "Aluminum Concentration (mg/kg)",
    y = "Shannon Diversity Index (H')"
  ) +
  theme_minimal() # Makes the background clean and white instead of gray
	
#Run a linear regression model
#Syntax: lm(Dependent_Variable ~ Independent_Variable, data = Your_Dataframe)
eco_model <- lm(Biodiversity ~ Al, data = varechem)

#Printing the statistical summary
summary(eco_model)
