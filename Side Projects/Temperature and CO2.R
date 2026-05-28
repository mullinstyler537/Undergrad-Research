#Project: Tracking 50 Years of Climate Change (Temperature and CO2 Trends)

library(tidyverse)
library(lubridate)

#Import the offical NOAA Mauna Loa CO2 monthly data
co2_data <- read_csv("https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_mlo.csv", comment = "#")

#Looking at the column names
head(co2_data)

#Print only the column names of the dataset
colnames(co2_data)

#Clean up the columns
co2_clean <- co2_data %>%
	select(Year = year, Month = month, CO2 = average)

#View the cleaned structure
head(co2_clean)

# Filter rows where the Year is greater than or equal to 1970
co2_modern <- co2_clean %>% 
  filter(Year >= 1970)

# Check the summary to ensure our minimum year is now 1970
summary(co2_modern)

#Create a single date column
co2_final <- co2_modern %>%
  mutate(Date = make_date(year = Year, month = Month, day = 1))

#Look at the final wrangled dataset
head(co2_final)

#Plotting the continuous CO2 line
ggplot(data = co2_final, aes(x = Date, y = CO2)) +
  geom_line(color = "darkblue", linewidth = 0.7) +
  theme_minimal() +
  labs(
    title = "The Keeling Curve: Atmospheric CO2 Concentration",
    subtitle = "Mauna Loa Observatory, Hawaii (1970 - Present)",
    x = "Year",
    y = "Carbon Dioxide (parts per million)"
   )

#Zoomed in plot to see seasonal fluctuations from 2020-2023
ggplot(data = co2_final, aes(x = Date, y = CO2)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_point(color = "darkgreen", size = 1.5, alpha = 0.5) + # Add points to see monthly data
  coord_cartesian(xlim = as.Date(c("2020-01-01", "2023-12-01"))) + # Zoom to 2020-2023
  theme_light() +
  labs(
    title = "Seasonal Atmospheric CO2 Fluctuations",
    subtitle = "Zoomed view showing the Earth 'breathing' over a 3-year cycle",
    x = "Month / Year",
    y = "Carbon Dioxide (parts per million)"
  )

#Group by year and calculate the average CO2 level for each year
co2_annual <- co2_final %>%
  group_by(Year) %>%
  summarize(Mean_CO2 = mean(CO2))

#Looking at the new annual summary table
print(co2_annual)

#Calculate the difference between the current year and the previous year
co2_rate <- co2_annual %>%
  mutate(Annual_Increase = Mean_CO2 - lag(Mean_CO2))

#View the table with the new growth rate column
print(co2_rate)

#Overall average increase per year across the whole dataset
mean(co2_rate$Annual_Increase, na.rm = TRUE)

#Look at the rate of increase in just the last 5 years to compare
co2_rate %>% filter(Year >= 2020)

#Fit a quadratic regression model using the annual data
#I(Year^2) tells R to treat Year-squared as a literal math operation
climate_model <- lm(Mean_CO2 ~ Year + I(Year^2), data = co2_annual)

#View the model summary
summary(climate_model)

#Create a small dataframe containing the year 2050
future_year <- data.frame(Year = 2050)

#Predict CO2 for the year 2050
predicted_co2 <- predict(climate_model, newdata = future_year)

#Print the result to the console
print(predicted_co2)

#Part 2: Predicting the CO2 Levels Until 2050

#Create the sequence of years from 2026 to 2050
future_years <- data.frame(Year = 2026:2050)

#Use the quadratic model from earlier to predict CO2 for ALL of those years
future_prediction <- predict(climate_model, newdata = future_years)

#Combine the years and predictions into a clean dataframe
future_df <- future_years %>%
  mutate(
    Mean_CO2 = future_prediction, 
    Data_Type = "Predicted"
  )

#Prepare historical data for merging
historical_df <- co2_annual %>%
  mutate(Data_Type = "Historical")

#Combine both dataframes into one
full_timeline <- bind_rows(historical_df, future_df)

#Plot the full timeline up to 2050
ggplot(data = full_timeline, aes(x = Year, y = Mean_CO2, color = Data_Type)) +
  geom_line(linewidth = 1.2) +
  # Add a dashed vertical line at the current year (2026) to show where prediction starts
  geom_vline(xintercept = 2026, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Historical" = "darkblue", "Predicted" = "darkred")) +
  theme_minimal() +
  labs(
    title = "Atmospheric CO2 Projections through 2050",
    subtitle = "Historical NOAA data modeled with a quadratic forecast into the future",
    x = "Year",
    y = "Annual Mean CO2 (ppm)",
    color = "Data Type"
    )

#Save the plot with the timeline to 2050
ggsave("predicted_plot.png", width = 6, height = 4, dpi = 300)
ggsave("predicted_plot.pdf", width = 7, height = 5)