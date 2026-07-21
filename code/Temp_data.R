# Creating function to incorporate Greenhouse data in R ----
## Be sure to delete additional information in column names from excel data sheet
## ex: ,C5 or ,COM

# Packages ----
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("tidyr")
#install.packages("tidyverse")
#install.packages("lubridate")
#install.packages("readxl")
#install.packages("degday")
library(degday)
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(lubridate)

#Exporting data ----

temperature_df_raw <- read_excel("data/data_modif.xlsx", col_names = TRUE)
colnames(temperature_df_raw)

#Creating function ----

replace_commas <- function(x, positions) {
  commas <- gregexpr(",", x)[[1]]
  if (commas[1] == -1) return(x)
  for (pos in positions) {
    if (pos <= length(commas)) {
      substr(x, commas[pos], commas[pos]) <- "."
    }
  }
  return(x)
}

# Choosing which separator to transform in dot ----

comma_positions <- c(2, 4, 6, 8, 14, 18)

# Apply function to temperature_df ----

original_colnames <- colnames(temperature_df_raw)

temperature_df_raw[[1]] <- sapply(temperature_df_raw[[1]], replace_commas, positions = comma_positions)

temperature_df <- read_csv(I(temperature_df_raw[[1]]), col_names = FALSE)

colnames(temperature_df) <- original_colnames

temperature_df <- temperature_df %>% mutate(across(where(is.character), ~ parse_guess(.)))

# Verification ----

glimpse(temperature_df)
summary(temperature_df)
temperature_df

# Data management for degday function ----
# info needed per day: Tmin, Tmax

# Convert Celcius into Farenheit

temperature_df <- temperature_df %>%
  mutate(temp_f = `Climate Temperature(°C)` * 9/5 + 32)

# Extract Tmin and Tmax

temp_journaliere <- temperature_df %>%
  mutate(date = as.Date(`Date/Time`)) %>%
  group_by(date) %>%
  summarise(
    temp_min_f = min(temp_f, na.rm = TRUE),
    temp_max_f = max(temp_f, na.rm = TRUE)
  )
# Add a variable for the minimal temp of next week

temp_journaliere <- temp_journaliere %>%
  arrange(date) %>%             
  mutate(
    tmin_next = lead(temp_min_f, 1) 
  )
