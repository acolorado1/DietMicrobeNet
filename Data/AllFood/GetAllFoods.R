# Get all foods from FooDB
library(readr)
library(dplyr)

food <- read_csv("Data/Food.csv")

# create file formatted like streamlit web app output 
all_food <- food %>% 
  dplyr::select(id, name, name_scientific) %>% 
  dplyr::mutate(food_frequency = 1)

# this file can be found in the Data/AllFood file 
# write_csv(x = all_food, file = "foodb_allfoods_dataframe.csv")
