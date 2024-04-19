# -----------------------------------------------------------------------------------------------------------------------
# Code to read in climate and dengue data for Singapore from 2000 to 2022
# Author: Emilie Finch
# -----------------------------------------------------------------------------------------------------------------------

# Load packages and fns for analysis ------------------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  here, dplyr, janitor, data.table, openxlsx, lubridate, ggplot2, ggpubr, tidyr, nnet, splines, zoo, corrplot,
  INLA, tibble, stringr, showtext, sysfonts, cowplot, scoringutils, purrr, RColorBrewer, pROC, ggpubr, scales,
  MetBrewer, tidyquant, zoo, hydroGOF, showtext
)


# Load data -------------------------------------------------------------------------------------------------------------

dengue_singapore <- read.csv(here("data", "dengue-cases-climate.csv"))
dengue_singapore <- dengue_singapore |> 
  mutate(date = as.Date(date, format = "%d/%m/%Y"))