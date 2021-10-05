library(tidyverse)

dat <- read_csv(
    file = "~/Documents/Projects/Cascading Effects/data/161025_RegimeShiftsDatabase.CSV"
)


dat %>% pull(3)
