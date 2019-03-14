

raw_com <- read.csv("../Data/SmMammRawCom.csv")
env <- read.csv("../Data/Nocturnal only site data.csv")
all <- left_join(raw_com, env)
write.csv(all, "Data/large_mammal.csv", row.names = F)

# Load community data
raw_com <- read.csv("../Data/Lg_mammal_com_dom_incl.csv")

# Load environmental data
env <- read.csv("../Data/Large mammal data.csv")
