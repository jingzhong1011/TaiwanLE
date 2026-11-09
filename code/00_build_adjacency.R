########################
### Adjacency matrix ###
########################
# Queen contiguity of townships for the ICAR terms (data/NB_mat.csv).
# Offshore counties (Kinmen, Lienchiang; code "09") are excluded and townships
# are sorted by MOI to match the row order of DeathCount_v5.csv.
# Run from the repository root.

###################
### Environment ###
###################
library(sf)
library(spdep)

#################
### Adjacency ###
#################
TWbord <- TWmap::bord_town
TWbord <- TWbord[substr(TWbord$town, 1, 2) != "09", ]
TWbord <- TWbord[order(as.numeric(TWbord$town)), ]

NB <- nb2mat(poly2nb(TWbord, queen = T), style = "B", zero.policy = T)
colnames(NB) <- paste0("V", 1: ncol(NB))

write.csv(NB, "data/NB_mat.csv", row.names = F)
