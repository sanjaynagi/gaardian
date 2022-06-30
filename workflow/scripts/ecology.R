library(tidyverse)
library(data.table)
rm(list=ls())

load("Obuasi.Rdata")
Obuasi = Obuasi %>% rename(latitude=y, longitude=x)

metadata = fread("config/metadata.tsv")

# Get nearest points in Obuasi data
library(rgeos)
set1sp <- SpatialPoints(metadata[, c('longitude', 'latitude')])
set2sp <- SpatialPoints(Obuasi[, c('longitude', 'latitude')])

metadata$ObuasiRow = apply(gDistance(set2sp, set1sp, byid=TRUE), 1, which.min)
Obuasi$ObuasiRow = seq(1, nrow(Obuasi))
Obuasi = Obuasi %>% select(-c('longitude', 'latitude'))

metadata = left_join(metadata, Obuasi)
env_columns = colnames(Obuasi)[1:13]

## check correlatins between environmental variables
library(pheatmap)
env_cor = cor(Obuasi[, env_columns])
pheatmap(env_cor, legend=TRUE,)

env_columns = env_columns[!env_columns %in% c('aetV', 'LC2', 'aetM', 'temM')]

env_cor2 = cor(Obuasi[, env_columns])
pheatmap(env_cor2, legend=TRUE,)

# center and scale environmental variables
metadata = metadata %>% select(-c('aetV', 'LC2', 'aetM', 'temM')) %>%  mutate_at(.vars = env_columns, 
          .funs = list("scaled" = scale), center=TRUE)


### map env variables ###


library(sf)

map = st_as_sf(metadata, coords = c('longitude', 'latitude'))
map = st_set_crs(map, crs = 4326)

ggplot(map) + geom_sf(aes(color=LC_scaled)) + theme_light()

ggplot(map) + geom_sf(aes(color=Elev_scaled)) + theme_light()




meta_ecol = cbind(metadata[, c("location", "latitude", "longitude", 'mining')], metadata[, grep("scaled", names(metadata)), with = FALSE]) %>% 
  distinct() %>% 
  data.frame()

meta_ecol = meta_ecol %>% group_by(location) %>% 
  mutate(duplication_id = seq(n())) %>% 
  ungroup() %>% filter(duplication_id == 1) %>% select(-duplication_id)

# Read Fst data 
fst = fread("results/spatial/Fst_fm_coluzzii.tsv") %>% rename(location=loc1, location2=loc2)
# Remove locations with less than n
n = metadata %>% group_by(location) %>% count()
n = n %>% filter(n > 10)
# Remove locs if there 
fst = fst %>% filter(location %in% n$location & location2 %in% n$location)

# join fst with meta ecology dataframe
fst = full_join(meta_ecol, fst, suffix=c("x", "y")) %>% select(-c("latitude", "longitude", "V1"))
meta_ecol = meta_ecol %>% mutate(location2=location) %>% select(-c("latitude", "longitude"))
# merge again 
meta_ecol = merge(meta_ecol, fst, by='location2') %>% select(-location2)

euclidean <- function(a, b) sqrt(sum((a - b)^2))

# calculate euclidean ecological distance for each comparison
ecol_dist = c()
for (i in 1:nrow(meta_ecol)){
    a = meta_ecol[i, 3:11] %>% as.vector() %>% t() %>% as.vector()
    b = meta_ecol[i, 14:22] %>% as.vector() %>% t() %>% as.vector()
    
    dist = euclidean(a,b)
    ecol_dist = c(ecol_dist, dist)
}

meta_ecol$ecol_dist = ecol_dist

meta_ecol$mining.x = as.numeric(as.factor(meta_ecol$mining.x))
meta_ecol$mining.y = as.numeric(as.factor(meta_ecol$mining.y))
meta_ecol$mining = abs(meta_ecol$mining.x - meta_ecol$mining.y)


# Pivot into separate matrices 
fst_matrix = meta_ecol %>% select(location.x, location.y, fst) %>% pivot_wider(., names_from = "location.x", values_from = 'fst') %>% column_to_rownames('location.y')
ecol_matrix = meta_ecol %>% select(location.x, location.y, ecol_dist) %>%  pivot_wider(., names_from = "location.x", values_from = 'ecol_dist') %>% column_to_rownames('location.y')
km_matrix = meta_ecol %>% select(location.x, location.y, km) %>%  pivot_wider(., names_from = "location.x", values_from = 'km') %>% column_to_rownames('location.y')
mining_matrix = meta_ecol %>% select(location.x, location.y, mining) %>%  pivot_wider(., names_from = "location.x", values_from = 'mining') %>% column_to_rownames('location.y')



tri_replace = function(m1){
  lt <- lower.tri(m1)
  ut <- upper.tri(m1)
  
  m1[lt] <- ifelse(is.na(lowerTriangle(m1, byrow=FALSE)), upperTriangle(m1, byrow=TRUE), lowerTriangle(m1, byrow=FALSE))
  upperTriangle(m1) <- lowerTriangle(m1, byrow=TRUE)
  return(m1)
}

library(gdata)
### sort and replace triangles with respective elemetns from other triangle 
fst_matrix = fst_matrix[sort(row.names(fst_matrix)),]  %>% 
  select(sort(current_vars()))
fst_matrix = tri_replace(fst_matrix)

km_matrix = km_matrix[sort(row.names(km_matrix)),] %>% 
  select(sort(current_vars()))
km_matrix = tri_replace(km_matrix)

ecol_matrix = ecol_matrix[sort(row.names(ecol_matrix)),] %>% 
  select(sort(current_vars()))
ecol_matrix = tri_replace(ecol_matrix)

mining_matrix= mining_matrix[sort(row.names(mining_matrix)),] %>% 
  select(sort(current_vars()))
mining_matrix = tri_replace(mining_matrix)

### Mantel tests ### 
library(vegan)


rda(X=fst_matrix, Y=ecol_matrix, Z=km_matrix, scale=FALSE)

meta_ecol$fst

rda(fst~)

?rda
### TODO 
# Control for km_matrix, using Fst and gold mine matrix...
# Add in gold mine into ecological distance matrix
fst_matrix = fst_matrix/(1-fst_matrix)


mantel(fst_matrix, km_matrix, permutations = 10000)
 
mantel(fst_matrix, ecol_matrix, permutations = 10000)

mantel(km_matrix, ecol_matrix, permutations = 10000)

mantel.partial(fst_matrix, km_matrix, ecol_matrix, permutations = 10000)

mantel.partial(fst_matrix, ecol_matrix, km_matrix, permutations = 10000)



mantel.partial(fst_matrix, km_matrix, mining_matrix, permutations = 10000)




mantel.partial(fst_matrix, mining_matrix, km_matrix, permutations = 10000, method='spearman')
?mantel.partial


mantel.partial(fst_matrix, mining_matrix, ecol_matrix, permutations = 10000)

mantel.partial(ecol_matrix, mining_matrix, km_matrix, permutations = 10000)
