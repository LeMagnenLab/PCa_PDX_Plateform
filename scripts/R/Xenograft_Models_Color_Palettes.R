library(paletteer)
library(colorRamp2)
library(colorRamps)

pal.PDOX.PDOXO.updated <- c(
  ########### P20.11 ####################
  "P20-11_Mouse_B_P1_Xenograft"           = "#8B0000FF",   # Dark Red
  "P20-11_Mouse_B_P1_Xenograft_Organoids" = "#CC5151FF",   # Indian Red
  "P20-11_Mouse_B_P3_Xenograft_Organoids" = "#E57E7EFF",   # Dark Orange
  "P20-11_Mouse_A_P3_Xenograft_Organoids" = "#99540FFF",   # Orange Red
  
  ########### P20.23 ####################
  "P20-23_Mouse_A_P1_Xenograft"           = "#00008BFF",   # Dark Blue
  "P20-23_Mouse_A_P1_Xenograft_Organoids" = "#4169E1FF",   # Royal Blue
  "P20-23_Mouse_A_P3_Xenograft"           = "#0F6B99FF",   # Slate Blue
  "P20-23_Mouse_A_P3_Xenograft_Organoids" = "#7EC3E5FF",   # Light Blue
  "P20-23_Mouse_A_P9_Xenograft_Organoids" = "#9370DBFF"    # Medium Purple
)

pal.PDOX.PDOXO.updated.short.id <- c(
  ########### P20.11 ####################
  "P20-11_BP1_PDX" = "#8B0000FF",   # Dark Red
  "P20-11_BP1_Org" = "#CC5151FF",   # Indian Red
  "P20-11_BP3_Org" = "#E57E7EFF",   # Dark Orange
  "P20-11_AP3_Org" = "#99540FFF",   # Orange Red
  
  ########### P20.23 ####################
  "P20-23_AP1_PDX" = "#00008BFF",   # Dark Blue
  "P20-23_AP1_Org" = "#4169E1FF",   # Royal Blue
  "P20-23_AP3_PDX" = "#0F6B99FF",   # Slate Blue
  "P20-23_AP3_Org" = "#7EC3E5FF",   # Light Blue
  "P20-23_AP9_Org" = "#9370DBFF"    # Medium Purple
)

pal.patient.id = c(
 "P20-23" = "#00008BFF",   # Dark Blue
  "P20-11" = "#8B0000FF"   # Dark Red
)

# Teatment paelette
pal.treatment = c("noDHT_noENZA" = "#F3C300", 
                  "noDHT_ENZA" = "#00FFFF", 
                  "DHT_noENZA" ="#875692", 
                  "DHT_noENZ" ="#875692", 
                  "DHT_ENZA" = "#00FF66")

# Teatment paelette
pal.dht.status = c("noDHT" = "#F3C300", 
            "DHT" = "#875692")

#  Clustering
pal_clustering = c(
  "1" = "#F38400",
  "2" = "#A1CAF1", 
  "3" = "#BE0032", 
  "4" = "#DCD300", 
  "5" = "#008856", 
  "6" = "#4B0082", 
  "7" = "#E68FAC", 
  "8" = "#0067A5", 
  "9" = "#F99379", 
  "10" = "#604E97", 
  "11"  =  "#F6A600", 
  "12" = "#B3446C", 
  "13" = "#848482", 
  "14" = "#882D17", 
  "15" =  "#8DB600", 
  "16" = "#654522", 
  "17" = "#E25822",
  "18" = "#2B3D26",
  "19" = "#C2B280")

pal_clustering_bis = c(
  "1" = "#FED439FF",
  "2" = "#709AE1FF", 
  "3" = "#8A9197FF", 
  "4" = "#D2AF81FF", 
  "5" = "#D5E4A2FF", 
  "6" = "#4B0082")



pal_zscore <- c(
  "#26456EFF", "#244C7CFF", "#21538BFF", "#1C5A99FF", "#1C63A1FF", "#1C6CAAFF", "#1F74B1FF", 
  "#2B7BB4FF", "#3482B6FF", "#3F8BBAFF", "#4F98C4FF", "#5EA5CEFF", "#78B1D3FF",
  "#9CBBCFFF", "#FFFFFF", "#D6BFBBFF", "#E9A79DFF", "#F78E80FF", "#F6796AFF", 
  "#EC6857FF", "#E25644FF", "#DC4636FF", "#D73529FF", "#D21E1CFF", "#CB1618FF", 
  "#C51517FF", "#BE1316FF", "#B3101BFF", "#A70C20FF", "#9C0824FF"
)

pal.log2 = color_function = c("lightgrey",rev(paletteer_c("viridis::magma", 30)))

# cell_cycle_phase
pal_cell_cycle_phase = c(
  "S" ="#1F77B4FF",
  "G1" = "#F8766D",  
  "G2M" = "#2CA02CFF")

# Model System
pal.system = c(
  "Xenograft_Organoids" = "#ADD8E6", 
  "Xenograft" = "#FFD1DC")
