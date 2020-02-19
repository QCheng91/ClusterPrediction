###################Get the raw data###########################
library(flowCore)
library(tidyverse)
library(stringr)

#2020 data
get_pro_data <- function() {
  
  #K1
  pro1 <- read.FCS("Raw Data/cytof-erythroid/2020/FCS-Files-Ungated-for-Ed-Qian-Kinetic1/experiment_276252_illustration_344261_spill_applied_fcs_files/c01_20191126_Helios1-1_Barcodedx20_01_0_concatenated_Ungated_Day_0_Ungated.fcs", transformation=FALSE)
  pro2 <- read.FCS("Raw Data/cytof-erythroid/2020/FCS-Files-Ungated-for-Ed-Qian-Kinetic1/experiment_276252_illustration_344261_spill_applied_fcs_files/c02_20191126_Helios1-1_Barcodedx20_01_0_concatenated_Ungated_Day_2_Ungated.fcs", transformation=FALSE)
  pro3 <- read.FCS("Raw Data/cytof-erythroid/2020/FCS-Files-Ungated-for-Ed-Qian-Kinetic1/experiment_276252_illustration_344261_spill_applied_fcs_files/c03_20191126_Helios1-1_Barcodedx20_01_0_concatenated_Ungated_Day_4_Ungated.fcs", transformation=FALSE)
  pro4 <- read.FCS("Raw Data/cytof-erythroid/2020/FCS-Files-Ungated-for-Ed-Qian-Kinetic1/experiment_276252_illustration_344261_spill_applied_fcs_files/c04_20191126_Helios1-1_Barcodedx20_01_0_concatenated_Ungated_Day_6_Ungated.fcs", transformation=FALSE)
  pro5 <- read.FCS("Raw Data/cytof-erythroid/2020/FCS-Files-Ungated-for-Ed-Qian-Kinetic1/experiment_276252_illustration_344261_spill_applied_fcs_files/c05_20191126_Helios1-1_Barcodedx20_01_0_concatenated_Ungated_Day_8_Ungated.fcs", transformation=FALSE)
  pro6 <- read.FCS("Raw Data/cytof-erythroid/2020/FCS-Files-Ungated-for-Ed-Qian-Kinetic1/experiment_276252_illustration_344261_spill_applied_fcs_files/c06_20191126_Helios1-1_Barcodedx20_01_0_concatenated_Ungated_Day_10_Ungated.fcs", transformation=FALSE)
  pro7 <- read.FCS("Raw Data/cytof-erythroid/2020/FCS-Files-Ungated-for-Ed-Qian-Kinetic1/experiment_276252_illustration_344261_spill_applied_fcs_files/c07_20191126_Helios1-1_Barcodedx20_01_0_concatenated_Ungated_Day_11_Ungated.fcs", transformation=FALSE)
  pro8 <- read.FCS("Raw Data/cytof-erythroid/2020/FCS-Files-Ungated-for-Ed-Qian-Kinetic1/experiment_276252_illustration_344261_spill_applied_fcs_files/c08_20191126_Helios1-1_Barcodedx20_01_0_concatenated_Ungated_Day_12_Ungated.fcs", transformation=FALSE)
  pro9 <- read.FCS("Raw Data/cytof-erythroid/2020/FCS-Files-Ungated-for-Ed-Qian-Kinetic1/experiment_276252_illustration_344261_spill_applied_fcs_files/c09_20191126_Helios1-1_Barcodedx20_01_0_concatenated_Ungated_Day_14_Ungated.fcs", transformation=FALSE)
  pro10 <- read.FCS("Raw Data/cytof-erythroid/2020/FCS-Files-Ungated-for-Ed-Qian-Kinetic1/experiment_276252_illustration_344261_spill_applied_fcs_files/c10_20191126_Helios1-1_Barcodedx20_01_0_concatenated_Ungated_Day_16_Ungated.fcs", transformation=FALSE)
  
  pro1_df = flowCore::exprs(pro1)
  pro2_df = flowCore::exprs(pro2)
  pro3_df = flowCore::exprs(pro3)
  pro4_df = flowCore::exprs(pro4)
  pro5_df = flowCore::exprs(pro5)
  pro6_df = flowCore::exprs(pro6)
  pro7_df = flowCore::exprs(pro7)
  pro8_df = flowCore::exprs(pro8)
  pro9_df = flowCore::exprs(pro9)
  pro10_df = flowCore::exprs(pro10)
  
  g1 = data.frame(pro1_df) 
  g2 = data.frame(pro2_df) 
  g3 = data.frame(pro3_df) 
  g4 = data.frame(pro4_df) 
  g5 = data.frame(pro5_df) 
  g6 = data.frame(pro6_df) 
  g7 = data.frame(pro7_df) 
  g8 = data.frame(pro8_df) 
  g9 = data.frame(pro9_df) 
  g10 = data.frame(pro10_df)
  
  ###############################################################################
  sur_marker = c("CD123", "CD135", "CD235ab", "CD45RA", "CD34", "CD36", "CD49f",
                 "CD44",  "CD71",   "CD38",    "CD90", "CD203c", "CD41")
  
  other_marker = c("HBA")

  tran_factor = c("TBX15", "NFE2p45", "GATA1-PE", "KLF1_2", "ProMBP1", "MAFG",
                  "PU1",   "CEBPa",    "RUNX1",   "IKZF1",   "FLI1_2", "BACH1", 
                  "Ki67",  "CyclinB1", "IdU", "GATA2")

  DNA = c("Ir191", "Ir193", "Ce140", "Pt195")
  
  
  isotops = c("Eu151Di", "Lu175Di", "Pr141Di", "Nd143Di", "Sm149Di", "Gd155Di", "Dy164Di",
              "Er166Di", "Er168Di", "Yb172Di", "Dy161Di", "Nd146Di", "Y89Di",
              
              "Er170Di",
              
              "Yb171Di", "Sm154Di", "Gd156Di", "Ho165Di", "Nd148Di", "Sm152Di",
              "Er167Di", "Nd145Di", "Gd158Di", "Tb159Di", "Tm169Di", "Nd150Di",
              "Dy162Di", "Eu153Di", "I127Di", "Dy163Di",
              
              "Ir191Di",  "Ir193Di", "Ce140Di", "Pt195Di")

  pro_marker = c(sur_marker, other_marker, tran_factor, DNA)
  
  seq=c(1:34)
  
  day_levels = c(0, 2, 4, 6, 8, 10, 11, 12, 14, 16)
  
  prot_list = list()
  
  for(i in seq) {
    
    marker = isotops[i]
    
    #K1
    protein1 = bind_rows( g1[marker] %>% mutate(obstime = 0, col_id = 1:nrow(g1)),
                          g2[marker]%>% mutate(obstime = 2, col_id = 1:nrow(g2)),
                          g3[marker]%>% mutate(obstime = 4, col_id = 1:nrow(g3)),
                          g4[marker]%>% mutate(obstime = 6, col_id = 1:nrow(g4)),
                          g5[marker]%>% mutate(obstime = 8, col_id = 1:nrow(g5)),
                          g6[marker]%>% mutate(obstime = 10, col_id = 1:nrow(g6)),
                          g7[marker]%>% mutate(obstime = 11, col_id = 1:nrow(g7)),
                          g8[marker]%>% mutate(obstime = 12, col_id = 1:nrow(g8)),
                          g9[marker]%>% mutate(obstime = 14, col_id = 1:nrow(g9)),
                          g10[marker]%>% mutate(obstime = 16, col_id = 1:nrow(g10))) %>%
      
      mutate(obstime = factor(obstime, levels = day_levels))
    
    protein1 = protein1 %>% mutate(prot_name = pro_marker[i])
    colnames(protein1) = c("Intensity", "obstime", "cell_id", "cell.type")
    prot_list = c(prot_list, list(protein1))
  }   
  return(prot_list)
}

all_data = as_tibble(bind_rows(get_pro_data()))

all_data = all_data %>% mutate(cell = factor(paste(cell_id, obstime)))
all_data = all_data %>% mutate(capture = obstime)

all_data %>% ggplot(aes(x=obstime, y=asinh(Intensity), fill = obstime)) +
  geom_violin(size=0.1, scale = "width", bw = 0.25) +
  facet_wrap(~cell.type)

gene_data0 = all_data %>% spread(cell.type, Intensity) 

mean_irid = function(ir1, ir2){
  reg = lm(ir1~ir2)
  # reg = lm(data = gated_data, Ir191~Ir193)
  (predict(reg) + ir1)/2
}

gene_data = gene_data0 %>% filter(Ce140 <= sinh(5.0)) %>% 
  mutate(mean_Ir = mean_irid(Ir191, Ir193)) %>% 
  filter(Pt195 <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5)) #%>% 
  filter(mean_Ir >= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5),
         mean_Ir <= sinh(8.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5))

gene_data %>% group_by(obstime) %>% summarise(count = n())

gene_data %>% gather(protein, Intensity, BACH1:TBX15) %>% 
    ggplot(aes(x=obstime, y=asinh(Intensity), fill = obstime)) +
    geom_violin(size=0.1, scale = "width", bw = 0.25) +
    facet_wrap(~protein)

gene_data %>% ggplot(aes(x=asinh(ProMBP1), y=asinh(GATA2), colour = asinh(CD203c))) +
  geom_point(size=1.5, alpha = 0.5) +
  facet_wrap(~obstime)


gene_data %>% filter(ProMBP1>= sinh(7.5), GATA2 <= sinh(6.5), CD123 >= sinh(5.5)) %>% 
  ggplot(aes(x=asinh(ProMBP1), y=asinh(GATA2), colour = asinh(CD123))) +
  geom_point(size=1.5, alpha = 0.5) +
  facet_wrap(~obstime)



