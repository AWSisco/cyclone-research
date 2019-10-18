#Get text files for the matplotlib viridis colormaps
library(viridis)
#A: magma; B: inferno; C: plasma; D: viridis; E: cividis)

magma <- subset(viridis.map, viridis.map$opt=="A")[,1:3]
inferno <- subset(viridis.map, viridis.map$opt=="B")[,1:3]
plasma <- subset(viridis.map, viridis.map$opt=="C")[,1:3]
viridis <- subset(viridis.map, viridis.map$opt=="D")[,1:3]
cividis <- subset(viridis.map, viridis.map$opt=="E")[,1:3]

write.table(magma, file = "mpl_magma.txt", col.names = FALSE, row.names = FALSE, sep = " ")
write.table(inferno, file = "mpl_inferno.txt", col.names = FALSE, row.names = FALSE, sep = " ")
write.table(plasma, file = "mpl_plasma.txt", col.names = FALSE, row.names = FALSE, sep = " ")
write.table(viridis, file = "mpl_viridis.txt", col.names = FALSE, row.names = FALSE, sep = " ")
write.table(cividis, file = "mpl_cividis.txt", col.names = FALSE, row.names = FALSE, sep = " ")

spectral_hex <- c("#9e0142","#a00343","#a20643","#a40844","#a70b44","#a90d45","#ab0f45",
                  "#ad1245","#af1446","#b11646","#b31947","#b51b47","#b71d48","#ba2048",
                  "#bc2248","#be2449","#c02749","#c12949","#c32b4a","#c52d4a","#c7304a",
                  "#c9324a","#cb344b","#cd364b","#ce384b","#d03b4b","#d23d4b","#d33f4b",
                  "#d5414b","#d7434b","#d8454b","#da474a","#db494a","#dd4b4a","#de4d4a",
                  "#df4f4a","#e1514a","#e2534a","#e35549","#e45749","#e65949","#e75b49",
                  "#e85d49","#e95f49","#ea6149","#eb6349","#ec6549","#ed6749","#ee6a49",
                  "#ef6c49","#f06e4a","#f0704a","#f1724a","#f2744b","#f3774b","#f3794c",
                  "#f47b4d","#f47e4d","#f5804e","#f6824f","#f68550","#f78750","#f78951",
                  "#f88c52","#f88e53","#f89154","#f99356","#f99557","#f99858","#fa9a59",
                  "#fa9c5a","#fa9f5c","#fba15d","#fba35e","#fba660","#fba861","#fcaa62",
                  "#fcad64","#fcaf65","#fcb167","#fcb368","#fcb56a","#fdb86b","#fdba6d",
                  "#fdbc6e","#fdbe70","#fdc071","#fdc273","#fdc474","#fdc676","#fdc878",
                  "#fdca79","#fecc7b","#fecd7d","#fecf7e","#fed180","#fed382","#fed584",
                  "#fed685","#fed887","#feda89","#fedb8b","#fedd8d","#fede8f","#fee090",
                  "#fee192","#fee394","#fee496","#fee698","#fee79a","#fee89b","#feea9d",
                  "#feeb9f","#feeca1","#feeda2","#feefa4","#fef0a5","#fef1a7","#fef2a8",
                  "#fdf3a9","#fdf3aa","#fdf4ab","#fdf5ac","#fcf6ad","#fcf6ae","#fcf7af",
                  "#fbf7af","#fbf8b0","#faf8b0","#faf9b0","#f9f9b0","#f9f9b0","#f8f9b0",
                  "#f7faaf","#f7faaf","#f6faae","#f5faae","#f4f9ad","#f3f9ac","#f2f9ac",
                  "#f2f9ab","#f0f9aa","#eff8a9","#eef8a8","#edf8a7","#ecf7a7","#ebf7a6",
                  "#e9f6a5","#e8f6a4","#e7f5a3","#e5f5a2","#e4f4a2","#e2f3a1","#e0f3a1",
                  "#dff2a0","#ddf1a0","#dbf19f","#d9f09f","#d7ef9f","#d6ee9f","#d4ee9f",
                  "#d2ed9e","#d0ec9e","#cdeb9f","#cbea9f","#c9e99f","#c7e89f","#c5e89f",
                  "#c3e79f","#c0e6a0","#bee5a0","#bce4a0","#b9e3a0","#b7e2a1","#b4e1a1",
                  "#b2e0a1","#b0dfa1","#addea2","#abdda2","#a8dca2","#a6dba3","#a3daa3",
                  "#a0d9a3","#9ed8a3","#9bd7a3","#99d6a4","#96d5a4","#94d4a4","#91d3a4",
                  "#8ed1a4","#8cd0a4","#89cfa5","#87cea5","#84cda5","#82cba5","#7fcaa6",
                  "#7dc9a6","#7ac7a6","#77c6a6","#75c5a7","#73c3a7","#70c2a8","#6ec0a8",
                  "#6bbea8","#69bda9","#66bba9","#64b9aa","#62b8aa","#60b6ab","#5db4ac",
                  "#5bb2ac","#59b0ad","#57aeae","#55acae","#53aaaf","#51a8af","#50a6b0",
                  "#4ea4b1","#4ca2b1","#4ba0b2","#499db2","#489bb3","#4799b3","#4697b3",
                  "#4595b4","#4492b4","#4390b4","#438eb4","#428cb5","#4289b5","#4287b4",
                  "#4285b4","#4283b4","#4280b4","#437eb3","#437cb3","#447ab3","#4577b2",
                  "#4575b1","#4673b1","#4771b0","#486eaf","#4a6caf","#4b6aae","#4c68ad",
                  "#4e65ac","#4f63ab","#5161aa","#525fa9","#545ca8","#555aa7","#5758a6",
                  "#5956a5","#5b53a4","#5c51a3","#5e4fa2")

spectral_rgb <- data.frame(t(col2rgb(spectral_hex))/255.0)
write.table(spectral_rgb, file = "cbr_spectral.txt", 
            col.names = FALSE, row.names = FALSE, sep = " ")
