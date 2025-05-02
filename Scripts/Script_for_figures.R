#     _          _                 _   _                       _                                         _          
#    / \   _ __ | |_ __ _ _ __ ___| |_(_) ___   _ __ ___   ___| |_ __ _  __ _  ___ _ __   ___  _ __ ___ (_) ___ ___ 
#   / _ \ | '_ \| __/ _` | '__/ __| __| |/ __| | '_ ` _ \ / _ \ __/ _` |/ _` |/ _ \ '_ \ / _ \| '_ ` _ \| |/ __/ __|
#  / ___ \| | | | || (_| | | | (__| |_| | (__  | | | | | |  __/ || (_| | (_| |  __/ | | | (_) | | | | | | | (__\__ \
# /_/   \_\_| |_|\__\__,_|_|  \___|\__|_|\___| |_| |_| |_|\___|\__\__,_|\__, |\___|_| |_|\___/|_| |_| |_|_|\___|___/
#                                                                       |___/                                       

# Jiarui Sun (jiarui.sun@uq.edu.au) and Paul Dennis

library(gsveasyr) # To install this package, contact Jiarui
library(vegan)
source("Functions.R")

load("../Data/ANT_CARD_w_env.RData")

# Data frame  |
# ------------------------------------------------------------------------------
# ARO_Refs    | Information of antimicrobial resistance ontology (ARO) terms
# card.a      | Number of detected ORFs for each ARO in samples
# rel.card.a  | 

################################################################################
# Step 0 - Prepare data set
################################################################################
# Temperature (MASAT) colors:
masat_col.df <- read.csv("../Data/masat_color_code.csv", header = T, row.names = 1) 
all(rownames(masat_col.df) == rownames(env.df)) # [1] TRUE
env.df$masat_col <- c(masat_col.df$Colour)

# Rename AROs by gene name:
all(colnames(rel.card.a) == ARO_Refs$ARO_ID) # [1] TRUE
colnames(rel.card.a) <- ARO_Refs$GeneName
# Re-order samples by Latitude (instead of temperature), and re-order genes by
# gene names:
rel.card.a <- rel.card.a[order(env.df$Lat),order(colnames(rel.card.a))]
env.df <- env.df[order(env.df$Lat),]
# ARGs present in over half samples:
sel.rel.card.a <- rel.card.a[,which(colSums(decostand(rel.card.a, "pa"))>14)]
# 14 ARGs present in > 5 samples;  13 ARGs present in > 10 samples; 12 ARGs present in > 11 samples


# Alpha diversity indices:
tot <- rowSums(rel.card.a) # Total abundances
sobs <- rowSums(decostand(rel.card.a, method='pa')) # Sobs
shan <- diversity(rel.card.a) # Shannon index

# Composition / Beta diversity as Primary and secondary principle components:
pca <- rda(sqrt(rel.card.a))
PC1 <- pca$CA$u[,1] # scores(pca)$sites[,1]
PC2 <- pca$CA$u[,2] # scores(pca)$sites[,2]


################################################################################
# Fig 1 - Correlations between environmental factors
################################################################################
# 1) Correlation plot (diagonal):
tmp.df <- env.df[,c('Lat','Long','Alt','Temp',sort(env_params[2:21]))]
#   // Add total + shan
tmp.df$Total <- tot
tmp.df$Sobs <- sobs
tmp.df$Shannon <- shan
#  // Add PC1, PC2 
envfit.res <- envfit(pca, tmp.df[,1:24], na.rm = T)
tmp.df$PC1 <-  pca$CA$u[,1]
tmp.df$PC2 <-  pca$CA$u[,2]
tmp.df$HFP <- env.df$HFP 
#  // As NH4 in Tub.58 is NA, we need re-calculate correlations of NH4 without this sample:
tmp.df_NH4 <- tmp.df[-which(is.na(env.df$NH4)),] # This is to get testRes.R and testRes.P for NH4!
#  // Correlation plot:
testRes.R <- as.data.frame(cor(tmp.df))
testRes.R.NH4 <- as.data.frame(cor(tmp.df_NH4))
testRes.R$NH4 <- testRes.R.NH4$NH4
testRes.R[15,] <- as.matrix(testRes.R.NH4[15,])
testRes.P <- cor.mtest(tmp.df, conf.level = 0.95)
for (i in c(1:ncol(testRes.P$p))) {
  testRes.P$p[i,i] = 1
}
corrplot(as.matrix(testRes.R), p.mat = testRes.P$p, method = 'circle', type = 'lower', 
         insig = 'blank', pch.cex = 0.8, tl.col='black', diag = T, cl.pos = 'n',
         col = colorRampPalette(c("#0072B2", "white", "#D55E00"))(200))


# 2) lm plots:
require(gridExtra)
gp1 <- ggplot(env.df, aes(x=Lat, y=CN)) + geom_point(colour = "white") + geom_smooth(method="lm", col="white") +
  theme(panel.border = element_rect(fill='transparent', linewidth = 1, linetype = 'solid',colour = "white"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_line(linetype="solid", linewidth=0.3, colour = "white"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'))

gp2 <- ggplot(env.df, aes(x=Lat, y=HFP)) + geom_point(colour = "white") + geom_smooth(method="lm", col="white") +
  theme(panel.border = element_rect(fill='transparent', linewidth = 1, linetype = 'solid',colour = "white"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_line(linetype="solid", linewidth=0.3, colour = "white"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'))

gp3 <- ggplot(env.df, aes(x=Lat, y=Temp)) + geom_point(colour = "white") + geom_smooth(method="lm", col="white") +
  theme(panel.border = element_rect(fill='transparent', linewidth = 1, linetype = 'solid',colour = "white"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_line(linetype="solid", linewidth=0.3, colour = "white"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'))

gp4 <- ggplot(env.df, aes(x=Temp, y=HFP)) + geom_point(colour = "white") + geom_smooth(method="lm", col="white") +
  theme(panel.border = element_rect(fill='transparent', linewidth = 1, linetype = 'solid',colour = "white"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_line(linetype="solid", linewidth=0.3, colour = "white"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'))

grid.arrange(gp1, gp2, gp3, gp4, nrow=2, ncol=3)


################################################################################
# Fig 2 - Heatmap and bar graph
################################################################################
rel.card.a_fig2 <- rel.card.a
# Add Total abundance:
rel.card.a_fig2$Total <- tot
pheatmap(sqrt(rel.card.a_fig2), cluster_cols = F, cluster_rows = F, gaps_col = 29,
         cellheight = 11, color = colorRampPalette(c("white","grey15"))(100))

# Add annotations of Diversity indices:
#all(rownames(env.df) == names(tot)) # [1] TRUE
div.df <- data.frame(Total = tot, Sobs = sobs, Shannon=shan, 
                     PC1 = pca$CA$u[,1], PC2 = pca$CA$u[,2],
                     row.names = rownames(env.df))
par(mfrow=c(3,1))
barplot(div.df$PC2, las=2, main="PC2")
barplot(div.df$PC1, las=2, main="PC1")
barplot(div.df$Shannon, las=2, main="shannon")
#barplot(div.df$Total, names.arg = rownames(div.df), las=2, main="total")
par(mfrow=c(1,1))


################################################################################
# Fig 2 - Taxonomic annotations on Heatmap
################################################################################
taxa.card$Tub <- factor(taxa.card$Tub, levels = rownames(env.df)[order(env.df$Lat)])
# Set up drug classes:
taxa.card$DrugClass <- ARO_Refs[taxa.card$ARO_id, "Drug_class"]
drugclasses <- c("aminoglycoside", "fluoroquinolone", "glycopeptide", "rifamycin", "tetracycline")
for (c in drugclasses) {
  taxa.card[[c]] <- 0
  taxa.card[grep(c, taxa.card$DrugClass), c] <- 1
}
taxa.card$others <- 0
taxa.card[which(rowSums(taxa.card[,drugclasses]) == 0), "others"] <- 1
# Set up phyla (the most abundant 10 phyla):
top10.p <- c("Acidobacteriota","Actinomycetota","Armatimonadota",
             "Bacteroidota","Chloroflexota","Cyanobacteriota",#"Deinococcota",
             "Myxococcota","Patescibacteria","Pseudomonadota","Verrucomicrobiota")
taxa.card[which(!(taxa.card$phylum%in%top10.p)),"phylum"]="Others"
taxa.card$phylum <- factor(taxa.card$phylum, levels = c("Others",rev(top10.p)))
top10.p.col <- rev(c("#466791","#dc4555","#953ada","#4fbe6c","#ce49d3",
                     "#a7b43d","#d49f36","#552095","#507f2d","#5a51dc","grey"))
#------------------
# // Phylum-by-ARG (Top annotation panel) --
ggplot(taxa.card, aes(x = Best_Hit_ARO_short, fill = phylum)) + 
  geom_bar(position = 'fill') + scale_fill_manual(values = top10.p.col) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
# // Phylum-by-Site (Left annotation panel) --
ggplot(taxa.card, aes(x = Tub, fill = phylum)) + 
  geom_bar(position="fill") + scale_fill_manual(values = top10.p.col) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# // Phylum-by-DrugClass (Bottom left annotation panel)--
gp1 <- ggplot(taxa.card, aes(x = aminoglycoside, fill = phylum)) + 
  geom_bar(position="fill",show.legend = F) + scale_fill_manual(values = top10.p.col) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gp2 <- ggplot(taxa.card, aes(x = fluoroquinolone, fill = phylum)) + 
  geom_bar(position="fill",show.legend = F) + scale_fill_manual(values = top10.p.col) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gp3 <- ggplot(taxa.card, aes(x = glycopeptide, fill = phylum)) + 
  geom_bar(position="fill",show.legend = F) + scale_fill_manual(values = top10.p.col) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gp4 <- ggplot(taxa.card, aes(x = rifamycin, fill = phylum)) + 
  geom_bar(position="fill",show.legend = F) + scale_fill_manual(values = top10.p.col) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gp5 <- ggplot(taxa.card, aes(x = tetracycline, fill = phylum)) + 
  geom_bar(position="fill",show.legend = F) + scale_fill_manual(values = top10.p.col) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gp6 <- ggplot(taxa.card, aes(x = others, fill = phylum)) + 
  geom_bar(position="fill",show.legend = F) + scale_fill_manual(values = top10.p.col) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
require(gridExtra)
grid.arrange(gp1, gp2, gp3, gp4, gp5, gp6, ncol=6, nrow=1) # Then clean up in Inkscape


################################################################################
# Fig S1 - Phylum-level taxonomic distribution of ARGs
################################################################################
ggplot(taxa.card, aes(x = Best_Hit_ARO_short, fill = phylum)) + 
  geom_bar() + scale_fill_manual(values = top10.p.col) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


################################################################################
# Fig S2 - Phylum by site
################################################################################
aro10 <- c("adeF", "rsmA", "vanYA", "vanYB", "vanYF", "vanYG","vanYM", "vanWG", "vanTG", "vanWI") # colnames(sel.rel.card.a)
# Draw stacked barplot:
require(gridExtra)
top10.p.col.z <- c(top10.p.col,'black')
gps <- list()
for (aro in aro10) {
  tmp.tmp <- data.frame(Tub=rownames(env.df), phylum=rep("Zzz",29))
  tmp.df1 <- taxa.card[which(taxa.card$Best_Hit_ARO_short==aro),c("Tub", "phylum")]
  tmp.df1 <- rbind(tmp.tmp, tmp.df1) # Add pseudo phylum 'Zzz' a pseodo count 1 to all samples!
  tmp.df1$phylum <- factor(tmp.df1$phylum, levels = c("Others",rev(top10.p),'Zzz'))
  tmp.df1$Tub <- factor(tmp.df1$Tub, levels = rownames(env.df))
  tmp.col <- top10.p.col.z[which(summary(tmp.df1$phylum)>0)]
  gps[[aro]] <- ggplot(tmp.df1, aes(x = Tub, fill = phylum)) + ggtitle(aro) + 
    geom_bar() + scale_fill_manual(values = tmp.col) +
    labs(y = "Count") + guides(fill="none") + scale_y_continuous(limits = c(1, NA)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
grid.arrange(gps[[1]],gps[[2]],gps[[3]],gps[[4]],gps[[5]],
             gps[[6]],gps[[7]],gps[[8]],gps[[9]],gps[[10]], ncol=2, nrow=5)
# As a pseudo phylum with count 1 was added to the dataset, so that every sample
# is displayed on the plot, remember to adjust the y-axis value in Inkscape


################################################################################
# Fig 3 - Heatmap Family annotations
################################################################################
# card.f.df <- summarise_table(taxa.card, 'd_p_c_o_f', 'Best_Hit_ARO_short', summarise_by_column = NULL,
#                              method = 'n') # This function is in 8_BGC_analysis.R
# write.csv(card.f.df,"CARD-mmseq2taxonomy_TableS2.csv") # Then add col for OTU ids!
load_singlem_cmm("../Data/CARD-mmseq2taxonomy_TableS2.csv") # You can put any number for rarefaction level here - we are not using rarefied dataset
rm(otu.arc, otu.bac, otu.bac.r, rel.otu.all, rel.otu.arc, rel.otu.bac)
card.f.df <- otu.all/100 # As load_singlem_cmm() function multiplies 100 for the need of singleM data
std.card.f.df <- decostand(card.f.df, method= 'standardize') # Scaled to mean and std of each ARG!
# Select the Taxa that have total abundance over 0.25%
sel.std.card.f.df <- std.card.f.df[,which(colSums(card.f.df/18.31)>0.25)]
draw_cmm_heatmap(sel.std.card.f.df, taxonomy_df = taxonomy.all) 

# Decorate row/col heatmap:
arg_tot <- data.frame(ARG_Total_68=rowSums(card.f.df[,colnames(sel.std.card.f.df)])/18.31,
                      ARG_Total_All=rowSums(card.f.df)/18.31,
                      row.names = rownames(card.f.df))
fam_tot <- data.frame(Family_Total=colSums(card.f.df[,colnames(sel.std.card.f.df)])/18.31,
                      Family_Total2=colSums(card.f.df[,colnames(sel.std.card.f.df)])/18.31,
                      row.names = colnames(sel.std.card.f.df))
pheatmap(arg_tot, cluster_rows = F, cluster_cols = F, 
         color = colorRampPalette(c("white","darkgreen"))(100))
pheatmap(fam_tot, cluster_rows = F, cluster_cols = F, 
         color = colorRampPalette(c("white","orange"))(100))

# Then merger the above plot in Inkscape.


################################################################################
# Fig 5a - Stepwise linear correlations
################################################################################
# Fig 5a - Corrplot
fig4_df <- data.frame(Total=tot, Shan=shan, PC1 = pca$CA$u[,1], PC2 = pca$CA$u[,2], 
                      row.names = names(tot))
sel.rel.card.a <- rel.card.a[,which(colSums(decostand(rel.card.a, "pa"))>10)]
fig4_df <- cbind(sel.rel.card.a, fig4_df)
res.df <- auto_lm_forward_new(fig4_df, c("HFP", env_params), env.df, hm_title = "Linear (fwd)")
corrplot(as.matrix(res.df$t_values), 
         method = 'circle', type = 'full', insig = 'blank', pch.cex = 0.8, 
         tl.col='black', diag = T, is.corr = FALSE, 
         col = colorRampPalette(c("#0072B2", "white", "#D55E00"))(200)) # Also scale the color to t values

#### Then test if Latitude is significant when adding to parsimonous model
#### (Table S7):
summary(lm(fig4_df$Total ~ Temp + CN + pH + Lat, data = env.df))
summary(lm(fig4_df$Shan ~ Mg + Lat, data = env.df))
summary(lm(fig4_df$`otr(A)` ~ CN + pcMC + Lat, data = env.df))
summary(lm(fig4_df$vanYF ~ Temp + pH + Lat, data = env.df))
summary(lm(fig4_df$vanYM ~ Temp + Mg + Lat, data = env.df))
summary(lm(fig4_df$vanTG ~ Temp + Lat, data = env.df))


################################################################################
# Fig 5b - PERMANOVA & RDA
################################################################################
# To get the best model of PERMANOVA (with different data transformation):
res <- auto_adonis_forward(list(raw=rel.card.a, pa=decostand(rel.card.a,"pa"),
                                hell=decostand(rel.card.a, "hellinger"),
                                sqrt=sqrt(rel.card.a), std=decostand(rel.card.a,"standardize")), 
                           env_params, env.df, trans.method = NULL)
#### Best model is ~ pH + Temp + DOC + CN :
adonis2(sqrt(rel.card.a) ~ pH + Temp + DOC + CN, data = env.df, na.action = na.omit, method = "euc") 

#### Test if Latitude contribute to Beta:
adonis2(sqrt(rel.card.a) ~ pH + Temp + CN + DOC + Lat, data = env.df, na.action = na.omit, method = "euc") # P(lat) = 0.0014
adonis2(sqrt(rel.card.a) ~ Lat, data = env.df, na.action = na.omit, method = "euc") # p = 0.017
#         RDA:
permutest(rda(sqrt(rel.card.a)~ pH + Temp + CN + DOC + Lat, env.df), by = "terms")
# Partial RDA:
permutest(rda(sqrt(rel.card.a) ~ Lat + Condition(pH + Temp + CN + DOC), env.df), by = "terms") # p =0.01

# Fig 5b - RDA constraint to PERMANOVA fwd selected model
datapca <- rda(decostand(rel.card.a,'hellinger') ~ pH + Temp + DOC + CN, env.df, na.action=na.omit) # This model is selected from data that NH4 for Tub 58 is NA!!!
# Set up colours and shape/sizes of points for the ordination plot:
masat_ph.df <- read.csv("../Data/masat_color_code.csv", header = T, row.names = 1) #all(rownames(masat_col.df) == rownames(env.df)) # [1] TRUE
masat_col <- masat_ph.df$masat_col  # This will be the color for samples
ph_grp <- as.factor(masat_ph.df$pH_group)
ph_grp_cex <- c(0.8, 1.2, 1.6, 2, 2.4, 2.8)   # 'ph_grp_cex[ph_grp]' will be the sampel point sizes

draw_rda_plot(datapca, env.df = env.df[,c("Temp","pH","CN","DOC")], sd.val = 2, 
              text.cex = 0.8, site.col = masat_col, site.cex = ph_grp_cex[ph_grp])
legend("topleft", legend = levels(ph_grp), pch = 1, pt.bg = "grey15", 
       pt.cex = ph_grp_cex)
legend('topleft', legend = factor(c(-12:0)), pch = 15, pt.cex = 1, 
       col = colorRampPalette(c("#005CE6","#8FE5FF"))(13))
# Then adjust the legend and labels in Inkscape.


################################################################################
# Fig S4 - Linear correlations
################################################################################
# Fig S4a,b - linear plot
df.tmp <- env.df[,c('Lat','Long','Alt','Temp','pH', 'CN', 'Mg')]
#   // Add total + shan
df.tmp$Total <- tot
df.tmp$Shannon <- shan
require(gridExtra)
gp1 <- ggplot(df.tmp, aes(x=Temp, y=Total)) + geom_point() + geom_smooth(method="lm", col="blue")
gp2 <- ggplot(df.tmp, aes(x=pH, y=Total)) + geom_point() + geom_smooth(method="lm", col="blue")
gp3 <- ggplot(df.tmp, aes(x=CN, y=Total)) + geom_point() + geom_smooth(method="lm", col="blue")
#gp4 <- ggplot(df.tmp, aes(x=Temp, y=Sobs)) + geom_point() + geom_smooth(method="lm", col="blue")
#gp5 <- ggplot(df.tmp, aes(x=Mg, y=Sobs)) + geom_point() + geom_smooth(method="lm", col="blue")
gp6 <- ggplot(df.tmp, aes(x=Mg, y=Shannon)) + geom_point() + geom_smooth(method="lm", col="blue")
grid.arrange(gp1, gp2, gp3, gp6, nrow=2, ncol=2)

# Fig S4d Pairwise linear regression
env_params1 <- c("Lat","Long","Alt","Temp","Ca","Cl","CN","Cu","DOC","EC","Fe","K","Mg","Mn","NH4","NO3","P","pcMC","pH","PO4","SO4","TC","TN","Zn","HFP")
figS2a_df <- fig4_df
colnames(figS2a_df)[2] <- "otrA"
figS2a.res.df <- auto_lm_forward(figS2a_df, env_params1, env.df, corr.test = T,
                                 hm_title = "Pearson Correlation (t value)")
corrplot(as.matrix(figS2a.res.df$Corr.Summary), 
         method = 'circle', type = 'full', insig = 'blank', pch.cex = 0.8, 
         tl.col='black', diag = T, is.corr = FALSE, 
         col = colorRampPalette(c("#0072B2", "white", "#D55E00"))(200))


################################################################################
# Fig S4c - PCA
################################################################################
datapca <- rda(sqrt(rel.card.a))
ev <- envfit(datapca, env.df[,env_params], choices = 1:2, scaling = "symmetric",permutations = 5000, na.rm = TRUE)
plot(datapca, scaling = "symmetric", display = "sites",cex=1,xlim=c(-1.5,3), 
     xlab=paste("PC1 (",axis.percent(datapca)[[1]],"%)",sep=""),
     ylab=paste("PC2 (",axis.percent(datapca)[[2]],"%)",sep=""))
points(datapca,dis='sp',pch=4,col='grey',cex=0.6,scaling=3)
points(datapca, dis='sites',pch=19,col='#c90076', cex=1.5,scaling=3)
sd.val=2
points(scores(datapca,scaling=3)$sp[which(scores(datapca,scaling=3)$sp[,1] > sd.val * sd(scores(datapca)$sp[,1])),],pch=4,col='black',cex=0.6)	
points(scores(datapca,scaling=3)$sp[which(scores(datapca,scaling=3)$sp[,1] < 0 - (sd.val * sd(scores(datapca)$sp[,1]))),],pch=4,col='black',cex=0.6) 
points(scores(datapca,scaling=3)$sp[which(scores(datapca,scaling=3)$sp[,2] > sd.val * sd(scores(datapca)$sp[,2])),],pch=4,col='black',cex=0.6)
points(scores(datapca,scaling=3)$sp[which(scores(datapca,scaling=3)$sp[,2] < 0 - (sd.val * sd(scores(datapca)$sp[,2]))),],pch=4,col='black',cex=0.6) 
ev2 <- envfit(datapca, env.df[,c("Temp","pH","Mg","DOC", "NH4")], choices = 1:2, scaling = "symmetric",permutations = 5000, na.rm = TRUE)
plot(ev2,add=TRUE,cex=.8)
# plot(datapca, scaling = "symmetric", type='t') 
# Then edit in Inkscape

# sd.val=3
# highlight_otus <- c(which(scores(datapca,scaling=3)$sp[,1] > sd.val * sd(scores(datapca)$sp[,1])),
#                     which(scores(datapca,scaling=3)$sp[,1] < 0 - (sd.val * sd(scores(datapca)$sp[,1]))),
#                     which(scores(datapca,scaling=3)$sp[,2] > sd.val * sd(scores(datapca)$sp[,2])),
#                     which(scores(datapca,scaling=3)$sp[,2] < 0 - (sd.val * sd(scores(datapca)$sp[,2]))))
# scores(datapca,scaling=3)$sp[highlight_otus,]


################################################################################
# Fig S5 - DeepARG
################################################################################
# Load data:
deeparg.all <- read.csv('../Data/deeparg_assembly_all.tsv', sep = '\t', header = T)
# **** Pick the dataset you wanna use and uncomment **** #
# // (1) 'True' hits (probablity > 0.8): 
# arg.df <- summarise_table(deeparg.all[which(deeparg.all$FileType=="Probable"),], group_for_row = 'Sample', group_for_col = 'predicted_ARG.class', 
#                           summarise_by_column = 'counts', method = 'sum', env.df = env.df)
# // (2) Potential hits (probablity < 0.8): 
# arg.df <- summarise_table(deeparg.all[which(deeparg.all$FileType=="Potential"),], group_for_row = 'Sample', group_for_col = 'predicted_ARG.class', 
#                           summarise_by_column = 'counts', method = 'sum', env.df = env.df)
# // (3) All hits (both 'True' and 'Potential'): 
# arg.df <- summarise_table(deeparg.all,group_for_row = 'Sample', group_for_col = 'predicted_ARG.class', 
#                           summarise_by_column = 'counts', method = 'sum', env.df = env.df)
# // (4) All hits, dereplicated, and excl CARD hits: 
# deeparg.all.derep <- deeparg.all[which(deeparg.all$FileType!="Potential_loser"),]
# arg.df <- summarise_table(deeparg.all.derep[which(is.na(deeparg.all.derep$In_CARD)),], group_for_row = 'Sample', group_for_col = 'predicted_ARG.class', 
#                           summarise_by_column = 'counts', method = 'sum', env.df = env.df)
# // (5) Final: hits > 0.8 probablity, dereplicated, and excl CARD hits:
deeparg.all.derep <- deeparg.all[which(deeparg.all$FileType=="Probable"),]
arg.df <- summarise_table(deeparg.all.derep[which(is.na(deeparg.all.derep$In_CARD)),], group_for_row = 'Sample', group_for_col = 'predicted_ARG.class', 
                          summarise_by_column = 'counts', method = 'sum', env.df = env.df)
arg.df <- arg.df[rownames(env.df),]
# Number of ARGs too low. Skip rarefaction!
rel.arg.df <- 1000000*arg.df/env.df$assembly_gene_count # ARG per million genes
tot.deeparg <- rowSums(rel.arg.df)

# 1) Overview of Deeparg results:
# pheatmap(t(rel.arg.df), cluster_rows = F, cluster_cols = F)


# 2) Plot: Pairwise linear regression (Fig S5a)
env_params1 <- c("Lat","Long","Alt","Temp","pH","EC","TC","TN","CN",
                 "Ca","Cu","Fe","K","Mg","Mn","P","Zn",
                 "DOC","NO3","NH4","SO4","PO4","Cl","pcMC",
                 "HFP")
fig_df <- rel.arg.df
fig_df$Total <- tot.deeparg
fig.res.df <- auto_lm_forward(fig_df, env_params = env_params1, env_df = env.df, 
                              test_method = "F", corr.test = T, 
                              hm_title = "Pearson Correlation (t value)")
corrplot(as.matrix(fig.res.df$Corr.Summary), 
         method = 'circle', type = 'full', insig = 'blank', pch.cex = 0.8, 
         tl.col='black', diag = T, is.corr = FALSE, 
         col = colorRampPalette(c("#0072B2", "white", "#D55E00"))(200))

# 3) Plot: linear regression of factors significant with 'Total' (Fig S5b-e):
require(gridExtra)
fig_df <- data.frame(Total=fig_df$Total, MASAT = env.df$Temp, DOC = env.df$DOC, Mg = env.df$Mg, NH4 = env.df$NH4)
gp1 <- ggplot(fig_df, aes(x=MASAT, y=Total)) + geom_point(colour = "black") + geom_smooth(method="lm", col="blue")
gp2 <- ggplot(fig_df, aes(x=DOC, y=Total)) + geom_point(colour = "black") + geom_smooth(method="lm", col="blue")
gp3 <- ggplot(fig_df, aes(x=Mg, y=Total)) + geom_point(colour = "black") + geom_smooth(method="lm", col="blue")
gp4 <- ggplot(fig_df, aes(x=NH4, y=Total)) + geom_point(colour = "black") + geom_smooth(method="lm", col="blue")
grid.arrange(gp1, gp2, gp3, gp4, nrow=2, ncol=2)
