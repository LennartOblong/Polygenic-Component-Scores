mydata <- as.data.frame(read.table("/projects/0/einf2700/oblongl/melodic-mix"))
IDPs <- read.table("/projects/0/einf2700/oblongl/all-idp-2240-code.txt")



#fig_heat <- mydata %>% 
#  ggplot( aes(fill=coefficient, x=comp, y=interaction(sub_scale, scale))) +
#  geom_tile(aes(fill = coefficient, width=0.9, height=0.9)) + theme(text=element_text(family="Arial", size =9),
#                                                             axis.text.x = element_text(angle = 45, hjust= 0.2), 
#                                                             axis.ticks = element_blank(),
#                                                             legend.text = element_text(size = 8), 
#                                                             legend.key.size = unit(0.3, "cm")) +
#  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, limits =c(-0.4, 0.4), na.value = "#CCCCCC",oob = scales::squish) +
#  labs(x="", y="") + scale_x_discrete(position = "top") +
#scale_y_discrete(limits=rev, guide = "axis_nested") + theme(ggh4x.axis.nestline = element_line(linetype = 1, colour = "red")) +
#  geom_text(aes(label=plevel)) 




fig_heat_file <- file.path(outdir, "sig_unadj_heatmap.png")
ggsave(fig_heat_file, fig_heat, dpi=300, width=12, height = 8, units = "cm")

fig_heat

