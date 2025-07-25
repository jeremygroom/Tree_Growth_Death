


# Finding domain summaries
domain.summaries <- domain.index %>% 
  map(\(d) domain.sum.fcn(bootstrap_results, d)) %>%
  do.call(rbind, .) %>%
  arrange(Species, Domain)

saveRDS(list(bootstrap_results = bootstrap_results, domain.summaries = domain.summaries), file = paste0(save.loc, "Domain_Analysis_Output.RDS"))



 # - BREAK: all domain types included and saved above, only 6 quantiles continued below.





## -- Figure code 



# Read the data files
#dvs_data <- read.csv("DvS.csv")
#d_lmh_data <- read.csv("D_LMH.csv")
#s_lmh_data <- read.csv("S_LMH.csv")



diff.panel.fcn <- function(diff.dat, remove.y, fig.title, lab.right) {

ggplot(diff.dat, aes(y = as.numeric(species_label) + y.offset)) + 
  geom_vline(xintercept = 0, linewidth = 0.5) + 
  geom_segment(aes(x = LCI.95, xend = UCI.95, 
                   y = as.numeric(species_label) + y.offset, 
                   yend = as.numeric(species_label) + y.offset,
                   color = Domain), 
               linewidth = 0.8) + 
  geom_point(aes(x = Means, y = as.numeric(species_label) + y.offset,
                 fill = ifelse(significant, "black", "white"),
                 color = Domain),
             size = 2, shape = 21, stroke = 0.8) +
  #scale_y_continuous(breaks = 1:length(levels(DvS2$species_label)),
  #                    labels = levels(DvS2$species_label)) +
  scale_y_continuous(breaks = 1:length(levels(diff.dat$species_label)),
                     labels = levels(diff.dat$species_label)) +
  scale_color_viridis(discrete = TRUE, name = "Domain", option = "B", begin = 0.2, end = 0.6) +
  scale_fill_identity() + 
  labs(x = NULL, y = NULL, title = fig.title) +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10, face = "italic"),
        legend.position = "inside",
        legend.position.inside = c(if(lab.right) 0.8 else 0.2, 0.95),
        legend.background = element_rect(fill = "white", color = "gray80"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"))  + 
    {if(remove.y) {
      theme(axis.text.y = element_blank())
    }}
}

p1 <- diff.panel.fcn(DvS2, remove.y = FALSE, fig.title = "Dry vs. Stable", lab.right = FALSE)
p2 <- diff.panel.fcn(D_LMH2, remove.y = TRUE, fig.title = "Dry, High/Med/Low", lab.right = FALSE)
p3 <- diff.panel.fcn(S_LMH2, remove.y = TRUE, fig.title = "Stable, High/Med/Low", lab.right = TRUE)

grand.x.lab <- ggdraw() + draw_label("Annual Mortality Rate", x = 0.6, y = 0.5) + theme_bw() + theme(rect = element_blank())

mort.plt <- plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(0.9, 0.5, 0.5)) 

mort.plt2 <- plot_grid(mort.plt, grand.x.lab, 
                      ncol = 1, 
                      rel_heights = c(1, 0.03)) 

ggsave("MortPlot.png", plot = mort.plt2, device = "png", width = 10, height = 10, units = 'in')




