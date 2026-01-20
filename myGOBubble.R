myGOBubble <- function (data, display, title, colour, labels, ID = T, table.legend = T, 
                        table.col = T, bg.col = F) 
{
  zscore <- adj_pval <- category <- count <- id <- term <- NULL
  if (missing(display)) 
    display <- "single"
  if (missing(title)) 
    title <- ""
  if (missing(colour)) 
    cols <- c("chartreuse4", "brown2", "cornflowerblue")
  else cols <- colour
  if (missing(labels)) 
    labels <- 5
  if (bg.col == T & display == "single") 
    cat("Parameter bg.col will be ignored. To use the parameter change display to 'multiple'")
  colnames(data) <- tolower(colnames(data))
  if (!"count" %in% colnames(data)) {
    rang <- c(5, 5)
    data$count <- rep(1, dim(data)[1])
  }
  else {
    rang <- c(1, 10)
  }
  data$adj_pval <- -log(data$adj_pval, 10)
  sub <- data[!duplicated(data$term), ]
  g <- ggplot(sub, aes(zscore, adj_pval, fill = category,size = count)) + 
    labs(title = title, x = "z-score", y = "-log (adj p-value)") + 
    geom_point(shape = 21, col = "black", alpha = 1/2) + 
    geom_hline(yintercept = 1.3, col = "orange") + scale_size(range = rang)
  if (!is.character(labels)) 
    sub2 <- subset(sub, subset = sub$adj_pval >= labels)
  else sub2 <- subset(sub, sub$id %in% labels | sub$term %in% 
                        labels)
  if (display == "single") {
    g <- g + scale_fill_manual("Category", values = cols#, labels = c("BP", "CC", "KEGG", "MF")
    ) + 
      theme(legend.position = "bottom") + 
      annotate("text", x = min(sub$zscore) + 0.2, y = 1.4, 
               label = "Threshold", colour = "orange", size = 4)
    if (ID) 
      g <- g + geom_text_repel(data = sub2, aes(x = zscore, y = adj_pval, 
                                                label = id), size = 4)
    else g <- g + geom_text_repel(data = sub2, aes(x = zscore, 
                                                   y = adj_pval, label = term), size = 4)
    if (table.legend) {
      if (table.col) 
        table <- draw_table(sub2, col = cols)
      else table <- draw_table(sub2)
      g <- g + theme(axis.text = element_text(size = 10), 
                     axis.line = element_line(colour = "black"), 
                     axis.ticks = element_line(colour = "black"), 
                     axis.title = element_text(size = 10, face = "plain"), 
                     panel.background = element_blank(), panel.grid.minor = element_blank(), 
                     #panel.grid.major = element_line(colour = "grey80"),
                     legend.key = element_rect(fill = "white"),
                     legend.direction = "vertical",legend.position = "right",
                     plot.background = element_blank())
      graphics::par(mar = c(0.1, 0.1, 0.1, 0.1))
      grid.arrange(g, table, ncol = 2)
    }
    else {
      g + theme(axis.text = element_text(size = 10), axis.line = element_line(colour = "black"), 
                axis.ticks = element_line(colour = "black"), 
                axis.title = element_text(size = 10, face = "plain"), 
                panel.background = element_blank(), panel.grid.minor = element_blank(), 
                #panel.grid.major = element_line(colour = "grey80"), 
                legend.key = element_rect(fill = "white"),
                legend.direction = "vertical",legend.position = "right",
                plot.background = element_blank())
    }
  }
  else {
    if (bg.col) {
      dummy_col <- data.frame(category = c("BP", "CC", 
                                           "MF"), adj_pval = sub$adj_pval[1:3], zscore = sub$zscore[1:3], 
                              size = 1:3, count = 1:3)
      g <- g + geom_rect(data = dummy_col, aes(fill = category), 
                         xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
                         alpha = 0.1) + facet_grid(. ~ category, space = "free_x", 
                                                   scales = "free_x") + scale_fill_manual(values = cols, 
                                                                                          guide = "none")
    }
    else {
      g <- g + facet_grid(. ~ category, space = "free_x", 
                          scales = "free_x") + scale_fill_manual(values = cols, 
                                                                 guide = "none")
    }
    if (ID) {
      g + geom_text(data = sub2, aes(x = zscore, y = adj_pval, 
                                     label = id), size = 5) + 
        theme(axis.title = element_text(size = 14, 
                                        face = "plain"), axis.text = element_text(size = 14),
              axis.line = element_line(colour = "black"),
              axis.ticks = element_line(colour = "black"),
              panel.border = element_rect(fill = "transparent",colour = "black"), 
              panel.background = element_blank(),
              panel.grid = element_blank(), 
              plot.background = element_blank())
    }
    else {
      g + geom_text(data = sub2, aes(x = zscore, y = adj_pval, 
                                     label = term), size = 5) + 
        theme(axis.title = element_text(size = 14,face = "plain"), 
              axis.text = element_text(size = 14), 
              axis.line = element_line(colour = "black"),
              axis.ticks = element_line(colour = "black"),
              panel.border = element_rect(fill = "transparent", colour = "black"), 
              panel.background = element_blank(),
              panel.grid = element_blank(), 
              plot.background = element_blank())
    }
  }
}