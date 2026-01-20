draw_table1 <- function(data, col){
  id <- term <- NULL
  colnames(data) <- tolower(colnames(data))
  if (missing(col)){
    tt1 <- ttheme_default(base_size = 6)
  }else{
    text.col <- c(rep(col[1], sum(data$category == 'BP')), rep(col[2], sum(data$category == 'CC')), rep(col[3], sum(data$category == 'MF')))
    tt1 <- ttheme_minimal(base_size = 6,
                          core = list(bg_params = list(fill = text.col, col=NA, alpha= 1/3)), 
                          colhead = list(fg_params = list(col = "black")))
  }
  table <- tableGrob(subset(data, select = c(id, term)), cols = c('ID', 'Description'), rows = NULL, theme = tt1)
  return(table)
}

myGOCircle <- function (data, title, nsub, rad1, rad2, table.legend = T, zsc.col, 
                        lfc.col, label.size, label.fontface) 
{
  xmax <- y1 <- zscore <- y2 <- ID <- logx <- logy2 <- logy <- logFC <- NULL
  if (missing(title)) 
    title <- ""
  if (missing(nsub)) 
    if (dim(data)[1] > 10) 
      nsub <- 10
  else nsub <- dim(data)[1]
  if (missing(rad1)) 
    rad1 <- 2
  if (missing(rad2)) 
    rad2 <- 3
  if (missing(zsc.col)) 
    zsc.col <- c("red", "white", "blue")
  if (missing(lfc.col)) 
    lfc.col <- c("firebrick1","cornflowerblue")
  else lfc.col <- rev(lfc.col)
  if (missing(label.size)) 
    label.size = 5
  if (missing(label.fontface)) 
    label.fontface = "bold"
  data$adj_pval <- -log(data$adj_pval, 10)
  suby <- data[!duplicated(data$term), ]
  if (is.numeric(nsub) == T) {
    suby <- suby[1:nsub, ]
  }
  else {
    if (strsplit(nsub[1], ":")[[1]][1] == "GO") {
      suby <- suby[suby$ID %in% nsub, ]
    }
    else {
      suby <- suby[suby$term %in% nsub, ]
    }
    nsub <- length(nsub)
  }
  N <- dim(suby)[1]
  r_pval <- round(range(suby$adj_pval), 0) + c(-2, 2)
  ymax <- c()
  for (i in 1:length(suby$adj_pval)) {
    val <- (suby$adj_pval[i] - r_pval[1])/(r_pval[2] - r_pval[1])
    ymax <- c(ymax, val)
  }
  df <- data.frame(x = seq(0, 10 - (10/N), length = N), xmax = rep(10/N - 0.2, N), 
                   y1 = rep(rad1, N), y2 = rep(rad2, N), ymax = ymax, 
                   zscore = suby$zscore, ID = suby$ID)
  # scount <- data[!duplicated(data$term), which(colnames(data) == 
  #                                                "count")][1:nsub]
  #idx_term <- which(!duplicated(data$term) == T)
  scount <- suby$count
  idx_term <- as.numeric(rownames(suby))
  xm <- c()
  logs <- c()
  # for (sc in 1:length(scount)) {
  #   idx <- c(idx_term[sc], idx_term[sc] + scount[sc] - 1)
  #   val <- stats::runif(scount[sc], df$x[sc] + 0.06, (df$x[sc] + df$xmax[sc] - 0.06))
  #   xm <- c(xm, val)
  #   r_logFC <- round(range(data$logFC[idx[1]:idx[2]]), 0) + 
  #     c(-1, 1)
  #   for (lfc in idx[1]:idx[2]) {
  #     val <- (data$logFC[lfc] - r_logFC[1])/(r_logFC[2] - 
  #                                              r_logFC[1])
  #     logs <- c(logs, val)
  #   }
  # }
  # cols <- c()
  # for (ys in 1:length(logs)) cols <- c(cols, ifelse(data$logFC[ys] > 0, "upregulated", "downregulated"))
  cols <- c()
  for (sc in 1:length(scount)) {
    idx <- c(idx_term[sc], idx_term[sc] + scount[sc] - 1)
    val <- stats::runif(scount[sc], df$x[sc] + 0.06, (df$x[sc] + df$xmax[sc] - 0.06))
    xm <- c(xm, val)
    r_logFC <- round(range(data$logFC[idx[1]:idx[2]]), 0) + 
      c(-1, 1)
    for (lfc in idx[1]:idx[2]) {
      val <- (data$logFC[lfc] - r_logFC[1])/(r_logFC[2] - 
                                               r_logFC[1])
      logs <- c(logs, val)
      cols <- c(cols, ifelse(data$logFC[lfc] > 0, "upregulated", "downregulated"))  
    }
  }

  
  dfp <- data.frame(logx = xm, logy = logs, logFC = factor(cols,levels = c("upregulated", "downregulated")), 
                    logy2 = rep(rad2, length(logs)))
  c <- ggplot() + geom_rect(data = df, aes(xmin = x, xmax = x + 
                                             xmax, ymin = y1, ymax = y1 + ymax, fill = zscore), 
                            colour = "white") +                 #中心柱状图
    geom_rect(data = df, aes(xmin = x, xmax = x + xmax, ymin = y2, 
                             ymax = y2 + 1), fill = "gray80") + #灰色扇形
    geom_rect(data = df, 
              aes(xmin = x, xmax = x + xmax, ymin = y2 + 0.5, ymax = y2 + 0.5), colour = "white") + #灰色扇形中间的白线
    geom_rect(data = df, aes(xmin = x,
                             xmax = x + xmax, ymin = y2 + 0.25, ymax = y2 + 0.25), colour = "white") + #灰色扇形四分之一处白线
    geom_rect(data = df, aes(xmin = x,
                             xmax = x + xmax, ymin = y2 + 0.75, ymax = y2 + 0.75), colour = "white") + #灰色扇形四分之三处白线
    geom_text(data = df, aes(x = x +(xmax/2), y = y2 + 1.3, label = ID, angle = 360 - (x = x + (xmax/2))/(10/360)),
              size = label.size, fontface = label.fontface) +    #圈图外围标签
    coord_polar() + labs(title = title) + ylim(1, rad2 + 1.6) + xlim(0, 10) + 
    theme_blank + scale_fill_gradient2("z-score",space = "Lab", low = zsc.col[3], mid = zsc.col[2], 
                                       high = zsc.col[1],guide = guide_colourbar(title.position = "top", title.hjust = 0.5),
                                       breaks = c(min(df$zscore), max(df$zscore)), labels = c("decreasing","increasing")) + 
    theme(panel.spacing = unit(0,"cm"),plot.margin = margin(0,0,0,0),
          legend.position = "bottom",
          legend.background = element_blank(),
          legend.key.size = unit(6, "pt"),
          legend.title = element_text(color = "black", size = 6),
          legend.text = element_text(color = "black", size = 6),
          legend.margin = margin(0,0,0,0),legend.box.margin = margin(0,0,0,0),legend.box.spacing = unit(0,"cm"), 
          legend.spacing = unit(0,"cm"),
          legend.box = "vertical", legend.direction = "horizontal") + 
    
    # geom_point(data = dfp, aes(x = logx, y = logy2 + logy), 
    #            pch = 21, fill = "transparent", colour = "transparent", 
    #            size = 3) + 
    geom_point(data = dfp, aes(x = logx,y = logy2 + logy, colour = logFC), 
               size = 1, alpha = 0.6) + #灰色扇形上的气泡
    scale_colour_manual(values = lfc.col,
                        guide = guide_legend(title.position = "top", title.hjust = 0.5))
  if (table.legend) {
    table <- draw_table1(suby) #绘制表格
    graphics::par(mar = c(0.01, 0.01, 0.01, 0.01))
    grid.arrange(c, table, ncol = 2,widths= c(1,1)) #圈图与表格合并
  }
  else {
    c + theme(plot.background = element_rect(fill = "white"), 
              panel.background = element_rect(fill = "white"))
  }
}