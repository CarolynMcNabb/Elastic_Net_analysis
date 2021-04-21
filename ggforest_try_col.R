ggforest <- function (x, intercept = F, labels = NULL, hetero = "none", palette = wesanderson::wes_palettes$Darjeeling1) 
{
  studies <- x$metas[[1]]$slab
  x <- x$metas
  require("ggplot2")
  rma2df = function(x) {
    df <- as.data.frame(data.table::rbindlist(lapply(x, function(x) {
      rbind(data.frame(Study = " ", LogFC = x$b, 
                       CILB = x$ci.lb, CIUB = x$ci.ub, p = x$pval, group = "Meta-\nanalysis", 
                       I2 = x$I2, Qp = x$QEp, stringsAsFactors = FALSE), 
            data.frame(Study = x$slab, LogFC = x$yi, CILB = x$yi - 
                         2 * sqrt(x$vi), CIUB = x$yi + 2 * sqrt(x$vi), 
                       p = x$pval, group = "Elastic Net Regression models", I2 = NA, 
                       Qp = NA, stringsAsFactors = FALSE))
    })))
    df$predictor <- rep(names(x), each = nrow(df)/length(x))
    df
  }
  remresdf = rma2df(x)
  remresdf <- transform(remresdf, interval = CIUB - CILB)
  remresdf <- transform(remresdf, RelConf = 1/interval)
  if (!intercept) {
    remresdf <- remresdf[remresdf$predictor != "(Intercept)", 
                         ]
    remresdf$predictor <- factor(remresdf$predictor, levels = unique(remresdf$predictor))
    if (!is.null(labels)) {
      remresdf$predictor <- factor(remresdf$predictor, 
                                   labels = labels)
    }
  }
  p = ggplot(remresdf, aes(LogFC, factor(Study, c(rev(studies), 
                                                  " ")),
                           colour = factor(Study, c(studies, 
                                                    " ")),
                                                    xmax = CIUB, xmin = CILB, shape = group 
                           )) + geom_vline(xintercept = 0, linetype = 2, 
                                                             alpha = 0.75) + geom_errorbarh(alpha = 0.5, color = "black", 
                                                                                            height = 0.1) + geom_point(aes(size = RelConf)) + scale_size(range = c(2, 
                                                                                                                                                                   5), guide = FALSE) + geom_point(data = subset(remresdf, 
                                                                                                                                                                                                                 Study == " "), size = 3) + theme_pubclean() + theme(text = element_text(size = 11)) + 
    facet_grid(group ~ ., scales = "free", space = "free_y") + 
    labs(x = "Beta weight (95% CI) of \npredicted distance on \nobserved distance", y = NULL) + scale_colour_manual(values = palette) + 
    guides(colour = F, shape = F) + scale_x_continuous(breaks = scales::pretty_breaks(3), 
                                                       labels = function(x) round(as.numeric(x), digits = 2)) + 
    expand_limits(x = c(-0.02, 0.02))+
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
          axis.ticks.y = element_blank())
  if (hetero == "I2") {
    p <- p + geom_label(aes(x = Inf, y = -Inf, hjust = 1, 
                            vjust = 0, label = ifelse(is.na(I2), I2, paste0("I² = ", 
                                                                            round(I2, 2), "%"))), colour = "black", label.size = 0)
  }
  else if (hetero == "Qp") {
    p <- p + geom_label(aes(label = ifelse(is.na(Qp), NA, 
                                           ifelse(Qp < 0.001, "***", ifelse(Qp < 0.01, "**", 
                                                                            ifelse(Qp < 0.05, "*", ""))))), position = position_nudge(y = 0.2), 
                        colour = "black", label.size = 0)
  }
  return(p)
}
  
  
