plot_data_one_iso_custom = function (mix, source, discr, filename, plot_save_pdf, plot_save_png) 
{
  if (length(grep("C", mix$iso_names)) == 1) 
    x_label <- expression(paste(delta^13, "C (‰)", sep = ""))
  if (length(grep("N", mix$iso_names)) == 1) 
    x_label <- expression(paste(delta^15, "N (‰)", sep = ""))
  if (length(grep("S", mix$iso_names)) == 1) 
    x_label <- expression(paste(delta^34, "S (‰)", sep = ""))
  if (length(grep("O", mix$iso_names)) == 1) 
    x_label <- expression(paste(delta^18, "O (‰)", sep = ""))
  if (length(grep("SP", mix$iso_names)) == 1) 
    x_label <- expression(paste(delta^15, "N-SP (‰)", 
                                sep = ""))
  if (!exists("x_label")) 
    x_label <- mix$iso_names
  y_data <- 0.5
  y <- rep(y_data, mix$N)
  df <- data.frame(x = mix$data_iso, y = y)
  spacing <- 0.1
  if (!is.na(source$by_factor)) {
    source_linetype <- sort(rep(1:source$n.sources, source$S_factor_levels))
    source_color <- factor(as.numeric(source$S_factor1))
    index <- seq(from = 1, to = 1 + (source$n.sources - 
                                       1) * source$S_factor_levels, by = source$S_factor_levels)
    discr_mu_plot <- array(NA, dim = c(length(source$S_MU[, 
                                                          1]), mix$n.iso))
    discr_sig2_plot <- array(NA, dim = c(length(source$S_MU[, 
                                                            1]), mix$n.iso))
    for (i in 1:source$n.sources) {
      discr_mu_plot[index[i]:(index[i] + source$S_factor_levels - 
                                1), ] <- matrix(rep(discr$mu[i], source$S_factor_levels), 
                                                nrow = source$S_factor_levels, ncol = mix$n.iso, 
                                                byrow = T)
      discr_sig2_plot[index[i]:(index[i] + source$S_factor_levels - 
                                  1), ] <- matrix(rep(discr$sig2[i], source$S_factor_levels), 
                                                  nrow = source$S_factor_levels, ncol = mix$n.iso, 
                                                  byrow = T)
    }
    y_sources <- seq(y_data + 0.2, (source$S_factor_levels * 
                                      source$n.sources * spacing) - spacing + y_data + 
                       0.2, by = spacing)
    MU_plot <- source$S_MU[, mix$iso_names] + discr_mu_plot
    SIG_plot <- sqrt(source$S_SIG[, mix$iso_names]^2 + discr_sig2_plot)
  }
  else {
    source_linetype <- 1:source$n.sources
    source_color <- factor(rep("black", source$n.sources))
    index <- 1:source$n.sources
    discr_mu_plot <- discr$mu
    discr_sig2_plot <- discr$sig2
    y_sources <- seq(y_data + 0.2, (source$n.sources * spacing) - 
                       spacing + y_data + 0.2, by = spacing)
    source$S_factor_levels <- 0.5
    MU_plot <- source$S_MU + discr_mu_plot
    SIG_plot <- sqrt(source$S_SIG^2 + discr_sig2_plot)
  }
  MU_plot <- as.vector(MU_plot)
  SIG_plot <- as.vector(SIG_plot)
  df_sources <- data.frame(x = MU_plot, xmin = MU_plot - SIG_plot, 
                           xmax = MU_plot + SIG_plot, y = y_sources, linetype = source_linetype, 
                           scolour = source_color)
  source.labels <- data.frame(x = MU_plot[index], y = y_sources[index] + 
                                spacing * source$S_factor_levels, label = source$source_names)
  .e <- environment()
  dev.new()
  if (mix$n.effects == 2) {
    shapes <- c(16, 17, 15, 3, 7, 8, 1, 6, 35, 36, 37, 4, 
                18, 14, 11, 9, 13)
    shapes <- shapes[1:mix$FAC[[2]]$levels]
    if (!is.na(source$by_factor)) {
      g <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, 
                                                   y = y), environment = .e) + ggplot2::geom_jitter(ggplot2::aes(colour = factor(mix$FAC[[1]]$values), 
                                                                                                                shape = factor(mix$FAC[[2]]$values)),width = 0,  height = 0.1, show.legend = T) + ggplot2::scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)), 
                                                                                                                                                                                                                                                  labels = mix$FAC[[1]]$labels) + ggplot2::scale_shape_manual(values = shapes, 
                                                                                                                                                                                                                                                                                                              labels = mix$FAC[[2]]$labels) + ggplot2::geom_point(data = df_sources, 
                                                                                                                                                                                                                                                                                                                                                                  ggplot2::aes(x = x, y = y, colour = scolour), 
                                                                                                                                                                                                                                                                                                                                                                  size = 2, show.legend = F) + ggplot2::geom_errorbarh(data = df_sources, 
                                                                                                                                                                                                                                                                                                                                                                                                                       ggplot2::aes(xmin = xmin, xmax = xmax, colour = scolour), 
                                                                                                                                                                                                                                                                                                                                                                                                                       size = 1, height = 0, linetype = source_linetype, 
                                                                                                                                                                                                                                                                                                                                                                                                                       show.legend = F) + ggplot2::geom_text(data = source.labels, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                             ggplot2::aes(x = x, y = y, label = label), show.legend = F) + 
        ggplot2::scale_y_continuous(breaks = NULL) + 
        ggplot2::ylab("") + ggplot2::xlab(x_label) + 
        ggplot2::theme_bw() + ggplot2::theme(legend.position = c(0, 
                                                                 1), legend.justification = c(0, 1), legend.title = ggplot2::element_blank())
      print(g)
    }
    else {
      g <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, 
                                                   y = y), environment = .e) + ggplot2::geom_jitter(ggplot2::aes(colour = factor(mix$FAC[[1]]$values), 
                                                                                                                shape = factor(mix$FAC[[2]]$values)), width = 0, height = 0.1, show.legend = T) + ggplot2::scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)), 
                                                                                                                                                                                                                                                  labels = mix$FAC[[1]]$labels) + ggplot2::scale_shape_manual(values = shapes, 
                                                                                                                                                                                                                                                                                                              labels = mix$FAC[[2]]$labels) + ggplot2::geom_point(data = df_sources, 
                                                                                                                                                                                                                                                                                                                                                                  ggplot2::aes(x = x, y = y), size = 2, show.legend = F) + 
        ggplot2::geom_errorbarh(data = df_sources, ggplot2::aes(xmin = xmin, 
                                                                xmax = xmax), size = 1, height = 0, linetype = source_linetype, 
                                show.legend = F) + ggplot2::geom_text(data = source.labels, 
                                                                      ggplot2::aes(x = x, y = y, label = label), show.legend = F) + 
        ggplot2::scale_y_continuous(breaks = NULL) + 
        ggplot2::ylab("") + ggplot2::xlab(x_label) + 
        ggplot2::theme_bw() + ggplot2::theme(legend.position = c(0, 
                                                                 1), legend.justification = c(0, 1), legend.title = ggplot2::element_blank())
      print(g)
    }
  }
  if (mix$n.effects == 1) {
    if (!is.na(source$by_factor)) {
      g <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, 
                                                   y = y), environment = .e) + ggplot2::geom_jitter(ggplot2::aes(colour = factor(mix$FAC[[1]]$values)),width = 0, height = 0.1, 
                                                                                                   show.legend = T) + ggplot2::scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)), 
                                                                                                                                                     labels = mix$FAC[[1]]$labels) + ggplot2::geom_point(data = df_sources, 
                                                                                                                                                                                                         ggplot2::aes(x = x, y = y, colour = scolour), 
                                                                                                                                                                                                         size = 2, show.legend = F) + ggplot2::geom_errorbarh(data = df_sources, 
                                                                                                                                                                                                                                                              ggplot2::aes(xmin = xmin, xmax = xmax, colour = scolour), 
                                                                                                                                                                                                                                                              size = 1, height = 0, linetype = source_linetype, 
                                                                                                                                                                                                                                                              show.legend = F) + ggplot2::geom_text(data = source.labels, 
                                                                                                                                                                                                                                                                                                    ggplot2::aes(x = x, y = y, label = label), show.legend = F) + 
        ggplot2::scale_y_continuous(breaks = NULL) + 
        ggplot2::ylab("") + ggplot2::xlab(x_label) + 
        ggplot2::theme_bw() + ggplot2::theme(legend.position = c(0, 
                                                                 1), legend.justification = c(0, 1), legend.title = ggplot2::element_blank())
      print(g)
    }
    else {
      g <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, 
                                                   y = y), environment = .e) + ggplot2::geom_jitter(ggplot2::aes(colour = factor(mix$FAC[[1]]$values)), width = 0, height = 0.1, 
                                                                                                   show.legend = T) + ggplot2::scale_colour_discrete(breaks = levels(factor(mix$FAC[[1]]$values)), 
                                                                                                                                                     labels = mix$FAC[[1]]$labels) + ggplot2::geom_point(data = df_sources, 
                                                                                                                                                                                                         ggplot2::aes(x = x, y = y), size = 2, show.legend = F) + 
        ggplot2::geom_errorbarh(data = df_sources, ggplot2::aes(xmin = xmin, 
                                                                xmax = xmax), size = 1, height = 0, linetype = source_linetype, 
                                show.legend = F) + ggplot2::geom_text(data = source.labels, 
                                                                      ggplot2::aes(x = x, y = y, label = label), show.legend = F) + 
        ggplot2::scale_y_continuous(breaks = NULL) + 
        ggplot2::ylab("") + ggplot2::xlab(x_label) + 
        ggplot2::theme_bw() + ggplot2::theme(legend.position = c(0, 
                                                                 1), legend.justification = c(0, 1), legend.title = ggplot2::element_blank())
      print(g)
    }
  }
  if (mix$n.effects == 0) {
    g <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, 
                                                 y = y)) + ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.2, 
                                                                                                                   height = 0.1)) + ggplot2::geom_point(data = df_sources, 
                                                                                                                                                        ggplot2::aes(x = x, y = y), size = 2, show.legend = F) + 
      ggplot2::geom_errorbarh(data = df_sources, ggplot2::aes(xmin = xmin, 
                                                              xmax = xmax), size = 1, height = 0, linetype = source_linetype, 
                              show.legend = F) + ggplot2::geom_text(data = source.labels, 
                                                                    ggplot2::aes(x = x, y = y, label = label), show.legend = F) + 
      ggplot2::scale_y_continuous(breaks = NULL) + ggplot2::ylab("") + 
      ggplot2::xlab(x_label) + ggplot2::theme_bw() + ggplot2::theme(legend.position = c(0, 
                                                                                        1), legend.justification = c(0, 1), legend.title = ggplot2::element_blank())
    print(g)
  }
  if (plot_save_pdf == TRUE) {
    mypath <- file.path(paste(getwd(), "/", filename, ".pdf", 
                              sep = ""))
    cairo_pdf(filename = mypath, width = 7, height = 7)
    print(g)
    dev.off()
  }
  if (plot_save_png == TRUE) {
    mypath <- file.path(paste(getwd(), "/", filename, ".png", 
                              sep = ""))
    png(filename = mypath)
    print(g)
    dev.off()
  }
}
