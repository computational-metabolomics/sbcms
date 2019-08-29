#' @import ggplot2
#' @import reshape2
#' @importFrom grDevices pdf dev.off
NULL

##' Plot the output from signal batch correction for the first 100 features
##'
##' @param df A data frame containing the original data before correction
##'  (samples in columns, features in rows).
##' @param corrected_df A data frame containing the corrected data,
##'  as output by QCRSC
##' @param batch A vector indicating the batch each sample was measured in.
##'  If only one batch was measured then all values should be set to 1
##' @param classes A factor or character vector of sample classes.
##'  All QC samples should be labelled "QC".
##' @param order A numeric vector  indicating the order in which samples
##'  were measured.
##' @param output Filename of the output pdf file. Can include the path.
##' 
##' @return Pdf file showing data before and after signal correction
##' 
##' @examples 
##' classes <- sbcdata$class
##' batch <- sbcdata$batch
##' order <- c(1:nrow(sbcdata$data))
##' data <- t(sbcdata$data[, 1:20])
##' out <- QCRSC(df = data, order = order, batch = batch, classes = classes,
##'  spar = 0, minQC = 4)
##' sbcmsPlot (df = data, corrected_df = out, classes, batch, order,
##'  output="sbcms_plots.pdf") 
##'  
##' @export

sbcmsPlot <- function (df, corrected_df, classes, batch, order,
  output="sbcms_plots.pdf")
{
  shapes <- rep(1,length(classes))
  shapes[classes == "QC"] <- 3

  manual_color <- c("#386cb0", "#ef3b2c", "#7fc97f", "#fdb462", "#984ea3",
    "#a6cee3", "#778899", "#fb9a99", "#ffff33")

  gg_THEME <- theme(
    panel.background=element_blank(),
    panel.grid.major=element_line(color="gray80", size=0.3),
    axis.line=element_line(color="black"),
    axis.text=element_text(color="black"),
    axis.title=element_text(color="black"),
    panel.grid.minor.x=element_line(color="gray80", size=0.3,
      linetype="dashed"),
    panel.grid.minor.y=element_line(color="gray80", size=0.3)

  )

  plots <- list()

  plotLim <- 100
  if (nrow(df) < 100) plotLim <- nrow(df)

  for (peakn in seq_len(plotLim)){
    cat (peakn, "\n")
    A <- data.frame (x=c(seq_len(ncol(df))), original=log(df[peakn, ],10),
                     fitted=log(corrected_df[peakn, ],10),
                     batch=as.factor(batch), shapes=shapes)
    A <- melt(A, id.vars=c("x","batch","shapes"))

    plots[[peakn]] <- ggplot (A, aes_(~x, ~value, col=~batch, shape=~shapes))+
      facet_grid(variable ~ .)+
      geom_point() +scale_shape_identity()+
      scale_colour_manual(values=manual_color)+
      ggtitle(row.names(df)[peakn])+
      ylab("log10(intensity)")+
      xlab("injection order")+
      gg_THEME
  }

  pdf (output)
  invisible(lapply(plots, print))
  dev.off()

}
