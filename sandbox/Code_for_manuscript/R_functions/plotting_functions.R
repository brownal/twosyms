# Function to set transparency of a color
transparent <- function(col, alpha) {
  tmp <- col2rgb(col)
  rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255, alpha)
}

# A function for plotting histograms colored by outcomes
# data should be a list of the form (data1=data to plot 1st, data2=data to plot 2nd, ...)
stackedHist <- function(data, colors=NULL, legend=NULL, legendLocation=NULL, plotLegend=FALSE, breaks="Sturges", ...) {
  plotargs <- list(...)

  # Assign colors as shades of gray if no colors are given
  if (is.null(colors)) {
    colors <- lapply(seq(0, 0.9, length.out=length(data)), function(x) {rgb(x, x, x)})
  }

  # Create but don't plot a histogram of all the data
  # This will be used get the breaks for the histograms we "stack" and to normalize their heights
  histAll <- hist(unlist(data), breaks=breaks, plot=FALSE)
  breaks <- histAll$breaks # get the breaks to use for the other histograms
  normalizeBy <- sum(histAll$counts) # Get the total number of simulations to normalize by

  # Plot the data as a stacked histogram.
  # We go through the data backwards, first plotting a histogram of everything, then
  # removing one dataset at a time and adding new histogram on top the other histograms
  # This creates a "stacked" effect, with the first dataset at the bottom of the stack
  histAll$density <- histAll$counts / normalizeBy * 100
  do.call(plot, c(list(histAll, freq=FALSE, col=colors[[length(data)]]), plotargs))

  for (i in (length(data)-1):1) {
    histi <- hist(unlist(data[1:i]), breaks=breaks, plot=FALSE)
    histi$density <- histi$counts / normalizeBy * 100
    plot(histi, freq=FALSE, add=TRUE, col=colors[[i]])
  }

  # Add a legend
  if (plotLegend) {
    if (is.null(legend)) {
      legend <- lapply(seq(1, length(data)), function(x) {paste0("data ", x)})
    }

    if (is.null(legendLocation)) {
      legendLocation <- "topright"
    }

    legend(legendLocation, legend, pt.bg=colors, col="black", pch=22)
  }
}
