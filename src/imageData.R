library(R6)
library(ggplot2)

imageData <- R6Class("imageData",
  public = list(
    initialize = function(dataMatrix) {
      if(!is.matrix(dataMatrix)){
        stop("dataMatrix must be a matrix.")
      }
      private$.dataMatrix <- dataMatrix
      private$.dim <- dim(dataMatrix)
      private$.valueSet <- unique(as.vector(dataMatrix))
      private$.ncolors <- length(private$.valueSet)
    },
    
    plot = function(){
      df <- data.frame(x = as.vector(row(private$.dataMatrix)),
                       y = as.vector(col(private$.dataMatrix)),
                       value = as.factor(as.vector(private$.dataMatrix)))
      
      ggplot(df, aes(x = x, y = y, fill = value)) +
        geom_tile() +
        scale_fill_brewer(palette = "Dark2") +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.title = element_blank(),
              legend.text = element_text(margin = margin(l = 5)))
    }
  ),

  active = list(
    matrix = function() {
      private$.dataMatrix
    },
    dim = function() {
      private$.dim
    },
    valueSet = function() {
      private$.valueSet
    },
    ncolors = function() {
      private$.ncolors
    }
    
  ),

  private = list(
    .dataMatrix = NULL,
    .dim = NULL,
    .valueSet = NULL,
    .ncolors = NULL
  )
)

#n <- 100
#m <- 120
#a <- matrix(rbinom(n*m, size = 1, prob = 0.3), nrow = n, ncol = m)

#x <- imageData$new(a)
#x$plot()
