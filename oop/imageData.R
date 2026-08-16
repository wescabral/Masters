library(R6)
library(ggplot2)

imageData <- R6Class("imageData",
  public = list(
    initialize = function(dataMatrix, continuous = NULL) {
      if(!is.matrix(dataMatrix)){
        stop("dataMatrix must be a matrix.")
      }
      private$.dataMatrix <- dataMatrix
      private$.dim <- dim(dataMatrix)

      if(is.null(continuous)){
        if(is.numeric(dataMatrix)){
          vals <- as.vector(dataMatrix)
          epsilon <- .Machine$double.eps^0.5
          if(all(abs(vals - round(vals)) < epsilon, na.rm = TRUE)){
            private$.continuous <- FALSE
          } else {
            private$.continuous <- TRUE
          }
        } else {
          private$.continuous <- FALSE
        }
      } else {
        private$.continuous <- isTRUE(continuous)
      }

      if(private$.continuous){
        private$.valueRange <- range(as.vector(dataMatrix), na.rm = TRUE)
        private$.valueSet <- NULL
        private$.ncolors <- NA
      } else {
        private$.valueSet <- unique(as.vector(dataMatrix))
        private$.ncolors <- length(private$.valueSet)
        private$.valueRange <- range(as.vector(dataMatrix), na.rm = TRUE)
      }
    },

    setMatrix = function(newdataMatrix, continuous = NULL) {
      if(is.vector(newdataMatrix) && length(newdataMatrix) == prod(private$.dim)){
        newdataMatrix <- matrix(newdataMatrix, nrow = private$.dim[1], ncol = private$.dim[2])
      }
      if(!is.matrix(newdataMatrix)){
        stop("dataMatrix must be a matrix or a vector with compatible length as the prior matrix.")
      }
      private$.dataMatrix <- newdataMatrix
      private$.dim <- dim(newdataMatrix)

      if(is.null(continuous)){
        if(is.numeric(private$.dataMatrix)){
          vals <- as.vector(private$.dataMatrix)
          epsilon <- .Machine$double.eps^0.5
          if(all(abs(vals - round(vals)) < epsilon, na.rm = TRUE)){
            private$.continuous <- FALSE
          } else {
            private$.continuous <- TRUE
          }
        } else {
          private$.continuous <- FALSE
        }
      } else {
        private$.continuous <- isTRUE(continuous)
      }

      if(private$.continuous){
        private$.valueRange <- range(as.vector(private$.dataMatrix), na.rm = TRUE)
        private$.valueSet <- NULL
        private$.ncolors <- NA
      } else {
        private$.valueSet <- unique(as.vector(private$.dataMatrix))
        private$.ncolors <- length(private$.valueSet)
        private$.valueRange <- range(as.vector(private$.dataMatrix), na.rm = TRUE)
      }
      invisible(self)
    },
    
    plot = function(){
      df <- data.frame(x = as.vector(row(private$.dataMatrix)),
                       y = as.vector(col(private$.dataMatrix)),
                       value = as.factor(as.vector(private$.dataMatrix)))
      
      base_plot <- ggplot(df, aes(x = x, y = y))

      if(private$.continuous){
        base_plot <- base_plot +
          geom_tile(aes(fill = value)) +
          scale_fill_gradient(low = "white", high = "steelblue") +
          labs(fill = NULL)
      } else {
        base_plot <- base_plot +
          geom_tile(aes(fill = as.factor(value))) +
          scale_fill_brewer(palette = "Dark2") +
          labs(fill = NULL)
      }

      base_plot +
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
    },
    continuous = function() {
      private$.continuous
    },
    valueRange = function() {
      private$.valueRange
    }
  ),

  private = list(
    .dataMatrix = NULL,
    .dim = NULL,
    .valueSet = NULL,
    .ncolors = NULL,
    .continuous = FALSE,
    .valueRange = NULL
  )
)