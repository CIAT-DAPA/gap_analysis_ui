###############################################################################################
####Creating arguments to run MaxEnt in the SDM Packages (Modified of Jorge Velasquez Script###
###############################################################################################

CreateMXArgs <- function(calibration, use.maxnet = TRUE){
  mxnt.args <- c("linear")
if(use.maxnet){
  
  if(!is.null(calibration)){
    best.ind <- which.min(calibration$deltaAICc)
    args <- calibration[best.ind, ]
    features <- args$classes
    betamultiplier <- args$regMult
  
  if(grepl("q", features)){
    mxnt.args <- c(mxnt.args, "quadratic")
  } else {
    mxnt.args <- c(mxnt.args, "")
  }
  if(grepl("h", features)){
    mxnt.args <- c(mxnt.args, "hinge")
  } else {
    mxnt.args <- c(mxnt.args, "")
  }
  if(grepl("p", features)){
    mxnt.args <- c(mxnt.args, "product")
  } else {
    mxnt.args <- c(mxnt.args, "")
  }
  if(grepl("t", features)){
    mxnt.args <- c(mxnt.args, "threshold")
  } else {
    mxnt.args <- c(mxnt.args, "")
  }
  mxnt.args <- c(mxnt.args, paste0("betamultiplier=", betamultiplier))
  }else{
   mxnt.args <- c(mxnt.args, "quadratic", "hinge", "product", "", "betamultiplier=1.0")
 }
  
  mxnt.args <- mxnt.args[which(mxnt.args != "")]
  
  
}else{
  if(!is.null(calibration)){
    
    best.ind <- which.min(calibration$deltaAICc)
    features <-as.data.frame(cbind(calibration$linear, calibration$quadratic, calibration$product, calibration$hinge, calibration$threshold))
    names(features) <- c("l", "q", "p", "h", "t")
    features <- features[best.ind,]
    features <- unlist(lapply(1:ncol(features),function(i){
      if(features[,i] == TRUE){
        x <- colnames(features[i])
      } else {
        x <- NULL
      }
      return(x)
    }))
    features <- paste(features, collapse = "")
    betamultiplier <- calibration$regMult[best.ind]
    
    if(grepl("q", features)){
      mxnt.args <- c(mxnt.args, "quadratic")
    } else {
      mxnt.args <- c(mxnt.args, "")
    }
    if(grepl("h", features)){
      mxnt.args <- c(mxnt.args, "hinge")
    } else {
      mxnt.args <- c(mxnt.args, "")
    }
    if(grepl("p", features)){
      mxnt.args <- c(mxnt.args, "product")
    } else {
      mxnt.args <- c(mxnt.args, "")
    }
    if(grepl("t", features)){
      mxnt.args <- c(mxnt.args, "threshold")
    } else {
      mxnt.args <- c(mxnt.args, "")
    }
    mxnt.args <- c(mxnt.args, paste0("betamultiplier=", betamultiplier))
  } else {
     mxnt.args <- c(mxnt.args, "quadratic", "hinge", "product", "", "betamultiplier=1.0")
   }
  
  mxnt.args <- mxnt.args[which(mxnt.args != "")]
 
}  
  return(mxnt.args)
  
}

###############################################################################################
####Calibration function using wright et al., 2014 approach###
###############################################################################################

Calibration_function <- function(spData, bg_data, occName, out_path, ommit, use.maxnet = TRUE){
  cat("Initializing calibration step \n")
  #suppressMessages(if(!require(pacman)){install.packages("pacman");library(pacman)}else{library(pacman)})
  #pacman::p_load(devtools, maxnet)
  # if(!require(enmSdm)){
  #   devtools::install_github('adamlilith/omnibus')
  #   devtools::install_github('adamlilith/statisfactory')
  #   devtools::install_github('adamlilith/enmSdm')
  #   library(omnibus)
  #   library(enmSdm)
  # } else {
  #   library(omnibus)
  #   library(statisfactory)
  #   library(enmSdm)
  #   require(maxnet)
  # }
  
  if(ommit == F){
   
    # Calibration using MaxEnt instead of Maxnet R package.
    
      if(use.maxnet){ 
        #use maxnet
        if(!file.exists(out_path)){
          
          
          if("status" %in% names(spData)){
            spData$status <- NULL
          }
          if("predicted" %in% names(spData)){
            spData$predicted <- NULL
          }
          if("database_id" %in% names(spData)){
            spData$database_id <- NULL
          }
          if("source_db" %in% names(spData)){
            spData$source_db <- NULL
          }
          
          
          
          spData <- spData[complete.cases(spData),]
          
          
          spData <- dplyr::bind_rows(spData %>% 
                                       dplyr::mutate(Y = 1), 
                                     bg_data %>% 
                                       dplyr::mutate(Y = 0)) %>% 
            dplyr::rename(!!occName := Y) %>% 
            dplyr::mutate(across(everything(.), as.numeric))
          
          rm(bg_data)
          
          
          
        cat("Calculating best parameters for maxNet \n")
        cat("This process will take several minutes, please be patient. \n")
        
        data_train <- spData %>% 
          dplyr::select(-any_of(c("drymonths_2_5_min", "ethnicity", "monthCountByTemp10"))) %>% 
          as.matrix(.)
        
        #adding all presence points to background
        pres_to_bg <- data_train[which(data_train[, 1] == 1), ]
        pres_to_bg[,1] <- rep(0, length(pres_to_bg [, 1]))
        
        p <- c(data_train[, 1], pres_to_bg[,1])#adding all presence points to background
        data <- rbind(data_train[, -c(1,2,3)], pres_to_bg[, -c(1,2,3)])#adding all presence points to background
        
        data <- data.frame(occName = p, data)
        
        
          calibration <- trainMaxNet(data = data, regMult = seq(0.5, 6, 0.5), out = 'tuning', verbose = FALSE)
          
          write.csv(calibration, out_path, quote = F, row.names = F)
        }else{
          cat("Calibration File already created, importing it \n")
          
          calibration <- read.csv(out_path)
        }
        
       
        
      }
    
      feat <- CreateMXArgs(calibration, use.maxnet)
      beta <- feat[(grepl("betamultiplier=", feat))]
      beta <- as.numeric(gsub("betamultiplier=", "", beta))
      feat <- feat[(!grepl("betamultiplier=", feat))]
   
    
  } else {
    cat("Ommiting calibration step\n")
    calibration <- NULL
    cat("Setting default parameters","\n")
    beta <- 1
    feat <- c("linear", "quadratic", "hinge", "product")
  
  }
  
  
  
  return(list(features = feat, beta = beta))
  cat("Process done... \n")
}
