
# Create dataset for RNN
# the size of the data is 30, other signatures will be padded with zeros.
get_rnn_dataset <- function()
{
  data_dir = "rnn_input/"
  suppressWarnings(dir.create(data_dir))
  max_signatures <- 30
  rescaled = F
  
  for (dir_ in list.files(DIR_RESULTS, full.names = T))
  {
    for (example in list.files(dir_)) {
      example_path <- paste0(dir_, "/", example, "/")
      
      if (!file.exists(example_path))
        next
      
      file_summary <- paste0(example_path, "mixtures")
      if (rescaled) {
        file_summary <- paste0(file_summary, ".rescaled")
      }
      file_summary <- paste0(file_summary, ".csv")
      
      if (!file.exists(file_summary))
      {
        print(paste0("No overall summary for ", dir_))
        next
      }
      
      mixtures <- read.csv(file_summary)
      mixtures.data <- mixtures[,-1]
      
      if (ncol(mixtures) < 5)
        next
      
      signatures <- sapply(mixtures[,1], toString)
      signatures <- as.numeric(gsub("S([\\d]*)", "\\1", signatures))
      missing_sigs <-  setdiff(1:max_signatures, signatures)
      initial_data <- cbind(signatures, mixtures.data)
      padded_data <- cbind(missing_sigs, data.matrix(matrix(0, ncol=ncol(mixtures.data), nrow=length(missing_sigs))))
      colnames(padded_data) <- colnames(initial_data)
      padded_data <- rbind(initial_data, padded_data)
      padded_data <- padded_data[order(padded_data[,1]),]
      padded_data <- padded_data[,-1]
      padded_data <-  round( padded_data, digits = 3)
      options(scipen=999)
      write.table(padded_data, file = paste0(data_dir, example, ".txt"), 
                  sep="\t", quote = F, row.names = F, col.names = F)
    }
  }
}