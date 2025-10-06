# ============================================================
# Correlation-Based Redundancy Filter (Pairwise Pearson | abs >= 0.90)
# ------------------------------------------------------------
# Reads a feature table, marks highly correlated (redundant) columns,
# removes them, writes a filtered file, then reloads it and prints the
# Output column. Functionality is intentionally unchanged.
#
# Notes:
#  • The last column is assumed to be the label (Output); only the first
#    (col - 1) columns are considered for correlation pruning.
#  • Threshold is |r| >= 0.90 (Pearson).
#  • Bare variable names (i, j, col1, col2) are left as-is to print in console.
#  • All absolute Windows paths kept exactly as provided.
# ============================================================

## ============= read from the original feature set ========================= ##
# Original (commented) path retained:
# file <- read.csv("D:\\Research Work\\Disertation Project 2\\Raw Data\\featureExtraction\\dataSplit\\addNewSeq\\validate_train_data.csv", header = TRUE)

file <- read.csv(
  "D:\\Research_Work\\Disertation_Project_2\\Raw Data\\featureExtraction\\reviewer_response/test_result/dataSplit\\validate_train_data.csv",
  header = TRUE
)

# Optional preprocessing ideas from original code are kept commented:
# file <- file %>% select(where(~ n_distinct(.) > 1))
# file$Output[file$Output == -1] <- 0
# file <- as.numeric(file$Output)
# normalize <- function(x) { (x - min(x)) / (max(x) - min(x)) }
# file <- as.data.frame(lapply(file, normalize))

# ---------------------------
# Bookkeeping
# ---------------------------
col <- ncol(file)         # total number of columns (features + Output)
col                        # print for visibility (unchanged behavior)

# flag: 0 = keep (default), 1 = remove
flag <- c(1:col)
flag[1:col] <- 0

# Rlist holds indices of columns to remove due to high correlation
Rlist <- c()

# Exclude the last column from correlation pruning (assumed label)
col <- col - 1
col                         # print (unchanged)

# ============================================================
# Build list of redundant columns via pairwise abs(Pearson r) >= 0.90
# Nested loop is preserved. If flag[i] is already 1, we skip marking j.
# ============================================================
for (i in 1:col) {
  for (j in 1:col) {
    
    if (flag[i] == 0) {
      # Grab the two candidate columns
      col2 <- file[[i]]
      col1 <- file[[j]]
      
      # The following bare variables are intentionally preserved.
      # In R scripts, this causes their values to print to console.
      i
      j
      col2
      col1
      
      # Absolute Pearson correlation (unchanged)
      corre <- abs(cor(col1, col2, method = "pearson"))
      
      # res is built but not used later; kept for parity
      res <- c(corre, i, j)
      
      # If highly correlated (and not the same column), mark j for removal
      if (corre >= 0.90 && i != j) {
        flag[j] <- 1
        Rlist <- c(Rlist, j)
        # (index-tracking code from original left commented)
        # Rlist[index] <- j
        # index <- index + 1
      }
    }
  }
}

## ================== De-duplicate and sort removal list =================== ##
Rlist <- unique(Rlist)
Rlist
Rlist <- sort(Rlist)
Rlist

## ================== Remove marked columns from original file ============= ##
corrFST <- file
corrFST <- corrFST[, -Rlist]   # column-wise drop by indices

length(corrFST)                # prints number of remaining columns
corrFST                        # prints resulting data frame (unchanged)

# Original alternative removal-by-loop kept commented:
# row <- length(Rlist)
# corrFST <- file
# count <- 0
# x <- 0
# for (i in 1:row) {
#   x <- as.numeric(Rlist[i]) - count
#   corrFST <- corrFST[-x]
#   count <- count + 1
# }

# ---------------------------
# Write filtered feature set
# ---------------------------
# (Other historical output paths kept commented)
# write.csv(corrFST, file = "D:\\Bacteriocin\\data\\featureExtraction\\CorrFeatureSet.csv", row.names = FALSE)
# write.csv(corrFST, file = "D:\\Research Work\\Disertation Project 2\\Raw Data\\featureExtraction\\pairsoncorrelationSVCL1.csv", row.names = FALSE)

write.csv(
  corrFST,
  file = "D:\\Research_Work\\Disertation_Project_2\\Raw Data\\featureExtraction\\reviewer_response/test_result\\pairsoncorrelation.csv",
  row.names = FALSE
)

# ---------------------------
# Reload the filtered file and show Output column
# ---------------------------
file_read <- read.csv(
  "D:\\Research_Work\\Disertation_Project_2\\ Raw Data\\featureExtraction\\reviewer_response/test_result/pairsoncorrelation.csv",
  header = TRUE, sep = ","
)
length(file_read)     # print the column count
file_read$Output      # print label column (unchanged)
