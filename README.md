# Flavonols-peak-areas
### MAKE SURE TO MAKE THE WORKING DIRECTORY THE DIRECTORY WHERE ALL THE MZXML 
### OR MZML FILES ARE STORED. As a simplest practice, store them in the same file
### as this code. 

## for xcms see
## http://bioconductor.org/packages/devel/bioc/html/xcms.html
library(xcms)

# Note: use / not \ in your file path. If copied from Windows, you can use the find and replace method to quickly correct this issue.
setwd("F:/DANA/Cody/mzxml_041816")

# This is the fuction that needs to be called to run everything. 
# All parameters are lists.
# All lists need to be same length except files.
allPeakAreas <- function(files, mz, mz.tol, rt.min, rt.max, compound_names)
{
  # assigns ranges all of the mz value ranges and rt ranges that need to be evaluated. 
  mzmin <- mz - mz.tol
  mzmax <- mz + mz.tol
  
  peakranges <- data.frame(mzmin = mzmin, mzmax = mzmax, rtmin = rt.min*60, rtmax = rt.max*60)
  pr.matrix <- as.matrix(peakranges)
  
  # Loops through all files to extract peak areas from each file for each mz/rt
  all_files_data <- lapply(files, function(file) {

    # Creates and xcmsRaw object from the file
    xraw <- xcmsRaw(filename = file)
    
    # Finds peaks and returns data about each
    # peaks <- getPeaks(xraw, ranges)
    
    # Changed step size
    peaks <- getPeaks(xraw, pr.matrix, step = .05)
    
    # Attaches the compound names onto the peak data
    peaks <- cbind(compound_names, peaks)
    
    # Trims out unneeded data and adds file names
    peaks <- apply(peaks, 1, function(peak){
      
      c(file, peak[1], peak[2], peak[5], peak[8])  
        
      })
    
    # returns each set of peaks to the outter loop
    print(paste(file, " complete."))
    peaks
  })
  
  # Transposes returned data for readability and converts to a dataframe
  all_files_data <- data.frame(all_files_data)
  all_files_data <- t(all_files_data)
  
  # Changes the row names for readability
  return_df <- all_files_data[,-1]
  rownames(return_df) <- all_files_data[,1]
  
  # Passes back a dataframe with files as row names and columns:
  # Compound_names | mz | rt | into    Note: into = peak area
  return_df
}

### Input data goes down here. Files can be either mzML or mzXML files. 

# Make sure each filename is unique and that they are all in the same directory! 
# Here is where you can swap out what files are being checked.
# filenames.txt should be a list of the names of the mzML or mzXML files you wish to analyze. The should each be on their own line in a .txt file.
files <- scan("filenames.txt", what = "", sep="\n")

#### INPUT OPTION 1 #####

# These are the mz values and retention times that you will look at.
# Make sure mz, mz.tol, rt.min, rt.max, and compound_names are the same length.
#mz <- c(90.0555,175.1195,134.0453,241.0311,148.061,76.03985,156.0773,132.1025,147.1134,150.0589,166.0868,116.0712,106.0504,120.661,182.0817,118.0868)
#mz.tol <- c(.001,.001,.001,.001,.01,.001,.001,.001,.001,.001,.001,.01,.001,.01,.001,.001)

# Retention Times (for future use these need to be rt.min and rt.max values)
#rt.min <- c(11,19,15.4,20.2,13.68,12.42,19.24,8.1,19.2,9.52,8.74,9.8,14.25,12.5,11.55,9.2)
#rt.max <- c(11.34,19.6,16.54,20.64,14.12,12.74,19.64,8.54,19.66,9.8,9.08,10.2,14.65,12.8,11.85,9.45)

# List of strings (every individual name needs to be in quotes)
#compound_names <- c("A","R","D","Cystine","E","G","H","I/L","K","M","F","P","S","T","Y","V")



#### INPUT OPTION 2 #####

# Can also run this command to get data from csv template. CSV columns need to have same headers as template.
filename <- "050916_Codyfinal.csv"

input_data <- read.csv(filename, header = TRUE)

compound_names <- as.character(input_data$name) 
mz <- input_data$mass
mz.tol <- input_data$mz.tol
rt.min <- input_data$rt.min
rt.max <- input_data$rt.max

#########################

# Command to run function once all variables are set and functions have been entered into the environment. 
data <- allPeakAreas(files, mz, mz.tol, rt.min, rt.max, compound_names)

# If you want to generate a csv
write.csv(data, "outputcodyfinal.csv")
