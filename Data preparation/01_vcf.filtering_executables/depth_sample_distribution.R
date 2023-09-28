# With this script I will study the distribution of the depth per coverage
# in each lynx pardinus sample, in order to set limits.
#install.packages("tidyverse")
#install.packages("haven")
library(tidyverse)

#Set the working directory we are going to work on
setwd("/Users/lorenalorenzo/Documents/depth_filter/samples/")
# Create an object with the working directory where there are the depth files.
wd<-"/Users/lorenalorenzo/Documents/depth_filter/samples/"

# Create a list with the names of the files
list.files(wd)
sample_files <- list.files(wd, pattern="*masked.depth")

# Create an empty table to fill with samples info in each for loop
depth_per_sample <- data.frame()

################  DEPTH CALCULATIONS ######################

for (i in 1:length(sample_files)){

  #Import the sample file table. Use paste0 to tell R to read next command
  #combined (in this case, wd followed by sample_file)
  input.depth <- read_delim(paste0(wd, sample_files[i]), col_names=F, delim='\t')

  # Create a frequency table of the depth values
  depth_freq_table <- as.data.frame(table(input.depth$X3))

  # Define the functions for mean and standard deviation,
  # and define maximum and minimum depth based on different criteria
  depth_mean <- mean(input.depth$X3)
  depth_sd <- sd(input.depth$X3)
  depth_mode<- which.max(depth_freq_table$Freq)

  maxDepth= 1.8 * depth_mode
  minDepth= 3

  # Calculate percentage of variants left after filtering for min Depth and max Depth
  PercentLeft <- length(which(input.depth$X3 > minDepth & input.depth$X3 < maxDepth)) / nrow(input.depth) *100

  # Save the name of each sample an add all info to the table
  # previously created empty. Then, save this info in our laptop.
  sample <- substr(sample_files[i], start = 1, stop = 12)
  depth_per_sample <- rbind(depth_per_sample, data.frame(sample, depth_mean, depth_sd,
                                depth_mode, minDepth, maxDepth,  PercentLeft))
  write.csv(depth_per_sample, file="samples_depth.csv")

  # Draw and save a graph of the distribution of depth values,
  # with upper and lower depth limits
  p<- ggplot(depth_freq_table, aes(x = as.numeric(Var1), y = Freq)) +
      geom_bar(stat = "identity", color = "black") +
      scale_x_continuous(breaks = seq(0, 250, by=5), limits = c(0, 2*maxDepth)) +
      scale_y_continuous(expand=c(0,0),limits = c(0, 1.2*max(depth_freq_table$Freq))) +
      geom_vline(xintercept=maxDepth,linetype="dashed", size=0.5) +
      geom_vline(xintercept=minDepth,linetype="dashed", size=0.5) +
      theme_classic() +
      theme(text = element_text(size=10))
  p + ggtitle(sample) + xlab("depth distribution")
  ggsave (filename = (paste0(sample,"_depth_distribution.pdf")), path = wd)

}
