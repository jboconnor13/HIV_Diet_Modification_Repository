###This function provides an simplified field in which rarer items are grouped into "Other

## Inputs:

## data=data set
## name_field=the field that you want to be categorized into more abundant names or "Other"
## type=how you want to define the abundance 
##    "Average"= All names with average abundance below threshold categorized as "Other" 
##    "Median"= All names with median abundance below threshold categorized as "Other"
##    "Max"= All names with max abundance below threshold categorized as "Other"
##    "Max per Sample"= All names with max abundance in all samples below threshold categorized as "Other"
##    "Max Relative per Sample"= All names with max relative abundance in all samples below threshold categorized as "Other"
##    "Top Average"= All names not in the top n average abundance values (specified by threshold) categorized as "Other" 
##    "Top Median"= All names not in the top n median abundance values (specified by threshold) categorized as "Other"
##    "Top Max"= All names not in the top n max abundance values (specified by threshold) categorized as "Other"
##    "Frequency"=  All names detected with a frequency in samples set less than specified threshold are categorized as "Other"
##    "Percentage"= All names detected with a percentage in samples set less than specified threshold are categorized as "Other"
##threshold=The threshold for cutoffs of what is determined to be "Other" if it is a top type it is the threshold is a list, otherwise it is a cutoff in which any amount under is categorized as Other  
##abundance_field=the quantifiable field (default: "Abundance")
##id_var=the id variable for samples, needed for frequency and percentage fields (default: "Sample")

## Outputs:

## A vector of the values if the name_field with names either the same (if included by threshold) or "other" (if cut off by the threshold)

other_generator <- function(data, name_field, type, threshold, abundance_field="Abundance", id_var="Sample") {
library(dplyr)


#First we have to pull all the column numbers 
  
#We pull the column number for the name_field
nm_col <- which(colnames(data)==name_field)
#We pull the column number for the abundance_field
ab_col <- which(colnames(data)==abundance_field)
#We pull the column number for the id_field
id_col <- which(colnames(data)==id_var)

#Now we genearte a temporary data set 
data_temp <- data[,c(id_col,nm_col,ab_col)]

#For simplicity we rename the columns
colnames(data_temp) <- c("id","nm","ab")

###Continuous threshold types include names that reflect values above a continuous numerical threshold (i.e. only names with an average over 0.5) categorizing those below as "Other"

#If it is a "Average" type, we pull names with average values above the specified threshold 
if (type=="Average") {
  max_data_temp <- data_temp %>%
    group_by(nm) %>%
    summarise(avg = mean(as.numeric(ab))) %>%
    arrange(desc(avg))
  names_inc <- max_data_temp$nm[which(as.numeric(max_data_temp$avg)>threshold)]
}

#If it is a "Median" type, we pull names with median values above the specified threshold 
if (type=="Median") {
  max_data_temp <- data_temp %>%
    group_by(nm) %>%
    summarise(med = median(as.numeric(ab))) %>%
    arrange(desc(med))
  names_inc <- max_data_temp$nm[which(as.numeric(max_data_temp$med)>threshold)]
}

#If it is a "Max" type, we pull names with max values above the specified threshold 
if (type=="Max") {
  max_data_temp <- data_temp %>%
    group_by(nm) %>%
    summarise(max = max(as.numeric(ab))) %>%
    arrange(desc(max))
  names_inc <- max_data_temp$nm[which(as.numeric(max_data_temp$max)>threshold)]
}

#If it is a "Max per Sample" type, we pull the max of the abundance field grouped by sample (specified by the id_var)
  if (type=="Max per Sample") {
    max_data_temp <- data_temp %>%
      group_by(nm,id) %>%
      summarise(max = max(as.numeric(ab)))
    #We then pick names that have a max above the threshold
    names_inc <- max_data_temp$nm[which(as.numeric(max_data_temp$max)>threshold)]
  }
#If it is a "Max Relative per Sample" type, we pull the max of the relative abundance (abundance divided by total other abundance) field grouped by sample (specified by the id_var)
  if (type=="Max Relative per Sample") {
    max_data_temp <- data_temp
    #We create a "perc" field to include relative abundances "percent abundance"
    max_data_temp["perc"] <- NA
    #Then we go through each row 
    for (l in 1:nrow(max_data_temp)) {
      #We select the sample id (specified by id_var in column 1)
      stud_l <- max_data_temp[l,1]
      #We select the sum for that sample id
      sum_l <- sum(as.numeric(max_data_temp[which(as.vector(max_data_temp$id)==stud_l),3]))
      #We divide the abundance by that sum to get the relative abundance 
      max_data_temp$perc <- as.numeric(max_data_temp$ab[l])/as.numeric(sum_l)
    }
    max_data_temp_perc <- max_data_temp %>%
      group_by(nm,id) %>%
      summarise(max = max(as.numeric(perc)))
    #Now as done above with the max abundance type we group the max perc values and select names included above the threshold 
    names_inc <- max_data_temp_perc$nm[which(as.numeric(max_data_temp_perc$max)>threshold)]
  }



###Top threshold types include names that reflect the top n values specified by an integer threshold (i.e. only names with the top 10 averages) categorizing those not in the top n to be categorized as "Other"

#If it is a "Top Average" type, we pull the average amounts and select those above the top n names (specified by threshold) 
  if (type=="Top Average") {
    max_data_temp <- data_temp %>%
      group_by(nm) %>%
      summarise(avg = mean(as.numeric(ab))) %>%
      arrange(desc(avg))
    names_inc <- max_data_temp$nm[1:threshold]
  }
#If it is a "Top Median" type, we pull the median amounts and select those above the top n names (specified by threshold)
  if (type=="Top Median") {
    max_data_temp <- data_temp %>%
      group_by(nm) %>%
      summarise(med = median(as.numeric(ab))) %>%
      arrange(desc(med))
    names_inc <- max_data_temp$nm[1:threshold]
  }
#If it is a "Top Max" type, we pull the max amounts and select those above the top n names (specified by threshold)
  if (type=="Top Max") {
    max_data_temp <- data_temp %>%
      group_by(nm) %>%
      summarise(max = max(as.numeric(ab))) %>%
      arrange(desc(max))
    names_inc <- max_data_temp$nm[1:threshold]
  }
   
####Prevalence threshold types include names that reflect prevalence in the sample set (i.e. select names that appear in 50% of the samples)

#If it is "Frequency" type, we pull the frequency in which the name is detected in the sample and pull names above a sepcific frequency threshold    
  if (type=="Frequency") {
    data_temp <- data_temp[which(data_temp[,3]>0),]
    freq_data_temp <-as.data.frame(table(data_temp[which(colnames(data_temp)=="nm")]))
    names_inc <- freq_data_temp[which(as.numeric(freq_data_temp$Freq)>threshold),1]
  }
#If it is "Percentage" type, we pull the percentage in which the name is detected in the sample and pull names above a specific percentage threshold
  if (type=="Percentage") {
    data_temp <- data_temp[which(data_temp[,3]>0),]
    id_num <- length(unique(data[,which(colnames(data)==id_var)]))
    freq_data_temp <-as.data.frame(table(data_temp[which(colnames(data_temp)=="nm")]))
    names_inc <- freq_data_temp[which(as.numeric(freq_data_temp$Freq/id_num)>threshold),1]
  }
    
#Now we pull all the original names
  old_varnames <- as.vector(data[,nm_col])
#Create a new duplicate
  new_varnames <- old_varnames
#Now we go through the old_var_names and if they are in the names_inc, we include them otherwise the Other values remain
  for (r in 1:nrow(data)) {
    if (!(old_varnames[r] %in% names_inc)) {
      new_varnames[r] <- "Other"
    }
  }
  return(new_varnames)
} 
