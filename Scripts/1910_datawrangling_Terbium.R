# This script is to process the Python data
require(dplyr)
require(tidyr)
require(ggplot2)

# Load data
RTEoutputraw.df <- read.csv("C:/path/Python_file/NAME.csv",header = TRUE)

# Load configuration csv
config <- read.csv("C:/path/config_file/NAMEconfig.csv", header = TRUE)

# Join to get information into the data frame 


RTEoutput.m.df <- left_join(RTEoutputraw.df, config, by = "sample")
glimpse(RTEoutput.m.df)


RTEoutput.df <- RTEoutput.m.df 

# Rename column 2
colnames(RTEoutput.df)[2] <- c("target")

# Combine data for each nt by summarizing data from different primers:

RTEoutput.df <- RTEoutput.df %>%
  group_by(sample, target, nt, name, controlsample,
           reagent, replicate, ref, Experiment) %>%
  summarize_each(funs(sum), toA:RT)
glimpse(RTEoutput.df)


# Calculate Pterm, Pmut

RTEoutput.df <- RTEoutput.df %>%
  ungroup() %>%
  mutate(Pstop = stop/(RT + stop),
         Pin = in./depth,
         Pdel = del/depth,
         Pmut = mut/depth,
         PtoA = toA/depth,
         PtoG = toG/depth,
         PtoC = toC/depth,
         PtoT = toT/depth)
#glimpse(RTEoutput.df)


# make a dataframe of only the control samples:

cntlsamples <- levels(as.factor(RTEoutput.df$controlsample))
cntl.df <- RTEoutput.df %>%
  filter(sample %in% cntlsamples)


# Now give the columns of the cntl.df new names:

colnames(cntl.df) <- paste("c", colnames(cntl.df), sep = "")
#glimpse(cntl.df)



RTEoutput.df <- left_join(RTEoutput.df, cntl.df,
                      by = c("controlsample" = "csample",
                             "target" = "ctarget",
                             "nt" = "cnt"))
#glimpse(RTEoutput.df)  
rm(cntl.df)

# Make a function to calculate mean estimates for delta-Pmod:

normToControl <- function(df, col){
  ccol <- (paste("c", col, sep = ""))
  ncol <- (paste("n", col, sep = ""))
  mutate_call = lazyeval::interp(~ a - b, a = as.name(col), b = as.name(ccol))
  mutate_(df, .dots = setNames(list(mutate_call), ncol))
}


# Apply this fucntion to do all the relavent columns:

for(col in c("Pstop", "Pmut", "PtoA", "PtoC", "PtoG", "PtoT", "Pin", "Pdel")) {
  RTEoutput.df <- normToControl(RTEoutput.df, col)
}



# Make a dataframe with the correct structure and names:

RTEoutput.l <- RTEoutput.df %>%
  dplyr::select(target,
                nt,
                ref,
                sample,
                name,
                reagent,
                RT,
                cRT,
                depth,
                cdepth,
                Experiment
                ) %>%
  mutate(rstat = 2.1, cstat = 2.1, nstat = 2.1, stat = "a")
RTEoutput.l <- RTEoutput.l[0,] # Now the dataframe is empty.



for(col in c("Pstop", "Pmut", "PtoA", "PtoC", "PtoG", "PtoT", "Pin", "Pdel")) {
  df <- RTEoutput.df %>%
    dplyr::select(target,
                  nt,
                  ref,
                  sample,
                  name,
                  reagent,
                  RT,
                  cRT,
                  depth,
                  cdepth,
                  Experiment,
                  ends_with(col)) %>%
    mutate(stat = col)
  
  colnames(df)[12:15] <- c("rstat","cstat","nstat","stat")
  
  RTEoutput.l <- bind_rows(RTEoutput.l, df)
  
}
# Filter for Pstop
RTEoutput.l <- RTEoutput.l %>%
  filter(stat %in% c("Pstop"))

#Saving final dataframe
write.csv(RTEoutput.l, file = "Terbium_NAME.l.csv")

