# Config file for RTEventsCounter.py
# Copyright Yuk Pung Wang 2016
# Simon Lab, Yale University
# File name: conf.py
# Version: 1.0
# Python version: 2.7

# Set output file name (.csv)
outputFileName = "SP082220.csv"

# Overwrite or append when output file already exists
appendIfOutputFileExists = True

# Save log to "counterLog.txt"
saveLog = True

# Nucleotide numbering base in output (0-based or 1-based)
ntBase1 = True

# Output same-nucleotide toN as actual count or "NA"
sameNucAsNA = True

# Mapping qual threshold to include a SAM line
minMapQual = 0

# Long deletions (>1nts)
# "dels"/"stops"/"none"
longDelsAs = "dels"

# Include insertions counts
countInserts = True

# Paired reads behavior
alignPaired = True
ignoreUnpaired = True

# Choose to only use events within a certain range
# 0 = first nucleotide of start of aligned RT sequence (3 prime side in ref seq direction), excluding primer nucs
# [inclusive-start, exclusive-end]; in reverse to event strings (3 prime to 5 prime)
# Both must be non-negative integers
# Should be shorter than the total length of the ref seq
# [0, 0] to toggle off
eventsRange = [0, 1000]

# Primer mode
# Toggle random primer mode
randomPrimers = True

# Absolute number of bases to remove from the end for primers (For random primer mode only)
randPrimerMaskLen = 8

# Set primer input mode (For non-random primer mode only)
# True = sequences / False = coordinates
primerInputModeSeqs = True

# Toggle whether to use 5 prime nt number as primer name automatically (For non-random primer mode only)
startSitePrimerName = True

# Include misprimed reads as primer "NA" (For non-random primer mode only)
includeMisprimed = True

# Start site tolerance +/- for primer designation (For non-random primer mode only)
primerStartTolerance = 5
