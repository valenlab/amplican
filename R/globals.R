# Centralized global variables for data.table and other NSE contexts (including ggplot2)

utils::globalVariables(c(
  ".", "Forward", "ID", "Reverse", "Total", "consensus", "net_width", 
  "originally", "readType", "read_counts", "read_id", "replacement",
  "seqnames", "counts", "start", "end", "strand", "overlaps", "events", 
  "all_e_not_overlap", "unassignedTable", "Low_Score", "Reads_Filtered", 
  "Reads", "PRIMER_DIMER", "ID", "Barcode", "Forward_Reads_File", 
  "Reverse_Reads_File", "Group", "guideRNA", "Found_Guide", "Control", 
  "Forward_Primer", "Reverse_Primer", "Direction", "Amplicon", "Donor", 
  "fwdPrPosEnd", "rvePrPos", "Reads_Del", "Reads_In", "Reads_Edited", 
  "Reads_Frameshifted", "HDR", "fwd_idx", "rve_idx", "extra", "n", 
  "read_shares", "has_HDR", "has_Del", "has_In", "is_FS", "has_Edit",
  "i.HDR", "i.Reads_Del", "i.Reads_In", "i.Reads_Edited", "i.Reads_Frameshifted",
  "is_WT", "percentage", "variable", "value", "frequency", "num",
  "x", "y", "category", "group", "frequencyReal", "nucleotide", "upper",
  "xmin", "xmax", "ymin", "ymax", "codon", "position", "coverage"
))
