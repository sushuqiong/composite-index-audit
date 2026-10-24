# dump_registry_v33.R (PURE ASCII)
# Export the frozen registry as a CSV for the manuscript builder and the
# machine-readable results file.
source("index_registry_v33.R")
OUT <- Sys.getenv("CIA_OUT", file.path("..", "results"))
reg <- v33_registry()
cat("columns:", paste(names(reg), collapse = " | "), "\n")
print(head(as.data.frame(reg), 3))
write.csv(reg, file.path(OUT, "v33_registry_table.csv"), row.names = FALSE)
cat("wrote v33_registry_table.csv with", nrow(reg), "rows\n")
