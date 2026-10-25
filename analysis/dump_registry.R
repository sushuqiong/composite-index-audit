# dump_registry.R (PURE ASCII)
# Export the frozen registry as a CSV for the manuscript builder and the
# machine-readable results file.
source("index_registry.R")
OUT <- Sys.getenv("CIA_OUT", file.path("..", "results"))
reg <- registry_v34()
cat("columns:", paste(names(reg), collapse = " | "), "\n")
print(head(as.data.frame(reg), 3))
write.csv(reg, file.path(OUT, "v34_registry_table.csv"), row.names = FALSE)
cat("wrote v34_registry_table.csv with", nrow(reg), "rows\n")
