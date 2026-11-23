# Base-R check entry point. The test files remain compatible with testthat;
# see testthat-standard.R.txt for the optional standard test_check() entry.
library(CompositeIndexAudit)
source("base-test-harness.R", local = FALSE)
test_files <- sort(list.files("testthat", pattern = "^test.*[.]R$", full.names = TRUE))
stopifnot(length(test_files) == 9L)
for (file in test_files) {
  sys.source(file, envir = globalenv())
}
cat(sprintf("BASE_TEST_SUMMARY: blocks=%d PASS=%d FAIL=%d\n",
            .cia_harness$blocks, .cia_harness$pass, .cia_harness$fail))
stopifnot(.cia_harness$blocks == 49L,
          .cia_harness$pass == 447L,
          .cia_harness$fail == 0L)
