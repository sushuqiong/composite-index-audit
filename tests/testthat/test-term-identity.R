test_that("a term key carries all five fields", {
  k <- make_term_key("tg_mmol_L", "log10", "10", "mmol/L", "-")
  expect_equal(k, "tg_mmol_L|log10|10|mmol/L|-")
})

test_that("log base is part of the identity, not a decoration", {
  k_e  <- make_term_key("tg", "log",   "e",  "mmol/L", "-")
  k_10 <- make_term_key("tg", "log10", "10", "mmol/L", "-")
  expect_equal(unname(term_key_in(k_e, k_e)), "match")
  expect_equal(unname(term_key_in(k_10, k_e)), "mismatch")
})

test_that("declared units that differ never match", {
  k_mmol <- make_term_key("tg", "log10", "10", "mmol/L", "-")
  k_mg   <- make_term_key("tg", "log10", "10", "mg/dL", "-")
  expect_equal(unname(term_key_in(k_mg, k_mmol)), "mismatch")
})

test_that("undeclared units are reported as unkeyed, never as a match", {
  k_declared <- make_term_key("tg", "log10", "10", "mmol/L", "-")
  k_bare     <- make_term_key("tg", "log10", "10")
  expect_equal(unname(term_key_in(k_bare, k_declared)), "unkeyed")
})

test_that("a declared transform never matches identity", {
  k_log <- make_term_key("tg", "log10", "10", "mmol/L", "-")
  k_id  <- make_term_key("tg", "identity", "-", "mmol/L", "-")
  expect_equal(unname(term_key_in(k_id, k_log)), "mismatch")
})
