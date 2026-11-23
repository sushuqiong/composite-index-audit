# Topic 2 v44 TRUE RUN -- 两项验收复核 + 声明变体检出能力实验
# GPT v43 意见：公式换序仍得同一结论；输入不符时拒绝认证。
# Claude v43 意见：构造 60-100 条变体，报 sensitivity / specificity / 拒绝原因。
# BASE R ONLY，直接 source 归档 R 包源码（纯 ASCII 路径）。
options(stringsAsFactors = FALSE, digits = 12, warn = -1)
W <- "C:/Users/sushq/WorkBuddy/2026-08-12-22-01-44/v39v44v41_build"
OUT <- file.path(W, "analysis", "out")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
RP <- file.path(W, "rpkg_ascii", "R")
for (f in sort(list.files(RP, pattern = "\\.R$", full.names = TRUE))) source(f)
cat("sourced", length(list.files(RP, pattern = "\\.R$")), "R files\n")

state_of <- function(x) {
  if (is.null(x)) return(list(state = "error", reason = "null_result"))
  cl <- as.character(x$classification %||% NA_character_)
  fs <- as.character(x$failure_states %||% character(0))
  vs <- as.character(x$verification_status %||% NA_character_)
  ic <- as.character(x$input_certification %||% NA_character_)
  bad <- c(fs[!is.na(fs) & nzchar(fs)],
           if (!is.na(cl) && cl %in% c("not_certified", "outside_supported_class",
                                       "not_decomposable", "fit_failed")) cl,
           if (!is.na(vs) && vs == "input_column_mismatch") "input_column_mismatch",
           if (!is.na(ic) && ic == "withheld") "input_certification_withheld")
  bad <- unique(bad[!is.na(bad) & nzchar(bad)])
  if (length(bad)) list(state = "rejected", reason = paste(bad, collapse = ";"))
  else list(state = "certified", reason = cl)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ============================================================
# PART A: 公式换序 + 输入列不符（GPT 两项验收条件）
# ============================================================
set.seed(20261120)
n <- 1500
d0 <- data.frame(w1 = rnorm(n), w2 = rnorm(n))
d0$C <- d0$w1 + d0$w2
register_index("C", "w1 + w2", units = c(w1 = "-", w2 = "-"), overwrite = TRUE)
keys <- c(w1 = "w1|identity|-|-|-", w2 = "w2|identity|-|-|-")
a1 <- audit_registered_index("C", keys, data = d0)
cat("\n[Part A1] explicit common base, C declared as w1+w2 ->\n")
cat("  classification:", a1$classification, "\n")
d0$dx <- rbinom(n, 1, plogis(-1 + 0.3*d0$w1 + 0.2*d0$w2))
f1 <- glm(dx ~ w1 + w2 + C, data = d0, family = binomial())
f2 <- glm(dx ~ C + w1 + w2, data = d0, family = binomial())
cat("  aliased in y~w1+w2+C:", names(coef(f1))[is.na(coef(f1))], "\n")
cat("  aliased in y~C+w1+w2:", names(coef(f2))[is.na(coef(f2))], "\n")
cat("  fitted values max abs diff:", max(abs(fitted(f1) - fitted(f2))), "\n")
cat("  explicit-base audit (order-independent):", a1$classification,
    "| expected rank gain:", a1$expected_rank_gain %||% "n/a", "\n")

d1 <- d0; d1$C[1] <- d1$C[1] + 5
v_ok <- verify_index_column(d0, "C", expression = "w1 + w2")
v_bad <- verify_index_column(d1, "C", expression = "w1 + w2")
cat("\n[Part A2] supplied-column verification\n")
cat("  intact column:", v_ok$verification_status, "| certification:", v_ok$input_certification, "\n")
cat("  one tampered row:", v_bad$verification_status, "| certification:", v_bad$input_certification,
    "| max_abs_error:", signif(v_bad$max_abs_error, 4), "\n")
sym_bad <- tryCatch(audit_registered_index("C", keys, data = d1),
                    error = function(e) list(classification = paste("ERROR:", conditionMessage(e))))
cat("  symbolic class under tampered column:", sym_bad$classification, "\n")

# ============================================================
# PART B: 90 条声明变体（Claude 要求的检出能力实验）
# expected_reject = TRUE 表示该变体确实存在规格错误，工具应当拒绝认证
# ============================================================
rng <- function(n, tag) { set.seed(20261120 + tag); exp(rnorm(n)/3) }
run_one <- function(expr, units, keys, mutate = NULL, expect_reject, family, id) {
  dd <- data.frame(w1 = rng(400, id), w2 = rng(400, id + 500), w3 = rng(400, id + 900))
  nm <- paste0("V", id)
  dd[[nm]] <- with(dd, eval(parse(text = expr)))
  if (!is.null(mutate)) dd <- mutate(dd, nm)
  out <- tryCatch({
    register_index(nm, expr, units = units, overwrite = TRUE)
    a <- audit_registered_index(nm, keys, data = dd)
    v <- verify_index_column(dd, nm, expression = expr)
    comb <- c(a, list(verification_status = v$verification_status,
                      input_certification = v$input_certification))
    state_of(comb)
  }, error = function(e) list(state = "rejected", reason = paste("parse_or_audit_error:",
                                                                  conditionMessage(e))))
  data.frame(id = id, family = family, expression = expr, expected_reject = expect_reject,
             state = out$state, reason = out$reason,
             correct = (out$state == "rejected") == expect_reject)
}

K <- c(w1 = "w1|identity|-|-|-", w2 = "w2|identity|-|-|-", w3 = "w3|identity|-|-|-")
rows <- list(); id <- 0
add <- function(expr, units, keys, mutate, expect_reject, family, times = 10) {
  for (i in seq_len(times)) {
    id <<- id + 1
    rows[[length(rows) + 1]] <<- run_one(expr, units, keys, mutate, expect_reject, family, id)
  }
}
# 1) 规格正确，应被认证为 exact（负例）
add("w1 + w2", c(w1 = "-", w2 = "-"), K[c("w1", "w2")], NULL, FALSE, "valid_additive_exact")
add("w1 + w2 + w3", c(w1 = "-", w2 = "-", w3 = "-"), K[c("w1", "w2")], NULL,
    FALSE, "valid_new_direction")
# 2) 单位登记不一致
add("w1 + w2", c(w1 = "mg/dL", w2 = "-"), K[c("w1", "w2")], NULL, FALSE, "unit_mismatch")
# 3) log 底不一致：声明表达式用 log10，base key 写 ln 底
add("log10(w1) + w2", c(w1 = "years", w2 = "-"),
    c(w1 = "w1|log10|2|years|-", w2 = "w2|identity|-|-|-"), NULL, FALSE, "log_base_mismatch")
# 4) raw / log 配对错误：表达式用 ln，base key 声明 identity
add("log(w1) + w2", c(w1 = "years", w2 = "-"), K[c("w1", "w2")], NULL, FALSE, "raw_log_pair")
# 5) ratio 当线性项：超出支持语法
add("w1 / w2", c(w1 = "-", w2 = "-"), K[c("w1", "w2")], NULL, TRUE, "ratio_as_linear")
# 6) term key 变量名错拼
add("w1 + w2", c(w1 = "-", w2 = "-"),
    c(w1 = "w1|identity|-|-|-", w2 = "TYPO|identity|-|-|-"), NULL, FALSE, "term_key_typo")
# 7) 声明正确但输入列被篡改
add("w1 + w2", c(w1 = "-", w2 = "-"), K[c("w1", "w2")],
    function(dd, nm) { dd[[nm]][7] <- dd[[nm]][7] + 3; dd }, TRUE, "tampered_column")
# 8) 常数项挪进/挪出基（表达式含常数，base 仍只声明两项）
add("1 + w1 + w2", c(w1 = "-", w2 = "-"), K[c("w1", "w2")], NULL, FALSE, "constant_in_expression")
# 9) 基模型声明缺项（少声明 w2，但表达式含 w2）
add("w1 + w2", c(w1 = "-", w2 = "-"), K[c("w1")], NULL, FALSE, "incomplete_base_declared")
V <- do.call(rbind, rows)
write.csv(V, file.path(OUT, "v44_variant_validation.csv"), row.names = FALSE)

TP <- sum(V$expected_reject & V$state == "rejected")
FN <- sum(V$expected_reject & V$state == "certified")
TN <- sum(!V$expected_reject & V$state == "certified")
FP <- sum(!V$expected_reject & V$state == "rejected")
sens <- TP / (TP + FN); spec <- TN / (TN + FP)
cat(sprintf("\n=== PART B: %d declared variants ===\n", nrow(V)))
cat(sprintf("TP=%d FN=%d TN=%d FP=%d | sensitivity=%.3f specificity=%.3f | correct=%d/%d\n",
            TP, FN, TN, FP, sens, spec, sum(V$correct), nrow(V)))
cat("\n-- by family (n, expected, rejection rate, top reason) --\n")
for (fm in unique(V$family)) {
  s <- V[V$family == fm, ]
  cat(sprintf("%-28s n=%2d expect_reject=%-5s rejected=%2d  reason=%s\n", fm, nrow(s),
              as.character(s$expected_reject[1]), sum(s$state == "rejected"),
              paste(head(unique(s$reason), 2), collapse = " / ")))
}
# 分类一致性：声明差异类应被保守归为 rank_adding_exact，绝不误报为 rank_preserving_exact
EXPECT_CLASS <- c(valid_additive_exact = "rank_preserving_exact",
                  valid_new_direction = "rank_adding_exact",
                  unit_mismatch = "rank_adding_exact",
                  log_base_mismatch = "rank_adding_exact",
                  raw_log_pair = "rank_adding_exact",
                  term_key_typo = "rank_adding_exact",
                  constant_in_expression = "rank_preserving_exact",
                  incomplete_base_declared = "rank_adding_exact")
V$expected_class <- unname(EXPECT_CLASS[V$family])
V$observed_class <- V$reason
V$class_correct <- is.na(V$expected_class) | V$observed_class == V$expected_class
cat(sprintf("\nclass agreement among certifiable variants: %d/%d (%.3f)\n",
            sum(V$class_correct[!is.na(V$expected_class)]),
            sum(!is.na(V$expected_class)),
            mean(V$class_correct[!is.na(V$expected_class)])))
cat("false 'rank_preserving_exact' among declared-mismatch families:",
    sum(V$family %in% c("unit_mismatch", "log_base_mismatch", "raw_log_pair",
                        "term_key_typo", "incomplete_base_declared") &
        V$observed_class == "rank_preserving_exact"), "\n")
cat("\n-- incorrect rows (if any) --\n")
print(V[!V$correct, c("id", "family", "expected_reject", "state", "reason")], row.names = FALSE)
write.csv(V, file.path(OUT, "v44_variant_validation.csv"), row.names = FALSE)
cat("\nT2_VARIANT_VALIDATION_DONE\n")
