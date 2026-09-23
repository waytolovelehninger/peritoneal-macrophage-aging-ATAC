###############################################################################
# engulfment_ATAC_promoter_plots.R
#
# 목적: young vs old mouse peritoneal macrophage ATAC-seq에서
#       engulfment 유전자 promoter peak의 chromatin accessibility 변화 정리
# 데이터: Lu et al., bioRxiv 2025.06.11.659200 (Benayoun lab), PRJNA1095555
#         bioRxiv supplementary media-3.xlsx (= Supplementary Table S2)
#         sheet B = ATAC-seq peak 단위 DESeq2 결과
#         _F = female, _M = male (성별별 DESeq2)
#         log2FoldChange = log2 fold change per month of age (양수 = old에서 증가)
#         나이: README 참고. 논문 본문 young 4-5 / old 20-21 mo, Table S1 계산값
#               young 6.54 / old 22.57 mo (불일치는 확인 필요). 두 기준 모두 차이 약 16개월
#
# 분석 방침
#   - 성별 분리 (female old vs young, male old vs young). 성별 통합 없음
#   - 유전자별 promoter peak 1개: Distance.to.TSS in [-1000, +1000] (메인; 인자로 변경 가능),
#     여러 개면 |Distance.to.TSS| 최소 (동률이면 Start가 작은 것). 유의성 기준 선택 안 함
#   - per-month log2FC, lfcSE x 16 -> old vs young 근사값 (p, padj는 그대로)
#   - padj < 0.05 (* <0.05, ** <0.01, *** <0.001)
#
# 실행 (저장소 최상위 폴더에서):
#   메인 (+-1kb, results/):
#     Rscript scripts/engulfment_ATAC_promoter_plots.R
#   민감도 (+-2kb, results/sensitivity_2kb/, 메인 결과와 비교표 생성):
#     Rscript scripts/engulfment_ATAC_promoter_plots.R 2000 results/sensitivity_2kb results 1000
#   인자 (모두 선택, 순서대로):
#     1) window_bp   : promoter 범위 +-bp (기본 1000)
#     2) out_dir     : 출력 폴더 (기본 results)
#     3) compare_dir : 비교할 기존 결과 폴더 (intermediate/plot_table.rds 필요). 주면 비교표를
#                      compare_dir/tables/promoter_window_comparison.csv 로 저장
#     4) compare_window_bp : compare_dir 결과의 promoter 범위 +-bp (기본 1000)
#
# 출력 (out_dir 아래):
#   figures/      barplot, dotplot (PNG, PDF)
#   tables/       promoter_peaks_selected.csv, plot_table.csv, Tim4_summary.txt
#   intermediate/ .rds, promoter_candidates_all.csv, 00_run_summary.txt
#   sessionInfo.txt는 메인 실행이면 저장소 최상위, 그 외에는 out_dir/intermediate/
###############################################################################


# ============================================================================
# 0. 패키지 설치 / 로드
# ============================================================================
required_pkgs <- c("readxl", "dplyr", "tidyr", "ggplot2", "readr")

# TRUE로 바꾸면 누락 패키지를 CRAN에서 설치. 기본은 설치하지 않고 멈춤
auto_install <- FALSE

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  if (auto_install) {
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
  } else {
    stop("누락된 패키지: ", paste(missing_pkgs, collapse = ", "),
         "\n설치: install.packages(c(",
         paste0('"', missing_pkgs, '"', collapse = ", "), "))")
  }
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})


# ============================================================================
# 1. 조정 가능한 설정값
# ============================================================================
# 모든 경로는 저장소 최상위 폴더 기준 상대경로
if (!file.exists(file.path("scripts", "engulfment_ATAC_promoter_plots.R"))) {
  stop("저장소 최상위 폴더에서 실행하세요: Rscript scripts/engulfment_ATAC_promoter_plots.R")
}

input_xlsx  <- file.path("data", "media-3.xlsx")     # bioRxiv Supplementary Table S2
input_sheet <- "B"                                   # ATAC-seq peak DESeq2 결과
check_xlsx  <- file.path("data", "ATAC.xlsx")        # (선택) sheet B만 따로 저장한 파일과 동일성 확인
check_against_atac_xlsx <- FALSE

if (!file.exists(input_xlsx)) {
  stop("입력 파일 없음: ", input_xlsx,
       "\n받는 방법은 data/README_data.md 참고")
}

main_out_dir    <- "results"
main_window_bp  <- 1000

# 명령행 인자 (없으면 메인 설정). 사용법은 파일 상단 참고
cli_args <- commandArgs(trailingOnly = TRUE)
arg_or <- function(i, default) {
  if (length(cli_args) >= i && nzchar(cli_args[i])) cli_args[i] else default
}
window_bp         <- as.numeric(arg_or(1, main_window_bp))
out_dir           <- arg_or(2, main_out_dir)
compare_dir       <- arg_or(3, NA_character_)
compare_window_bp <- as.numeric(arg_or(4, main_window_bp))
stopifnot(!is.na(window_bp), window_bp > 0, !is.na(compare_window_bp))

# 메인 결과 보호: 메인 폴더에는 메인 범위로만 쓸 수 있음
same_dir <- function(a, b) {
  normalizePath(a, winslash = "/", mustWork = FALSE) ==
    normalizePath(b, winslash = "/", mustWork = FALSE)
}
is_main_run <- same_dir(out_dir, main_out_dir)
if (is_main_run && window_bp != main_window_bp) {
  stop("메인 결과 덮어쓰기 방지: window_bp = ", window_bp,
       " 이면 out_dir을 메인 폴더가 아닌 곳으로 지정하세요.")
}

# promoter 정의 (bp, Distance.to.TSS 기준, 양 끝 포함)
tss_window <- c(-window_bp, window_bp)

# 값 스케일: "approx_total" = per-month 계수 x months_diff (old vs young 근사)
#            "per_month"    = 원본 per-month 계수 그대로
scale_mode  <- "approx_total"
months_diff <- 16                                    # 두 그룹 나이 차 (약 16개월, README 참고)
scale_factor <- if (scale_mode == "approx_total") months_diff else 1

# 유의성
padj_cutoff <- 0.05
star_breaks <- c(0.001, 0.01, 0.05)                  # ***, **, *

# 색: old에서 감소 = 파랑, 증가 = 빨강
col_down  <- "#2166AC"
col_up    <- "#B2182B"
col_nopk  <- "grey50"

# 그림 크기 (inch) / 해상도
bar_width_in  <- 9;   bar_height_in <- 8.5
dot_width_in  <- 5.8; dot_height_in <- 8.5
png_dpi       <- 300

# 유전자 패널 (표시 이름, 파일의 공식 심볼). 이 순서대로 위에서 아래로 표시
gene_panel <- data.frame(
  display = c("Tim4", "Tim1", "Tim3", "Mertk", "Axl", "Tyro3", "Gas6", "Pros1",
              "Mfge8", "Itgav", "Itgb3", "Itgb5", "Bai1", "Bai3", "Stab1", "Stab2",
              "Cd36", "Cd14", "Cd300lb", "Cd300lf", "Lrp1", "Megf10", "Scarf1",
              "Rage", "DD1a"),
  symbol  = c("Timd4", "Havcr1", "Havcr2", "Mertk", "Axl", "Tyro3", "Gas6", "Pros1",
              "Mfge8", "Itgav", "Itgb3", "Itgb5", "Adgrb1", "Adgrb3", "Stab1", "Stab2",
              "Cd36", "Cd14", "Cd300lb", "Cd300lf", "Lrp1", "Megf10", "Scarf1",
              "Ager", "Vsir"),
  stringsAsFactors = FALSE
)
gene_panel$label <- ifelse(gene_panel$display == gene_panel$symbol,
                           gene_panel$symbol,
                           paste0(gene_panel$display, " (", gene_panel$symbol, ")"))

sexes <- c(F = "Female", M = "Male")


value_caption <- if (scale_mode == "approx_total") {
  paste0("log2FC = DESeq2 per-month coefficient x ", months_diff,
         " (approximate old vs young difference)")
} else {
  "log2FC = DESeq2 per-month coefficient (log2 fold change per month of age)"
}
# error bar 설명은 barplot에만 붙임 (dotplot에는 error bar 없음)
errbar_caption <- if (scale_mode == "approx_total") {
  paste0("error bars = lfcSE x ", months_diff)
} else {
  "error bars = lfcSE"
}
sig_caption <- "* padj < 0.05, ** padj < 0.01, *** padj < 0.001 (DESeq2 padj, unscaled)"
promoter_caption <- paste0("Promoter peak = peak closest to TSS within ",
                           tss_window[1], " to +", tss_window[2],
                           " bp (Distance.to.TSS); one peak per gene")


# ============================================================================
# 2. 출력 폴더
# ============================================================================
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
int_dir <- file.path(out_dir, "intermediate")
for (d in c(fig_dir, tab_dir, int_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
out_fig <- function(f) file.path(fig_dir, f)
out_tab <- function(f) file.path(tab_dir, f)
out_int <- function(f) file.path(int_dir, f)
session_info_path <- if (is_main_run) "sessionInfo.txt" else out_int("sessionInfo.txt")

write_utf8 <- function(lines, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con))
  writeLines(lines, con)
}

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  message(msg)
  run_log <<- c(run_log, msg)
}

logmsg("== Settings ==")
logmsg("Promoter window: ", tss_window[1], " to +", tss_window[2], " bp")
logmsg("Output dir: ", out_dir)
if (!is.na(compare_dir)) {
  logmsg("Compare with: ", compare_dir, " (+-", compare_window_bp, " bp)")
}
logmsg("")


# ============================================================================
# 3. 입력 읽기 및 확인
# ============================================================================
# 메모리 부족 시 대안: 전체 시트(약 7.2만 행 x 20열)를 text로 읽는 단계가 가장 무겁다.
#   (a) read_excel(..., range = cell_cols("A:T"))로 필요한 열만 지정하거나,
#   (b) 한 번 CSV로 변환해 두고 readr::read_csv(col_select = ...)로 필요한 열만 읽거나,
#   (c) check_against_atac_xlsx <- FALSE로 두 번째 파일 읽기를 생략 (기본값).
logmsg("== Input ==")
logmsg("File: ", input_xlsx, " / sheet: ", input_sheet)
logmsg("Sheets in file: ", paste(excel_sheets(input_xlsx), collapse = ", "))

raw_txt <- read_excel(input_xlsx, sheet = input_sheet, col_types = "text", na = "")
logmsg("Rows: ", nrow(raw_txt), " / Columns: ", ncol(raw_txt))
logmsg("Column names: ", paste(names(raw_txt), collapse = ", "))

na_string_counts <- vapply(raw_txt, function(x) sum(x == "NA", na.rm = TRUE), integer(1))
logmsg("Literal \"NA\" strings: ",
       paste0(names(na_string_counts)[na_string_counts > 0], "=",
              na_string_counts[na_string_counts > 0], collapse = ", "))

if (check_against_atac_xlsx && file.exists(check_xlsx)) {
  chk <- read_excel(check_xlsx, sheet = excel_sheets(check_xlsx)[1],
                    col_types = "text", na = "")
  same <- identical(as.data.frame(chk), as.data.frame(raw_txt))
  logmsg("Identical to ", basename(check_xlsx), " (all cells as text): ", same,
         " (rows ", nrow(chk), ", cols ", ncol(chk), ")")
  rm(chk); invisible(gc())
}

num_cols <- c("Start", "End", "Distance.to.TSS",
              as.vector(outer(c("baseMean", "log2FoldChange", "lfcSE", "stat",
                                "pvalue", "padj"), names(sexes), paste, sep = "_")))
atac <- raw_txt %>%
  mutate(across(everything(), ~ ifelse(.x == "NA", NA_character_, .x))) %>%
  mutate(across(all_of(num_cols), as.numeric))

# 숫자 변환 실패 확인 (원래 값이 있었는데 NA가 된 셀)
conv_fail <- vapply(num_cols, function(cn) {
  sum(!is.na(raw_txt[[cn]]) & raw_txt[[cn]] != "NA" & is.na(atac[[cn]]))
}, integer(1))
if (any(conv_fail > 0)) {
  logmsg("WARNING numeric conversion failures: ",
         paste0(names(conv_fail)[conv_fail > 0], "=", conv_fail[conv_fail > 0],
                collapse = ", "))
} else {
  logmsg("Numeric conversion: no failures")
}
rm(raw_txt); invisible(gc())

saveRDS(atac, out_int("atac_sheetB_parsed.rds"))


# ============================================================================
# 4. 패널 유전자의 promoter 후보 peak
# ============================================================================
logmsg("")
logmsg("== Gene panel ==")
absent <- setdiff(gene_panel$symbol, atac$Gene.Name)
logmsg("Panel symbols with no peak at all in file: ",
       if (length(absent)) paste(absent, collapse = ", ") else "none")

peaks_per_gene <- atac %>%
  filter(Gene.Name %in% gene_panel$symbol) %>%
  count(Gene.Name, name = "n_peaks_all")

candidates <- atac %>%
  filter(Gene.Name %in% gene_panel$symbol,
         !is.na(Distance.to.TSS),
         Distance.to.TSS >= tss_window[1],
         Distance.to.TSS <= tss_window[2]) %>%
  mutate(abs_dist_to_TSS = abs(Distance.to.TSS)) %>%
  group_by(Gene.Name) %>%
  arrange(abs_dist_to_TSS, Start, .by_group = TRUE) %>%   # TSS 최근접, 동률이면 Start 작은 것
  mutate(n_promoter_candidates = n(),
         candidate_rank = row_number()) %>%
  ungroup()

write_csv(candidates, out_int("promoter_candidates_all.csv"), na = "")


# ============================================================================
# 5. 유전자별 promoter peak 1개 선택 (성별 공통: 거리 기준이므로 동일 peak)
# ============================================================================
selected <- candidates %>% filter(candidate_rank == 1)

for (s in names(sexes)) {
  selected[[paste0("log2FC_x", scale_factor, "_", s)]] <-
    selected[[paste0("log2FoldChange_", s)]] * scale_factor
  selected[[paste0("lfcSE_x", scale_factor, "_", s)]] <-
    selected[[paste0("lfcSE_", s)]] * scale_factor
}

selected_out <- gene_panel %>%
  select(display, symbol) %>%
  inner_join(selected, by = c("symbol" = "Gene.Name")) %>%
  left_join(peaks_per_gene, by = c("symbol" = "Gene.Name"))
write_csv(selected_out, out_tab("promoter_peaks_selected.csv"), na = "")
saveRDS(selected_out, out_int("promoter_peaks_selected.rds"))

coverage <- gene_panel %>%
  left_join(peaks_per_gene, by = c("symbol" = "Gene.Name")) %>%
  left_join(distinct(candidates, Gene.Name, n_promoter_candidates),
            by = c("symbol" = "Gene.Name")) %>%
  mutate(n_peaks_all = coalesce(n_peaks_all, 0L),
         n_promoter_candidates = coalesce(n_promoter_candidates, 0L))

cov_one  <- coverage$symbol[coverage$n_promoter_candidates == 1]
cov_multi <- coverage$symbol[coverage$n_promoter_candidates > 1]
cov_none <- coverage$symbol[coverage$n_promoter_candidates == 0]
logmsg("Panel genes: ", nrow(gene_panel))
logmsg("  promoter peak 1개: ", length(cov_one), " genes (", paste(cov_one, collapse = ", "), ")")
logmsg("  promoter peak 2개 이상 (TSS 최근접 선택): ", length(cov_multi), " genes",
       if (length(cov_multi)) paste0(" (", paste(cov_multi, collapse = ", "), ")") else "")
logmsg("  promoter peak 없음: ", length(cov_none), " genes (", paste(cov_none, collapse = ", "), ")")


# ============================================================================
# 6. plot용 long table (유전자 x 성별)
# ============================================================================
star_of <- function(p) {
  case_when(is.na(p)             ~ "",
            p < star_breaks[1]   ~ "***",
            p < star_breaks[2]   ~ "**",
            p < star_breaks[3]   ~ "*",
            TRUE                 ~ "")
}

plot_tab <- bind_rows(lapply(names(sexes), function(s) {
  gene_panel %>%
    left_join(selected, by = c("symbol" = "Gene.Name")) %>%
    transmute(
      display, symbol, label,
      sex = sexes[[s]],
      PeakName, Chr, Start, End, Annotation, Distance.to.TSS,
      baseMean        = .data[[paste0("baseMean_", s)]],
      log2FC_per_month = .data[[paste0("log2FoldChange_", s)]],
      lfcSE_per_month  = .data[[paste0("lfcSE_", s)]],
      stat            = .data[[paste0("stat_", s)]],
      pvalue          = .data[[paste0("pvalue_", s)]],
      padj            = .data[[paste0("padj_", s)]]
    )
})) %>%
  mutate(
    scale_factor = scale_factor,
    log2FC_plot  = log2FC_per_month * scale_factor,
    lfcSE_plot   = lfcSE_per_month * scale_factor,
    status       = ifelse(is.na(PeakName), "no promoter peak", "promoter peak"),
    significant  = !is.na(padj) & padj < padj_cutoff,
    sig_label    = ifelse(status == "promoter peak" & is.na(padj), "padj NA", star_of(padj)),
    direction    = case_when(is.na(log2FC_plot) ~ NA_character_,
                             log2FC_plot < 0    ~ "Decreased in old",
                             TRUE               ~ "Increased in old"),
    label        = factor(label, levels = rev(gene_panel$label)),
    sex          = factor(sex, levels = unname(sexes))
  )

write_csv(plot_tab, out_tab("plot_table.csv"), na = "")
saveRDS(plot_tab, out_int("plot_table.rds"))


# ============================================================================
# 7. Barplot
# ============================================================================
pk   <- filter(plot_tab, status == "promoter peak")
nopk <- filter(plot_tab, status == "no promoter peak")

x_lim <- max(abs(c(pk$log2FC_plot - pk$lfcSE_plot, pk$log2FC_plot + pk$lfcSE_plot)),
             na.rm = TRUE) * 1.3
star_off <- x_lim * 0.03
pk <- pk %>%
  mutate(star_x = ifelse(log2FC_plot >= 0,
                         log2FC_plot + lfcSE_plot + star_off,
                         log2FC_plot - lfcSE_plot - star_off),
         star_hjust = ifelse(log2FC_plot >= 0, 0, 1))

x_title <- if (scale_mode == "approx_total") {
  "log2 fold change, old vs young (approx.)"
} else {
  "log2 fold change per month of age"
}
dir_fill <- c("Decreased in old" = col_down, "Increased in old" = col_up)

p_bar <- ggplot() +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_col(data = pk, aes(x = log2FC_plot, y = label, fill = direction), width = 0.7) +
  geom_errorbar(data = pk, aes(y = label, xmin = log2FC_plot - lfcSE_plot,
                               xmax = log2FC_plot + lfcSE_plot),
                width = 0.3, linewidth = 0.35, colour = "grey20") +
  geom_text(data = pk, aes(x = star_x, y = label, label = sig_label, hjust = star_hjust),
            size = 3.6, vjust = 0.75) +
  # 흰 배경 박스로 0 기준선을 가림 (기준선보다 나중에 그려야 함)
  geom_label(data = nopk, aes(x = 0, y = label, label = "no promoter peak"),
             colour = col_nopk, fontface = "italic", size = 2.8,
             fill = "white", border.colour = NA, label.padding = unit(0.12, "lines")) +
  facet_wrap(~ sex, nrow = 1, drop = FALSE) +
  scale_fill_manual(values = dir_fill, drop = FALSE, name = NULL) +
  scale_x_continuous(limits = c(-x_lim, x_lim)) +
  scale_y_discrete(drop = FALSE) +
  labs(title = "Promoter chromatin accessibility of engulfment genes",
       subtitle = "Peritoneal macrophages, old vs young, ATAC-seq, by sex",
       x = x_title, y = NULL,
       caption = paste(paste0(value_caption, "; ", errbar_caption),
                       sig_caption, promoter_caption,
                       "Source: Lu et al., bioRxiv 2025.06.11.659200, Table S2 sheet B",
                       sep = "\n")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold"),
        legend.position = "top",
        plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey25"))

ggsave(out_fig("ATAC_promoter_barplot.pdf"), p_bar, width = bar_width_in, height = bar_height_in)
ggsave(out_fig("ATAC_promoter_barplot.png"), p_bar, width = bar_width_in, height = bar_height_in,
       dpi = png_dpi)


# ============================================================================
# 8. Dotplot
# ============================================================================
pk <- pk %>% mutate(neglog10_padj = -log10(padj))
fill_lim <- max(abs(pk$log2FC_plot), na.rm = TRUE)

p_dot <- ggplot() +
  geom_point(data = pk, aes(x = sex, y = label, fill = log2FC_plot, size = neglog10_padj),
             shape = 21, colour = "grey30", stroke = 0.3) +
  geom_text(data = pk, aes(x = sex, y = label, label = sig_label),
            position = position_nudge(x = 0.28), size = 3.4, vjust = 0.75) +
  geom_text(data = nopk, aes(x = sex, y = label, label = "no promoter peak"),
            colour = col_nopk, fontface = "italic", size = 2.4) +
  scale_fill_gradient2(low = col_down, mid = "white", high = col_up, midpoint = 0,
                       limits = c(-fill_lim, fill_lim), name = x_title) +
  scale_size_continuous(range = c(1.5, 7), name = "-log10(padj)") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(title = "Promoter ATAC, old vs young",
       subtitle = "Peritoneal macrophages, by sex",
       x = NULL, y = NULL,
       caption = paste(value_caption, sig_caption, promoter_caption, sep = "\n")) +
  guides(fill = guide_colourbar(order = 1, barheight = unit(3.5, "cm")),
         size = guide_legend(order = 2)) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.title = element_text(size = 8.5),
        plot.caption = element_text(hjust = 0, size = 6.5, colour = "grey25"))

ggsave(out_fig("ATAC_promoter_dotplot.pdf"), p_dot, width = dot_width_in, height = dot_height_in)
ggsave(out_fig("ATAC_promoter_dotplot.png"), p_dot, width = dot_width_in, height = dot_height_in,
       dpi = png_dpi)


# ============================================================================
# 9. Tim4 (Timd4) 결과 요약
# ============================================================================
fmt <- function(x, d = 4) ifelse(is.na(x), "NA", formatC(x, format = "g", digits = d))

tim4 <- filter(plot_tab, symbol == "Timd4")
tim4_all_n <- coverage$n_peaks_all[coverage$symbol == "Timd4"]
tim4_lines <- c(
  "Tim4 (Timd4) promoter ATAC-seq: old vs young peritoneal macrophages",
  "Source: Lu et al., bioRxiv 2025.06.11.659200, Supplementary Table S2 sheet B (DESeq2, sex-stratified)",
  paste0("log2FC x", scale_factor, " = DESeq2 per-month coefficient x ", scale_factor,
         " (approximate old vs young difference). p, padj = original values."),
  paste0("Timd4 assigned peaks in file: ", tim4_all_n,
         "; promoter peaks (", tss_window[1], " to +", tss_window[2], " bp): ",
         coverage$n_promoter_candidates[coverage$symbol == "Timd4"]),
  ""
)
if (all(tim4$status == "no promoter peak")) {
  tim4_lines <- c(tim4_lines, "No promoter peak for Timd4.")
} else {
  r1 <- tim4[1, ]
  tim4_lines <- c(tim4_lines,
    paste0("Peak: ", r1$PeakName, "  (", r1$Chr, ":", r1$Start, "-", r1$End, ")"),
    paste0("Distance to TSS: ", r1$Distance.to.TSS, " bp;  Annotation: ", r1$Annotation),
    "")
  for (i in seq_len(nrow(tim4))) {
    r <- tim4[i, ]
    tim4_lines <- c(tim4_lines,
      paste0("[", r$sex, "]"),
      paste0("  log2FC x", scale_factor, " = ", fmt(r$log2FC_plot),
             "   (per-month = ", fmt(r$log2FC_per_month), ")"),
      paste0("  lfcSE  x", scale_factor, " = ", fmt(r$lfcSE_plot),
             "   (per-month = ", fmt(r$lfcSE_per_month), ")"),
      paste0("  pvalue = ", fmt(r$pvalue), ";  padj = ", fmt(r$padj),
             "  ", ifelse(r$significant, paste0("significant (", r$sig_label, ")"),
                          paste0("not significant (padj >= ", padj_cutoff, ")"))),
      paste0("  direction: ", r$direction),
      "")
  }
}
write_utf8(tim4_lines, out_tab("Tim4_summary.txt"))
cat(tim4_lines, sep = "\n")


# ============================================================================
# 9b. (선택) 다른 promoter 범위 결과와 비교 -- compare_dir 인자가 있을 때만
# ============================================================================
if (!is.na(compare_dir)) {
  ref_path <- file.path(compare_dir, "intermediate", "plot_table.rds")
  if (!file.exists(ref_path)) stop("비교 대상 없음: ", ref_path)
  if (same_dir(compare_dir, out_dir)) stop("compare_dir와 out_dir가 같음")
  ref_tab <- readRDS(ref_path)
  if (!all(ref_tab$scale_factor == scale_factor)) {
    stop("비교 대상의 scale_factor가 다름: ", paste(unique(ref_tab$scale_factor), collapse = ","))
  }

  tag_ref <- paste0(compare_window_bp / 1000, "kb")
  tag_cur <- paste0(window_bp / 1000, "kb")
  lfc_name <- paste0("log2FC_x", scale_factor)

  # 유전자당 1행: peak, TSS 거리, 성별 log2FC xN / padj / 별표 (peak는 성별 공통)
  widen <- function(pt, tag) {
    w <- pt %>%
      mutate(sex_s = substr(as.character(sex), 1, 1)) %>%
      select(symbol, sex_s, PeakName, Distance.to.TSS, log2FC_plot, padj, sig_label) %>%
      pivot_wider(id_cols = c(symbol, PeakName, Distance.to.TSS),
                  names_from = sex_s,
                  values_from = c(log2FC_plot, padj, sig_label),
                  names_glue = "{.value}_{sex_s}")
    names(w) <- sub("^log2FC_plot", lfc_name, names(w))
    rename_with(w, ~ paste0(.x, "_", tag), -symbol)
  }

  cmp <- gene_panel %>%
    select(display, symbol) %>%
    left_join(widen(ref_tab, tag_ref), by = "symbol") %>%
    left_join(widen(plot_tab, tag_cur), by = "symbol")

  p_ref <- cmp[[paste0("PeakName_", tag_ref)]]
  p_cur <- cmp[[paste0("PeakName_", tag_cur)]]
  cmp$peak_changed <- !((is.na(p_ref) & is.na(p_cur)) |
                          (!is.na(p_ref) & !is.na(p_cur) & p_ref == p_cur))
  cmp$change_type <- case_when(
    is.na(p_ref) & is.na(p_cur) ~ "no promoter peak in both",
    is.na(p_ref)                ~ paste0("peak added (none at +-", tag_ref, ")"),
    is.na(p_cur)                ~ paste0("peak lost (none at +-", tag_cur, ")"),
    p_ref == p_cur              ~ "same peak",
    TRUE                        ~ "different peak"
  )

  col_order <- c("display", "symbol",
                 paste0("PeakName_", c(tag_ref, tag_cur)),
                 paste0("Distance.to.TSS_", c(tag_ref, tag_cur)),
                 "peak_changed", "change_type",
                 as.vector(t(outer(
                   c(paste0(lfc_name, "_F"), "padj_F", "sig_label_F",
                     paste0(lfc_name, "_M"), "padj_M", "sig_label_M"),
                   c(tag_ref, tag_cur), paste, sep = "_"))))
  cmp <- cmp[, col_order]
  cmp_path <- file.path(compare_dir, "tables", "promoter_window_comparison.csv")
  write_csv(cmp, cmp_path, na = "")

  # 유의 유전자(성별 구분) 목록 비교
  sig_keys <- function(pt) {
    pt %>%
      filter(significant) %>%
      arrange(match(symbol, gene_panel$symbol), sex) %>%
      transmute(k = paste0(symbol, "(", substr(as.character(sex), 1, 1), ")")) %>%
      pull(k)
  }
  s_ref <- sig_keys(ref_tab)
  s_cur <- sig_keys(plot_tab)
  fmt_keys <- function(k) if (length(k)) paste(k, collapse = ", ") else "없음"
  added   <- setdiff(s_cur, s_ref)
  removed <- setdiff(s_ref, s_cur)

  logmsg("")
  logmsg("== Promoter window sensitivity (+-", tag_ref, " vs +-", tag_cur, ") ==")
  logmsg("Comparison table: ", cmp_path)
  logmsg("선택 peak 변경: ", sum(cmp$peak_changed), " genes",
         if (any(cmp$peak_changed)) {
           paste0(" (", paste0(cmp$symbol[cmp$peak_changed], ": ",
                               cmp$change_type[cmp$peak_changed], collapse = "; "), ")")
         } else "")
  logmsg(if (length(added) + length(removed) == 0) {
    paste0("기준을 +-", tag_cur, "로 넓혀도 유의한 유전자 목록(padj < ", padj_cutoff,
           ")은 바뀌지 않음: ", fmt_keys(s_cur))
  } else {
    paste0("기준을 +-", tag_cur, "로 넓히면 유의한 유전자 목록(padj < ", padj_cutoff,
           ")이 바뀜: 추가 = ", fmt_keys(added), "; 제외 = ", fmt_keys(removed),
           " (+-", tag_ref, ": ", fmt_keys(s_ref), ")")
  })
}


# ============================================================================
# 10. 실행 요약 + 확인 필요 사항
# ============================================================================
caveats <- c(
  "",
  "== 확인 필요 ==",
  "- Distance.to.TSS 부호(음수 = upstream 여부)는 확인되지 않음. 단, promoter 범위가",
  paste0("  ", tss_window[1], " ~ +", tss_window[2],
         " bp 대칭이고 |거리| 최소로 고르므로 peak 선택 결과에는 영향 없음."),
  paste0("- peak는 최근접 유전자 1개에만 할당되어 있음. 'no promoter peak'는 '해당 유전자에 할당된 +-",
         window_bp, " bp 내"),
  "  peak 없음'이지 promoter가 닫혀 있다는 뜻이 아님 (이웃 유전자에 할당되었거나 대체 TSS 사용 가능).",
  "- 나이: 논문 본문(young 4-5, old 20-21개월)과 Table S1 계산값(young 6.54, old 22.57개월)이",
  "  불일치 (확인 필요). 두 기준 모두 차이가 약 16개월이라 x16 환산과 p/padj에는 영향 없음.",
  paste0("  DESeq2에서 나이를 어떤 값으로 넣었는지도 확인 필요. x", months_diff, "은 근사값이며 ",
         "log2FC/lfcSE x", months_diff, "에만 해당."),
  "- DESeq2 모델 세부(나이를 연속 변수로 사용, 공변량 등)는 원 논문 methods로 확인 필요."
)
write_utf8(c(run_log, caveats), out_int("00_run_summary.txt"))


# ============================================================================
# 11. sessionInfo
# ============================================================================
si <- capture.output(sessionInfo())
cat(si, sep = "\n")
write_utf8(si, session_info_path)
