
# 02_appendix_indicator_table.R
# Builds the appendix table of ESG indicators from ind_class.rds,
# which 01_data_preparation.R creates. Every value in the table comes from
# that file:
#   type, direction, unit -> indicator classification (01, section 7)
#   imputation            -> rule that follows from type and direction (01, section 7)
#   report_rate           -> share of training firm-years FY2015-2020 with a
#                            reported value (01, section 6)
#   make_M                -> TRUE if report_rate < 0.99, i.e. a non-disclosure
#                            dummy M_k is created (01, section 7)
# Output: appendix_indicators.tex  (use \input{appendix_indicators} or paste)
# =============================================================================

library(dplyr)

setwd("~/Desktop/Thesis_ESG_data/Data")
ind <- readRDS("ind_class.rds")
stopifnot(nrow(ind) == 27)

# Summary used in the text (Sections 4.2.1 and 4.4.3)
print(count(ind, type))            # input / outcome
print(count(ind, type, make_M))    # which indicators get a dummy

unit_lab <- c(binary = "Binary", percent = "Percent", count = "Count",
              absolute = "Absolute$^{\\dagger}$", intensity = "Intensity",
              rate = "Rate")
imp_lab  <- c(median = "Median", p25 = "25th pct.", p75 = "75th pct.")

row_tex <- function(d) {
  sprintf("%s & %s & %s & %s & %.1f & %s \\\\",
          gsub("&", "\\\\&", d$indicator),
          unit_lab[d$unit],
          tools::toTitleCase(d$direction),
          imp_lab[d$imputation],
          100 * d$report_rate,
          ifelse(d$make_M, "Yes", "No"))
}

block <- function(t, title) {
  d <- ind |> filter(type == t) |> arrange(desc(report_rate))
  c(sprintf("\\multicolumn{6}{@{}l}{\\textit{%s}} \\\\", title), row_tex(d))
}

head_tex <- r"---(\begin{table}[htbp]
\centering
\footnotesize
\caption{ESG indicators used in the analysis}
\label{tab:indicators}
\renewcommand{\arraystretch}{1.12}
\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}X l l l r c@{}}
\toprule
\textbf{Indicator} & \textbf{Unit} & \textbf{Higher is} & \textbf{Imputed at} & \textbf{Reported (\%)} & \textbf{$M_k$} \\
\midrule)---"

foot_tex <- r"---(\bottomrule
\end{tabularx}

\vspace{4pt}
{\footnotesize\itshape\raggedright
Notes: The 27 indicators reported for at least 20\% of training firm-years (FY2015--FY2020), sorted by reporting rate. Input metrics capture policies, commitments or actions; outcome metrics capture realised performance (Christensen et al., 2022). Imputed at: value of the industry-year distribution used for missing values (Section~\ref{sec:imputation}). Reported: percentage of training firm-years in which the indicator is reported. $M_k$: whether a non-disclosure indicator is created (indicators reported in less than 99\% of training firm-years). $^{\dagger}$Transformed using $\log(1+x)$.\par}
\end{table})---"

body <- c(block("input",   "Input metrics"),
          "\\addlinespace[3pt]",
          block("outcome", "Outcome metrics"))

writeLines(c(head_tex, body, foot_tex), "appendix_indicators.tex")
cat("Written:", file.path(getwd(), "appendix_indicators.tex"), "\n")

# =============================================================================
# 03_appendix_indicator_selection.R
# Appendix C: how 65 raw indicators become the K = 27 used in the models.
# Uses indicator_rates.rds (01_data_preparation_revised.R, section 6):
# one row per indicator after the 7 removals, with its training reporting rate.
# =============================================================================

library(dplyr); library(tibble)
setwd("~/Desktop/Thesis_ESG_data/Data")

rates <- readRDS("indicator_rates.rds")     # 58 indicators, report_rate = training years

# Same seven measures as the `removed` list in section 6, with the reason
removed <- tribble(
  ~indicator,                                       ~reason,
  "Number of Employees",                            "Measures firm size, not ESG performance",
  "Current CEO Chairman Duality",                   "Duplicate of CEO Chairman Duality",
  "CO2 Equivalent Emissions Total Market-based",    "Duplicate of CO2 Equivalent Emissions Total",
  "CO2 Equivalent Emissions Total Location-based",  "Duplicate of CO2 Equivalent Emissions Total",
  "Energy Use Total",                               "Intensity version available",
  "Training Hours Total",                           "Per-employee version available",
  "Water Pollutant Emissions",                      "Intensity version available")

dropped <- rates |> filter(report_rate <  0.20) |> arrange(desc(report_rate))
kept    <- rates |> filter(report_rate >= 0.20)

# The numbers for Section 4.2.1
c(raw = nrow(rates) + nrow(removed), after_removal = nrow(rates),
  dropped_by_threshold = nrow(dropped), kept = nrow(kept))   # 65, 58, 31, 27

# LaTeX output
esc <- function(x) gsub("([&%_#$])", "\\\\\\1", x)   # escape special characters

rows_A <- sprintf("%s & %s \\\\", esc(removed$indicator), removed$reason)
rows_B <- sprintf("%s & %.1f \\\\", esc(dropped$indicator), 100 * dropped$report_rate)

tex <- c(
  r"---(\begin{table}[htbp]
\centering
\footnotesize
\caption{Selection of ESG indicators}
\label{tab:selection}
\renewcommand{\arraystretch}{1.12}
\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}X >{\raggedright\arraybackslash}p{6cm}@{}}
\toprule
\multicolumn{2}{@{}l}{\textit{Panel A: Removed before selection}} \\
\midrule
\textbf{Indicator} & \textbf{Reason} \\
\midrule)---",
rows_A,
r"---(\midrule
\multicolumn{2}{@{}l}{\textit{Panel B: Excluded by the reporting threshold}} \\
\midrule
\textbf{Indicator} & \textbf{Reported in training years (\%)} \\
\midrule)---",
rows_B,
r"---(\bottomrule
\end{tabularx}

\vspace{4pt}
{\footnotesize\itshape\raggedright
Notes: The raw dataset contains 65 individual ESG indicators. The seven measures in Panel A are removed first, leaving 58. Of these, the indicators in Panel B are reported in less than 20\% of training firm-years (FY2015--FY2020) and are excluded. The remaining 27 indicators are listed in Appendix B.\par}
\end{table})---")

writeLines(tex, "appendix_indicator_selection.tex")
