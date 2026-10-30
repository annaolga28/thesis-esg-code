# =============================================================================
# 01_data_preparation.R
# From the LSEG Excel export to the model-ready datasets (Phase 1).
#
# Input : Data_thesis.xlsx  (9 tabs, LSEG Workspace export of 20 March 2026)
# Output: esg_panel.rds        full firm-year panel incl. sector (all firms with ESG data)
#         esg_est.rds          estimation sample (common sample, before imputation)
#         sample_flow.rds      numbers for Figure fig:sample
#         appendix_sectors.csv sector table (Appendix tab:sectors)
#         indicator_rates.rds  reporting rate per indicator (train = selection; test = info)
#         ind_class.rds        indicator classification (Appendix indicator table)
#         data_nan.rds         transformed data, missing values left in (robustness R3)
#         data_main.rds        main imputation (median / p25 / p75)
#         data_r4a.rds         robustness: all values imputed at the median
#         data_r4b.rds         robustness: outcome metrics at p10 / p90
#         imputation_log.rds   where each imputed value came from
#
# (LSEG name <-> short name: ind_class$indicator <-> ind_class$var):
#   ids / outcome : RIC, year, split, Y
#   controls      : Size, Trend, PriorControversy, sector25
#   scores        : esg_score, e_pillar, s_pillar, g_pillar
#   indicators    : ind_class$var            dummies: M_<var>
#   industry_group: TRBC Industry Group, used for imputation only (not a predictor)
# =============================================================================

library(readxl); library(dplyr); library(tidyr); library(purrr); library(tibble)

setwd("~/Desktop/Thesis_ESG_data/Data")
path <- "Data_thesis.xlsx"


# 1. Read all tabs of the Excel export
# -----------------------------------------------------------------------------
# Each tab has two header rows: row 1 = metric name (merged across its year
# columns), row 2 = fiscal year (FY0, FY-1, ... FY-9). Columns are renamed to
# "metric__FYx" so metric and year can be separated later. The first four
# columns are the firm identifiers.
# The metric name is cut at the first line break (\r or \n) and trimmed
# BEFORE "__FYx" is added. Otherwise a hidden \r stays inside the name and
# columns such as `ESG Controversies Score` cannot be found later.

read_one <- function(sheet) {
  x <- read_excel(path, sheet = sheet, col_names = FALSE, col_types = "text",
                  .name_repair = "minimal")
  metric <- tibble(m = trimws(sub("[\r\n].*", "", unlist(x[1, ])))) |>
    fill(m) |> pull(m)
  fy     <- trimws(unlist(x[2, ]))
  nm     <- ifelse(is.na(fy), metric, paste(metric, fy, sep = "__"))
  nm[1:4] <- c("RIC", "Company", "Country", "MktCap_EURm")
  dat <- x[-(1:2), ]
  names(dat) <- nm
  dup <- unique(nm[duplicated(nm)])
  if (length(dup)) message(sheet, ": duplicate columns, first kept: ",
                           paste(dup, collapse = " | "))
  dat[, !duplicated(names(dat))] |> filter(!is.na(RIC))
}

sheets  <- excel_sheets(path)
esg_raw <- map(sheets, read_one) |> set_names(sheets)
map_int(esg_raw, nrow)                           # check data: 4,064 firms per tab



# 2. Build the panel: one row per firm per fiscal year
# -----------------------------------------------------------------------------
# Firm-level fields (industry, market cap at FY0) come from the main tab.
# Yearly metrics from all tabs are stacked, FY labels turned into calendar
# years (FY0 = 2024 ... FY-9 = 2015), then spread to one column per metric.
# FY0 = FY2024 for every firm (confirmed in LSEG, 18 Oct 2026).

static <- esg_raw[["ESG Europe 50mln"]] |>
  select(RIC, Company, Country, MktCap_EURm, !contains("__")) |>
  select(-starts_with("ESG Score Grade")) |>     # letter grades are not used
  mutate(MktCap_EURm = as.numeric(MktCap_EURm))  # stored as text in the export

long <- map_dfr(esg_raw, \(d) d |>
                  select(RIC, contains("__FY")) |>
                  pivot_longer(-RIC, names_to = c("metric", "fy"),
                               names_sep = "__", values_to = "value")) |>
  filter(!is.na(value)) |>
  mutate(year = 2024 + as.integer(sub("FY", "", fy))) |>
  select(-fy)

# Metrics repeated across tabs must carry the same value in every tab
long |> distinct(RIC, year, metric, value) |>
  count(RIC, year, metric) |> filter(n > 1) |> nrow()   # check: 0

# type.convert() only recognises TRUE/FALSE in capitals; "True"/"true" would
# stay text (and as.integer() would later turn them into NA). Upper-case first.
upper_tf <- function(x) {
  tf <- tolower(x) %in% c("true", "false")
  x[tf] <- toupper(x[tf]); x
}

panel <- long |>
  distinct(RIC, year, metric, .keep_all = TRUE) |>
  pivot_wider(names_from = metric, values_from = value) |>
  mutate(across(-c(RIC, year), \(x) type.convert(upper_tf(x), as.is = TRUE))) |>
  select(-starts_with("ESG Score Grade")) |>
  left_join(static, by = "RIC") |>
  arrange(RIC, year)

stopifnot(!any(grepl("[\r\n]", names(panel))))   # no hidden line breaks left
range(panel$year)                                # check: 2015 2024
dim(panel); n_distinct(panel$RIC)                # check: 16,125 rows, 2,214 firms


# 3. Next year's controversy score (lead of one year, within firm)
# -----------------------------------------------------------------------------
# cont_t1 is only filled if the firm's next row is exactly year + 1;
# gaps are not bridged. FY2024 rows have no next year and get NA.

panel <- panel |>
  group_by(RIC) |>
  arrange(year, .by_group = TRUE) |>
  mutate(cont_t1 = if_else(lead(year) == year + 1,
                           as.numeric(lead(`ESG Controversies Score`)), NA_real_)) |>
  ungroup()


# 4. Industry sector (25 groups, each >= 20 firms)  -> Ch. 4.3
# -----------------------------------------------------------------------------
# The TRBC Industry Groups are mapped to their TRBC Business Sector.
# Sectors with fewer than 20 firms are merged with the closest related sector
# (marked "merged"). Uranium and Schools are in the export but have no firm
# in the estimation sample. Industry is time-invariant (current TRBC group).
# Used for the industry controls and as the middle level of the imputation.

sector_map <- tribble(
  ~group,                                              ~sector25,
  "Machinery, Tools, Heavy Vehicles, Trains & Ships",  "Industrial Goods",
  "Aerospace & Defense",                               "Industrial Goods",
  "Professional & Commercial Services",                "Industrial & Commercial Services",
  "Construction & Engineering",                        "Industrial & Commercial Services",
  "Diversified Industrial Goods Wholesale",            "Industrial & Commercial Services",
  "Professional & Business Education",                 "Industrial & Commercial Services",  # merged
  "Miscellaneous Educational Service Providers",       "Industrial & Commercial Services",  # merged
  "Schools, Colleges & Universities",                  "Industrial & Commercial Services",  # merged
  "Freight & Logistics Services",                      "Transportation",
  "Passenger Transportation Services",                 "Transportation",
  "Transport Infrastructure",                          "Transportation",
  "Software & IT Services",                            "Software & IT Services",
  "Financial Technology (Fintech) & Infrastructure",   "Software & IT Services",            # merged
  "Semiconductors & Semiconductor Equipment",          "Technology Equipment",
  "Electronic Equipment & Parts",                      "Technology Equipment",
  "Communications & Networking",                       "Technology Equipment",
  "Computers, Phones & Household Electronics",         "Technology Equipment",
  "Office Equipment",                                  "Technology Equipment",
  "Integrated Hardware & Software",                    "Technology Equipment",
  "Telecommunications Services",                       "Telecommunications Services",
  "Banking Services",                                  "Banking & Investment Services",
  "Investment Banking & Investment Services",          "Banking & Investment Services",
  "Insurance",                                         "Insurance",
  "Collective Investments",                            "Collective Investments",
  "Investment Holding Companies",                      "Investment Holding Companies",
  "Real Estate Operations",                            "Real Estate",
  "Residential & Commercial REITs",                    "Real Estate",
  "Metals & Mining",                                   "Mineral Resources",
  "Construction Materials",                            "Mineral Resources",
  "Chemicals",                                         "Chemicals",
  "Paper & Forest Products",                           "Applied Resources",
  "Containers & Packaging",                            "Applied Resources",
  "Food & Tobacco",                                    "Food & Beverages",
  "Beverages",                                         "Food & Beverages",
  "Food & Drug Retailing",                             "Food & Drug Retailing",
  "Personal & Household Products & Services",          "Personal & Household Products",     # merged
  "Consumer Goods Conglomerates",                      "Personal & Household Products",     # merged
  "Pharmaceuticals",                                   "Pharmaceuticals & Medical Research",
  "Biotechnology & Medical Research",                  "Pharmaceuticals & Medical Research",
  "Healthcare Equipment & Supplies",                   "Healthcare Services & Equipment",
  "Healthcare Providers & Services",                   "Healthcare Services & Equipment",
  "Hotels & Entertainment Services",                   "Cyclical Consumer Services",
  "Media & Publishing",                                "Cyclical Consumer Services",
  "Automobiles & Auto Parts",                          "Automobiles & Auto Parts",
  "Specialty Retailers",                               "Retailers",
  "Diversified Retail",                                "Retailers",
  "Textiles & Apparel",                                "Cyclical Consumer Products",
  "Homebuilding & Construction Supplies",              "Cyclical Consumer Products",
  "Household Goods",                                   "Cyclical Consumer Products",
  "Leisure Products",                                  "Cyclical Consumer Products",
  "Electric Utilities & IPPs",                         "Utilities",
  "Multiline Utilities",                               "Utilities",
  "Natural Gas Utilities",                             "Utilities",
  "Water & Related Utilities",                         "Utilities",
  "Oil & Gas",                                         "Energy",
  "Oil & Gas Related Equipment and Services",          "Energy",
  "Coal",                                              "Energy",
  "Renewable Energy",                                  "Energy",                            # merged
  "Uranium",                                           "Energy"                             # merged
)

nrow(sector_map); n_distinct(sector_map$sector25)   # check: 59 groups, 25 sectors
stopifnot(!anyDuplicated(sector_map$group))

panel <- panel |> left_join(sector_map, by = c("TRBC Industry Group Name" = "group"))
stopifnot(!anyNA(panel$sector25))                # every group matched

saveRDS(panel, "esg_panel.rds")                  # saved AFTER the sector join



# 5. Estimation sample: outcome, controls, common sample, train/test split
# -----------------------------------------------------------------------------
# Y                = 1 if a controversy is recorded in t+1 (score < 100)
# PriorControversy = 1 if a controversy is recorded in t
# Common sample    = ESG Score and all three pillar scores available in t
# Split            = train FY2015-2020, test FY2021-2023 (predictor year t)

# Intermediate step for the sample flow: ESG Score in t and outcome in t+1
with_t1 <- panel |> filter(!is.na(`ESG Score`), !is.na(cont_t1))

est <- panel |>
  mutate(Y                = as.integer(cont_t1 < 100),
         PriorControversy = as.integer(`ESG Controversies Score` < 100),
         Trend            = year - 2015,
         Size             = log(MktCap_EURm),
         split            = if_else(year <= 2020, "train", "test")) |>
  filter(!is.na(`ESG Score`), !is.na(`Environmental Pillar Score`),
         !is.na(`Social Pillar Score`), !is.na(`Governance Pillar Score`),
         !is.na(Y))

# Checks
table(est$Y, useNA = "ifany")                    # only 0/1; 1s total 2,146
range(est$Trend)                                 # 0 8
min(est$MktCap_EURm)                             # > 50 (sample filter)
if (anyNA(est$PriorControversy) || anyNA(est$Size))
  stop("NA in PriorControversy or Size: decide how to treat these rows")

# The observation(s) lost between with_t1 and the common sample (footnote 4.1)
with_t1 |> anti_join(est, by = c("RIC", "year")) |>
  select(RIC, year, `ESG Score`, `Environmental Pillar Score`,
         `Social Pillar Score`, `Governance Pillar Score`)

# Ch. 4.1: sample flow (Figure fig:sample)
sample_flow <- tribble(
  ~step,                                      ~firms,                   ~firm_years,
  "LSEG export (market cap > EUR 50m)",       nrow(static),             NA_integer_,
  "With ESG data, FY2015-FY2024",             n_distinct(panel$RIC),    nrow(panel),
  "ESG Score in t and outcome in t+1",        n_distinct(with_t1$RIC),  nrow(with_t1),
  "Common sample (ESG Score + 3 pillars)",    n_distinct(est$RIC),      nrow(est))
sample_flow                     # expect 4,064 / 2,214 + 16,125 / 2,174 + 13,866 / 2,174 + 13,865
sum(est$Y); round(mean(est$Y), 3)                # 2,146 / 0.155

# Ch. 4.1: firm-years per predictor year (was 849-2,160 on the 13,866 sample)
est |> count(year, name = "firm_years") |>
  left_join(est |> group_by(year) |> summarise(Y1 = sum(Y)), by = "year")

# Ch. 4.5: size and class balance per period
est |> group_by(split) |>
  summarise(firm_years = n(), firms = n_distinct(RIC),
            Y1 = sum(Y), share_Y1 = round(mean(Y), 3))

# Sanity check (size-leakage limitation): Size at FY0 vs. PriorControversy and Y
round(c(cor_Size_Prior = cor(est$Size, est$PriorControversy),
        cor_Size_Y     = cor(est$Size, est$Y)), 3)

saveRDS(est, "esg_est.rds")
saveRDS(sample_flow, "sample_flow.rds")

# Ch. 4.3 / Appendix tab:sectors, built on the final common sample
n_distinct(est$`TRBC Industry Group Name`)                           # 57 groups
n_distinct(substr(as.character(est$`TRBC Industry Code`), 1, 4))     # 29 business sectors
n_distinct(est$sector25)                                             # 25

appendix_sectors <- est |>
  distinct(RIC, sector25, `TRBC Industry Group Name`) |>
  group_by(Sector = sector25) |>
  summarise(`Industry groups included` =
              paste(sort(unique(`TRBC Industry Group Name`)), collapse = "; "),
            Firms = n_distinct(RIC)) |>
  arrange(desc(Firms))

print(appendix_sectors, n = Inf)                 # smallest = 20 firms
stopifnot(sum(appendix_sectors$Firms) == n_distinct(est$RIC))   # sums to 2,174
write.csv(appendix_sectors, "appendix_sectors.csv", row.names = FALSE)


# 6. Indicator selection (training period only)  -> Chapter 4.2.1
# -----------------------------------------------------------------------------
# Raw indicators = all columns that are not identifiers, controls, LSEG
# derived scores, letter grades or controversy variables (65 columns).
# Seven are removed: duplicates, absolute measures for which an intensity
# version exists, and one size measure (58 remain). An indicator is kept if
# it is reported for >= 20% of TRAINING firm-years (K = 27).
# The test-period rate is stored for information only (disclosure trend).

not_indicators <- c("RIC", "year", "Company", "Country", "MktCap_EURm", "IPO Date",
                    "Company Market Capitalization", "Revenue 5 YR Hist Growth",
                    "Revenue 10 YR Hist Growth", "TRBC Industry Name",
                    "TRBC Industry Group Name", "TRBC Industry Code", "sector25",
                    "cont_t1", "Y", "PriorControversy", "Trend", "Size", "split")

removed <- c("Number of Employees",                            # size measure
             "Current CEO Chairman Duality",                   # duplicate
             "CO2 Equivalent Emissions Total Market-based",    # variant of CO2 total
             "CO2 Equivalent Emissions Total Location-based",  # variant of CO2 total
             "Energy Use Total",                               # intensity version exists
             "Training Hours Total",                           # per-employee version exists
             "Water Pollutant Emissions")                      # intensity version exists

candidates <- setdiff(names(est), not_indicators)
candidates <- candidates[!grepl("Score|Grade|Controvers", candidates)]
stopifnot(all(removed %in% candidates))          # no typos in `removed`
indicators <- setdiff(candidates, removed)
length(candidates); length(indicators)           # check data: 65 and 58

# every indicator must be numeric or logical after type.convert()
non_numeric <- names(Filter(\(x) !(is.numeric(x) || is.logical(x)), est[indicators]))
if (length(non_numeric)) stop("Not numeric/logical: ", paste(non_numeric, collapse = " | "))

rates <- est |>
  group_by(split) |>
  summarise(across(all_of(indicators), \(x) mean(!is.na(x)))) |>
  pivot_longer(-split, names_to = "indicator", values_to = "rate") |>
  pivot_wider(names_from = split, values_from = rate, names_prefix = "rate_") |>
  rename(report_rate = rate_train, report_rate_test = rate_test) |>
  arrange(desc(report_rate))

kept <- rates |> filter(report_rate >= 0.20) |> pull(indicator)
length(kept)                                     # K = 27
saveRDS(rates, "indicator_rates.rds")



# 7. Indicator classification (Appendix indicator table)
# -----------------------------------------------------------------------------
# type      input (policy, commitment, action) / outcome (realised result)
# direction higher value = better / worse / neutral (used for H4 polarity)
# unit      binary, percent, rate, count, intensity, absolute (-> log1p)
# var       short syntactic name used in the model datasets
# The imputation rule follows from type and direction. A missingness dummy
# M_k is created only if the indicator is reported for < 99% of training
# firm-years (otherwise the dummy is (almost) always 0).

clean_name <- function(x) gsub("^_|_$", "", tolower(gsub("[^A-Za-z0-9]+", "_", x)))

ind_class <- tribble(
  ~indicator,                                    ~type,     ~direction, ~unit,
  "Policy Board Size",                           "input",   "better",   "binary",
  "CEO Chairman Duality",                        "input",   "worse",    "binary",
  "Policy Employee Health & Safety",             "input",   "better",   "binary",
  "Policy Board Independence",                   "input",   "better",   "binary",
  "Board Gender Diversity, Percent",             "outcome", "better",   "percent",
  "Board Size",                                  "input",   "neutral",  "count",
  "Policy Business Ethics",                      "input",   "better",   "binary",
  "Policy Environmental Supply Chain",           "input",   "better",   "binary",
  "Environmental Supply Chain Management",       "input",   "better",   "binary",
  "Policy Human Rights",                         "input",   "better",   "binary",
  "Policy Child Labor",                          "input",   "better",   "binary",
  "Policy Water Efficiency",                     "input",   "better",   "binary",
  "Environmental Investments Initiatives",       "input",   "better",   "binary",
  "Environmental Restoration Initiatives",       "input",   "better",   "binary",
  "Health & Safety Policy",                      "input",   "better",   "binary",
  "Environmental Expenditures Investments",      "input",   "better",   "binary",
  "Women Employees",                             "outcome", "better",   "percent",
  "CO2 Equivalent Emissions Total",              "outcome", "worse",    "absolute",
  "Total Energy Use To Revenues USD in million", "outcome", "worse",    "intensity",
  "Turnover of Employees",                       "outcome", "worse",    "percent",
  "Water Use To Revenues USD in million",        "outcome", "worse",    "intensity",
  "Average Training Hours",                      "input",   "better",   "intensity",
  "Employee Fatalities",                         "outcome", "worse",    "absolute",
  "Total Injury Rate Total",                     "outcome", "worse",    "rate",
  "Employee Accidents",                          "outcome", "worse",    "absolute",
  "Total Injury Rate Employees",                 "outcome", "worse",    "rate",
  "Renewable Energy Use Ratio",                  "outcome", "better",   "percent"
) |>
  left_join(rates, by = "indicator") |>
  mutate(var        = clean_name(indicator),
         imputation = case_when(type == "input" | direction == "neutral" ~ "median",
                                direction == "better" ~ "p25",
                                direction == "worse"  ~ "p75"),
         make_M     = report_rate < 0.99)

stopifnot(setequal(ind_class$indicator, kept),   # all 27, no typos
          !anyDuplicated(ind_class$var))
count(ind_class, type)                           # 16 input / 11 outcome
count(filter(ind_class, make_M), type)           # dummies by type (H5b comparison)
saveRDS(ind_class, "ind_class.rds")


# 8. Model columns, transformation, missingness dummies  -> Ch. 4.3.3
# -----------------------------------------------------------------------------
# Only the modelling columns are kept, so nothing else (cont_t1, controversy
# counts, unused indicators) can slip into a model by accident. Indicators get
# their short names from ind_class$var.
# Order matters: dummies are created BEFORE imputation, otherwise the
# information on which values were unreported is lost.

keep <- c(RIC = "RIC", year = "year", split = "split", Y = "Y",
          Size = "Size", Trend = "Trend", PriorControversy = "PriorControversy",
          sector25 = "sector25", industry_group = "TRBC Industry Group Name",
          esg_score = "ESG Score",            e_pillar = "Environmental Pillar Score",
          s_pillar  = "Social Pillar Score",  g_pillar = "Governance Pillar Score",
          setNames(ind_class$indicator, ind_class$var))

bin_vars <- ind_class |> filter(unit == "binary")   |> pull(var)
abs_vars <- ind_class |> filter(unit == "absolute") |> pull(var)
m_vars   <- ind_class |> filter(make_M)             |> pull(var)

dat <- est |>
  select(all_of(keep)) |>
  mutate(across(all_of(bin_vars), as.integer),   # TRUE/FALSE -> 1/0
         across(all_of(abs_vars), log1p))        # absolute amounts -> log(1 + x)

stopifnot(all(unlist(dat[bin_vars]) %in% c(0L, 1L, NA)))   # binaries really 0/1

for (v in m_vars) dat[[paste0("M_", v)]] <- as.integer(is.na(dat[[v]]))
length(m_vars)                                   # number of M_k dummies (report in 4.3.3)

saveRDS(dat, "data_nan.rds")                     # missing values left in (robustness R3)


# 9. Imputation (three versions)  -> Chapter 4.2.2
# -----------------------------------------------------------------------------
# A missing value gets a percentile of the values REPORTED by other firms in
# the same TRBC Industry Group and year. If fewer than 5 firms reported, the
# 25-sector group of that year is used, then all firms of that year.
# Same-year statistics are used in test years too (cross-sectional values
# available at t, no look-ahead). Binary indicators always receive an
# observed value (0 or 1).
#
# Percentile per scheme:
#   main : input/neutral = median, outcome better = p25, outcome worse = p75
#   r4a  : everything at the median
#   r4b  : input/neutral = median, outcome better = p10, outcome worse = p90

pick_p <- function(type, direction, scheme) {
  if (scheme == "r4a" || type == "input" || direction == "neutral") return(0.50)
  lo <- if (scheme == "main") 0.25 else 0.10
  if (direction == "better") lo else 1 - lo
}

impute_one <- function(d, v, p, min_n = 5) {
  qtype <- if (v %in% bin_vars) 1 else 7
  q  <- function(x) quantile(x, p, na.rm = TRUE, type = qtype, names = FALSE)
  qn <- function(x) if (sum(!is.na(x)) >= min_n) q(x) else NA_real_
  d[[v]] <- as.numeric(d[[v]])                   # same type for reported and imputed
  d |>
    group_by(industry_group, year) |> mutate(.g = qn(.data[[v]])) |>
    group_by(sector25, year)       |> mutate(.s = qn(.data[[v]])) |>
    group_by(year)                 |> mutate(.y = q(.data[[v]]))  |>
    ungroup() |>
    mutate(.lvl = case_when(!is.na(.data[[v]]) ~ "reported",
                            !is.na(.g) ~ "group", !is.na(.s) ~ "sector", TRUE ~ "year"),
           !!v := coalesce(.data[[v]], .g, .s, .y))
}

make_imputed <- function(scheme) {
  d <- dat; log <- list()
  for (i in seq_len(nrow(ind_class))) {
    r <- ind_class[i, ]
    d <- impute_one(d, r$var, pick_p(r$type, r$direction, scheme))
    log[[r$var]] <- count(d, .lvl)
    d <- select(d, -.g, -.s, -.y, -.lvl)
  }
  list(data = d, log = bind_rows(log, .id = "var"))
}

main <- make_imputed("main")
r4a  <- make_imputed("r4a")
r4b  <- make_imputed("r4b")



# 10. Checks
# -----------------------------------------------------------------------------

imputed <- list(main = main$data, r4a = r4a$data, r4b = r4b$data)

# No missing values may remain in the imputed datasets (all FALSE)
sapply(imputed, \(d) anyNA(d[ind_class$var]))

# Reported values and dummies are untouched by the imputation (all TRUE)
sapply(imputed, \(d) all(map_lgl(ind_class$var, \(v) {
  rep <- !is.na(dat[[v]]); isTRUE(all.equal(d[[v]][rep], dat[[v]][rep]))
})) && isTRUE(all.equal(as.data.frame(d[paste0("M_", m_vars)]),
                         as.data.frame(dat[paste0("M_", m_vars)]))))

# Chapter 4.2.2: where imputed values came from (78.2% group, 10.4% sector, 11.4% year)
main$log |>
  filter(.lvl != "reported") |>
  group_by(.lvl) |> summarise(n = sum(n)) |>
  mutate(share = round(n / sum(n), 3))

saveRDS(main$data, "data_main.rds")
saveRDS(r4a$data,  "data_r4a.rds")
saveRDS(r4b$data,  "data_r4b.rds")
saveRDS(main$log,  "imputation_log.rds")
