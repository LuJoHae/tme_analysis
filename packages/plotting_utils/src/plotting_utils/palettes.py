"""Nature Methods wireframe and Okabe-Ito Colorblind-Safe scientific palettes.
"""

from __future__ import annotations

from typing import Final

# Nature Methods Minimalist Wireframe Palette
COLOR_BG_WHITE: Final[str] = "#FFFFFF"
COLOR_CANVAS_BG: Final[str] = "#FFFFFF"
COLOR_CARD_BG: Final[str] = "#FFFFFF"
COLOR_PANEL_BG: Final[str] = "#FFFFFF"

COLOR_HAIRLINE: Final[str] = "#CBD5E1"
COLOR_BORDER_HAIRLINE: Final[str] = "#CBD5E1"
COLOR_DIVIDER_RULE: Final[str] = "#E2E8F0"
COLOR_SUBTLE_FILL: Final[str] = "#F8FAFC"

COLOR_TEXT_PRIMARY: Final[str] = "#0F172A"
COLOR_TEXT_SECONDARY: Final[str] = "#334155"
COLOR_TEXT_MUTED: Final[str] = "#64748B"
COLOR_TEXT_HAIRLINE: Final[str] = "#94A3B8"
COLOR_LIGHT_GREY: Final[str] = "#E2E8F0"
COLOR_NEUTRAL_GREY: Final[str] = "#94A3B8"

THRESHOLD_COLOR: Final[str] = "#DC2626"

# Publication Typography
FONT_SANS: Final[str] = (
    "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif"
)
FONT_SERIF_MATH: Final[str] = "'Times New Roman', Times, Georgia, serif"

# Okabe-Ito Colorblind-Safe Scientific Palette
OKABE_BLACK: Final[str] = "#000000"
OKABE_BLUE: Final[str] = "#0072B2"           # Enriched in Non-Responders
OKABE_VERMILION: Final[str] = "#D55E00"      # Enriched in Responders
OKABE_ORANGE: Final[str] = "#E69F00"         # Private Clonal Spike / Intermediate
OKABE_SKY_BLUE: Final[str] = "#56B4E9"
OKABE_BLUISH_GREEN: Final[str] = "#009E73"
OKABE_YELLOW: Final[str] = "#F0E442"
OKABE_REDDISH_PURPLE: Final[str] = "#CC79A7"
OKABE_GREY: Final[str] = "#94A3B8"

OKABE_PALETTE: tuple[str, ...] = (
    OKABE_BLUE,
    OKABE_ORANGE,
    OKABE_BLUISH_GREEN,
    OKABE_YELLOW,
    OKABE_SKY_BLUE,
    OKABE_VERMILION,
    OKABE_REDDISH_PURPLE,
    OKABE_GREY,
)

# 13-Fitter Scientific Palette (Okabe-Ito + Paul Tol Colorblind-Safe)
COLOR_FITTER_LASSO: Final[str] = "#0072B2"          # Dark Blue (Okabe-Ito)
COLOR_FITTER_ELASTIC_NET: Final[str] = "#56B4E9"     # Sky Blue (Okabe-Ito)
COLOR_FITTER_LOGISTIC: Final[str] = "#009E73"        # Bluish Green (Okabe-Ito)
COLOR_FITTER_RF: Final[str] = "#E69F00"              # Orange (Okabe-Ito)

COLOR_FITTER_COHORT_ADJ: Final[str] = "#CC79A7"      # Reddish Purple (Okabe-Ito)
COLOR_FITTER_GROUP_LASSO: Final[str] = "#D55E00"     # Vermilion (Okabe-Ito)
COLOR_FITTER_MERF: Final[str] = "#882255"            # Wine / Deep Purple (Paul Tol)
COLOR_FITTER_MULTITASK: Final[str] = "#117733"       # Forest Green (Paul Tol)
COLOR_FITTER_META: Final[str] = "#4363D8"            # Royal Blue
COLOR_FITTER_INVARIANT: Final[str] = "#999933"       # Olive Sand (Paul Tol)
COLOR_FITTER_GLMM: Final[str] = "#44AA99"            # Teal (Paul Tol)
COLOR_FITTER_OSCAR: Final[str] = "#88CCEE"           # Cyan / Sky Blue (Paul Tol)
COLOR_FITTER_SLOPE: Final[str] = "#332288"           # Indigo (Paul Tol)

FITTER_COLOR_MAP: Final[dict[str, str]] = {
    "lasso": COLOR_FITTER_LASSO,
    "elastic_net": COLOR_FITTER_ELASTIC_NET,
    "logistic": COLOR_FITTER_LOGISTIC,
    "rf": COLOR_FITTER_RF,
    "oscar": COLOR_FITTER_OSCAR,
    "slope": COLOR_FITTER_SLOPE,
    "cohort_adjusted": COLOR_FITTER_COHORT_ADJ,
    "group_lasso": COLOR_FITTER_GROUP_LASSO,
    "merf": COLOR_FITTER_MERF,
    "multitask_logistic": COLOR_FITTER_MULTITASK,
    "meta_analysis": COLOR_FITTER_META,
    "multistudy_invariant": COLOR_FITTER_INVARIANT,
    "glmm_lasso": COLOR_FITTER_GLMM,
}

FITTER_DISPLAY_NAMES: Final[dict[str, str]] = {
    "lasso": "Lasso",
    "elastic_net": "Elastic Net",
    "logistic": "Logistic",
    "rf": "Random Forest",
    "oscar": "OSCAR (Octagonal)",
    "slope": "SLOPE (Sorted L1)",
    "cohort_adjusted": "Cohort-Adjusted",
    "group_lasso": "Group Lasso",
    "merf": "MERF (Mixed RF)",
    "multitask_logistic": "Multi-Task Logistic",
    "meta_analysis": "Meta-Analysis",
    "multistudy_invariant": "Multi-Study Invariant",
    "glmm_lasso": "GLMM Lasso",
}

FITTER_BASELINE_ORDER: Final[tuple[str, ...]] = (
    "lasso",
    "elastic_net",
    "logistic",
    "rf",
    "oscar",
    "slope",
)

FITTER_MULTI_COHORT_ORDER: Final[tuple[str, ...]] = (
    "cohort_adjusted",
    "group_lasso",
    "merf",
    "multitask_logistic",
    "meta_analysis",
    "multistudy_invariant",
    "glmm_lasso",
)

# Deconvolution Method Scientific Palette (Okabe-Ito Colorblind-Safe)
DECONV_METHOD_ORDER: Final[tuple[str, ...]] = (
    "Unregularized (NNLS)",
    "RegDeconv (Graph Lap)",
    "Rectangle (DWLS-QP)",
    "CIBERSORT (reimpl., nu-SVR)",
    "CIBERSORTx (Docker)",
    "InstaPrism",
    "BayesPrism (Gibbs)",
)

DECONV_COLOR_MAP: Final[dict[str, str]] = {
    "Unregularized (NNLS)": OKABE_VERMILION,        # #D55E00
    "RegDeconv (Graph Lap)": OKABE_BLUE,            # #0072B2
    "Rectangle (DWLS-QP)": OKABE_REDDISH_PURPLE,    # #CC79A7
    "CIBERSORT (reimpl., nu-SVR)": OKABE_ORANGE,    # #E69F00
    "CIBERSORTx (Docker)": OKABE_BLACK,             # #000000
    "InstaPrism": OKABE_SKY_BLUE,                   # #56B4E9
    "BayesPrism (Gibbs)": OKABE_BLUISH_GREEN,       # #009E73
}
