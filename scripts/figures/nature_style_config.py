#!/usr/bin/env python3
"""
Shared Nature Methods Design System Configuration.

Enforces:
- Okabe-Ito Colorblind-Safe scientific palette
- Hairline wireframe borders (0.5 pt / 0.75 pt)
- Zero drop shadows, square corners (rx=0)
- Formal typography hierarchy (Arial/Helvetica for text, Times/Georgia for math)
- Standard Inkscape layer annotations
"""

# Canvas Dimensions (Double-column landscape, 180 mm journal width)
WIDTH = 1400
HEIGHT = 630

# Nature Methods Pure Background & Wireframe
COLOR_CANVAS_BG = "#FFFFFF"       # Pure white canvas
COLOR_PANEL_BG = "#FFFFFF"        # Pure white panel backgrounds
COLOR_CONTAINER_BG = "#FFFFFF"    # No colored pastel washes

# Hairlines and Dividers
COLOR_BORDER_HAIRLINE = "#CBD5E1" # 0.5 pt - 0.75 pt Slate 300
COLOR_DIVIDER_RULE = "#E2E8F0"    # Hairline dividing rule
COLOR_SUBTLE_FILL = "#F8FAFC"     # Ultra-subtle 1% background for inner equation boxes

# High-Contrast Publication Typography
COLOR_TEXT_PRIMARY = "#0F172A"    # Near black / Slate 900
COLOR_TEXT_SECONDARY = "#334155"  # Slate 700
COLOR_TEXT_MUTED = "#64748B"      # Slate 500
COLOR_TEXT_HAIRLINE = "#94A3B8"   # Slate 400

# Okabe-Ito Colorblind-Safe Palette (Strictly for data marks, curves, DAG nodes)
OKABE_BLACK = "#000000"
OKABE_ORANGE = "#E69F00"          # Priors & Likelihood parameters
OKABE_SKY_BLUE = "#56B4E9"        # Single-cell reference profiles
OKABE_BLUISH_GREEN = "#009E73"    # Fractions, Convergence, Positive biomarkers
OKABE_YELLOW = "#F0E442"
OKABE_BLUE = "#0072B2"            # Bulk observed data & primary lineage
OKABE_VERMILION = "#D55E00"       # Latent counts Z, Tumor programs, Negative biomarkers
OKABE_REDDISH_PURPLE = "#CC79A7"  # MCMC & Bayesian inference

# Font Family Definitions
FONT_SANS = "Arial, Helvetica, -apple-system, sans-serif"
FONT_SERIF_MATH = "'Times New Roman', Times, 'STIX Two Text', Georgia, serif"
