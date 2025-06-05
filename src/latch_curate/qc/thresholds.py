from typing import TypedDict

class CountsThresholds(TypedDict):
    min_genes: int
    max_total_counts: int | None
    max_pct_mt: int

TECH_THRESHOLDS: dict[str, CountsThresholds] = {
    # droplet‑based
    "10x Chromium v2/v3 (3′/5′)":     {"min_genes": 200,  "max_total_counts": 25_000, "max_pct_mt": 15},
    "10x Chromium Flex (3′)":         {"min_genes": 100,  "max_total_counts": 15_000, "max_pct_mt": 20},
    "Drop‑seq / Seq‑Well":            {"min_genes": 200,  "max_total_counts": 20_000, "max_pct_mt": 20},
    "inDrops":                        {"min_genes": 200,  "max_total_counts": 15_000, "max_pct_mt": 20},
    "SPLiT‑seq / sci‑RNA‑seq":        {"min_genes": 100,  "max_total_counts": 10_000, "max_pct_mt": 10},

    # nucleus droplet
    "10x Single‑nucleus (snRNA‑seq)": {"min_genes": 100,  "max_total_counts": 10_000, "max_pct_mt":  5},

    # full‑length plate
    "Smart‑seq2":                     {"min_genes": 1_000,"max_total_counts": None,   "max_pct_mt": 30},
    "Smart‑seq3":                     {"min_genes": 800,  "max_total_counts": None,   "max_pct_mt": 30},
    "FLASH‑seq":                      {"min_genes": 800,  "max_total_counts": None,   "max_pct_mt": 30},
    "CEL‑Seq2":                       {"min_genes": 500,  "max_total_counts": 50_000, "max_pct_mt": 25},

    # spatial
    "10x Visium":                     {"min_genes": 200,  "max_total_counts": 50_000, "max_pct_mt": 20},
    "Slide‑seq v2":                   {"min_genes": 100,  "max_total_counts": 12_000, "max_pct_mt": 15},
}

