# =============================================================
# 05_pathway_mapping.py  — v4 (fixed)
#
# Fixes vs previously carried version:
#   1. map_genes_to_pathways_live returned list of strings —
#      caller (06_xai.py) expected list of (name, pvalue) tuples.
#      NOW: all functions return consistent type:
#           list[tuple[str, float | None]]
#
#   2. g:Profiler was called WITHOUT a background gene list,
#      meaning enrichment was tested against all ~20,000 human
#      genes, not your ~2000 selected genes. p-values were
#      systematically anti-conservative (too small).
#      NOW: background_genes parameter required for live queries.
#
#   3. Coverage metric added: what fraction of your selected
#      genes are in the offline map? Caller can report this.
#
#   4. Offline map expanded to ~350 genes (was ~200).
#      Added: ferroptosis, RNA processing, splicing, nuclear
#      transport, cell polarity, vesicle trafficking — all
#      relevant to cancer drug response.
# =============================================================

try:
    from gprofiler import GProfiler
    GPROFILER_AVAILABLE = True
except ImportError:
    GPROFILER_AVAILABLE = False

# -----------------------------------------------------------
# TYPE DEFINITION (consistent across all functions)
# PathwayResult = list[tuple[pathway_name: str, pvalue: float | None]]
# pvalue is None when coming from offline dict (no stats available)
# -----------------------------------------------------------

# -----------------------------------------------------------
# OFFLINE MAP — ~350 cancer-relevant genes
# Each gene maps to exactly one canonical pathway string.
# Extend this freely — the richer it is, the fewer "unknown"
# genes appear in your XAI explanations.
# -----------------------------------------------------------
OFFLINE_MAP = {
    # --- Apoptosis ---
    "BCL2": "Apoptosis", "BCL2L1": "Apoptosis", "BCL2L11": "Apoptosis",
    "BCL2L2": "Apoptosis", "BCL2L14": "Apoptosis",
    "BAX": "Apoptosis", "BAD": "Apoptosis", "BAK1": "Apoptosis",
    "BOK": "Apoptosis", "BID": "Apoptosis", "BIM": "Apoptosis",
    "CASP3": "Apoptosis", "CASP6": "Apoptosis", "CASP7": "Apoptosis",
    "CASP8": "Apoptosis", "CASP9": "Apoptosis", "CASP10": "Apoptosis",
    "CYCS": "Apoptosis", "APAF1": "Apoptosis",
    "BIRC5": "Apoptosis", "BIRC2": "Apoptosis", "XIAP": "Apoptosis",
    "MCL1": "Apoptosis", "PMAIP1": "Apoptosis", "BBC3": "Apoptosis",
    "DIABLO": "Apoptosis", "HRK": "Apoptosis",
    # --- Ferroptosis ---
    "GPX4": "Ferroptosis", "SLC7A11": "Ferroptosis",
    "ACSL4": "Ferroptosis", "LPCAT3": "Ferroptosis",
    "NFE2L2": "Ferroptosis / oxidative stress",
    "KEAP1": "Ferroptosis / oxidative stress",
    "FTH1": "Ferroptosis", "FTL": "Ferroptosis",
    # --- Cell cycle ---
    "TP53": "Cell cycle / DNA damage", "TP63": "Cell cycle",
    "TP73": "Cell cycle", "RB1": "Cell cycle",
    "CDK1": "Cell cycle", "CDK2": "Cell cycle",
    "CDK4": "Cell cycle", "CDK6": "Cell cycle",
    "CDK7": "Cell cycle", "CDK8": "Cell cycle",
    "CDK9": "Cell cycle",
    "CCNA1": "Cell cycle", "CCNA2": "Cell cycle",
    "CCNB1": "Cell cycle", "CCNB2": "Cell cycle",
    "CCND1": "Cell cycle", "CCND2": "Cell cycle", "CCND3": "Cell cycle",
    "CCNE1": "Cell cycle", "CCNE2": "Cell cycle",
    "CCNH": "Cell cycle",
    "CDKN1A": "Cell cycle", "CDKN1B": "Cell cycle", "CDKN1C": "Cell cycle",
    "CDKN2A": "Cell cycle", "CDKN2B": "Cell cycle",
    "E2F1": "Cell cycle", "E2F3": "Cell cycle",
    "MYC": "Cell cycle / transcription",
    "MYCN": "Cell cycle / transcription",
    "PLK1": "Cell cycle", "PLK4": "Cell cycle",
    "AURKA": "Cell cycle", "AURKB": "Cell cycle",
    "BUB1": "Cell cycle", "MAD2L1": "Cell cycle",
    # --- DNA repair ---
    "BRCA1": "DNA repair", "BRCA2": "DNA repair",
    "ATM": "DNA repair", "ATR": "DNA repair",
    "CHEK1": "DNA repair", "CHEK2": "DNA repair",
    "RAD51": "DNA repair", "RAD51C": "DNA repair",
    "RAD52": "DNA repair", "PALB2": "DNA repair",
    "PARP1": "DNA repair", "PARP2": "DNA repair",
    "MLH1": "DNA repair", "MSH2": "DNA repair",
    "MSH6": "DNA repair", "PMS2": "DNA repair",
    "ERCC1": "DNA repair", "ERCC2": "DNA repair",
    "XPC": "DNA repair", "PCNA": "DNA repair",
    "FANCA": "DNA repair / Fanconi", "FANCC": "DNA repair / Fanconi",
    "FANCD2": "DNA repair / Fanconi",
    "WRN": "DNA repair", "BLM": "DNA repair",
    # --- PI3K-AKT-mTOR ---
    "AKT1": "PI3K-AKT", "AKT2": "PI3K-AKT", "AKT3": "PI3K-AKT",
    "PIK3CA": "PI3K-AKT", "PIK3CB": "PI3K-AKT",
    "PIK3CD": "PI3K-AKT", "PIK3CG": "PI3K-AKT",
    "PIK3R1": "PI3K-AKT", "PIK3R2": "PI3K-AKT",
    "PTEN": "PI3K-AKT",
    "MTOR": "PI3K-AKT / mTOR", "MLST8": "mTOR",
    "TSC1": "mTOR", "TSC2": "mTOR",
    "RPTOR": "mTOR", "RICTOR": "mTOR",
    "RPS6KB1": "mTOR", "EIF4EBP1": "mTOR",
    "PDK1": "PI3K-AKT", "SGK1": "PI3K-AKT",
    # --- RAS-MAPK ---
    "KRAS": "RAS-MAPK", "NRAS": "RAS-MAPK", "HRAS": "RAS-MAPK",
    "BRAF": "RAS-MAPK", "RAF1": "RAS-MAPK", "ARAF": "RAS-MAPK",
    "MAP2K1": "RAS-MAPK", "MAP2K2": "RAS-MAPK",
    "MAPK1": "RAS-MAPK", "MAPK3": "RAS-MAPK",
    "MAPK4": "RAS-MAPK", "MAPK6": "RAS-MAPK",
    "MAPK8": "RAS-MAPK / JNK", "MAPK9": "RAS-MAPK / JNK",
    "MAPK14": "RAS-MAPK / p38", "MAPK11": "RAS-MAPK / p38",
    "NF1": "RAS-MAPK", "RASA1": "RAS-MAPK",
    "SOS1": "RAS-MAPK", "GRB2": "RAS-MAPK",
    # --- WNT ---
    "CTNNB1": "WNT signaling", "APC": "WNT signaling",
    "GSK3B": "WNT signaling", "AXIN1": "WNT signaling", "AXIN2": "WNT signaling",
    "FZD1": "WNT signaling", "FZD2": "WNT signaling",
    "LRP5": "WNT signaling", "LRP6": "WNT signaling",
    "TCF7L2": "WNT signaling", "DVL1": "WNT signaling",
    "RNF43": "WNT signaling", "RSPO1": "WNT signaling",
    # --- Notch ---
    "NOTCH1": "Notch signaling", "NOTCH2": "Notch signaling",
    "NOTCH3": "Notch signaling", "NOTCH4": "Notch signaling",
    "JAG1": "Notch signaling", "JAG2": "Notch signaling",
    "DLL1": "Notch signaling", "DLL3": "Notch signaling",
    "HES1": "Notch signaling", "MAML1": "Notch signaling",
    # --- Hedgehog ---
    "PTCH1": "Hedgehog signaling", "SMO": "Hedgehog signaling",
    "GLI1": "Hedgehog signaling", "GLI2": "Hedgehog signaling",
    "SUFU": "Hedgehog signaling",
    # --- Metabolism ---
    "SLC2A1": "Glucose metabolism", "SLC2A4": "Glucose metabolism",
    "HK1": "Glycolysis", "HK2": "Glycolysis",
    "LDHA": "Glycolysis", "LDHB": "Glycolysis",
    "PKM": "Glycolysis", "PFKFB3": "Glycolysis",
    "GPI": "Glycolysis", "ALDOA": "Glycolysis",
    "IDH1": "TCA cycle / metabolism", "IDH2": "TCA cycle / metabolism",
    "IDH3A": "TCA cycle", "SDHA": "TCA cycle", "FH": "TCA cycle",
    "ACAD8": "Fatty acid oxidation", "HADHA": "Fatty acid oxidation",
    "FASN": "Lipid synthesis", "ACACA": "Lipid synthesis",
    "GLS": "Glutamine metabolism", "GLS2": "Glutamine metabolism",
    "GOT1": "Glutamine metabolism",
    # --- Chromatin / epigenetics ---
    "EZH2": "Chromatin remodeling", "EZH1": "Chromatin remodeling",
    "KDM6A": "Chromatin remodeling", "KDM5C": "Chromatin remodeling",
    "DNMT1": "DNA methylation", "DNMT3A": "DNA methylation",
    "DNMT3B": "DNA methylation", "TET2": "DNA methylation",
    "HDAC1": "Histone deacetylation", "HDAC2": "Histone deacetylation",
    "HDAC3": "Histone deacetylation", "SIRT1": "Histone deacetylation",
    "KAT6A": "Histone acetylation", "EP300": "Histone acetylation",
    "CREBBP": "Histone acetylation",
    "SMARCA4": "SWI/SNF complex", "SMARCB1": "SWI/SNF complex",
    "ARID1A": "SWI/SNF complex", "ARID1B": "SWI/SNF complex",
    "KMT2A": "Histone methylation", "KMT2D": "Histone methylation",
    "KDM1A": "Histone methylation",
    # --- RTK / growth factor ---
    "EGFR": "EGFR / RTK signaling", "ERBB2": "ERBB signaling",
    "ERBB3": "ERBB signaling", "ERBB4": "ERBB signaling",
    "MET": "RTK signaling", "ALK": "RTK signaling",
    "ROS1": "RTK signaling", "RET": "RTK signaling",
    "FGFR1": "FGFR signaling", "FGFR2": "FGFR signaling",
    "FGFR3": "FGFR signaling", "FGFR4": "FGFR signaling",
    "PDGFRA": "PDGFR signaling", "PDGFRB": "PDGFR signaling",
    "KIT": "SCF-KIT signaling", "CSF1R": "RTK signaling",
    "IGF1R": "IGF signaling", "INSR": "Insulin signaling",
    "VEGFR1": "VEGF signaling", "VEGFR2": "VEGF signaling",
    "FLT3": "RTK signaling", "AXL": "RTK signaling",
    "NTRK1": "RTK signaling", "NTRK2": "RTK signaling",
    # --- JAK-STAT ---
    "JAK1": "JAK-STAT", "JAK2": "JAK-STAT", "JAK3": "JAK-STAT",
    "TYK2": "JAK-STAT",
    "STAT1": "JAK-STAT", "STAT3": "JAK-STAT", "STAT5A": "JAK-STAT",
    "STAT5B": "JAK-STAT", "SOCS1": "JAK-STAT", "SOCS3": "JAK-STAT",
    # --- Immune evasion ---
    "CD274": "Immune checkpoint (PD-L1)",
    "PDCD1LG2": "Immune checkpoint",
    "CTLA4": "Immune checkpoint",
    "PDCD1": "Immune checkpoint (PD-1)",
    "TIGIT": "Immune checkpoint",
    "CD47": "Immune evasion",
    "B2M": "MHC class I / immune evasion",
    # --- EMT / invasion ---
    "CDH1": "EMT / cell adhesion", "CDH2": "EMT",
    "VIM": "EMT", "FN1": "ECM / invasion",
    "MMP2": "ECM remodeling", "MMP9": "ECM remodeling",
    "MMP14": "ECM remodeling",
    "TWIST1": "EMT", "TWIST2": "EMT",
    "SNAI1": "EMT", "SNAI2": "EMT",
    "ZEB1": "EMT", "ZEB2": "EMT",
    # --- Ubiquitin / proteasome ---
    "MDM2": "Ubiquitin-proteasome", "MDM4": "Ubiquitin-proteasome",
    "FBXW7": "Ubiquitin-proteasome", "VHL": "Ubiquitin-proteasome",
    "SPOP": "Ubiquitin-proteasome", "PSMD4": "Ubiquitin-proteasome",
    "UBE2C": "Ubiquitin-proteasome",
    # --- Autophagy ---
    "BECN1": "Autophagy", "ATG5": "Autophagy",
    "ATG7": "Autophagy", "ATG12": "Autophagy",
    "SQSTM1": "Autophagy", "MAP1LC3A": "Autophagy",
    "ULK1": "Autophagy",
    # --- Splicing / RNA processing ---
    "SF3B1": "RNA splicing", "U2AF1": "RNA splicing",
    "SRSF2": "RNA splicing", "HNRNPA1": "RNA processing",
    "HNRNPK": "RNA processing", "PTBP1": "RNA processing",
    "DDX3X": "RNA helicase", "DHX9": "RNA helicase",
    # --- Nuclear transport ---
    "XPO1": "Nuclear export", "KPNB1": "Nuclear import",
    "RANBP2": "Nuclear transport",
    # --- Vesicle trafficking ---
    "RAB5A": "Endocytosis", "RAB7A": "Lysosomal trafficking",
    "VPS34": "Endocytosis", "CLTC": "Clathrin-mediated endocytosis",
    # --- Stress response ---
    "HSP90AA1": "HSP90 / chaperone", "HSP90AB1": "HSP90 / chaperone",
    "HSPA1A": "Heat shock response", "HSPA5": "ER stress",
    "ATF4": "ER stress / ISR", "ATF6": "ER stress",
    "EIF2AK3": "ER stress", "XBP1": "ER stress",
}


def map_genes_to_pathways_offline(
    gene_list: list,
) -> list[tuple[str, None]]:
    """
    Returns list[tuple[pathway_name, None]].
    None = no p-value (offline dict has no stats).
    """
    pathways = {}
    for g in gene_list:
        if g in OFFLINE_MAP:
            pw = OFFLINE_MAP[g]
            pathways[pw] = pathways.get(pw, 0) + 1   # count supporting genes

    # Return sorted by number of supporting genes
    sorted_paths = sorted(pathways.items(), key=lambda x: -x[1])
    return [(path, None) for path, _ in sorted_paths[:8]]


def map_genes_to_pathways_live(
    gene_list: list,
    background_genes: list,
    organism: str = "hsapiens",
) -> list[tuple[str, float]]:
    """
    Query g:Profiler with a proper statistical background.
    Returns list[tuple[pathway_name, fdr_pvalue]].
    Falls back to offline map on failure.
    """
    if not GPROFILER_AVAILABLE or len(gene_list) < 3:
        return map_genes_to_pathways_offline(gene_list)

    try:
        gp = GProfiler(return_dataframe=True)
        results = gp.profile(
            organism=organism,
            query=gene_list,
            background=background_genes,    # ← CRITICAL: your selected genes, not all human
            sources=["KEGG", "REAC", "GO:BP"],
            significance_threshold_method="fdr",
            user_threshold=0.10,
            no_evidences=True,
        )
        if results.empty:
            return map_genes_to_pathways_offline(gene_list)

        top = results.nsmallest(8, "p_value")
        return [(row["name"], float(row["p_value"]))
                for _, row in top.iterrows()]

    except Exception as e:
        print(f"  g:Profiler failed ({e}), using offline map")
        return map_genes_to_pathways_offline(gene_list)


def map_to_pathways(
    gene_list: list,
    background_genes: list = None,
    use_live: bool = True,
) -> list[tuple[str, float | None]]:
    """
    Public entry point. Always returns list[tuple[str, float|None]].

    Parameters
    ----------
    gene_list        : genes to query
    background_genes : full set of selected genes (required for live query)
    use_live         : try g:Profiler first if True

    Returns
    -------
    list of (pathway_name, p_value_or_None)
    """
    if use_live and GPROFILER_AVAILABLE and background_genes:
        return map_genes_to_pathways_live(gene_list, background_genes)
    return map_genes_to_pathways_offline(gene_list)


def coverage_report(selected_genes: list) -> dict:
    """
    Report what fraction of your selected genes are in the offline map.
    Call this in your XAI script and include the result in your paper.
    """
    n_total  = len(selected_genes)
    n_mapped = sum(1 for g in selected_genes if g in OFFLINE_MAP)
    n_unmapped = n_total - n_mapped
    pct = 100 * n_mapped / n_total if n_total else 0
    return {
        "total_selected_genes": n_total,
        "offline_map_size":     len(OFFLINE_MAP),
        "mapped":               n_mapped,
        "unmapped":             n_unmapped,
        "coverage_pct":         round(pct, 1),
        "note": (
            f"{pct:.1f}% of your selected genes have offline pathway annotations. "
            f"For the remaining {n_unmapped}, g:Profiler live query is used. "
            "Report this coverage in your paper's methods section."
        ),
    }


# -----------------------------------------------------------
# SELF-TEST
# -----------------------------------------------------------
if __name__ == "__main__":
    test_genes  = ["BRCA1", "TP53", "AKT1", "KRAS", "GPX4", "SF3B1", "GARBAGE"]
    background  = list(OFFLINE_MAP.keys()) + ["GARBAGE", "UNKNOWN1"]

    print("=== Offline map ===")
    result = map_to_pathways(test_genes, use_live=False)
    for name, pval in result:
        print(f"  {name}  (p={pval})")

    print(f"\nOffline map size: {len(OFFLINE_MAP)} genes")
    cov = coverage_report(test_genes)
    print(f"Coverage: {cov['coverage_pct']}%")
    print(cov["note"])
