shared_prompt = """
Use these NHANES raw files (merge on SEQN):
DEMO_J.xpt: demographics + survey design vars + weights (WTMEC2YR, SDMVPSU, SDMVSTRA).
DPQ_J.xpt: PHQ-9 depression screener items (DPQ010–DPQ090).
VID_J.xpt: serum 25-hydroxyvitamin D; key variable LBXVIDMS (total 25(OH)D, sum of D2 + D3, excluding epi-D3).
Vitamin D cross-cycle comparability has known assay issues historically; NHANES provides an analytical note (relevant if you later pool cycles).
"""