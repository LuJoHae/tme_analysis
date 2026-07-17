#!/usr/bin/env python3
"""
Orchestration Script for Re-running Linear Deconvolution and Predictions.
Overwrites cached deconvolution files in Lair and runs all prediction and analysis scripts.
"""

import sys
import subprocess
from pathlib import Path

# Ensure package is imported
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))
import datalair
from single_cell_datasets import SingleCellDeconvolution, SingleCellBagaevDeconvolution

def main():
    lair = datalair.Lair("/storage/halu/lair")
    python_bin = sys.executable
    
    print("====================================================")
    print("STEP 1: Reusing SingleCellDeconvolution from Lair...")
    print("====================================================")
    ds_deconv = SingleCellDeconvolution()
    lair.safe_derive(ds_deconv, overwrite=False)
    
    print("\n====================================================")
    print("STEP 2: Reusing SingleCellBagaevDeconvolution from Lair...")
    print("====================================================")
    ds_bagaev = SingleCellBagaevDeconvolution()
    lair.safe_derive(ds_bagaev, overwrite=False)
    
    print("\n====================================================")
    print("STEP 3: Rerunning Standard Deconvolution Predictions...")
    print("====================================================")
    subprocess.run([python_bin, "scripts/iAtlas-predict-response-deconv.py"], check=True)
    
    print("\n====================================================")
    print("STEP 4: Rerunning Advanced Predictions (Multi-Seed + Error Bars)...")
    print("====================================================")
    subprocess.run([python_bin, "scripts/iAtlas-predict-response-advanced.py"], check=True)
    
    print("\n====================================================")
    print("STEP 5: Rerunning Bagaev Constrained Comparative Predictions...")
    print("====================================================")
    subprocess.run([python_bin, "scripts/iAtlas-bagaev-constrained-prediction.py"], check=True)
    
    print("\n====================================================")
    print("STEP 6: Rerunning Cohort Investigation & Diagnostics...")
    print("====================================================")
    subprocess.run([python_bin, "scripts/investigate-cohort-clustering.py"], check=True)
    
    print("\n====================================================")
    print("STEP 7: Rerunning Response Correlation Analysis (Self-Constructed Reference)...")
    print("====================================================")
    subprocess.run([python_bin, "scripts/iAtlas-self-ref-correlation.py"], check=True)
    
    print("\n====================================================")
    print("Pipeline Rerun Complete!")
    print("====================================================")

if __name__ == "__main__":
    main()
