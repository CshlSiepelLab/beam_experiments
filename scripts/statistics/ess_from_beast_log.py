#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

def calculate_ess(trace, burn_in=0.1, max_lag=2000):
    # Remove burn-in
    n_burn_in = int(len(trace) * burn_in)
    trace = trace[n_burn_in:]
    
    # Input validation
    if len(trace) == 0:
        raise ValueError("Trace is empty after burn-in removal")
    
    # Calculate mean and auto-correlation
    mean = np.mean(trace)
    n = len(trace)
    max_lag = min(n - 1, max_lag)  # Ensure max_lag doesn't exceed n-1
    
    # Calculate auto-correlation for each lag
    acf = np.correlate(trace - mean, trace - mean, mode='full')[n-1:]
    acf = acf[:max_lag] / acf[0]
    
    # Calculate ESS using auto-correlation
    positive_acf = np.where(acf <= 0)[0]
    cutoff = len(acf) if len(positive_acf) == 0 else positive_acf[0]
    act = 1 + 2 * np.sum(acf[1:cutoff])
    
    # Handle edge cases where ACT might be too small
    if act < 1:
        act = 1
    
    ess = n / act
    return ess

def main():
    
    # File and parameter settings
    log_file = sys.argv[1]

    # Get parameters from command line
    parameters = sys.argv[2].split(',')
    
    # Read and process data
    data = pd.read_csv(log_file, sep='\t', comment='#')
    
    # Check if all parameters exist in the data
    missing_params = [p for p in parameters if p not in data.columns]
    if missing_params:
        raise ValueError(f"Parameters not found in log file: {', '.join(missing_params)}")
    
    # Calculate and print ESS for each parameter
    for parameter in parameters:
        ess = calculate_ess(data[parameter].values)
        print(f"ESS for {parameter}: {ess}")

if __name__ == "__main__":
    main()
