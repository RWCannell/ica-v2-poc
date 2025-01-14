import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csv_file = "python/dragen_runtime_with_approximate_cost.csv"
    
df = pd.read_csv(csv_file)  
   
dragen_runtime_values = df['dragen_runtime'].values
sample_size_values = df['sample_size'].values
approximate_cost_values = df['approximate_cost'].values

print(df.describe())
