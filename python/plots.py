import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csv_file = "python/runtimes_and_costs_for_SGDVP_data.csv"
    
df = pd.read_csv(csv_file)  
   
dragen_runtime_values = df['dragen_runtime'].values
workflow_runtime_values = df['workflow_runtime'].values
approximate_cost_values = df['approximate_cost'].values
actual_cost_values = df['actual_cost'].values
sample_size_values = df['sample_size'].values
approximate_cost_values = df['approximate_cost'].values

plt.scatter(dragen_runtime_values, approximate_cost_values)
plt.show()

print(df.describe())
