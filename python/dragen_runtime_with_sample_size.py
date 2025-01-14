import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csv_file = "python/dragen_runtime_with_output_size.csv"
    
df = pd.read_csv(csv_file)  
   
dragen_runtime_values = df['dragen_runtime'].values
sample_size_values = df['sample_size'].values
output_size_values = df['output_size'].values
total_workflow_runtime_values = df['total_workflow_runtime'].values

print(df.describe())
