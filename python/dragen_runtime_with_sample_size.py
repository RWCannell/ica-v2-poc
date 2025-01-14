import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csv_file = "python/dragen_runtime_with_output_size.csv"
    
df = pd.read_csv(csv_file)  
   
# print(df)

dragen_runtime_values = df['dragen_runtime'].values
print(dragen_runtime_values)