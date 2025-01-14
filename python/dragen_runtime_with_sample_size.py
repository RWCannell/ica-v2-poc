import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

file_name = "python/dragen_runtime_with_output_size.csv"

# columns = []
# rows = []

# # reading csv file
# with open(file_name, 'r') as csv_file:
#     csv_reader = csv.reader(csv_file)

#     columns = next(csv_reader)

#     for row in csv_reader:
#         rows.append(row)

#     print("Total number of rows: %d" % (csv_reader.line_num))

# print('Column names are:' + ', '.join(column for column in columns))

# print('\nFirst 10 rows are:\n')
# for row in rows[:10]:
#     for col in row:
#         print("%10s" % col, end=" "),
#     print('\n')
    
df = pd.read_csv(file_name)  
   
# output the dataframe 
print(df)

# dragen_runtime_values = df['dragen_runtime'].values
# print(dragen_runtime_values)