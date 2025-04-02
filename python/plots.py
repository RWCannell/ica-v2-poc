import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csv_file_SGDVP = "runtimes_in_minutes_and_costs_for_SGDVP_data.csv"
csv_file_1000G = "runtimes_in_minutes_and_costs_for_1000G_data.csv"
    
df_SGDVP = pd.read_csv(csv_file_SGDVP)
df_1000G = pd.read_csv(csv_file_1000G)  
   
dragen_runtime_values_SGDVP = df_SGDVP['dragen_runtime'].values
workflow_runtime_values_SGDVP = df_SGDVP['workflow_runtime'].values
approximate_cost_values_SGDVP = df_SGDVP['approximate_cost'].values
actual_cost_values_SGDVP = df_SGDVP['actual_cost'].values
sample_size_values_SGDVP = df_SGDVP['sample_size'].values
approximate_cost_values_SGDVP = df_SGDVP['approximate_cost'].values

dragen_runtime_values_1000G = df_1000G['dragen_runtime'].values
workflow_runtime_values_1000G = df_1000G['workflow_runtime'].values
actual_cost_values_1000G = df_1000G['actual_cost'].values
sample_size_values_1000G = df_1000G['sample_size'].values

# plt.scatter(sample_size_values_SGDVP, actual_cost_values_SGDVP, label='cost for SGDVP samples')
# plt.title('Cost vs Sample Size for SGVDP Samples')
# plt.xlabel('Sample Size (GB)')
# plt.ylabel('Cost (USD)')
# plt.show()

# plt.scatter(sample_size_values_1000G, actual_cost_values_1000G, label='cost for 1000G samples')
# plt.title('Cost vs Sample Size for 1000G Samples')
# plt.xlabel('Sample Size (GB)')
# plt.ylabel('Cost (USD)')
# plt.show()

# plt.scatter(dragen_runtime_values_SGDVP, actual_cost_values_SGDVP)
# plt.title('DRAGEN Runtime vs Sample Size for SGDVP Samples')
# plt.xlabel('DRAGEN Runtime (seconds)')
# plt.ylabel('Cost (USD)')
# plt.show()

# plt.scatter(dragen_runtime_values_1000G, actual_cost_values_1000G)
# plt.title('Cost vs DRAGEN Runtime for 1000G Samples')
# plt.xlabel('DRAGEN Runtime (seconds)')
# plt.ylabel('Cost (USD)')
# plt.show()

# plt.scatter(dragen_runtime_values_SGDVP, actual_cost_values_SGDVP)
# plt.title('Cost vs DRAGEN Runtime for SGDVP Samples')
# plt.xlabel('DRAGEN Runtime (seconds)')
# plt.ylabel('Cost (USD)')
# plt.show()

# plt.scatter(sample_size_values_SGDVP, dragen_runtime_values_SGDVP)
# plt.title('DRAGEN Runtime vs Sample Size for SGDVP Samples')
# plt.xlabel('Sample Size (GB)')
# plt.ylabel('DRAGEN Runtime (s)')
# plt.show()

# plt.scatter(sample_size_values_1000G, dragen_runtime_values_1000G)
# plt.title('DRAGEN Runtime vs Sample Size for 1000G Samples')
# plt.xlabel('Sample Size (GB)')
# plt.ylabel('DRAGEN Runtime (s)')
# plt.show()

# def generate_scatter_plot(x_values, y_values, x_label, y_label, plot_title):
#     plt.scatter(x_values, y_values)
#     plt.title(plot_title)
#     plt.xlabel(x_label)
#     plt.ylabel(y_label)
#     plt.show()


dragen_runtime_values_SGDVP = df_SGDVP['dragen_runtime'].values
workflow_runtime_values_SGDVP = df_SGDVP['workflow_runtime'].values

dragen_runtime_values_1000G = df_1000G['dragen_runtime'].values
workflow_runtime_values_1000G = df_1000G['workflow_runtime'].values
    
if __name__ == '__main__':
    # generate_scatter_plot(
    #     sample_size_values_1000G, 
    #     dragen_runtime_values_1000G, 
    #     'Sample Size (GB)', 
    #     'DRAGEN Runtime (s)', 
    #     'DRAGEN Runtime vs Sample Size for 1000G Samples'
    # )
    
    # df_SGDVP['dragen_runtime'] = df_SGDVP['dragen_runtime'] / 60
    # df_SGDVP['dragen_runtime'] = np.round(df_SGDVP['dragen_runtime'], decimals=2)
    # df_SGDVP['workflow_runtime'] = df_SGDVP['workflow_runtime'] / 60
    # df_SGDVP['workflow_runtime'] = np.round(df_SGDVP['workflow_runtime'], decimals=2)
    
    # df_1000G['dragen_runtime'] = df_1000G['dragen_runtime'] / 60
    # df_1000G['dragen_runtime'] = np.round(df_1000G['dragen_runtime'], decimals=2) 
    # df_1000G['workflow_runtime'] = df_1000G['workflow_runtime'] / 60
    # df_1000G['workflow_runtime'] = np.round(df_1000G['workflow_runtime'], decimals=2)
    
    df_SGDVP.to_csv('runtimes_in_seconds_and_costs_for_SGDVP_data.csv')
    df_1000G.to_csv('runtimes_in_seconds_and_costs_for_1000G_data.csv')
    print(df_SGDVP.head())
    print(df_1000G.head())
    