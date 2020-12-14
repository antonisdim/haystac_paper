import pandas as pd 
import glob 
from sys import argv 

script, folder = argv  

folder = folder.rstrip('/')

df_list = []  

for file_path in glob.glob(folder + '/*time.log'): 
    stats_dict = {} 
    with open(file_path, 'r') as fin: 
        for line in fin: 
            line = line.rstrip().lstrip('\t').split(': ') 
            stats_dict[line[0]] = line[1] 
    stats_df = pd.DataFrame(stats_dict, index=[0])     
    stats_df['File'] = file_path 
    df_list.append(stats_df) 
#     print(df_list)
#     print(file_path) 
    
total_df = pd.concat(df_list)


total_df.to_csv(folder + 'stats.tsv', index = False)
