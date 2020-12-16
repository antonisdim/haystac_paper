import datetime
import pandas as pd 
import numpy as np
from sys import argv 

script, input_file = argv 

test_results = pd.read_csv(input_file, sep =',') 

test_small = test_results[['Maximum resident set size (kbytes)', 'Elapsed (wall clock) time (h:mm:ss or m:ss)', 'File']]

test_small.rename(columns={'Maximum resident set size (kbytes)': 'RSS', 'Elapsed (wall clock) time (h:mm:ss or m:ss)': 'Runtime'}, inplace=True) 

test_small['Runtime_new'] = test_small['Runtime'].apply(lambda x: datetime.datetime.strptime(x,'%M:%S.%f')) 

test_small['Timedelta'] = test_small['Runtime_new'] - datetime.datetime.strptime('00:00.0','%M:%S.%f') 

test_small['Runtime_s'] = test_small['Timedelta'].apply(lambda x: x / np.timedelta64(1, 's')) 

test_small_format = test_small[['File', 'RSS', 'Runtime', 'Runtime_s']] 

test_small_format['Runtime_m'] = test_small_format['Runtime_s']/60

test_small_format['Method'] = '' 

for idx, row in test_small_format.iterrows(): 
	if 'malt' in row['File']: 
		test_small_format.at[idx, 'Method'] = 'malt' 
	elif 'haystack' in row['File']: 
		test_small_format.at[idx, 'Method'] = 'haystac' 
	elif 'kraken' in row['File']: 
		test_small_format.at[idx, 'Method'] = 'kraken' 
	elif 'sigma' in row['File']: 
		test_small_format.at[idx, 'Method'] = 'sigma' 
	
test_small_format['Number of Reads'] = '' 


for idx, row in test_small_format.iterrows(): 
	if '10K' in row['File']: 
		test_small_format.at[idx, 'Number of Reads'] = 10000 
	elif '100K' in row['File']: 
		test_small_format.at[idx, 'Number of Reads'] = 100000  
	elif '1M' in row['File']: 
		test_small_format.at[idx, 'Number of Reads'] = 1000000 

test_small_format['Database'] = '' 

for idx, row in test_small_format.iterrows(): 
	if '10_species' in row['File']: 
		test_small_format.at[idx, 'Database'] = 10 
	elif '100_species' in row['File']: 
		test_small_format.at[idx, 'Database'] = 100 
	elif '500_species' in row['File']: 
		test_small_format.at[idx, 'Database'] = 500 


test_small_format[test_small_format['Number of Reads'] != ""].to_csv('Sample_benchmark_new.txt', sep='\t', index=False) 

test_small_format[(test_small_format['Number of Reads'] == "") & (~test_small_format['File'].str.contains('8000'))].to_csv('DB_benchmark_new.txt', sep='\t', index=False) 

test_small_format.to_csv('Summary_Table.tsv', sep='\t', index=False) 

