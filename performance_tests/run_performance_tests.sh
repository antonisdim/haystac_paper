#!/usr/bin/env bash

# Author:    Evangelos A. Dimopoulos, Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Oxford
# Email:     antonisdim41@gmail.com
# License:   MIT

bash run_kraken_db.sh &> kraken_db.log
bash run_kraken_samples.sh &> kraken_samples.log
bash run_kraken_1M_vs_dbs.sh &> kraken_1M_vs_dbs.log

bash run_rip_db_mem.sh &> rip_db_mem.log
bash run_rip_samples.sh &> rip_samples.log
bash run_rip_1M_vs_dbs.sh &> rip_1M_vs_dbs.log
bash run_rip_db_no_mem.sh &> rip_db_no_mem.log 

bash run_sigma_db.sh &> sigma_db.log
bash run_sigma_samples.sh &> sigma_samples.log
bash run_sigma_1M_vs_dbs.sh &> sigma_1M_vs_dbs.log

bash run_malt_db.sh &> malt_db.log
bash run_malt_samples.sh &> malt_samples.log
bash run_malt_1M_vs_dbs.sh &> malt_1M_vs_dbs.log
