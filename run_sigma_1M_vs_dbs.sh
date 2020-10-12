conda activate rip

mkdir sigma_outputs/input_1M_10sp
mkdir sigma_outputs/input_1M_100sp
mkdir sigma_outputs/input_1M_1000sp
mkdir sigma_outputs/input_1M_500sp

time -v ./Sigma/bin/sigma-align-reads -c sigma_config_input_1M_10sp.cfg -w sigma_outputs/input_1M_10sp -p `grep -c ^processor /proc/cpuinfo`; ./Sigma/bin/sigma-build-model -c sigma_config_input_1M_10sp.cfg -w sigma_outputs/input_1M_10sp/; ./Sigma/bin/sigma-solve-model -c sigma_config_input_1M_10sp.cfg -w sigma_outputs/input_1M_10sp/ -t `grep -c ^processor /proc/cpuinfo`; mv sigma_out.* ./sigma_outputs/input_1M_10sp/

time -v ./Sigma/bin/sigma-align-reads -c sigma_config_input_1M_100sp.cfg -w sigma_outputs/input_1M_100sp -p `grep -c ^processor /proc/cpuinfo`; ./Sigma/bin/sigma-build-model -c sigma_config_input_1M_100sp.cfg -w sigma_outputs/input_1M_100sp/; ./Sigma/bin/sigma-solve-model -c sigma_config_input_1M_100sp.cfg -w sigma_outputs/input_1M_100sp/ -t `grep -c ^processor /proc/cpuinfo`; mv sigma_out.* ./sigma_outputs/input_1M_100sp/

time -v ./Sigma/bin/sigma-align-reads -c sigma_config_input_1M_1000sp.cfg -w sigma_outputs/input_1M_1000sp -p `grep -c ^processor /proc/cpuinfo`; ./Sigma/bin/sigma-build-model -c sigma_config_input_1M_1000sp.cfg -w sigma_outputs/input_1M_1000sp/; ./Sigma/bin/sigma-solve-model -c sigma_config_input_1M_1000sp.cfg -w sigma_outputs/input_1M_1000sp/ -t `grep -c ^processor /proc/cpuinfo`; mv sigma_out.* ./sigma_outputs/input_1M_1000sp/

time -v ./Sigma/bin/sigma-align-reads -c sigma_config_input_1M_500sp.cfg -w sigma_outputs/input_1M_500sp -p `grep -c ^processor /proc/cpuinfo`; ./Sigma/bin/sigma-build-model -c sigma_config_input_1M_500sp.cfg -w sigma_outputs/input_1M_500sp/; ./Sigma/bin/sigma-solve-model -c sigma_config_input_1M_500sp.cfg -w sigma_outputs/input_1M_500sp/ -t `grep -c ^processor /proc/cpuinfo`; mv sigma_out.* ./sigma_outputs/input_1M_500sp/

