#run sigma sample analysis performance test against databases of various sizes

MAX_CPU=$(grep -c ^processor /proc/cpuinfo)

mkdir sigma_outputs/input_1M_10sp
mkdir sigma_outputs/input_1M_100sp
mkdir sigma_outputs/input_1M_1000sp
mkdir sigma_outputs/input_1M_500sp

time -v ./Sigma/bin/sigma-align-reads \
	-c sigma_config_input_1M_10sp.cfg \
	-w sigma_outputs/input_1M_10sp \
	-p $MAX_CPU; \
	./Sigma/bin/sigma-build-model \
	-c sigma_config_input_1M_10sp.cfg \
	-w sigma_outputs/input_1M_10sp/; \
	./Sigma/bin/sigma-solve-model \
	-c sigma_config_input_1M_10sp.cfg \
	-w sigma_outputs/input_1M_10sp/ \
	-t $MAX_CPU; \
	mv sigma_out.* ./sigma_outputs/input_1M_10sp/

time -v ./Sigma/bin/sigma-align-reads \
	-c sigma_config_input_1M_100sp.cfg \
	-w sigma_outputs/input_1M_100sp \
	-p $MAX_CPU; \
	./Sigma/bin/sigma-build-model \
	-c sigma_config_input_1M_100sp.cfg \
	-w sigma_outputs/input_1M_100sp/; \
	./Sigma/bin/sigma-solve-model \
	-c sigma_config_input_1M_100sp.cfg \
	-w sigma_outputs/input_1M_100sp/ \
	-t $MAX_CPU; \
	mv sigma_out.* ./sigma_outputs/input_1M_100sp/

time -v ./Sigma/bin/sigma-align-reads \
	-c sigma_config_input_1M_1000sp.cfg \
	-w sigma_outputs/input_1M_1000sp \
	-p $MAX_CPU; \
	./Sigma/bin/sigma-build-model \
	-c sigma_config_input_1M_1000sp.cfg \
	-w sigma_outputs/input_1M_1000sp/; \
	./Sigma/bin/sigma-solve-model \
	-c sigma_config_input_1M_1000sp.cfg \
	-w sigma_outputs/input_1M_1000sp/ \
	-t $MAX_CPU; \
	mv sigma_out.* ./sigma_outputs/input_1M_1000sp/

time -v ./Sigma/bin/sigma-align-reads \
	-c sigma_config_input_1M_500sp.cfg \
	-w sigma_outputs/input_1M_500sp \
	-p $MAX_CPU; \
	./Sigma/bin/sigma-build-model \
	-c sigma_config_input_1M_500sp.cfg \
	-w sigma_outputs/input_1M_500sp/; \
	./Sigma/bin/sigma-solve-model \
	-c sigma_config_input_1M_500sp.cfg \
	-w sigma_outputs/input_1M_500sp/ \
	-t $MAX_CPU; \
	mv sigma_out.* ./sigma_outputs/input_1M_500sp/

