# todo 

# create conda env 

conda env create -f environment.yaml -n myenv

# activate conda env 

conda activate myenv 

# set the right options for ram usage for malt 

sed -i'' 's/-Xmx64G/-Xmx1024G/' ~/.conda/pkgs/malt-0.41-1/opt/malt-0.41/malt-run.vmoptions
sed -i'' 's/-Xmx64G/-Xmx1024G/' ~/.conda/pkgs/malt-0.41-1/opt/malt-0.41/malt-build.vmoptions

# install haystack 

python -m pip install git+https://github.com/antonisdim/Haystack

# run the performance tests

bash run_performance_test.sh &> performance_test.log 

