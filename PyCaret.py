
	# Install miniconda
	# create a conda environment: conda create --name pycaret python=3.8
	
	# activate conda environment: conda activate pycaret
	# To deactivate an active environment, use: conda deactivate
	# To remove: conda env remove --name pycaret
	

	# install pycaret: pip install pycaret
	# install full version: pip install 'pycaret[full]'
	# To avoid conflicts with function called _format_load_msg: pip install "joblib<1.4.0"
	# (configure): pip install ipykernel
	# (configure):  pip install notebook

	
	# Finally, initiate/create notebook kernel: 
	# python -m ipykernel install --user --name pycaret --display-name "pycaret"
	# jupyter kernelspec list
	# jupyter notebook 
