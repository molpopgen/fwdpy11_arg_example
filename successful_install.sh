#!/usr/bin/env bash
# WARNING this will overwrite an environment named fwdpy11_arg_test_env
conda create -yf -n fwdpy11_arg_test_env pip
source activate fwdpy11_arg_test_env
conda env update -q -n fwdpy11_arg_test_env -f environment.yml
pip install -e . && python setup.py install
