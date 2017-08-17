#!/usr/bin/env bash
conda create -y -n test_env pip
source activate test_env
pip install -e . || python setup.py install
