# !/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest
import os
import glob
import shutil
from mgktools.features_mol.features_generators import AVAILABLE_FEATURES_GENERATORS
from mgktools.exe.run import mgk_cache_data


CWD = os.path.dirname(os.path.abspath(__file__))
if os.path.exists(f"{CWD}/tmp"):
    shutil.rmtree(f"{CWD}/tmp", ignore_errors=True)
os.mkdir(f"{CWD}/tmp")
data_paths = glob.glob(f"{CWD}/data/*.csv")


def test_CacheData():
    arguments = [
        "--cache_path", f"{CWD}/tmp/cache.pkl",
        "--data_paths"] + data_paths + [
        "--smiles_columns", "smiles", "smiles1", "smiles2",
        "--features_generators_name"] + AVAILABLE_FEATURES_GENERATORS + [
        "--n_jobs", "4",
    ]
    mgk_cache_data(arguments)
