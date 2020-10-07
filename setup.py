#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

from setuptools import setup

name = "phospho_iLIR"
url = "https://github.com/ELELAB/phospho_iLIR"
author = "Valentina Sora, Matteo Tiberti, Elena Papaleo"
author_email = "sora.valentina1@gmail.com"
version = "0.1"
description = "phospho_iLIR pipeline"
package_data = {"phospho_iLIR" : ["config/*"]}
package_dir = {"phospho_iLIR" : "phospho_iLIR"}
packages = ["phospho_iLIR"]
entry_points = {"console_scripts" : \
                    ["phospho_iLIR = phospho_iLIR.phospho_iLIR:main"], \
               }
install_requires = ["dask", "distributed", "matplotlib", "pandas", "pyyaml"]


setup(name = name, \
      url = url, \
      author = author, \
      author_email = author_email, \
      version = version, \
      description = description, \
      package_data = package_data, \
      package_dir = package_dir, \
      packages = packages, \
      entry_points = entry_points, \
      install_requires = install_requires)