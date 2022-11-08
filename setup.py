#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    setup.py
#
#    Setup for the phospho-iLIR software.
#
#    Copyright (C) 2022 Valentina Sora 
#                       <sora.valentina1@gmail.com>
#                       Matteo Tiberti 
#                       <matteo.tiberti@gmail.com> 
#                       Elena Papaleo
#                       <elenap@cancer.dk>
#
#    This program is free software: you can redistribute it and/or
#    modify it under the terms of the GNU General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program. 
#    If not, see <http://www.gnu.org/licenses/>.



# Standard library
from setuptools import setup



# Name of the software
name = "phospho_iLIR"

# URL where to download the package
url = "https://github.com/ELELAB/phospho_iLIR"

# Author(s)
author = "Valentina Sora, Matteo Tiberti, Elena Papaleo"

# E-mail address of the main developer
author_email = "sora.valentina1@gmail.com"

# Version
version = "0.0.1"

# Software description
description = \
      "phospho_iLIR - a pipeline to produce all phospho-mimetic " \
      "variants of a putative phospho-regulated LIR, and to assess " \
      "the possible effect of phosphorylations on the variants' " \
      "secondary structure propensities."

# Directories containing the package(s)
package_dir = {"phospho_iLIR" : "phospho_iLIR"}

# Package(s) included
packages = ["phospho_iLIR"]

# Entry points
entry_points = \
      {"console_scripts" : \
            ["phospho_iLIR = phospho_iLIR.phospho_iLIR:main"], \
      }

# Dependencies
install_requires = \
      ["beautifulsoup4", "dask", "distributed", "html5lib",
       "matplotlib",  "pandas", "pyyaml", "requests", "six"]

# Launch the setup
setup(name = name,
      url = url,
      author = author,
      author_email = author_email,
      version = version,
      description = description,
      include_package_data = True,
      package_dir = package_dir,
      packages = packages,
      entry_points = entry_points,
      install_requires = install_requires)