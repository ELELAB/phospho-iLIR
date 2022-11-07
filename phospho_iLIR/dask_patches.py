#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    dask_patches.py
#
#    Patches to fix some undesired Dask behaviors.
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
import logging
# Third-party packages
from distributed.utils import DequeHandler



# To address a bug that resets the distributed.worker
# logger to WARNING level when a task is launched on
# the worker, no matter what the configuration was
def reset_worker_logger(log_file = None):
    """Utility function to reset a Dask logger handlers
    and level to desired values.
    """

    # Se the new level
    NEWLEVEL = logging.INFO
    
    # Get the logger
    logger = logging.getLogger("distributed.worker")

    # Define the handlers to keep
    h_to_keep = \
        [h for h in logger.handlers if type(h).__name__ == \
         DequeHandler.__name__]

    formatter = logging.Formatter("%(name)s:%(levelname)s:%(message)s")

    # Remove all the handlers from the logger
    for h in logger.handlers:
        logger.removeHandler(h)
    
    # For each of the handlers to keep
    for h in h_to_keep:
        
        # Set the new level for the handler
        h.setLevel(NEWLEVEL)
        h.setFormatter(formatter)
        # Add the handler to the logger
        logger.addHandler(h)

    if log_file:

        handler = logging.FileHandler(log_file, mode = "a")
        handler.setLevel(NEWLEVEL)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    # Reset the logger level to the new level
    logger.setLevel(NEWLEVEL)
    
    # Return the new logger
    return logger