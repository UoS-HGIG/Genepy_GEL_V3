#!/bin/python3.10
# -*- coding: utf-8 -*-
import sys
#sys.path.append("/drop/.local/lib/python3.10/site-packages")
from numba import njit, cuda, prange
import numpy as np
import pandas as pd

import math
import pyarrow.csv as pa_csv
from typing import Tuple
import argparse
import time
import os
import re
import gc

print("Done!")