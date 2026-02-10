
import os
import sys
import pandas as pd
import numpy as np

from graphposterior.matrix_utils import calculate_avg_informative_characters



input_file = sys.argv[1]

avg_informative = calculate_avg_informative_characters(input_file)

print(f"Result: {avg_informative:.2f}")
