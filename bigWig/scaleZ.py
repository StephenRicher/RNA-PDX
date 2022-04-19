#!/usr/bin/env python3

import sys
import pandas as pd
from scipy.stats import zscore

dtypes = {'chrom': str, 'start': int, 'end': int, 'score': float}
df = pd.read_csv(sys.stdin, dtype=dtypes, names=dtypes.keys(), sep='\t')
df['score'] = zscore(df['score'])
df.to_csv(sys.stdout, header=False, index=False, sep='\t')
