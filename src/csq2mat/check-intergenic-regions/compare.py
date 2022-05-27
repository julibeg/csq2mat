# %% ###################################################################
import numpy as np
import pandas as pd

# %% ###################################################################
# get the intergenic regions found with our own method
own_orig = pd.read_csv("intergenic_regions.csv")
own = own_orig[["start", "end"]]
# %% ###################################################################
# get the intergenic regions found by bedtools
bt = pd.read_csv(
    "intergenic_regions_from_bedtools.tsv", sep="\t", header=None, usecols=[1, 2]
)
bt.columns = ['start', 'end']
bt['start'] += 2
# %% ###################################################################
print(bt.head())
print('------------------')
print(own.head())
# %% ###################################################################
own.shape
bt.shape
# %% ###################################################################
own.equals(bt)
# --> looks like our intergenic regions are correct!