import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
from math import exp,sqrt, log
from scipy import optimize
import pandas as pd

# lines = open('lognor_distr_data_1.txt').read().splitlines()
# lines = lines[3:]
# df = pd.DataFrame(lines)
# print(df[:,1])
# df = df["2"].str.replace("   ","w")
# df.to_csv('test_1.csv')

df = pd.read_csv('lognor_distr_data_1.csv', sep=" ", header=None)
print(df[0])
# df = df.drop([0,1,2])
# df.to_csv('test_1.csv')
