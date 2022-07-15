import pandas as pd

name = 'HpolPA_180added'
open_file = name + '.dat'
close_file = name + '.csv'
df = pd.read_table(open_file, sep="\s+",header = None)

df.to_csv(close_file)

