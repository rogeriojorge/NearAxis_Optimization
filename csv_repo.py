#!/usr/bin/env python3

## Create CSV file to be read by pandas
## Easy to use data science/machine learning machinery on it
import pandas as pd

d = {
    'a': (1, 101),
    'b': (2, 202),
    'c': (3, 303)
}
df = pd.DataFrame.from_dict(d, orient="index")


df.to_csv("data.csv")


df = pd.read_csv("data.csv", index_col=0)


d = df.to_dict("split")
d = dict(zip(d["index"], d["data"]))
