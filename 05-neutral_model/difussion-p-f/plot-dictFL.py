"""plot.py"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random as rnd
from copy import copy

L = 1000
A = 4

p_values = [1, 0.9, 0.8, 0.7]  # I should automate this

plt.figure()

for p in p_values:
    df = pd.read_csv(f"./Data/DictFL_{p}_.csv", header=None)
    df.columns = ["mean", "variance"]
    plt.plot(
        np.arange(1, len(df["mean"]) + 1), df["mean"], markersize=1, label=f"f {p}"
    )

    # plt.plot(df_sim['m'], label='m simulation')
    plt.ylabel("<m>")
    plt.xlabel("t")
    plt.legend(loc="best")
    if p == 1:
        plt.plot(df["mean"], df["mean"], label="alpha=1")
plt.yscale("log")
plt.xscale("log")

plt.figure()

for p in p_values:
    df = pd.read_csv(f"./Data/DictFL_{p}_.csv", header=None)
    df.columns = ["mean", "variance"]
    plt.plot(
        np.arange(1, len(df["mean"]) + 1),
        df["variance"] / df["mean"],
        markersize=1,
        label=f"f {p}",
    )

    # plt.plot(df_sim['m'], label='m simulation')
    plt.ylabel("D")
    plt.xlabel("t")
    plt.legend(loc="best")


# plt.yscale("log")
# plt.xscale("log")
plt.show()
