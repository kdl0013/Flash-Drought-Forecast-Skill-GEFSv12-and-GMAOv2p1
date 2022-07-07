#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 09:19:12 2022

@author: kdl
"""
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np

np.random.seed(19680801)

x = np.arange(-0.5, 10,5)  # len = 11
x_len = len(x)
y = np.arange(4.5, 11, 2)  # len = 7
y_len = len(y)

Z = np.random.rand(y_len-1, x_len-1)

fig, ax = plt.subplots()
ax.pcolormesh(x, y, Z)