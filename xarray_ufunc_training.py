#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 17:43:30 2022

@author: kdl
"""
import xarray as xr
import numpy as np


def magnitude(a, b):
    func = lambda x, y: np.sqrt(x**2 + y**2)
    return xr.apply_ufunc(func, a, b)

array = xr.DataArray([1, 2, 3], coords=[("x", [0.1, 0.2, 0.3])])
magnitude(array, -array)

magnitude(3, 4)
magnitude(3, np.array([0, 4]))
magnitude(array, 0)
