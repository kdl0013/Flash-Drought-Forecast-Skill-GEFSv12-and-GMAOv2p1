#!/usr/bin/env python3
import time

cdef unsigned long long int total
cdef int k
cdef float t1, t2, t

t1 = time.time()

for k in range(1000000000):
    total = total + k
print "Total =", total

t2 = time.time()
t = t2-t1
print("%.100f" % t)
