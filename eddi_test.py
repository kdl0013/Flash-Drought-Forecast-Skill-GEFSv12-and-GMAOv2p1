#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 09:35:06 2022

@author: kdl
"""

import numpy as np
from scipy.stats import rankdata
from scipy import stats

#Weekly accumulated reference ETo from 40 weeks
data = np.array([14.420787811279297, 5.172389984130859, 6.041740417480469, 7.319391250610352, 5.513860702514648, 11.390810012817383, 1.5083518028259277, 8.953362464904785, 7.720783710479736, 9.219938278198242, 9.839365005493164, 5.6793599128723145, 5.595060348510742, 13.403731346130371, 11.687836647033691, 10.184643745422363, 5.020684242248535, 15.581283569335938, 13.49914264678955, 2.6485064029693604, 15.74343490600586, 10.484978675842285, 16.945966720581055, 12.788093566894531, 8.48973560333252, 5.729422092437744, 20.70046615600586, 12.117950439453125, 9.18446159362793, 11.720061302185059, 10.383907318115234, 11.404386520385742, 13.245698928833008, 6.936684608459473, 7.616708755493164, 13.68994426727295, 12.031225204467773, 4.708139896392822, 2.7529962062835693, 6.448971748352051])

print(f'Maximum ETo value: {np.max(data)}, minimum is: {np.min(data)}.')
print(data[:10])
#Apply EDDI formula
''' i = 1 for maximum ETo aggregation in time window'''

# We see that the first value is ranked 5, showing higher ETo summation which is good
ranking = rankdata([-1 * i for i in data]).astype(int)
print(ranking[0:10])

#Empirical probabilities look fine, showing lower probabilities for higher values
#and higher probabilities for lower values
tukey = (ranking - 0.33)/ (len(ranking) + 0.33)
print(tukey[0:10])

eddi = stats.norm.ppf(tukey)


w = []
reverse = [] #keep this to reverse the sign of EDDI if probability is greater than 0.5
'''Now we can update the w term in EDDI based on the probability'''
for probability in tukey:
    if probability <= 0.5:
        w.append(np.sqrt(-2*np.log(probability)))
        reverse.append(1)
    elif probability > 0.5:
        w.append((1-probability))  
        reverse.append(-1)
        
print(w[0:10])    

''''I'm not sure if the new probabilites are sufficient. I may have misinterpreted your
formula regarding when to actually change the value of EDDI. Either for all EDDI values
at the very end, or if it needed to be changed based on the new updated w value for
changing the empirical tukey plotting position?'''
  


#constants
c0 = 2.515517
c1 = 0.802853
c2 = 0.010328
d1 = 1.432788
d2 = 0.189269
d3 = 0.001308
    

#This eddi version does create actual negative values, but they don't correspond 
#with the actual values
#Must reverse the sign of EDDI,so multiply by -1
final_out_eddi = []
for idx, w_val in enumerate(w):
    final_out_eddi.append((w_val - ((c0 + c1*w_val + c2*w_val**2)/
                                    (1 + d1*w_val + d2*w_val**2 + d3*w_val**3)))* -1)

'''This interpretation does make negative values, but they are not meaningful 
because we would expect higher EDDI values for higher accumulated ETo. The first day
in the list should have a fairly high EDDI because it is ranked 5th for accumulated ETo'''

for idx, val in enumerate(final_out_eddi):
    if idx == 20:
        break
    print(f'ETo summation was {round(data[idx],3)}.  EDDI is {round(final_out_eddi[idx],3)}')



'''Here is the other interpretation where we only change the sign of EDDI if the 
changed empirical plotting position was greater than 0.5.
But these values are much worse. No negative values at all'''

final_out_eddi = []
for idx, w_val in enumerate(w):
    final_out_eddi.append((w_val - ((c0 + c1*w_val + c2*w_val**2)/
                                    (1 + d1*w_val + d2*w_val**2 + d3*w_val**3)))* reverse[idx]) #added the reverse sign here
    
for idx, val in enumerate(final_out_eddi):
    if idx == 20:
        break
    print(f'ETo summation was {round(data[idx],3)}.  EDDI is {round(final_out_eddi[idx],3)}')


