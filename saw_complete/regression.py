#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 01:18:04 2021

@author: samarth
"""


import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression

sizes = []

for i in range(40,205,20):
    sizes.append(i)
    
regr_df = pd.DataFrame(columns = ['N', 'Rg'])
    
for size in sizes:
    
    folder = 'outputs/'
    filename = folder+'tmp_'+str(size)+'.out'
    
    df = pd.read_csv(filename, sep = ' ', names = ['x', 'Rg'])
    rg = df.iloc[0]['Rg']
    regr_df.loc[len(regr_df)] = [size, rg]
    
regr_df['N']  = np.log(regr_df['N'])
regr_df['Rg'] = np.log(regr_df['Rg'])

model = LinearRegression()
x = regr_df['N'].values
x = x.reshape(-1,1)
y = regr_df['Rg'].values
model.fit(x,y)
y_preds = model.predict(x)
plt.figure(figsize = (6,4), dpi=300)
plt.scatter(x, y, color = 'red')

coeff = model.coef_
text  = 'Regression coefficient = '+str(coeff)
plt.text(3.8, 2.2, text)

plt.plot(x, y_preds)
plt.title('Rg vs N')
plt.xlabel('ln(N)')
plt.ylabel(r'ln($R_g$)')
plt.savefig('rg_vs_N.pdf', format='pdf')

regr_df.to_csv('rg_vs_N.csv', index=False)