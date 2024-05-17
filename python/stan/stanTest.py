# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 14:58:02 2024

@author: zhixi
""" 
import os

from cmdstanpy import CmdStanModel

stan_file = os.path.join('F:\stanModels', 'stan_qLearning_4params.stan')

model = CmdStanModel(stan_file=stan_file)


#%%