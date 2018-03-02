# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 14:01:48 2018

@author: cdsch
"""

from sympy import *
n = symbols('n')
rho = 90
MNEUTRON= 50
solve(Eq(236*n**(2.54)+n*MNEUTRON,90),n)