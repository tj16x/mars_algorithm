#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 11:30:13 2018

@author: borislav
"""

from scipy.stats import rv_continuous

import numpy as np

_norm_pdf_C = np.sqrt(2*np.pi)

def _norm_pdf(x):
    return np.exp(-x**2/2.0) / _norm_pdf_C


class johnsonsl_gen(rv_continuous):
    
    '''Class for generating Johnson SL distribution'''
    
    def _pdf(self, x, a, b):
        
        trm = _norm_pdf(a + b * np.log(x))
        return b*x*1.0*trm
    
    
    def rvs(self, *args, **kwds):
        
        discrete = kwds.pop('discrete', None)
        args, loc, scale, size = self._parse_args_rvs(*args, **kwds)
        
        if np.all(scale == 0):
            return loc*np.ones(size, 'd')

        # `size` should just be an argument to _rvs(), but for, um,
        # historical reasons, it is made an attribute that is read
        # by _rvs().
        self._size = size
        vals = self._rvs(*args)

        vals = vals * scale + loc

        # Cast to int if discrete
        if discrete:
            if size == ():
                vals = int(vals)
            else:
                vals = vals.astype(int)

        return vals
    
johnsonsl = johnsonsl_gen(name='johnsonsl')


class johnsonsn_gen(rv_continuous):
    
    '''Class for generating Johnson SN distribution'''
    
    def _pdf(self, x, a, b):
        
        trm = _norm_pdf(a + b*x)
        return b*1.0*trm


johnsonsn = johnsonsn_gen(name='johnsonsn')