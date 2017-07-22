import numpy as np
import tdm_module as tdm

a = np.random.rand(100)
b = np.random.rand(100)
c = np.random.rand(100)
d = np.random.rand(100)

print(tdm.TDMAsolver(a,b,c,d))