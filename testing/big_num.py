from time import time
import numpy as np

a = np.array([3,4,5],dtype=np.longdouble)
b = np.array([3,4,5],dtype=object)
a_ = np.zeros(3,dtype=np.longdouble)
b_ = np.zeros(3,dtype=object)

t0 = time()
a_[:] = 12345678**a[:]%a[:]
t1 = time()
print(a_)
print(t1-t0)

t2 = time()
b_[:] = 12345678**b[:]*b[:]
t3 = time()
print(b_)
print(t3-t2)



