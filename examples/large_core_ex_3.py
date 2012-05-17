# examples/large_core_ex_1.py
#
# In this example, we investigate "by hand" the
# large reactor example.

import yamamoto_benchmark
import numpy as np
import time
#def objective(k, p) :
#  return k

def objective(k, p):
    """  Default objective function.
    """
    delta = 0
    if k < 1.1 :
        delta = k - 1.1
    return 1.0 * (1.5 - p) + 50.0 * delta


# Get the reactor from the premade script.
reactor = yamamoto_benchmark.make_large_core()

# Evaluate the default pattern. 
k_best, p_best = reactor.evaluate()
f_best = objective(k_best, p_best)

print "initial k = ",k_best," p = ",p_best, " f = ",f_best


#reactor.plot_peaking()
stencil = reactor.core.stencil
pattern= reactor.core.pattern

kij = np.zeros((np.size(pattern), 2))

k = 0
for i in range(0, 9):
  for j in range(0, 9):
    if stencil[i][j] > 0 :
      kij[k][0] = i
      kij[k][1] = j
      k = k + 1

#t1 = time.time()
#for it in range(0, 10000):
#  reactor.evaluate()
#etime =  time.time()-t1
#print " elapsed time = ", etime, " seconds"
#print " time per eval = ", etime / 10000, " sec/eval"

evals = 0

# Single swap
#t1 = time.time()
#for it in range(0, 1) :
#  for k1 in range(1, len(pattern)) :
#    for k2 in range(1, len(pattern)) :
#      if k1 == k2 :
#        pass
#      else :
#        evals = evals + 1
#        i1 = kij[k1][0]
#        j1 = kij[k1][1]
#        i2 = kij[k2][0]
#        j2 = kij[k2][1]
#        reactor.swap([i1,j1],[i2,j2])
#        keff, peak = reactor.evaluate()
#        f = objective(keff, peak)
#        if f > f_best :
#          print "**** new k = ",keff," p = ",peak, " f = ",f, " evals = ", evals
#          k_best = keff
#          p_best = peak
#          f_best = f
#          pattern = reactor.core.pattern
#        else :
#          # didn't improve, so undo
#          reactor.swap([i1,j1],[i2,j2])
#etime =  time.time()-t1
#print " elapsed time = ", etime, " seconds"
#print " time per eval = ", etime / 10000, " sec/eval"


evals = 0

t1 = time.time()
done = False
# Dual swap
for it in range(0, 1) :
  for k1 in range(1, len(pattern)) :
    if done == True :
      break
    for k2 in range(1, len(pattern)) :
      if done == True :
        break
      if k1 == k2 :
        pass
      else :
        for k3 in range(1, len(pattern)) :
          if done == True :
            break
          for k4 in range(1, len(pattern)) :
            if k3 == k4 :
              pass
            else :
              evals = evals + 1
              #print evals
              if evals == 100000 :
                done = True
              i1 = kij[k1][0]
              j1 = kij[k1][1]
              i2 = kij[k2][0]
              j2 = kij[k2][1]
              i3 = kij[k3][0]
              j3 = kij[k3][1]
              i4 = kij[k4][0]
              j4 = kij[k4][1]
              reactor.swap([i1,j1],[i2,j2])
              reactor.swap([i3,j3],[i4,j4])
              keff, peak = reactor.evaluate()
              f = objective(keff, peak)
              
              if f > f_best :
                k_best = keff
                p_best = peak
                f_best = f
                pattern = reactor.core.pattern
                print "**** new k = ",keff," p = ",peak, " f = ",f, " evals = ", evals
              else :
                #print "     new k = ",keff," p = ",peak, " f = ",f
                # didn't improve, so undo
                reactor.swap([i3,j3],[i4,j4])
                reactor.swap([i1,j1],[i2,j2])
etime =  time.time()-t1
print " elapsed time = ", etime, " seconds"
print " time per eval = ", etime / evals, " sec/eval"

#print "final k = ",k_best," p = ",p_best, " f = ",f_best
#print " evals = ", evals

#reactor.plot_peaking()
