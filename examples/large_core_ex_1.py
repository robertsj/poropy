# examples/coretools/large_core_ex_1.py

import large_core
import time
import numpy as np
# Here, we'll investigate the small benchmark core.  First,
# build it using the premade script.
reactor = large_core.make_large_core()

#reactor.shuffle(np.array([48,38,19,2,31,15,34,37,20,9\
#                                 ,30,0,41,1,28,23,16,17,21,8,\
#                                 12,13,29,33,26,14,35,32,45,24,\
#                                 7,10,27,3,22,6,18,40,4,11,5,25,\
#                                 46,44,43,42,39,47,36]))

# FLARE best
#reactor.shuffle(np.array([48,46,35,30,34,25,12,42,7,32,3,43,11,38,33,1,39,24,13,44,8,27,19,31,21,29,4,16,23,47,10,41,22,17,14,20,15,40,2,37,36,5,18,28,6,45,0,9,26]))

# LABAN, order 4, best
#reactor.shuffle(np.array([48,27,24,28,18,0,6,41,30,11,46,34,9,43,37,25,20,12,21,15,4,45,16,17,19,26,13,31,38,33,29,32,3,22,47,44,7,1,8,2,39,14,10,23,5,36,42,35,40]))

# View all the diagnostics down the chain.
reactor.display()
print "len sten", len(reactor.core.stencil[1,:]), len(reactor.core.stencil[:,1])
# Evaluate the default pattern.  We can grab the eigenvalue
# and peaking as return values.
num = 1
t = time.time()
for i in range(0, num) :
    k, p = reactor.evaluate()
print "k = ",k," p = ",p
t2 = time.time()-t
print "elapsed time = ", t2
print " per / eval  = ", t2/num
print "k = ",k," p = ",p
# Alternatively, we can use print_params to display current
# values of all optimization parameters.  Currently only
# keff and the max peaking are retained.
reactor.print_params()
reactor.print_peaking()
## We can also print the power peaking.
#reactor.print_peaking()
#
## With this, we can try optimizing by hand a bit.  Peaking
## occurs at (0, 1).  Printing the pattern helps visualize this.
#reactor.print_pattern()

# Do a swap and evaluate.
reactor.swap([0,1],[0,2])
reactor.print_pattern()
reactor.evaluate()
reactor.print_params()
reactor.print_peaking()
reactor.plot_peaking()
#reactor.plot_pattern('burnup')

# That's a significant peaking reduction with just a slight decrease
# in keff.  However, there is a better pattern.  For this keff, 
# the tradeoff curve in the theory document suggests a peaking 
# on the order of 1.75 or maybe less.  But we won't find out
# by hand...
