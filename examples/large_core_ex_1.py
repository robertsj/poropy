# examples/large_core_ex_1.py
#
# In this example, we investigate "by hand" the
# large reactor example.

import large_core
import time
import numpy as np

# Get the reactor from the premade script.
reactor = large_core.make_large_core()

# View all the diagnostics down the chain.
reactor.display()

# Evaluate the default pattern.  We can grab the eigenvalue
# and peaking as return values.
k, p = reactor.evaluate()

print "k = ",k," p = ",p

# Alternatively, we can use print_params to display current
# values of all optimization parameters.  Currently only
# keff and the max peaking are retained.
reactor.print_params()

## We can also print the power peaking.
reactor.print_peaking()


# With this, we can try optimizing by hand a bit.  Peaking
# occurs at (0, 1).  Printing the pattern helps visualize this.
reactor.print_pattern()

# We can also see what fuel type is where by looking, for 
# example, at the burnup
reactor.print_map('burnup')
# or the enrichment
reactor.print_map('enrichment')

# Now, something that tends to work is to swap a peaking
# bundle with lower burnup with a lower peaking bundle
# with higher burnup.  Let's switch the peaker with
# its 15 GWd/MTU neighbor at [0,2].  Then print and
# and evaluate.
reactor.swap([0,1],[0,2])
reactor.print_pattern()
reactor.evaluate()
reactor.print_params()
reactor.print_peaking()

# That's a slight peaking reduction with ja slight increase
# in keff.  However, there is a better pattern.  Try the
# "ring of fire":
pattern = np.array([48,36,5,6,19,23,17,40,3,10,15,25,32,1,44,7,9,18,33,31,8,43,11,20,26,\
                    24,21,16,35,27,28,29,30,12,41,34,22,13,2,45,37,14,0,4,42,47,46,38,39])
reactor.shuffle(pattern)
reactor.evaluate()
reactor.print_params()
reactor.print_peaking()
reactor.plot_peaking()

