# examples/example01.py

import small_core

# Here, we'll investigate the small benchmark core.  First,
# build it using the premade script.
reactor = small_core.make_small_core()

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

# We can also print the power peaking.
reactor.print_peaking()

# With this, we can try optimizing by hand a bit.  Peaking
# occurs at (0, 1).  Printing the pattern helps visualize this.
reactor.print_pattern()

# Do a swap and evaluate.
reactor.swap([0,1],[0,2])
reactor.print_pattern()
reactor.evaluate()
reactor.print_params()
reactor.print_peaking()

# That's a significant peaking reduction with just a slight decrease
# in keff.  However, there is a better pattern.  For this keff, 
# the tradeoff curve in the theory document suggests a peaking 
# on the order of 1.75 or maybe less.  But we won't find out
# by hand...
