"""Print the preserved Chiou MATLAB u and v arrays for the R import fallback."""
import sys

import numpy as np
import scipy.io

state = scipy.io.loadmat(sys.argv[1], squeeze_me=True)
for value in np.concatenate((np.ravel(state["u"]), np.ravel(state["v"]))):
    print(repr(float(value)))
