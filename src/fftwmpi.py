#!/usr/bin/env python3

import h5py as h5
import mytools as mt
import matplotlib.pyplot as plt
import numpy as np

print("\n\n===========================\nPython here\n===========================\n\n")

filename = 'output.h5'
f = h5.File(filename, 'r')
this = mt.H5Drill(f)
tval = this.t.value
xval = this.x.value
yval = this.y.value
rval = this.r.value
for i, value in enumerate(tval[0:11]):
    print("t: {:0.3f}\tx: {:0.3f}\tr: {:0.3f}\ty: {:0.3f}".format(
            tval[i],
            xval[i],
            rval[i],
            yval[i],
            )
        )


f = 1000.0/2501.0*np.linspace(0, 1250, 1251)
plt.plot(f, this.out_final.value[0:1251])
plt.show()
