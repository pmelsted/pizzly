import h5py
import numpy as np
import sys


def get_cumulative_dist(fn):
    f = h5py.File(fn)
    x = np.asarray(f['aux']['fld'], dtype='float64')
    y = np.cumsum(x)/np.sum(x)
    f.close()
    return y

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python get_fragment_length.py H5FILE [cutoff]")
        print("")
        print("Prints 95 percentile fragment length")
        print("")
        print("       H5FILE:      abundance.h5 file produced by kallisto")
        print("       cutoff:      percentage cutoff to use (default .95)")
    else:
        fn = sys.argv[1]
        if len(sys.argv) >=3:
            cutoff = float(sys.argv[2])
            if cutoff <= 0 or cutoff >= 1.0:
                cutoff = 0.95
        else:
            cutoff = 0.95
        y = get_cumulative_dist(fn)
        fraglen = np.argmax(y > .95)
        print(fraglen)
