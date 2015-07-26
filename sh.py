#!/usr/bin/env python
# coding: utf-8

"""
   File Name: sh.py
      Author: Wan Ji
      E-mail: wanji@live.com
  Created on: Sun Jul 26 21:16:30 2015 CST
"""
DESCRIPTION = """
"""

import os
import argparse
import logging

import numpy as np


def eigs(X, npca):
    l, pc = np.linalg.eig(X)
    idx = l.argsort()[::-1][:npca]
    return pc[:, idx], l[idx]


def trainSH(X, nbits):
    #
    # Input
    #   X = features matrix [Nsamples, Nfeatures]
    #   param.nbits = number of bits (nbits do not need to be a multiple of 8)
    #
    #
    # Spectral Hashing
    # Y. Weiss, A. Torralba, R. Fergus.
    # Advances in Neural Information Processing Systems, 2008.

    [Nsamples, Ndim] = X.shape
    SHparam = {'nbits': nbits}

    # algo:
    # 1) PCA
    npca = min(nbits, Ndim)
    pc, l = eigs(np.cov(X.T), npca)
    pc[:, 0] *= -1
    X = X.dot(pc)   # no need to remove the mean

    # 2) fit uniform distribution
    eps = np.finfo(float).eps
    mn = np.percentile(X, 5)
    mx = np.percentile(X, 95)
    mn = X.min(0) - eps
    mx = X.max(0) + eps

    # 3) enumerate eigenfunctions
    R = mx - mn
    maxMode = np.ceil((nbits+1) * R / R.max())
    nModes = maxMode.sum() - maxMode.size + 1
    modes = np.ones((nModes, npca))
    m = 0
    for i in xrange(npca):
        modes[m+1:m+maxMode[i], i] = np.arange(1, maxMode[i]) + 1
        m = m + maxMode[i] - 1
    modes = modes - 1
    omega0 = np.pi / R
    omegas = modes * omega0.reshape(1, -1).repeat(nModes, 0)
    eigVal = -(omegas ** 2).sum(1)
    ii = (-eigVal).argsort()
    modes = modes[ii[1:nbits+1], :]

    SHparam['pc'] = pc
    SHparam['mn'] = mn
    SHparam['mx'] = mx
    SHparam['modes'] = modes
    return SHparam


def runcmd(cmd):
    """ Run command.
    """

    logging.info("%s" % cmd)
    os.system(cmd)


def getargs():
    """ Parse program arguments.
    """

    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     formatter_class=
                                     argparse.RawTextHelpFormatter)
    # Example:
    #   parser.add_argument('ipaddr', type=str, default="127.0.0.1",
    #       help='IP address of server')
    parser.add_argument("--log", type=str, default="INFO",
                        help="log level")

    return parser.parse_args()


def main(args):
    """ Main entry.
    """

    logging.info("Hello INFO")
    logging.warn("Hello WARN")
    logging.error("Hello ERROR")
    logging.fatal("Hello FATAL")

if __name__ == '__main__':
    args = getargs()
    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: " + args.log)
    logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s",
                        level=numeric_level)
    main(args)
