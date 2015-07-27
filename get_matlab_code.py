#!/usr/bin/env python
# coding: utf-8

"""
   File Name: get_matlab_code.py
      Author: Wan Ji
      E-mail: wanji@live.com
  Created on: Mon Jul 27 07:59:56 2015 CST
"""
SH_URL = "http://www.cs.huji.ac.il/~yweiss/SpectralHashing/sh.zip"
DESCRIPTION = """
This script is for downloading Y. Weiss's matlab code from:
    %s
""" % SH_URL

import os
import sys
import tempfile
import argparse
import logging
import requests


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
    parser.add_argument('target', type=str, nargs='?', default="matlab",
                        help='directory for the extracted matlab code')
    parser.add_argument("--log", type=str, default="INFO",
                        help="log level")

    return parser.parse_args()


def main(args):
    """ Main entry.
    """

    if os.path.exists(args.target):
        logging.error("Target directory already exists: %s" % args.target)
        sys.exit(1)

    resp = requests.get(SH_URL, timeout=10.0)
    if resp.status_code != 200:
        logging.error("Download failed!")
        sys.exit(1)

    tmp = tempfile.mktemp()
    with open(tmp, "wb") as tmpf:
        tmpf.write(resp.content)
    os.makedirs(args.target)
    runcmd("cd %s; unzip %s; cd -" % (args.target, tmp))
    runcmd("cd %s;" % args.target +
           "git init;" +
           "git config user.name 'Y. Weiss';" +
           "git config user.email yweiss@cs.huji.ac.il;" +
           "git add .; git commit -m \"Y. Weiss's matlab code\";" +
           "cd -")

if __name__ == '__main__':
    args = getargs()
    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: " + args.log)
    logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s",
                        level=numeric_level)
    main(args)
