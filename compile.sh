#!/bin/bash

#########################################################################
# File Name: compile.sh
# Author: Wan Ji
# mail: wanji@live.com
# Created Time: 2015年07月27日 星期一 08时29分31秒
#########################################################################

gcc -shared -o bitsop.so -O3 -fPIC bitsop.c
