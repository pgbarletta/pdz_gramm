#!/bin/bash

gfortran -O2 -shared -fPIC apc.f -o apc_wrapper.so

exit 0
