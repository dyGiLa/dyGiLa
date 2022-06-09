#!/bin/bash

# This is the shell script for managing the compling, running and vatulization of phase diagram code.
# pyhton3 version is 3.8+, gcc version 9.4+, numpy and matplotlib are nessary.
#
# the *csv files will be written into ./pd folder.

g++ -std=c++17 -g -o pd.app phases_diagram.cpp matep.cpp matep.hpp && ./pd.app

python3 ./pd/sfdata.py


