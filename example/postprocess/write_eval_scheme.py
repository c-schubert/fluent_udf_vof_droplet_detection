#!/usr/bin/env python3

from os import listdir
from os.path import isfile, join

mypath="F:\\myfolder"
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

for fn in onlyfiles:
    if fn.endswith(".dat"):
         print("file/read-data/"+fn)
         print("ok")
         print("")
         print("/define/user-defined/execute-on-demand","\"CPAD_oD::libudf\"")
