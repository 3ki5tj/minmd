#!/usr/bin/env python

import os

script_path = os.path.realpath(__file__)
root_path = os.path.dirname(script_path)
build_path = os.path.join(root_path, "build")
if not os.path.exists(build_path):
    os.path.makedirs(build_path)
os.chdir(build_path)
print(build_path, os.path.exists(build_path), os.getcwd())
curr_path = os.path.realpath(os.getcwd())
if curr_path != build_path:
    exit(1)

commands = [
        "rm -rf *",
        "cmake ..",
        "cmake --build ."
]

for cmd in commands:
    ret = os.system(cmd)
    print('build.py: "%s" returns %s' % (cmd, ret))

