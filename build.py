#!/usr/bin/env python

import os, sys


build_param_dict = {
    "debug": {
        "path": "debug",
        "cmake_options": " -DCMAKE_BUILD_TYPE=Debug -DMINMD_DOUBLE=1 ",
    },
    "release": {
        "path": "build",
        "cmake_options": " ",
    },
}


if len(sys.argv) > 1 and sys.argv[1] in ("debug", "dbg"):
    build_type = "debug"
else:
    build_type = "release"

build_param = build_param_dict[build_type]



print("Building %s..." % build_type)

script_path = os.path.realpath(__file__)
root_path = os.path.dirname(script_path)
build_path = os.path.join(root_path, build_param["path"])
if not os.path.exists(build_path):
    os.makedirs(build_path)
os.chdir(build_path)
print(build_path, os.path.exists(build_path), os.getcwd())
curr_path = os.path.realpath(os.getcwd())
if curr_path != build_path:
    exit(1)

commands = [
        "rm -rf *",
        "cmake %s .." % build_param["cmake_options"],
        "cmake --build ."
]

for cmd in commands:
    ret = os.system(cmd)
    print('build.py: "%s" returns %s' % (cmd, ret))

