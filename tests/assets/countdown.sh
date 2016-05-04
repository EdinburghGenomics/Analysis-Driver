#!/usr/bin/env bash
for x in 2 1 Done; do echo ${x}; done
if [ x$1 == "xdodgy" ]; then exit 13; else exit 0; fi
