#!/usr/bin/env bash

for x in 3 2 1
do
    echo ${x}
    sleep 1
done
echo Done


if [ x$1 == "xdodgy" ]
then
    exit 13
else
    exit 0
fi
