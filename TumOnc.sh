#!/bin/sh

cdir=`dirname $0`

if [ -f $cdir/TumOnc ]; then
    exec $cdir/TumOnc "$@"
else
    export PYTHONPATH=$cdir:$PYTHONPATH
    exec python3 -m TumOnc "$@"
fi
