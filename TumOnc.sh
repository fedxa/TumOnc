#!/bin/sh

export PYTHONPATH=`dirname $0`:$PYTHONPATH

python3 -m TumOnc $*
