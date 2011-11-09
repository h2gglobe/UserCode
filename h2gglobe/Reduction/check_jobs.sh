#!/bin/bash


if [[ -z $1 ]]; then
    echo "usage:  $0 <directory with log files>"
    exit 1
fi

tail -n 100 $(ls --color=none $1/*.log | sort ) | egrep -B 2 -A 3 '=>|rfcp|^[Rr]ed|Break' | egrep  '=>|rfcp|bytes in remote file|^[Rr]ed|Break'

