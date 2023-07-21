#!/bin/bash

awk '/^S/{print ">"$2"\n"$3}' $1 > `basename $1`.fasta
