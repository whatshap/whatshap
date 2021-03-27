#!/bin/bash

echo $#

if [ $1 -eq 1 ]; then
	echo parameter is true
fi

if [ $1 -eq 0 ]; then
	echo paramter is false
fi
