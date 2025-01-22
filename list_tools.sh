#!/bin/bash

data_dir=$HOME/lab-data
aqua_dir=$HOME/aqua_tools

echo "List of AQuA tools on this Tinkerbox"
echo
echo "---------------"
echo
ls -1 $aqua_dir/*.sh | tr '\n' '\0' | xargs -0 -n 1 basename
echo
echo
echo "Type the tool name from anywhere followed by '-h' to learn more"
