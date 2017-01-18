#!/bin/bash

for file in *.txt; do
  base=`echo "${file%.txt}"`
  mv -- "${file}" "${base}.png"
done
