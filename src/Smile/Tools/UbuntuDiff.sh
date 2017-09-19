#!/bin/bash
while true; do
  ./1001-rand > s.in
  ./1001-ans < s.in > s-ans.out
  ./1001-me < s.in > s.out
  diff s.out s-ans.out
  if [ $? -ne 0 ]; then 
    echo "=============" 
    break;
  else
    echo "No difference." 
  fi
done

