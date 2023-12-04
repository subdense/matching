#!/bin/sh
#PROCESS_COMMAND=$1 # pb quotes?
# works assuming single command running

rm profiling.log

while true
do
  time=`date +%s`
  #ps aux | grep "$PROCESS_COMMAND" | grep -v "grep" | awk -F" " '{printf "%.0f", $5}' >> profiling.log
  #ps aux | grep "python run.py" | grep -v "grep" | grep -v "sh"
  vsz=`ps aux | grep "python run.py" | grep -v "grep" | grep -v "sh" | awk -F" " '{print $5}'`
  rss=`ps aux | grep "python run.py" | grep -v "grep" | grep -v "sh" | awk -F" " '{print $6}'`
  np=`ps aux | grep "python run.py" | grep -v "grep" | grep -v "sh" | wc -l`
  echo "$time;$vsz;$rss;$np" >> profiling.log
  sleep 1
done
