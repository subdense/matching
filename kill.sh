ps aux | grep profiling | awk -F" " '{print "kill -9 "$2}'| sh

