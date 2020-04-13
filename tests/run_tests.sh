#!/usr/bin/env bash

# Select names from the currently available tests
TESTS="
starY16_dAu200
starY15_pp200a
starY15_pp200b
starY14_AuAu200a
starY14_AuAu200b
starY14_He3Au200
starY12_pp500
starY11_pp500"
#starY12_CuAu200
#starY10_AuAu62"

# Specific command line arguments for the tests
#
#     $ test_tpcrs <test name> [number of events] [electron cut off energy]
#
# Set the [number of events] to a negative number to process all available
# events in the corresponding test file.
ARGS_starY16_dAu200="  2       "
ARGS_starY15_pp200a="  2       "
ARGS_starY15_pp200b="  1       "
ARGS_starY14_AuAu200a="1       "
ARGS_starY14_AuAu200b="2 0.0001"
ARGS_starY14_He3Au200="1       "
ARGS_starY12_CuAu200=" 1       "
ARGS_starY12_pp500="   1       "
ARGS_starY11_pp500="   3 0.0001"
ARGS_starY10_AuAu62="  1 0.0001"

EXIT_CODE=0

handle_chld() {
    local still_running_pids=""
    for pid in $pids; do
        if [ ! -d /proc/$pid ]; then # job is not running
            wait $pid; pid_exit_code=$?; test_id=PID_$pid
            if [ $pid_exit_code -eq 0 ]
            then
                echo -n "$pid: ${!test_id} SUCCESS ($pid_exit_code) "
            else
                echo -n "$pid: ${!test_id} FAILURE ($pid_exit_code) "
            fi
            cat "$dir/$pid"
        else still_running_pids+=" $pid"
        fi
    done
    pids="$still_running_pids"
}

set -o monitor
trap handle_chld CHLD

dir=$(mktemp -d)
export TIME="real %E user %U sys %S"
pids=""

for test_id in $TESTS
do
    test_args=ARGS_$test_id
    (command time -o $dir/$BASHPID ./tests/test_tpcrs $test_id ${!test_args} &> log_${test_id} && diff ${test_id}_inp.log ${test_id}_out.log &> /dev/null) &
    pid=$!
    pids+=" $pid"
    declare PID_$pid=$test_id
    echo "$pid: $test_id"
    sleep 1
done

wait
rm -r "$dir"
exit $EXIT_CODE
