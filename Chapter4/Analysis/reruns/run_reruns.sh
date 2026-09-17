#!/bin/zsh
# Run ch4_rerun_all one stage per MATLAB session, retrying crashed sessions.
#
#   Chapter4/Analysis/reruns/run_reruns.sh                 every stage
#   Chapter4/Analysis/reruns/run_reruns.sh l1 load         only these
#   MATLAB=/path/to/matlab Chapter4/Analysis/reruns/run_reruns.sh
#
# Why a wrapper: on the author's Mac, MATLAB -batch aborts intermittently inside
# libcurl, at startup or minutes in, and concurrent sessions hang. So each stage
# gets its own session, waits until no MATLAB has run for 30 s, and is retried
# up to 10 times; every stage resumes from its saved runs. A stage that ends
# with a MATLAB error (not a crash) stops the script.
#
# Progress: Results/reruns/status.txt; one log per attempt beside it.

here=${0:A:h}
root=${here:h:h:h}
ML=${MATLAB:-matlab}
if ! command -v "$ML" >/dev/null 2>&1 && [[ -x /Applications/MATLAB_R2021b.app/bin/matlab ]]; then
  ML=/Applications/MATLAB_R2021b.app/bin/matlab
fi
out=$root/Results/reruns
mkdir -p $out
STATUS=$out/status.txt

stages=("$@")
if (( ${#stages} == 0 )); then
  stages=(robust case4 l1 load ch3tests ch3table long_robust long_sweep long_nine long_oos long_summary)
fi

wait_idle() {
  local quiet=0
  while (( quiet < 30 )); do
    if pgrep -f MATLAB >/dev/null; then quiet=0; else quiet=$((quiet+5)); fi
    sleep 5
  done
}

cd "$root"
for st in $stages; do
  done_ok=0
  for attempt in {1..10}; do
    wait_idle
    log=$out/${st}_try${attempt}.log
    rm -f $log
    echo "$(date +%H:%M:%S) $st attempt $attempt start" >> $STATUS
    "$ML" -batch "addpath(genpath(pwd)); diary('$log'); try, ch4_rerun_all('$st'); catch ME, disp(getReport(ME)); end; disp('JOB_END'); diary off" > $out/${st}_try${attempt}.stdout 2>&1
    if grep -q "RERUN_STAGE_DONE $st" $log 2>/dev/null; then
      echo "$(date +%H:%M:%S) $st OK" >> $STATUS; done_ok=1; break
    fi
    if grep -q "JOB_END" $log 2>/dev/null; then
      echo "$(date +%H:%M:%S) $st FAILED with a MATLAB error, see $log" >> $STATUS
      exit 1
    fi
    echo "$(date +%H:%M:%S) $st attempt $attempt crashed, resuming" >> $STATUS
  done
  if (( ! done_ok )); then
    echo "$(date +%H:%M:%S) $st gave up after 10 attempts" >> $STATUS
    exit 1
  fi
done
echo "$(date +%H:%M:%S) ALL_DONE" >> $STATUS
