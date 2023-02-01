#! /bin/bash

#################################
# Helper functions for all modes
#################################

function message() {
  local string=$1;
  local verbose=$2;
  if [ $verbose -gt 0 ]; then
    echo $string;
  fi
}

function report_status () {
  local status=$1;
  if [ $status -eq 0 ]; then
    echo "PASS";
  else
    echo "FAIL!";
  fi
}

function run_command() {

  # Get parameters
  local command=$1;
  local verbose=$2;
  local timeit=$3;

  # Run command
  #echo -e "\tPROG2:\t"$({ time eval $command; } 2>&1)
  if [ $verbose -gt 1 ]; then
    echo $command
  fi
  if [ ! -z "$timeit" ]; then
    (time eval $command) 2>/tmp/msamtools.time.out
  else
    eval $command
  fi
  local exit_code=$?

  # Check for error in exit code
  if [ $exit_code -ne 0 ]; then
    message "ERROR: Exit code $exit_code" "$verbose"
    return $exit_code
  else
    message "OK" "$verbose"
  fi

  # Print output
  if [ $verbose -gt 1 ]; then
    echo -e $(cat /tmp/msamtools.time.out)
  fi

  # Return success
  return 0
}

#################################
# Run valgrind for given params
#################################

function run_valgrind_check() {

  # Get parameters

  local program=$1;
  local in_params=$2;
  local verbose=$3;
  local tnum=$4;

  # Output file names
  local old_file=/tmp/msamtools.old.out
  local val_file=/tmp/msamtools.valgrind.out
  rm -rf $old_file $val_file

  ##########
  # Run program
  ##########

  # Form the command
  # We add 2>/dev/null as we catch the exit code.
  local params=$(echo $in_params | sed "s|__OUTFILE__|${old_file}|" | sed "s^__PROGRAM__^${program}^g");
  local command="valgrind --tool=memcheck --leak-check=full --show-reachable=yes --show-error-list=yes $params 2>$val_file";

  echo -n "Test $tnum: ";

  # Run command

  if [ $verbose -gt 1 ]; then
    echo -n "COMMAND: "
  fi
  run_command "$command" "$verbose" ""

  local status=1;
  if [ "$(grep ERROR $val_file | cut -f2- -d' ')" = "ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)" ]; then
    status=0
  fi
  message "Test $tnum: " "$verbose";
  report_status $status;
  message "$hline" "$verbose";

  return $exit_code
}

####################################
# Run two programs for given params
####################################

function compare_two_versions() {

  # Get parameters

  local program1=$1;
  local program2=$2;
  local in_params=$3;
  local verbose=$4;

  # Output file names
  local old_file=/tmp/msamtools.old.out
  local new_file=/tmp/msamtools.new.out;
  rm -rf $old_file $new_file; 

  ##########
  # Run old
  ##########

  # Form the command
  # We add 2>/dev/null as we catch the exit code.
  local params=$(echo $in_params | sed "s|__OUTFILE__|${old_file}|" | sed "s^__PROGRAM__^${program1}^g");
  local command="$params 2>/dev/null";

  # Run command

  if [ $verbose -gt 0 ]; then
    echo -n "PROG1: "
  fi
  run_command "$command" "$verbose" "time-me"
  local exit_code1=$?

  ##########
  # Run new
  ##########

  # Form the command
  params=$(echo $in_params | sed "s|__OUTFILE__|${new_file}|" | sed "s^__PROGRAM__^${program2}^g");
  command="$params 2>/dev/null";

  # Run command

  if [ $verbose -gt 0 ]; then
    echo -n "PROG2: "
  fi
  run_command "$command" "$verbose" "time-me"
  local exit_code2=$?

  # Compare
  # Text output
  if [ ! -z "$(echo $command | grep ' -bu ')" ]; then
    local viewer="samtools view";
    local digest1=$(eval "$viewer $old_file" | md5sum - | cut -f1 -d' ');
    local digest2=$(eval "$viewer $new_file" | md5sum - | cut -f1 -d' ');
  else
    local viewer="grep -Ev '^Unknown|^#'";
    if [ ! -z "$(file $old_file | grep 'gzip')" ]; then viewer="zgrep -Ev '^Unknown|^#'"; fi
    local digest1=$(eval "$viewer $old_file" | md5sum - | cut -f1 -d' ');
    local viewer="grep -Ev '^Unknown|^#'";
    if [ ! -z "$(file $new_file | grep 'gzip')" ]; then viewer="zgrep -Ev '^Unknown|^#'"; fi
    local digest2=$(eval "$viewer $new_file" | md5sum - | cut -f1 -d' ');
  fi

  if [ $verbose -gt 1 ]; then
    echo "$old_file - $digest1"
    echo "$new_file - $digest2"
  fi

  if [ "$digest1" != "$digest2" ]; then
    return 255
  else
    return $(expr $exit_code1 + $exit_code2)
  fi
}

function run_pairwise_test() {
  local prog1=$1;
  local prog2=$2;
  local command=$3;
  local verbose=$4;
  local tnum=$5;

  local hline="--------------------------------------------------------------------------------";

  echo -n "Test $tnum: ";
  local viewer="grep -Ev '^Unknown|^#'";
  if [ ! -z "$(echo $command | grep ' --gzip ')" ]; then viewer="zgrep -Ev '^Unknown|^#'"; fi
  compare_two_versions "$prog1" "$prog2" "$command" "$verbose";
  #echo "$prog1" "$prog2" "$command" "$viewer" "$verbose";
  local status=$?;
  message "Test $tnum: " "$verbose";
  report_status $status;
  message "$hline" "$verbose";
}

#################################
# Generate commands to be tested
#################################

function get_filter_commands() {
  local infile=$1;
  for l in 30 45; do
    if [ "$l" == "30" ]; then p=90; else p=95; fi
    for z in "" "-z 80" "-z 90"; do
      for special in "-b" "--besthit" "--uniqhit"; do
        command="__PROGRAM__ filter -l $l -p $p $z $special $infile > __OUTFILE__";
        echo "$command";
      done
    done
  done
}

function get_profile_commands() {
  local infile=$1;
  for total in "--total=60000" ""; do
      for multi in "all" "equal" "prop" "ignore"; do
        for unit in "" "--unit=rel" "--unit=ab" "--unit=tpm" "--unit=fpkm"; do
          for mincount in "" "--mincount=10"; do
            command="__PROGRAM__ profile --label test --multi=$multi $total $unit $mincount -o __OUTFILE__ $infile";
            echo "$command";
          done
        done
      done
  done
}

function get_coverage_commands() {
  local infile=$1;
  for summary in "" "--summary"; do
    command="__PROGRAM__ coverage $summary --gzip --skipuncovered -o __OUTFILE__ $infile";
    echo $command;
  done
}
  
function get_complex_commands() {
  local infile=$1;
  echo "__PROGRAM__ filter -l 45 -p 95 -z 90 -v $infile > __OUTFILE__";
  echo "__PROGRAM__ filter -l 10 -v $infile > __OUTFILE__";
  echo "__PROGRAM__ filter -l 45 --ppt 995 -z 90 -v $infile > __OUTFILE__";
  echo "__PROGRAM__ summary $infile > __OUTFILE__";
  echo "__PROGRAM__ filter -b -u -l 45 -p 95 -z 90 $infile | __PROGRAM__ profile --label test --multi=proportional --total=60000 -o __OUTFILE__ -";
  #echo "fixmate -i 200 input.bam > __OUTFILE__" \
}
