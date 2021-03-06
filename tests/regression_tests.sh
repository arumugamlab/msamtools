#! /bin/bash

# Sourced functions

script_location=$(dirname "$0")
source $script_location/functions.sh

#################################
# Get list of commands and run
# pairwise test on them
#################################

function run_pairwise_list() {
  local get_cmd_function=$1;
  local infile=$2;
  local tnum=0
  eval "$get_cmd_function $infile" | while read command
  do
        tnum=$(expr $tnum + 1);
        run_pairwise_test "$STABLE_EXE" "$DEVEL_EXE" "$command" "$verbose" "$tnum";
  done
}
  
#################################
# Get list of commands and run
# valgrind test on them
#################################

function run_valgrind_list() {
  local get_cmd_function=$1;
  local infile=$2;
  local tnum=0
  eval "$get_cmd_function $infile" | while read command
  do
        tnum=$(expr $tnum + 1);
        run_valgrind_check "$DEVEL_EXE" "$command" "$verbose";
  done
}

#for i in msam_*.c; do gcc -g -W -Wall -pedantic -DPROGRAM=\"msamtools\" -DBUILD=\"\" -DPACKAGE_NAME=\"maui\" -DPACKAGE_VERSION=\"2.0\" -c $i -o ${i%%.c}.o; done; gcc -g -W -Wall -static -pedantic -DPROGRAM=\"msamtools\" -DBUILD=\"\" -DPACKAGE_NAME=\"maui\" -DPACKAGE_VERSION=\"2.0\" -o test msam_*.o libbam.a libmaui.a libargtable2.a -lm -lz
#./test

function parse_cmdline() {
  local o
  while getopts "iacfpxv:q" o; do
    case "${o}" in
      i)
        echo "# Integrity tests = ON"
        integrity=1
      ;;
      a)
        echo "# all tests = ON"
        integrity=1
        coverage=1
        filter=1
        complex=1
        profile=1
      ;;
      c)
        echo "# coverage = ON"
        coverage=1
      ;;
      f)
        echo "# filter = ON"
        filter=1
      ;;
      p)
        echo "# profile = ON"
        profile=1
      ;;
      x)
        echo "# complex = ON"
        complex=1
      ;;
      v)
        verbose=${OPTARG}
        echo "# verbose = $verbose"
      ;;
      q)
        echo "# quick = ON"
        quick=1
      ;;
    esac
  done
}

#########
# Main
#########

# Set default options
verbose=0
integrity=0
coverage=0
filter=0
profile=0
complex=0
quick=0

# Parse command line
parse_cmdline $*

# Versions to compare
STABLE_EXE="$HOME/src/git/msamtools/msamtools.20200427";
DEVEL_EXE="$HOME/src/git/msamtools/msamtools";

# Files
small_file="input.uniq.bam"
tiny_file="tiny_aln.bam"
if [ "$quick" == "1" ]; then
  small_file=$tiny_file
fi

# Start the tests

#message "$hline" "$verbose"

for mode in filter coverage profile complex; 
do
  # Check if this mode should be run
  if [ "${!mode}" == "0" ]; then
    continue;
  fi
  echo "******************************"
  echo "* Testing $mode mode *"
  echo "******************************"

  # Run regression test
  echo "Regression tests:"
  echo "-----------------"
  run_pairwise_list "get_"${mode}"_commands" "$small_file"

  # Run integrity test if needed, but skip for complex mode
  if [ "$integrity" == "1" -a "$mode" != "complex" ]; then
    echo "Integrity tests:"
    echo "----------------"
    run_valgrind_list "get_"${mode}"_commands" "$tiny_file"
  fi

  echo ""
done
