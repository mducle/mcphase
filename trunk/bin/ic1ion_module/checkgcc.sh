#!/bin/sh

if test ! -x checkgcc
then
   $2 -static-libgcc checkgcc.cpp -o checkgcc 1>/dev/null 2>/dev/null
fi

./checkgcc $1
if [ "$?" -ne "0" ]
then
  if [ $# -gt 2 ]
  then
    echo Sorry, you version of GCC is older than $1.
    echo ic1ion produces inaccurate numerical results with these compilers.
    echo So you must upgrade your compiler to a newer version.
    exit 1
  fi
fi

