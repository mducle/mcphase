#!/usr/bin/perl
# batch to remove every nth line (n>2) in output.txt
# usage: nth_rm n < input.txt > output.txt 

# first~step
#               Match  every step'th line starting with line first.
#                 For example, ``sed -n 1~2p''  will  print  all  the
#               odd-numbered  lines  in  the  input stream, and the
#                 address 2~5 will match every fifth  line,  starting
#               with the second.  (This is a GNU extension.)
								                                   

$n=$ARGV[0];shift @ARGV;

while(<>){
$l=$_;
for($i=1;$i<$n;++$i){print $l;$l=<>;}
}

