#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

print "sleeping ".$ARGV[0]."s";

sleep $ARGV[0];



