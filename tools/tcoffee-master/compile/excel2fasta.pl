#!/usr/bin/env perl

open (F, $ARGV[0]);

while ( <F>)
  {
    @l=($_=~/(\S+)/g);
    
    $name=shift @l;
    
    print STDOUT "\n>$name\n";
    foreach $e (@l){$e=($e eq "0")?"O":"I";print "$e";}
  }
close (F);

		       
    
# Tue 13 Aug 2019 14:24:57 +0430
