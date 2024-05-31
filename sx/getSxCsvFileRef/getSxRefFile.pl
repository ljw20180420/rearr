#!/usr/bin/env -S perl -n

# Usage: getSxRefFile.pl <ref12.fa
# ref1 and ref2 appears alternatively in fasta form. This script transforms it a format ready to feed to rearrangement

BEGIN {
    $ext1up = $ARGV[0];
    $ext2up = $ARGV[1];
    $AG = $ARGV[2];
    for (my $i = 0; $i < 3; $i++){
        delete $ARGV[$i];
    }
}

if ($AG eq "A"){
    if ($. % 2){
        $spos = $ext1up + 4;
    }
    else{
        $spos = $ext2up + 4;
    }
    if (length($_) > $spos + 1){
        if (length($_) == $spos + 2){
            substr($_, $spos, 1) = "A";
        }
        else{
            substr($_, $spos, 2) = "AA";
        }
    }
}
if ($. % 2){
    chomp;
    print 0, "\t", $_, "\t", $ext1up;
}
else{
    chomp;
    print "\t", $ext2up, "\t", $_, "\t", length($_), "\n";
}