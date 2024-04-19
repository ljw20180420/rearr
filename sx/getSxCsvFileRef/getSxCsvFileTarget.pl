#!/usr/bin/env -S perl -anF,

# Usage: getSxCsvFileTarget.pl csvfile
# Extract target from sx csvfile

use feature 'say';

$rev = scalar reverse $F[1];
$rev =~ m/[acgt]/g;
$target = substr($F[1], length($F[1]) - pos($rev) + 1, 44);
substr($target, 16, 2) = "CC";
say $target;