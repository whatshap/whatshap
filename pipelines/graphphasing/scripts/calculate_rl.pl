#!/usr/bin/perl -w
#############################
###   Jennifer Meneghin   ###
###   August 27, 2009     ###
#############################
#-----------------------------------------------------------------------------------------------------------------------------------------
#Deal with passed parameters
#-----------------------------------------------------------------------------------------------------------------------------------------
if ($#ARGV == -1) {
    &usage;
}
$db_file = $ARGV[0];
unless ( open(DB, "$db_file") ) {
    print "\nCouldn't read fasta file: $db_file\n\n";
    &usage;
}
print "Parameters: database file = $db_file\n\n";
#-----------------------------------------------------------------------------------------------------------------------------------------
#The main event
#-----------------------------------------------------------------------------------------------------------------------------------------
%fasta_lengths = ();
$seq = "";
$sum = 0;
$count = 0;
$min = 9999;
$max = 0;
while (<DB>) {
    chomp;
    if (/^>/) {
	#finish up previous line.
	if (length($seq) > 0) {
	    $fasta_lengths{$id} = length($seq);
#	    print "$fasta_lengths{$id}...";
	    $sum = $sum + $fasta_lengths{$id};
	    $count++;    
	    if (length($seq) < $min) {
		$min = length($seq);
	    }
	    if (length($seq) > $max) {
		$max = length($seq);
	    }
	}
	$id = $_;
	$id =~ s/^>(.+?)\s.+$/$1/g;
	$seq = "";
    }
    else {
	$seq = $seq . $_;
    }
}
$fasta_lengths{$id} = length($seq);
#print "$fasta_lengths{$id}...";
$sum = $sum + $fasta_lengths{$id};
$count++;    
if ($count > 0) {
    $avg = $sum / $count;
    print "\n\nAVERAGE SEQUENCE LENGTH = $avg\n";
    print "TOTAL NUMBER OF BP = $sum\n";
    print "TOTAL NUMBER OF RECORDS = $count\n";
    print "MINIMUM SEQUENCE LENGTH = $min\n";
    print "MAXIMUM SEQUENCE LENGTH = $max\n\n";
}
else {
    print "ERROR: Couldn't find any sequences in the file.\n";
}
close(DB);
sub usage {
    print "\nUsage: seq_average_length.pl your_fasta_file\n";
    print "\nThis program reads a fasta file, then computes and reports the average sequence length,\n";
    print "the total number of bases, the total number of records, the minimum sequence length,\n";
    print "and the maximum sequence length.\n\n";
    print "Jennifer Meneghin\n";
    print "August 27, 2009\n";
    print "Updated October 28, 2010 JMeneghin\n\n";
    print "Usage: seq_average_length.pl <fasta file>\n\n";
    exit;
}

