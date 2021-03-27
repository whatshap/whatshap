use strict;

sub jobExtractReads($$$$$$) {
    my ($env, $cfg, $name, $reads, $extrReads, $baseSize) = @_;

    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};
    my $binPath = %$env{"BinPath"};
   
    my $threads = %$cfg{"THREADS"};
    my $genomeSize = %$cfg{"GENOME_SIZE"};
    my $coverage = %$cfg{"CNS_OUTPUT_COVERAGE"};

    # check exstrReads is gz file
    my $isGz = 0;
    my $extrReadsMid = $extrReads;
    if($extrReads =~ /.gz$/) {
        $isGz = 1;
        $extrReadsMid =~ s/.gz$//;
    }

    my @cmds = ();
    push @cmds, "$binPath/fsa_rd_tools longest --ifname=$reads --ofname=$extrReadsMid --base_size=$baseSize";

    if ($isGz) {
        push @cmds, "$binPath/pigz -f  -p $threads $extrReadsMid";
    }

    my $job = Job->new(
        name => "${name}_extract",
        ifiles => [$reads],
        ofiles => [$extrReads],
        gfiles => [$extrReads],
        mfiles => [],
        cmds => [@cmds],
        msg => "extracting reads for $name",
    );

    return $job;

}


1;