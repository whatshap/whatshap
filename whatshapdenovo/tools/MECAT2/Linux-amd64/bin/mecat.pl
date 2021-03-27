#!/usr/bin/env perl

use FindBin;
use lib $FindBin::RealBin;

use Plgd::Utils;
use Plgd::Script;
use Plgd::Project;

use strict;

sub defaultConfig() {
    return (
        PROJECT=>"",
        RAWREADS=>"",
        GENOME_SIZE=>"",
        THREADS=>4,
        MIN_READ_LENGTH=>500,
        CNS_OVLP_OPTIONS=>"",
        CNS_OPTIONS=>"-r 0.6 -a 1000 -c 4 -l 2000",
        CNS_OUTPUT_COVERAGE=>30,
        TRIM_OVLP_OPTIONS=>"-B",
        ASM_OVLP_OPTIONS=>"-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400",
        CLEANUP=>0,
        USE_GRID=>"false",
        GRID_NODE=>0,
        FSA_OL_FILTER_OPTIONS=>"--max_overhang=-1 --min_identity=-1",
        FSA_ASSEMBLE_OPTIONS=>"",
    );
}

sub loadMecatConfig($) {
    my ($fname) = @_;
    my %cfg = defaultConfig();

    loadConfig($fname, \%cfg);

    my @required = ("PROJECT", "GENOME_SIZE", "RAWREADS");
    foreach my $r (@required) {
        if (not exists($cfg{$r}) or $cfg{$r} eq "")  {
            plgdError("Not set config $r");
        }
    }


    return %cfg;
}

sub loadMecatEnv($) {
    my ($cfg) = @_;

    my %env = loadEnv($cfg);
    $env{"BinPath"} = $FindBin::RealBin;
    return %env;    
}

sub initializeMecatProject($) {
    my ($cfg) = @_;
    initializeProject($cfg);
}

sub runCorrectRawreads($$) {
    my ($env, $cfg) = @_;

    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};
    my $workDir = "$prjDir/1-consensus";
    mkdir $workDir;
    my $rawreads = %$cfg{"RAWREADS"};
    
    my $thread = %$cfg{"THREADS"};
    my $genomeSize = %$cfg{"GENOME_SIZE"};
    my $coverage = %$cfg{"CNS_OUTPUT_COVERAGE"};
    my $binPath = %$env{"BinPath"};
    my $cnsOvlpOptions = %$cfg{'CNS_OVLP_OPTIONS'};
    my $cnsOptions = %$cfg{'CNS_OPTIONS'};
   
    #my $job = Job->new(
    #    name => "cns_rawreads",
    #    ifiles => [$rawreads],
    #    ofiles => ["$workDir/cns_final.fasta"],
    #    gfiles => ["$workDir/cns_final.fasta", "$workDir/cns_reads.fasta", "$workDir/cns_pm*"],
    #    mfiles => ["$workDir/cns_pm*", "$workDir/cns_reads.fasta", "$workDir/cns_final.fasta.qual", 
    #               "$workDir/cns_final.fasta.qv", "$workDir/cns_final.frg"],
    #    cmds => ["$binPath/mecat2pw -j 0 -d $rawreads -o $workDir/cns_pm.can -w $workDir/cns_pm_dir -t $thread $cnsOvlpOptions",
    #             "$binPath/mecat2cns -i 0 -t $thread $cnsOptions $workDir/cns_pm.can $rawreads $workDir/cns_reads.fasta",
    #             "$binPath/extract_sequences $workDir/cns_reads.fasta $workDir/cns_final $genomeSize $coverage"],
    #    msg => "correcting rawreads",
    #);

    
    my $jobPw = Job->new(
        name => "cns_pw",
        ifiles => [$rawreads],
        ofiles => ["$workDir/cns_pm.can"],
        gfiles => ["$workDir/cns_pm*"],
        mfiles => [],
        cmds => ["$binPath/mecat2pw -j 0 -d $rawreads -o $workDir/cns_pm.can -w $workDir/cns_pm_dir -t $thread $cnsOvlpOptions"],
        msg => "correcting rawreads step 1 mecat2pw",
    );

    my $jobCns = Job->new(
        name => "cns_cns",
        ifiles => ["$workDir/cns_pm.can"],
        ofiles => ["$workDir/cns_reads.fasta"],
        gfiles => ["$workDir/cns_reads.fasta"],
        mfiles => [],
        cmds => ["$binPath/mecat2cns -i 0 -t $thread $cnsOptions $workDir/cns_pm.can $rawreads $workDir/cns_reads.fasta"],
        msg => "correcting rawreads step 2 mecat2cns",
    );

    my $jobExtr = Job->new(
        name => "cns_extract",
        ifiles => ["$workDir/cns_reads.fasta"],
        ofiles => ["$workDir/cns_final.fasta"],
        gfiles => ["$workDir/cns_final.fasta"],
        mfiles => [],
        #cmds => ["$binPath/extract_sequences $workDir/cns_reads.fasta $workDir/cns_final $genomeSize $coverage"],
        cmds => ["$binPath/mecat2elr $workDir/cns_reads.fasta $genomeSize $coverage $workDir/cns_final.fasta"],
        msg => "correcting rawreads step 3 extract_sequences",
    );

    my $job = Job->new (
        name => "cns_job",
        ifiles => [$rawreads],
        ofiles => ["$workDir/cns_final.fasta"],
        mfiles => ["$workDir/cns_pm*", "$workDir/cns_reads.fasta", "$workDir/cns_final.fasta.qual", 
                   "$workDir/cns_final.fasta.qv", "$workDir/cns_final.frg"],
        jobs => [$jobPw, $jobCns, $jobExtr],        
        msg => "correcting rawreads",
    );
    
    serialRunJobs($env, $cfg, $job);
}


sub runTrimReads($$) {
    my ($env, $cfg) = @_;

    my $prjDir = %$env{"WorkPath"} . "/" .%$cfg{"PROJECT"};
    my $workDir = "$prjDir/2-trim_bases";
    mkdir $workDir;
    my $volDir = "$workDir/trim_pm_dir";
    mkdir $volDir;

    my $cnsReads = "$prjDir/1-consensus/cns_final.fasta";
    my $trimReads = "$prjDir/2-trim_bases/trimReads.fasta"; 
    my $trimPm = "$volDir/trim_pm.m4";
    my $binPath = %$env{"BinPath"};
    my $options = %$cfg{"TRIM_OVLP_OPTIONS"};
    my $thread = %$cfg{"THREADS"};

    my $jobMkVol = Job->new(
        name => "tr_mk_vol",
        ifiles => [$cnsReads],
        ofiles => ["$volDir/num_volumes.txt"],
        gfiles => ["$volDir/num_volumes.txt"],
        mfiles => [],
        cmds => ["$binPath/v2mkvol $volDir $cnsReads"],
        msg => "making vol for trimming",
    );


    my $jobAlVol = Job->new(
        prefunc => sub($) {
            my ($job) = @_;

            my $count = getFileFirstItem("$volDir/num_volumes.txt", 0);

            for (my $i=0; $i<$count; $i=$i+1) {
                my $id = $i + 1;
                my $jobSub = Job->new(
                    name => "tr_al_vol_$i", 
                    ifiles => ["$volDir/num_volumes.txt"],
                    ofiles => ["$volDir/pm_$id.m4"],
                    gfiles => ["$volDir/pm_$id.m4"],
                    mfiles => ["$volDir/${id}_*.r"],
                    cmds => ["$binPath/v2asmpm -P$volDir -T$thread -S$id -E$count $options",
                            "cat $volDir/${id}_*.r > $volDir/pm_$id.m4"],
                    msg => "aligning volumn $id for trimming",
                );
                
                push @{$job->pjobs}, $jobSub;
                push @{$job->ofiles}, "$volDir/pm_$id.m4";
            }
        },
        name => "tr_al_vol",
        ifiles => ["$volDir/num_volumes.txt"],
        ofiles => [],   # prefunc
        mfiles => [],
        pjobs => [],    # prefunc
        msg => "aligning volumn for trimming",
    );

    my $jobCatVol = Job->new (
        prefunc => sub ($) {
            my ($job) = @_;

            my $subPmStr = join(" ", @{$jobAlVol->ofiles});
            $job->ifiles($jobAlVol->ofiles);
            $job->cmds(["cat $subPmStr  > $trimPm"]);

        },
        name => "tr_cat_vol",
        ifiles => [],   # prefunc
        ofiles => [$trimPm],
        gfiles => [$trimPm],
        mfiles => [],
        cmds => [],     # prefunc
        msg => "catenating pm for trimming",
    );


    my $lcrResult = "$workDir/lcr.txt";
    my $srResult = "$workDir/sr.txt";

    my $jobTrimCore = Job-> new (
        name => "tr_trim_read",
        ifiles => [$cnsReads, $trimPm],
        ofiles => [$trimReads],
        gfiles => [$trimReads],
        mfiles => [],
        cmds => ["$binPath/v2pm4 $volDir $trimPm 0.09 $thread",
                 "$binPath/v2lcr $trimPm $volDir 0.09 1 1 500 $lcrResult $thread",
                 "$binPath/v2sr $trimPm $volDir $lcrResult 500 $srResult $thread",
                 "$binPath/v2tb $cnsReads  $srResult $trimReads"],
        msg => "trimming reads for trimming",
    );

    my $job = Job->new (
        name => "tr_job",
        ifiles => [$cnsReads],
        ofiles => [$trimReads],
        mfiles => ["$volDir/trimReads_*.fasta", $volDir],
        jobs => [$jobMkVol, $jobAlVol, $jobCatVol, $jobTrimCore],
        msg => "trimming corrected reads",
    );
    
    serialRunJobs($env, $cfg, $job);
}

sub runAlignTReads($$) {
    my ($env, $cfg) = @_;

    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};
    my $workDir = "$prjDir/3-assembly";
    mkdir $workDir;

    my $volDir = "$workDir/asm_pm_dir";
    mkdir $volDir;

    my $binPath = %$env{"BinPath"};
    my $trimReads = "$prjDir/2-trim_bases/trimReads.fasta";

    my $asmPm = "$workDir/asm_pm.m4";
    my $options = %$cfg{"ASM_OVLP_OPTIONS"};
    my $thread = %$cfg{"THREADS"};

    my $jobMkVol = Job->new(
        name => "altr_mk_vol",
        ifiles => [$trimReads],
        ofiles => ["$volDir/num_volumes.txt"],
        gfiles => ["$volDir/num_volumes.txt"],
        mfiles => [],
        cmds => ["$binPath/v2mkvol $volDir $trimReads"],
        msg => "making vol for aligning trimmed reads",
    );

    my $jobAlVol = Job->new(
        prefunc => sub($) {
            my ($job) = @_;
            my $count = getFileFirstItem("$volDir/num_volumes.txt", 0);

            for (my $i=0; $i<$count; $i=$i+1) {
                my $id = $i + 1;
                my $jobSub = Job->new(
                    name => "altr_al_vol_$i", 
                    ifiles => ["$volDir/num_volumes.txt"],
                    ofiles => ["$volDir/pm_$id.m4"],
                    gfiles => ["$volDir/pm_$id.m4"],
                    mfiles => ["$volDir/${id}_*.r"],
                    cmds => ["$binPath/v2asmpm -P$volDir -T$thread -S$id -E$count $options",
                            "cat $volDir/${id}_*.r > $volDir/pm_$id.m4"],
                    msg => "aligning volumn $id for assembling",
                );
                push @{$job->pjobs}, $jobSub;
                push @{$job->ofiles}, "$volDir/pm_$id.m4";
            }

        },
        name => "altr_al_vol",
        ifiles => ["$volDir/num_volumes.txt"],
        ofiles => [],   # prefunc
        mfiles => [],
        pjobs => [],     # prefunc
        msg => "aligning volumn for assembling",
    );

    my $jobCatVol = Job->new (
        prefunc => sub ($) {
            my ($job) = @_;
            my $subPmStr = join(" ", @{$jobAlVol->ofiles});
            $job->ifiles($jobAlVol->ofiles);
            $job->cmds(["cat $subPmStr  > $asmPm"]);
        },

        name => "altr_cat_vol",
        ifiles => [],   # prefunc
        ofiles => [$asmPm],
        gfiles => [$asmPm],
        mfiles => [],
        cmds => [],     #prefunc
        msg => "catenating pm for assembling",
    );
    
    my $job = Job->new (
        name => "altr_job",
        ifiles => [$trimReads],
        ofiles => [$asmPm],
        mfiles => [$volDir],
        jobs => [$jobMkVol, $jobAlVol, $jobCatVol],
        msg => "aligning trimmed reads for assembling",
    );
    
    serialRunJobs($env, $cfg, $job);
}

sub runAssemble($$) {
    my ($env, $cfg) = @_;

    my $prjDir = %$env{"WorkPath"} . "/" .%$cfg{"PROJECT"};
    my $workDir = "$prjDir/4-fsa";
    mkdir $workDir;

    my $script = "$prjDir/scripts/assemble.sh";
    my $overlaps = "$prjDir/3-assembly/asm_pm.m4";
    my $reads = "$prjDir/2-trim_bases/trimReads.fasta";
    my $contigs = "$workDir/contigs.fasta";
    my $filtered_overlaps = "$workDir/filter.m4";

    my $binPath = %$env{"BinPath"}; 
    my $thread = %$cfg{"THREADS"};
    my $filterOptions = %$cfg{"FSA_OL_FILTER_OPTIONS"};
    if (%$cfg{"GENOME_SIZE"}) {
        $filterOptions = $filterOptions . " --genome_size=" . %$cfg{"GENOME_SIZE"};
    }
    my $assembleOptions = %$cfg{"FSA_ASSEMBLE_OPTIONS"};

    my $job = Job->new(
        name => "ass_job",
        ifiles => [$overlaps, $reads],
        ofiles => [$filtered_overlaps, $contigs],
        gfiles => [$filtered_overlaps, $contigs],
        mfiles => [],
        cmds => ["$binPath/fsa_ol_filter $overlaps $filtered_overlaps --thread_size=$thread --output_directory=$workDir $filterOptions", 
                 "$binPath/fsa_assemble $filtered_overlaps --read_file=$reads --thread_size=$thread --output_directory=$workDir $assembleOptions"],
        msg => "assembling",
    );

    serialRunJobs($env, $cfg, $job);
}


sub statCorrectedReads($$) {
    my ($env, $cfg) = @_;
    my $prjDir = %$env{"WorkPath"} . "/" .%$cfg{"PROJECT"};
 
    plgdInfo("N50 of corrected reads: $prjDir/1-consensus/cns_final.fasta");
    my $cmd = %$env{"BinPath"} . "/fsa_rd_stat $prjDir/1-consensus/cns_final.fasta";
    system($cmd);

}

sub statContigs($$) {
    my ($env, $cfg) = @_;
    my $prjDir = %$env{"WorkPath"} . "/" .%$cfg{"PROJECT"};
 
    my $cmd = %$env{"BinPath"} . "/fsa_rd_stat $prjDir/4-fsa/contigs.fasta";
    plgdInfo("N50 of contigs: $prjDir/4-fsa/contigs.fasta");
    system($cmd);
}

my %cfg = ();
my %env = ();

sub cmdCorrect($) {
    my ($fname) = @_;

    %cfg = loadMecatConfig($fname);
    %env = loadMecatEnv(\%cfg);
    initializeMecatProject(\%cfg);

    runCorrectRawreads(\%env, \%cfg);
    statCorrectedReads(\%env, \%cfg);
}

sub cmdAssemble($) {
    
    my ($fname) = @_;

    %cfg = loadMecatConfig($fname);
    %env = loadMecatEnv(\%cfg);
    initializeMecatProject(\%cfg);

    runCorrectRawreads(\%env, \%cfg);
    runTrimReads(\%env, \%cfg);
    runAlignTReads(\%env, \%cfg);
    runAssemble(\%env, \%cfg); 
    statContigs(\%env, \%cfg);
}

sub cmdConfig($) {
    my ($fname) = @_;

    my %cfg = defaultConfig();

    my @items = ("PROJECT", "RAWREADS", "GENOME_SIZE", "THREADS", "MIN_READ_LENGTH", 
                 "CNS_OVLP_OPTIONS", "CNS_OPTIONS", "CNS_OUTPUT_COVERAGE", "TRIM_OVLP_OPTIONS", "ASM_OVLP_OPTIONS", 
                 "FSA_OL_FILTER_OPTIONS", "FSA_ASSEMBLE_OPTIONS" );

    open(F, "> $fname") or die; 
    foreach my $k (@items) {
        if ($k =~ /OPTIONS/) {
            print F "$k=\"$cfg{$k}\"\n"
        } else {
            print F "$k=$cfg{$k}\n";
        }
    }
    
    foreach my $k (keys %cfg) {
        if (not grep /^$k$/, @items) {
            print F "$k=$cfg{$k}\n";
        }
    }
    close(F);

}



sub usage() {
    print "Usage: mecat.pl correct|assemble|config cfg_fname\n".
          "    correct:     correct rawreads\n" .
          "    assemble:    generate contigs\n" .
          "    config:      generate default config file\n"
}

sub main() {
    if (scalar @ARGV >= 2) {
        my $cmd = @ARGV[0];
        my $cfgfname = @ARGV[1];

        if ($cmd eq "correct") {
            cmdCorrect($cfgfname);
        } elsif ($cmd eq "assemble") {
            cmdAssemble($cfgfname);
        } elsif ($cmd eq "config") {
            cmdConfig($cfgfname);
        } else {
            usage();
        }
    } else {
        usage();
    }
}


$SIG{TERM}=$SIG{INT}=\& catchException;
sub catchException { 
    plgdInfo("Catch an Exception, and do cleanup");
    stopRunningScripts(\%env, \%cfg);
    exit -1; 
} 

#eval {
    main();
#};

if ($@) {
    catchException();
}

END {
    stopRunningScripts(\%env, \%cfg);
}
