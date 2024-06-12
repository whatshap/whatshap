#!/usr/bin/env perl

use FindBin;
use lib $FindBin::RealBin;
use Cwd;
use File::Path qw(make_path remove_tree);
use File::Basename;
use Carp;
use List::Util qw/max min/;

use Plgd::Utils;
use Plgd::Script;
use Plgd::Project;
use Necat;

use Env qw(PATH);

use strict;

my $defaultPairwiseMapingOptions = "-k 15 -z 10 -q 500 -b 2000 -s 3 -n 500 -a 500 -d 0.25 -e 0.5 -m 500 -t 1 -j 1 -u 0 -i 1";



my @defaultConfig = (
    ["PROJECT", "", ""],
    ["ONT_READ_LIST", "", ""],
    ["GENOME_SIZE", "", ""],
    ["THREADS", "4", ""],
    ["MIN_READ_LENGTH", "3000", ""],
    ["PREP_OUTPUT_COVERAGE", "40", ""],
    ["OVLP_FAST_OPTIONS", "-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000", ""],
    ["OVLP_SENSITIVE_OPTIONS", "-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000", ""],
    ["CNS_FAST_OPTIONS", "-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0", ""],
    ["CNS_SENSITIVE_OPTIONS", "-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0", ""],
    ["TRIM_OVLP_OPTIONS", "-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400", ""],
    ["ASM_OVLP_OPTIONS", "-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400", ""],
    ["NUM_ITER", "2", ""],
    ["CNS_OUTPUT_COVERAGE", "30", ""],
    ["CLEANUP", "1", ""],
    ["USE_GRID", "false", ""],
    ["GRID_NODE", "0", ""],
    ["GRID_OPTIONS", "", ""],
    ["SMALL_MEMORY", "0", ""],
    ["FSA_OL_FILTER_OPTIONS", "", ""],
    ["FSA_ASSEMBLE_OPTIONS", "", ""],
    ["FSA_CTG_BRIDGE_OPTIONS", "", ""],
    ["POLISH_CONTIGS", "true", ""],
#    ["COMPRESS", "1", ""]
);

sub defaultConfig() {
    my %cfg = ();
    for my $i (0 .. $#defaultConfig){
        $cfg{$defaultConfig[$i][0]} = $defaultConfig[$i][1];
    }
    return %cfg;
}


sub loadNecatConfig($) {
    my ($fname) = @_;
    my %cfg = ();
    loadConfig($fname, \%cfg);

    my @required = ("PROJECT", "GENOME_SIZE");
    foreach my $r (@required) {
        if (not exists($cfg{$r}) or $cfg{$r} eq "")  {
            plgdError("Not set config $r");
        }
    }

    if (not exists($cfg{"SMALL_MEMORY"}) or ($cfg{"SMALL_MEMORY"} eq "")) {
        $cfg{"SMALL_MEMORY"} = 0;
    }

    if ((not exists($cfg{"ONT_READ_LIST"}) or $cfg{"ONT_READ_LIST"} eq "") and
        (not exists($cfg{"CNS_READ_LIST"}) or $cfg{"CNS_READ_LIST"} eq "") )  {
        plgdError("Not set config ONT_READ_LIST or CNS_READ_LIST");
    }

    $cfg{"COMPRESS"} = 1;   # TODO

    return %cfg;
}

sub loadNecatEnv($) {
    my ($cfg) = @_;

    my %env = loadEnv($cfg);
    $env{"BinPath"} = $FindBin::RealBin;
    $env{"NumberOfScriptError"} = 0;
    return %env;    
}


sub initializeNecatProject($) {
    my ($cfg) = @_;

    initializeProject($cfg);
}


sub runCnsPrepare($$$) {
    my ($env, $cfg, $readList) = @_;

    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};
    my $workDir = "$prjDir/1-consensus/raw_reads";
    mkdir $workDir;

    #my $rawreadList = %$cfg{"ONT_READ_LIST"};
    my $filterdList = "$workDir/raw_read_list.txt";
    my $binPath = %$env{"BinPath"};
    my $threads = %$cfg{"THREADS"};
    my $minLength = %$cfg{"MIN_READ_LENGTH"};
    
    my $isGz = $cfg->{"COMPRESS"};
    my $filtered = $isGz ? "$workDir/filtered_raw_reads.fasta.gz" : "$workDir/filtered_raw_reads.fasta";
    my $filteredTmp =  "$workDir/filtered_raw_reads.fasta";
    my $baseSize = $cfg->{"GENOME_SIZE"}* $cfg->{"PREP_OUTPUT_COVERAGE"};

    my @cmds = ();
    push @cmds, "$binPath/fsa_rd_tools longest --ifname $readList --ofname $filteredTmp --min_length $minLength --base_size $baseSize --discard_illegal_read";
    if ($isGz) {
        push @cmds, "$binPath/pigz -f -p $threads $filteredTmp";
    }
    push @cmds, "echo $filtered > $filterdList";

    my $job = Job->new(
        prefunc => sub ($) { 
            my ($job) = @_;
            my @files = linesInFile($readList);
            for my $i (0..$#files) {
                if (not -e $files[$i] ) {  die "File not Exist: $files[$i]\n"; }
                push @{$job->ifiles}, $files[$i];
            }
        },
        name => "cns_pprr",
        ifiles => [$readList],   
        ofiles => [$filtered, $filterdList],   
        gfiles => [$filtered, $filterdList],   
        mfiles => [],
        cmds => [@cmds],
        msg => "filtering reads for consensus"
    );

    serialRunJobs($env, $cfg, $job);

}

sub runCnsAlign($$$) {
    my ($env, $cfg, $step) = @_;
    
    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};
    my $workDir = "$prjDir/1-consensus/cns_iter$step";
    my $volDir = "$workDir/PackedData";
    mkdir $volDir;
    
    my $candidates = "$workDir/cns_candidates.txt";
    my $binPath = %$env{"BinPath"};
    my $threads = %$cfg{"THREADS"};
    my $volFile = "$volDir/reads_info.txt";

    my $readList = "";
    my $alOptions = "";
    if ($step == 1) {
        $readList = "$prjDir/1-consensus/raw_reads/raw_read_list.txt";
        $alOptions = mergeOptionString($defaultPairwiseMapingOptions, %$cfg{"OVLP_SENSITIVE_OPTIONS"});

    } else {
        my $ii = $step -1;
        $readList = "$prjDir/1-consensus/cns_iter$ii/read_list.txt";
        $alOptions = mergeOptionString($defaultPairwiseMapingOptions, %$cfg{"OVLP_FAST_OPTIONS"});
    }

    my $jobMkVol = Job->new (
        name => "cns_${step}_mk_vol",
        ifiles => [$readList],
        ofiles => [$volFile],
        gfiles => [$volFile],
        mfiles => [],
        cmds => ["$binPath/oc2mkdb $volDir $readList"],
        msg => "making vol for cns $step",
    );

    my $jobAlVol = Job->new (
        prefunc => sub($) {
            my ($job) = @_;
            my $count = getFileFirstItem($volFile, 0);

            for (my $i=0; $i<$count; $i=$i+1) {
                my $jobSub = Job->new(
                    name => "cns_${step}_al_vol_$i", 
                    ifiles => [$volFile],
                    ofiles => ["$volDir/pm_result_$i"],
                    gfiles => ["$volDir/pm_result_$i"],
                    mfiles => [],
                    cmds => ["$binPath/oc2pmov $alOptions -t $threads $volDir $i $volDir/pm_result_$i"],
                    msg => "aligning vol $i for cns $step",
                );
                push @{$job->pjobs}, $jobSub;
                push @{$job->ofiles}, "$volDir/pm_result_$i";
            }
        },

        name => "cns_${step}_al_vol",
        ifiles => [$volFile],
        ofiles => [],   # prefunc
        mfiles => [],
        pjobs => [],    # prefunc
        msg => "aligning vol for cns $step",
    );

    my $jobCatVol = Job->new (
        prefunc => sub($) {
            my ($job) = @_;
            my $pms = join(" ", @{$jobAlVol->ofiles});
            $job->ifiles($jobAlVol->ofiles);
            $job->cmds(["cat $pms > $candidates"]);
        },
        name => "cns_${step}_cat_vol",
        ifiles => [], # prefunc
        ofiles => [$candidates],
        gfiles => [$candidates],
        mfiles => [],
        cmds => [], # prefunc
        msg => "catenating vol for cns $step",
    );
        
    serialRunJobs($env, $cfg, Job->new(
        name => "cns_${step}_align",
        ifiles => [$readList],
        ofiles => [$candidates],
        mfiles => [],
        jobs => [$jobMkVol, $jobAlVol, $jobCatVol],
        msg => "aliging for cns $step",
    ));
}


sub runCnsCorrect($$$) {
    my ($env, $cfg, $step) = @_;

    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};

    my $workDir = "$prjDir/1-consensus/cns_iter$step";
    mkdir $workDir;
    my $volDir = "$workDir/PackedData";
    mkdir $volDir;
    
    my $cnsOptions="";
    if ($step == 1) {
        $cnsOptions= %$cfg{"CNS_SENSITIVE_OPTIONS"} . " -r 0";
    } else {
        $cnsOptions= %$cfg{"CNS_FAST_OPTIONS"} . " -r 1";
    }

    if ($step == %$cfg{"NUM_ITER"}) {
        $cnsOptions= "$cnsOptions -f 0";
    } else {
        $cnsOptions= "$cnsOptions -f 1";
    }

    my $outReadList = "$workDir/read_list.txt";
    my $candidates = "$workDir/cns_candidates.txt";
    my $cnsReads = "$workDir/cns.fasta";
    my $uncnsReads = "$workDir/raw.fasta";
    my $cnsReadsTemp = "$workDir/tmp_cns.fasta";
    my $uncnsReadsTemp = "$workDir/tmp_raw.fasta";


    my $binPath = %$env{"BinPath"};
    my $threads = %$cfg{"THREADS"};
    my $smallMemory = %$cfg{"SMALL_MEMORY"};

    my $partFile = "${candidates}.partitions";

    my $jobPart = Job->new(
        name => "cns_${step}_part",
        ifiles => [$candidates],
        ofiles => [$partFile],
        gfiles => [$partFile],
        mfiles => [],
        cmds => ["$binPath/oc2pcan -t  $threads $volDir $candidates"],
        msg => "partitioning candidates",
    );

    my $jobCns = Job->new(
        prefunc => sub ($) {
            my ($job) = @_;

            my @subJobs = ();
            my $count = getFileFirstItem($partFile, 0);

            for (my $i=0; $i<$count; $i=$i+1) {
                my $jobSub = Job->new(
                    name => "cns_${step}_cns_part_$i", 
                    ifiles => [$partFile],
                    ofiles => ["${cnsReadsTemp}_$i", "${uncnsReadsTemp}_$i"],
                    gfiles => ["${cnsReadsTemp}_$i", "${uncnsReadsTemp}_$i"],
                    mfiles => [],
                    cmds => ["$binPath/oc2cns -s $smallMemory -t $threads $cnsOptions $volDir $candidates ${cnsReadsTemp}_$i ${uncnsReadsTemp}_$i -mn $i $count"],
                    msg => "aligning volumn $i for correcting $step",
                );

                push @{$job->pjobs}, $jobSub;
                push @{$job->ofiles}, "${cnsReadsTemp}_$i", "${uncnsReadsTemp}_$i";
            }

        },
        name => "cns_${step}_cns",
        ifiles => [$partFile],
        ofiles => [], # prefunc
        mfiles => [],
        pjobs => [], # prefunc
        msg => "correcting part for cns $step",

    );

    my $jobCat = Job->new(
        prefunc => sub($) {
            my ($job) = @_;

            my $cns = "";
            my $uncns = "";
            my $count = getFileFirstItem($partFile, 0);
            for (my $i=0; $i < $count; $i = $i + 1) {
                $cns = "$cns ${cnsReadsTemp}_$i";
                $uncns = "$uncns ${uncnsReadsTemp}_$i";
            }
            $job->ifiles($jobCns->ofiles);
            $job->cmds(["cat $cns > $cnsReadsTemp",
                        "cat $uncns > $uncnsReadsTemp"]);

            @{$job->mfiles} = @{$jobCns->ofiles};
        },
        name => "cns_${step}_cat_cns",
        ifiles => [],   # prefunc
        ofiles => [$cnsReadsTemp, $uncnsReadsTemp],
        gfiles => [$cnsReadsTemp, $uncnsReadsTemp],
        mfiles => [],
        cmds => [],     # prefunc
        msg => "catenating reads for cns $step.",
    );

    my $jobReorder = Job->new(
        name => "cns_${step}_reorder",
        ifiles => [$cnsReadsTemp, $uncnsReadsTemp],
        ofiles => ["$cnsReads.gz", "$uncnsReads.gz", $outReadList],
        gfiles => ["$cnsReads.gz", "$uncnsReads.gz", $outReadList],
        mfiles => [],
        cmds => ["$binPath/oc2ReorderCnsReads $cnsReadsTemp $uncnsReadsTemp $cnsReads $uncnsReads",
                 "$binPath/pigz -f  -p $threads $cnsReads",
                 "$binPath/pigz -f  -p $threads $uncnsReads",
                 "echo $cnsReads.gz > $outReadList",
                 "echo $uncnsReads.gz >> $outReadList" ],
        msg => "reordering for cns $step",
    );

    serialRunJobs($env, $cfg, Job->new(
        name => "cns_${step}_correct",
        ifiles => [$candidates, "$volDir/reads_info.txt"],
        ofiles => ["$cnsReads.gz", "$uncnsReads.gz", $outReadList],
        mfiles => [$candidates, "$candidates.p*", $cnsReadsTemp, $uncnsReadsTemp, "${cnsReadsTemp}_", "${uncnsReadsTemp}_*"],
        jobs => [$jobPart, $jobCns, $jobCat, $jobReorder],
        msg => "correcting for cns $step",
    ));
}

sub runCnsIterate($$$) {
    my ($env, $cfg, $step) = @_;

    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};
    my $workDir = "$prjDir/1-consensus/cns_iter$step";
    mkdir $workDir;
    my $volDir = "$workDir/PackedData";
    mkdir $volDir;
    
    my $readList = "$prjDir/1-consensus/raw_reads/raw_read_list.txt";
    if ($step > 1) {
        $readList = "$prjDir/1-consensus/cns_iter" . ($step -1) . "/read_list.txt";
    }
    my $outReadList = "$workDir/read_list.txt";
    my $cnsReads = "$workDir/cns.fasta.gz";
    my $uncnsReads = "$workDir/raw.fasta.gz";

    serialRunJobs($env, $cfg, Job->new(
        name => "cns_${step}",
        ifiles => [$readList],
        ofiles => [$cnsReads, $uncnsReads, $outReadList],
        mfiles => [$volDir],
        funcs => [sub($$) { runCnsAlign($env, $cfg, $step)} ,
                  sub($$) { runCnsCorrect($env, $cfg, $step);}],
        msg => "correcting iteratively for cns $step."
    ));
}

sub runCnsExtract($$$) {
    my ($env, $cfg, $lastCnsReads) = @_;

    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};
    my $workDir = "$prjDir/1-consensus";

    my $binPath = %$env{"BinPath"};
    my $threads = %$cfg{"THREADS"};
    my $count = %$cfg{"NUM_ITER"};
    #my $lastCnsReads = "$workDir/cns_iter$count/cns.fasta.gz"; 
    
    my $isGz = $cfg->{"COMPRESS"};
    my $finalCnsReads = $isGz ? "$workDir/cns_final.fasta.gz" : "$workDir/cns_final.fasta";
        
    my $base_size = %$cfg{"GENOME_SIZE"}*%$cfg{"CNS_OUTPUT_COVERAGE"};

    my $job = jobExtractReads($env, $cfg, "cns", $lastCnsReads, $finalCnsReads, $base_size);

    serialRunJobs($env, $cfg, $job);

}


sub runConsensusRaw($$) {
    my ($env, $cfg) = @_;

    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};
    my $workDir = "$prjDir/1-consensus";
    mkdir $workDir;
    
    my $isGz = $cfg->{"COMPRESS"};
    my $finalCnsReads = $isGz ? "$workDir/cns_final.fasta.gz" : "$workDir/cns_final.fasta";

    my $rawreadList = %$cfg{"ONT_READ_LIST"};
    my @funcs = (sub ($$) { runCnsIterate($env, $cfg, 1)},
                 sub ($$) { runCnsIterate($env, $cfg, 2)},
                 sub ($$) { runCnsIterate($env, $cfg, 3)},
                 sub ($$) { runCnsIterate($env, $cfg, 4)},
                 sub ($$) { runCnsIterate($env, $cfg, 5)},
                 sub ($$) { runCnsIterate($env, $cfg, 6)},
                 sub ($$) { runCnsIterate($env, $cfg, 7)},
                 sub ($$) { runCnsIterate($env, $cfg, 8)});

    
    my $count = %$cfg{"NUM_ITER"};
    plgdError("Parameter NUM_ITER must be in [1, " . scalar @funcs . "].") if ($count < 1 or $count > scalar @funcs) ;

    sub funcRawPrepare($$) {
        my ($env, $cfg) = @_;
        runCnsPrepare($env, $cfg, $rawreadList);
    };

    sub funcRawExtract($$) {
        my ($env, $cfg) = @_;
        my $count = %$cfg{"NUM_ITER"};
        my $lastCnsReads = "$workDir/cns_iter$count/cns.fasta.gz"; 
        runCnsExtract($env, $cfg, $lastCnsReads);

    }

    my $job = Job->new(
        prefunc => sub ($) {
            my ($job) = @_;
            my @files = linesInFile($rawreadList);
            for my $i (0..$#files) {
                if (not -e $files[$i] ) {  die "File not Exist: $files[$i]\n"; }
                push @{$job->ifiles}, $files[$i];

            }
        },

        name => "cns_job",
        ifiles => [$rawreadList],
        ofiles => [$finalCnsReads],
        mfiles => [],
        funcs => [\&funcRawPrepare, @funcs[0..($count-1)], \&funcRawExtract], #prefunc
        msg => "correcting rawreads",
    );

    serialRunJobs($env, $cfg, $job);
}


sub runConsensusCns($$) {
    my ($env, $cfg) = @_;

    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};
    my $workDir = "$prjDir/1-consensus";
    mkdir $workDir;

    my $readList = %$cfg{"CNS_READ_LIST"};
    
    sub funcCnsPrepare($$) {
        my ($env, $cfg) = @_;
        runCnsPrepare($env, $cfg, $readList);
    };

    sub funcCnsExtract($$) {
        my ($env, $cfg) = @_;

        runCnsExtract($env, $cfg, "$workDir/raw_reads/raw_read_list.txt");

    }
    
    my $job = Job->new(
        prefunc => sub ($) {
            my ($job) = @_;
            my @files = linesInFile($readList);
            for my $i (0..$#files) {
                if (not -e $files[$i] ) {  die "File not Exist: $files[$i]\n"; }
                push @{$job->ifiles}, $files[$i];

            }
        },

        name => "cns_job",
        ifiles => [$readList],  #prefunc
        ofiles => ["$workDir/cns_final.fasta.gz"],
        mfiles => [],
        funcs => [\&funcCnsPrepare, \&funcCnsExtract], 
        msg => "correcting cnsreads",
    );

    serialRunJobs($env, $cfg, $job);
}

sub runConsensus($$) {
    my ($env, $cfg) = @_;

    if (not exists($cfg->{"CNS_READ_LIST"}) or %$cfg{"CNS_READ_LIST"} eq "") {
        runConsensusRaw($env, $cfg);
    } else {
        runConsensusCns($env, $cfg);
    }
}

sub makePmJob($$$$$$$$) {
    my ($env, $cfg, $name, $desc, $readList, $pmm4, $options, $volDir) = @_;

    my $binPath = %$env{"BinPath"};
    my $threads = %$cfg{"THREADS"};
    
    my $volFile = "$volDir/reads_info.txt";

    my $jobMkVol = Job->new(
        prefunc => sub($) {
            my ($job) = @_;
            my @files = linesInFile($readList);
            for my $i (0..$#files) {
                if (not -e $files[$i] ) {  die "File not Exist: $files[$i]\n"; }
                push @{$job->ifiles}, $files[$i];
            }
        }, 

        name => "${name}_mk_vol",
        ifiles => [$readList],  # prefunc
        ofiles => [$volFile],
        gfiles => [$volFile],
        mfiles => [],
        cmds => ["$binPath/oc2mkdb $volDir $readList"],
        msg => "making vol for $desc",
    );

    my $jobAlVol = Job->new(
        prefunc => sub($) {
            my ($job) = @_;
            
            my $count = getFileFirstItem($volFile, 0);

            for (my $i=0; $i<$count; $i=$i+1) {
                my $jobSub = Job->new(
                    name => "${name}_al_vol_$i", 
                    ifiles => [$volFile],
                    ofiles => ["$volDir/pm_result_$i"],
                    gfiles => ["$volDir/pm_result_$i"],
                    mfiles => [],
                    cmds => ["$binPath/oc2asmpm $options -t $threads $volDir $i $volDir/pm_result_$i"],
                    msg => "aligning vol $i for $desc",
                );
                push @{$job->pjobs}, $jobSub;
                push @{$job->ofiles}, "$volDir/pm_result_$i";
            }

        },
        name => "${name}_al_vol",                    
        ifiles => [$volFile],
        ofiles => [], # prefunc
        mfiles => [],
        pjobs => [], # prefunc
        msg => "aligning vol for $desc",
    );

    my $jobCatVol = Job->new (
        prefunc => sub ($) {
            my ($job) = @_;
            my $subPmStr = join(" ", @{$jobAlVol->ofiles});
            $job->ifiles($jobAlVol->ofiles);
            if (scalar @{$job->ifiles} > 0) {
                $job->cmds(["cat $subPmStr  > $pmm4"]);
            } else {
                $job->cmds([":> $pmm4"]);

            }
        },

        name => "${name}_cat_vol",
        ifiles => [], # prefunc
        ofiles => [$pmm4],
        gfiles => [$pmm4],
        mfiles => [],
        cmds => [], # prefunc
        msg => "catenating vol for ${name}",
    );

    return Job->new (
        prefunc => sub($) {
            my ($job) = @_;
            my @files = linesInFile($readList);
            for my $i (0..$#files) {
                if (not -e $files[$i] ) {  die "File not Exist: $files[$i]\n"; }
                push @{$job->ifiles}, $files[$i];
            }
        }, 

        name => "${name}_align",
        ifiles => [$readList],  # prefunc
        ofiles => [$pmm4],
        mfiles => [],
        jobs => [$jobMkVol, $jobAlVol, $jobCatVol],
        msg => "aligning for $desc",
    );
}

sub makeRmJob($$$$$$$$$) {
    my ($env, $cfg, $name, $desc, $readList, $ref, $pmm4, $options, $volDir) = @_;


    my $binPath = %$env{"BinPath"};
    my $threads = %$cfg{"THREADS"};
    my $volFile = "$volDir/reads_info.txt";

    my $jobMkVol = Job->new (

        name => "${name}_mk_vol",
        ifiles => [$readList],
        ofiles => [$volFile],
        gfiles => [$volFile],
        mfiles => [],
        cmds => ["$binPath/oc2mkdb $volDir $readList"],
        msg => "making vol for $desc",
    );

    my $jobAlVol = Job->new (
        prefunc => sub($) {
            my ($job) = @_;
            my $count = getFileFirstItem($volFile, 0);
            
            for (my $i=0; $i<$count; $i=$i+1) {
                my $jobSub = Job->new(
                    name => "${name}_al_vol_$i", 
                    ifiles => [$volFile],
                    ofiles => ["$volDir/rm_result_$i"],
                    gfiles => ["$volDir/rm_result_$i"],
                    mfiles => [],
                    cmds => ["$binPath/oc2rm_worker -t $threads $volDir $ref $volDir/rm_result_$i -mn $i $count"],
                    msg => "aligning volumn $i for $desc",
                );
                push @{$job->pjobs}, $jobSub;
                push @{$job->ofiles}, "$volDir/rm_result_$i";
            }
        },

        name => "${name}_al_vol",
        ifiles => [$volFile],
        ofiles => [], # prefunc
        mfiles => [],
        pjobs => [], # prefunc
        msg => "aligning vol for $desc",
    );

    my $jobCatVol = Job->new (
        prefunc => sub ($) {
            my ($job) = @_;
            my $subPmStr = join(" ", @{$jobAlVol->ofiles});
            $job->ifiles($jobAlVol->ofiles);
            if (scalar @{$jobAlVol->ofiles} > 0) {
                $job->cmds(["cat $subPmStr  > $pmm4"]);
            } else {
                #$job->cmds(["rm -f $pmm4 && touch $pmm4"]);
                $job->cmds([":> $pmm4"]);
            }
        },

        name => "${name}_cat_vol",
        ifiles => [], # prefunc
        ofiles => [$pmm4],
        gfiles => [$pmm4],
        mfiles => [],   #prefunc
        cmds => [], # prefunc
        msg => "catenating vol for ${name}",
    );

}


sub runTrimBasesFast($$) {
    my ($env, $cfg) = @_;
    
    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};
    my $workDir = "$prjDir/2-trim_bases";
    mkdir $workDir;
    mkdir "$prjDir/3-assembly";   # 3-assembly 

    my $binPath = %$env{"BinPath"};
    my $threads = %$cfg{"THREADS"};

    my $cnsReads = "$prjDir/1-consensus/cns_final.fasta.gz";
    my $trimReads = "$prjDir/trimReads.fasta";
    my $renumReads = "$workDir/renum_reads.fasta";
    my $options = %$cfg{"ASM_OVLP_OPTIONS"};

    my $volDir = "$workDir/pac_in";
    my $u2uVolDir = "$workDir/uu_ovlps";
    my $c2uVolDir = "$workDir/cu_ovlps";

    my $trimReads = "$prjDir/trimReads.fasta";
    my $pmm4 = "$prjDir/3-assembly/pm.m4";
    my $trimPmm4 = "$workDir/pm.m4";
    my $tmpTrimReads = "$workDir/tmp_trimReads.fasta";
    my $tmpPmm4 = "$workDir/tmp_pm.m4";

    my $compReads = "$workDir/complete_reads.fasta";
    my $compReadList = "$workDir/complete_read_list.txt";
    my $uncompReads = "$workDir/uncomplete_reads.fasta";
    my $uncompReadList = "$workDir/uncomplete_read_list.txt";
    my $c2uPmm4 = "$workDir/cu_ovlps.rm.m4";
    my $u2uPmm4 = "$workDir/uu_ovlps.rm.m4";

    my $jobRenum = Job->new(
        name => "tr_renum",
        ifiles => [$cnsReads],
        ofiles => [$renumReads],
        gfiles => [$renumReads, "$workDir/read_list.txt"],
        mfiles => [],
        cmds => ["$binPath/oc2renumberSeqs $cnsReads $renumReads",
                 "echo $workDir/renum_reads.fasta > $workDir/read_list.txt"],
        msg => "renum reads for trimming reads",
    );

    my $jobAlign = makePmJob($env, $cfg, "tr", "trimmng reads", "$workDir/read_list.txt", $trimPmm4, $options . " -u 1", $volDir);

    my $errCut = 0.1;
    my $minOvlp = 1;
    my $minCov = 1;
    my $minSize = 1000;
    my $clippedRanges = "$workDir/clipped_ranges.txt";

    my $jobTrim1 = Job->new(
        name => "tr_trim",
        ifiles => [$trimPmm4, $renumReads], # prefunc
        ofiles => [$compReads, $uncompReads, $tmpPmm4, $uncompReadList],
        gfiles => [$compReads, $uncompReads, $tmpPmm4, $uncompReadList, $compReadList],
        mfiles => [$clippedRanges],
        cmds => ["$binPath/oc2pm4 $volDir $trimPmm4 $errCut $threads", 
                 "$binPath/oc2lcr $trimPmm4 $volDir $errCut $minOvlp $minCov $minSize $threads $clippedRanges",
                 "$binPath/oc2etr $clippedRanges $renumReads $trimPmm4 $compReads $uncompReads $tmpPmm4",
                 "echo $uncompReads > $uncompReadList",
                 "echo $compReads > $compReadList"], 

        msg => "trimming reads after aligning",
    );


    my $jobC2u = Job->new (

        prefunc => sub($) {
            my ($job) = @_;

            if ( (-z  $compReads) or (-z $uncompReads)) {
                push @{$job->cmds}, "touch $c2uPmm4";
            } else {
                push @{$job->jobs}, makeRmJob($env, $cfg, "tr_c2u", "algining complete to uncomplete for trimming",
                                             $compReadList, $uncompReads, $c2uPmm4, $options, $c2uVolDir);
            }

        },
        name => "tr_c2u",
        ifiles => [$compReads, $uncompReads], # prefunc
        ofiles => [$c2uPmm4],
        gfiles => [$c2uPmm4],
        mfiles => [],
        #cmds => ["if [[ ! -s $compReads ]] && [[ ! -s $uncompReads ]]; then \n ".
        #           "$binPath/oc2rm -k 13 -t $threads $compReads $uncompReads $c2uPmm4\n" .
        #         "else\n" .
        #           "touch $c2uPmm4\n".
        #         "fi"], 
        msg => "algining complete to uncomplete for trimming",
    );

    my $jobU2u = makePmJob($env, $cfg, "tr_u2u", "algining uncomplete", $uncompReadList, $u2uPmm4, $options, $u2uVolDir);

    my $jobOrder = Job->new (
        name => "tr_order",
        ifiles => [$compReads, $uncompReads, $c2uPmm4, $u2uPmm4], # prefunc
        ofiles => ["$trimReads.gz", "$pmm4.gz"],
        gfiles => ["$trimReads.gz", "$pmm4.gz"],
        mfiles => [$compReads, $uncompReads, $c2uPmm4, $u2uPmm4, $tmpTrimReads, $tmpPmm4],
        cmds => ["cat $compReads $uncompReads > $tmpTrimReads",
                  "cat $c2uPmm4 $u2uPmm4 >> $tmpPmm4",
                 "$binPath/oc2orderResults $tmpTrimReads $tmpPmm4 $trimReads $pmm4",
                 "$binPath/pigz  -f -p $threads $trimReads",
                 "$binPath/pigz  -f -p $threads $pmm4"], 

        msg => "ordering result for trimming",
    );

    serialRunJobs($env, $cfg, Job->new(
        name => "tr_job",
        ifiles => [$cnsReads],
        ofiles => ["$trimReads.gz", "$pmm4.gz"],
        mfiles => [$volDir, $pmm4, $trimPmm4, $renumReads, $u2uVolDir,"$workDir/pm.m4.p*", "$workDir/pm.m4.partitions", "$workDir/read_list.txt",
                   $uncompReadList, $compReadList],
        jobs => [$jobRenum, $jobAlign, $jobTrim1, $jobC2u, $jobU2u, $jobOrder],
        msg => "trimming reads",
    ));
}


sub runTrimAccurate0($$) {
    my ($env, $cfg) = @_;

    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};
    my $workDir = "$prjDir/2-trim_bases";
    mkdir $workDir;

    #my $cnsReads = "$prjDir/1-consensus/cns_iter" . %$cfg{"NUM_ITER"} . "/cns.fasta";
    my $cnsReads = "$prjDir/1-consensus/cns_final.fasta.gz";
    my $trimReads = "$prjDir/trimReads.fasta";
    my $renumReads = "$workDir/renum_reads.fasta";
    my $pmm4 = "$workDir/pm.m4";

    my $volDir = "$workDir/pac_in";

    my $binPath = %$env{"BinPath"};
    my $threads = %$cfg{"THREADS"};
    my $options = %$cfg{"TRIM_OVLP_OPTIONS"};
    my $volFile = "$volDir/reads_info.txt";
    my $trimBinPath = "trim_bases_accurate0";

    my $jobRenum = Job->new(
        name => "tr_renum",
        ifiles => [$cnsReads],
        ofiles => ["$workDir/renum_reads.fasta"],
        gfiles => ["$workDir/renum_reads.fasta"],
        mfiles => [],
        cmds => ["$binPath/oc2renumberSeqs $cnsReads $workDir/renum_reads.fasta"],
        msg => "renum reads for trimming reads",
    );

    my $jobMkVol = Job->new(
        name => "tr_mk_vol",
        ifiles => ["$workDir/renum_reads.fasta"],
        ofiles => [$volFile],
        gfiles => [$volFile],
        mfiles => [],
        cmds => ["echo $workDir/renum_reads.fasta > $workDir/read_list.txt",
                 "$binPath/oc2mkdb $volDir $workDir/read_list.txt"],
        msg => "making vol for trimming reads",
    );

    my $jobAlVol = Job->new(
        prefunc => sub($) {
            my ($job) = @_;
            
            my $count = getFileFirstItem("$volDir/reads_info.txt", 0);

            for (my $i=0; $i<$count; $i=$i+1) {
                my $jobSub = Job->new(
                    name => "tr_al_vol_$i", 
                    ifiles => ["$volDir/volume_names.txt"],
                    ofiles => ["$volDir/pm_result_$i"],
                    gfiles => ["$volDir/pm_result_$i"],
                    mfiles => [],
                    cmds => ["$binPath/oc2asmpm $options -t $threads $volDir $i $volDir/pm_result_$i"],
                    msg => "aligning vol $i for trimming reads",
                );
                push @{$job->pjobs}, $jobSub;
                push @{$job->ofiles}, "$volDir/pm_result_$i";
            }

        },
        name => "tr_al_vol",                    
        ifiles => ["$volDir/volume_names.txt"],
        ofiles => [], # prefunc
        mfiles => [],
        pjobs => [], # prefunc
        msg => "aligning vol for trimming reads",
    );

    my $jobCatVol = Job->new (
        prefunc => sub ($) {
            my ($job) = @_;
            my $subPmStr = join(" ", @{$jobAlVol->ofiles});
            $job->ifiles($jobAlVol->ofiles);
            $job->cmds(["cat $subPmStr  > $pmm4"]);
        },

        name => "tr_cat_vol",
        ifiles => [], # prefunc
        ofiles => [$pmm4],
        gfiles => [$pmm4],
        mfiles => [],
        cmds => [], # prefunc
        msg => "catenating vol for trimming reads",
    );

    my $errCut = 0.09;
    my $minOvlp = 1;
    my $minCov = 1;
    my $minSize = 500;
    #my $clippedRanges = "$workDir/clipped_ranges.txt";
    my $jobTrim = Job->new(

        name => "tr_trim",
        ifiles => [$pmm4, $renumReads], # prefunc
        ofiles => ["$trimReads.gz"],
        gfiles => ["$trimReads.gz"],
        mfiles => [],
        # chen modified trimming method to reduce overhangs. oc2trimBases was removed
        cmds => ["$binPath/$trimBinPath/oc2pm4 $workDir/pac_in $pmm4 $errCut $threads", 
                 "$binPath/$trimBinPath/oc2lcr $pmm4 $volDir $errCut $minOvlp $minCov $minSize $trimReads $threads",
                 "$binPath/pigz -f -p $threads $trimReads"], 
                #"$binPath/oc2lcr $pmm4 $packedDir $errCut $minOvlp $minCov $minSize $clippedRanges $thread", 
                #"$binPath/oc2trimBases $renumReads $clippedRanges $trimReads");

        msg => "trimming reads after aligning",
    );

    serialRunJobs($env, $cfg, Job->new(
        name => "tr_job",
        ifiles => [$cnsReads],
        ofiles => ["$trimReads.gz"],
        mfiles => [$volDir, $pmm4, $renumReads],
        jobs => [$jobRenum, $jobMkVol, $jobAlVol, $jobCatVol, $jobTrim],
        msg => "trimming reads",
    ));
}

sub runTrimAccurate($$) {
    my ($env, $cfg) = @_;

    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};
    my $workDir = "$prjDir/2-trim_bases";
    mkdir $workDir;

    my $cnsReads = "$prjDir/1-consensus/cns_final.fasta.gz";
    my $trimReads = "$prjDir/trimReads.fasta";
    my $renumReads = "$workDir/renum_reads.fasta";
    my $pmm4 = "$workDir/pm.m4";

    my $volDir = "$workDir/pac_in";

    my $binPath = %$env{"BinPath"};
    my $threads = %$cfg{"THREADS"};
    my $smallMemory = %$cfg{"SMALL_MEMORY"};
    my $options = %$cfg{"TRIM_OVLP_OPTIONS"};
    my $volFile = "$volDir/reads_info.txt";
    my $trimBinPath = "trim_bases_accurate";

    my $jobRenum = Job->new(
        name => "tr_renum",
        ifiles => [$cnsReads],
        ofiles => ["$workDir/renum_reads.fasta"],
        gfiles => ["$workDir/renum_reads.fasta"],
        mfiles => [],
        cmds => ["$binPath/oc2renumberSeqs $cnsReads $workDir/renum_reads.fasta"],
        msg => "renum reads for trimming reads",
    );

    my $jobMkVol = Job->new(
        name => "tr_mk_vol",
        ifiles => ["$workDir/renum_reads.fasta"],
        ofiles => [$volFile],
        gfiles => [$volFile],
        mfiles => ["$workDir/read_list.txt"],
        cmds => ["echo $workDir/renum_reads.fasta > $workDir/read_list.txt",
                 "$binPath/oc2mkdb $volDir $workDir/read_list.txt"],
        msg => "making vol for trimming reads",
    );

    my $jobAlVol = Job->new(
        prefunc => sub($) {
            my ($job) = @_;
            
            my $count = getFileFirstItem("$volDir/reads_info.txt", 0);

            for (my $i=0; $i<$count; $i=$i+1) {
                my $jobSub = Job->new(
                    name => "tr_al_vol_$i", 
                    ifiles => ["$volDir/volume_names.txt"],
                    ofiles => ["$volDir/pm_result_$i"],
                    gfiles => ["$volDir/pm_result_$i"],
                    mfiles => [],
                    cmds => ["$binPath/oc2asmpm $options -t $threads $volDir $i $volDir/pm_result_$i"],
                    msg => "aligning vol $i for trimming reads",
                );
                push @{$job->pjobs}, $jobSub;
                push @{$job->ofiles}, "$volDir/pm_result_$i";
            }

        },
        name => "tr_al_vol",                    
        ifiles => ["$volDir/volume_names.txt"],
        ofiles => [], # prefunc
        mfiles => [],
        pjobs => [], # prefunc
        msg => "aligning vol for trimming reads",
    );

    my $jobCatVol = Job->new (
        prefunc => sub ($) {
            my ($job) = @_;
            my $subPmStr = join(" ", @{$jobAlVol->ofiles});
            $job->ifiles($jobAlVol->ofiles);
            $job->cmds(["cat $subPmStr  > $pmm4"]);
        },

        name => "tr_cat_vol",
        ifiles => [], # prefunc
        ofiles => [$pmm4],
        gfiles => [$pmm4],
        mfiles => [],
        cmds => [], # prefunc
        msg => "catenating vol for trimming reads",
    );

    my $errCut = 0.09;
    my $minOvlp = 1;
    my $minCov = 1;
    my $minSize = 500;
    
    my $jobTrimBefore = Job->new(
        name => "tr_trim_b",
        ifiles => [$pmm4], # prefunc
        ofiles => ["$pmm4.partitions"],   # TODO
        gfiles => [],
        mfiles => [],
        
        cmds => ["$binPath/$trimBinPath/oc2pm4 $workDir/pac_in $pmm4 $errCut $threads"],     # generate pm.m4.p* pm.m4.partitions
                 
        msg => "trimming reads 1 after aligning",
    );

    
    my $jobTrimMid = Job->new(
        prefunc => sub($) {
            my ($job) = @_;
            
            my $count = getFileFirstItem("$volDir/reads_info.txt", 0);

            for (my $i=0; $i<$count; $i=$i+1) {
                my $jobSub = Job->new(
                    name => "tr_trim_m_$i", 
                    ifiles => ["$volDir/reads_info.txt", $pmm4],
                    ofiles => ["${trimReads}_$i"],
                    gfiles => ["${trimReads}_$i"],
                    mfiles => [],
                    #cmds => [ "$binPath/oc2lcr ${pmm4} $volDir $errCut $minOvlp $minCov $minSize $trimReads  $threads -mn $i $count"],
                    cmds => [ "$binPath/$trimBinPath/oc2lcr ${pmm4} $volDir $smallMemory $threads ${trimReads}_$i -mn $i $count"],
                    msg => "trimming reads 2 $i after aligning ",
                );
                push @{$job->pjobs}, $jobSub;
                push @{$job->ofiles}, "${trimReads}_$i";
            }

        },

        name => "tr_trim_m",
        ifiles => ["$pmm4.partitions", "$volDir/reads_info.txt"], # prefunc
        ofiles => [], # prefunc
        gfiles => [], # prefunc
        mfiles => [],
        pjobs => [], # prefunc
        msg => "trimming reads 2 after aligning",
    );

    my $jobTrimAfter = Job->new(
        prefunc => sub ($) {
            my ($job) = @_;
            my $subPmStr = join(" ", @{$jobTrimMid->ofiles});
            $job->ifiles($jobTrimMid->ofiles);
            $job->mfiles($jobTrimMid->ofiles);
            $job->cmds(["cat $subPmStr | $binPath/pigz -f -p $threads  - -c > $trimReads.gz"]);
        },

        name => "tr_trim_a",
        ifiles => [], # prefunc
        ofiles => ["$trimReads.gz"],
        gfiles => ["$trimReads.gz"],
        mfiles => [],
        cmds => [], #prefunc 
        msg => "trimming reads 3 after aligning",
    );

    serialRunJobs($env, $cfg, Job->new(
        name => "tr_job",
        ifiles => [$cnsReads],
        ofiles => ["$trimReads.gz"],
        mfiles => [$volDir, $pmm4, $renumReads, "$workDir/pm.m4.p*", "$workDir/pm.m4.partitions"],
        jobs => [$jobRenum, $jobMkVol, $jobAlVol, $jobCatVol, $jobTrimBefore, $jobTrimMid, $jobTrimAfter],
        msg => "trimming reads",
    ));
}

sub runAlignReads($$) {
    my ($env, $cfg) = @_;

    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};
    my $workDir = "$prjDir/3-assembly";
    mkdir $workDir;

    my $trimReads = "$prjDir/trimReads.fasta.gz";
    my $pmm4 = "$workDir/pm.m4.gz";
    my $volDir = "$workDir/pm_dir";

    my $options = %$cfg{"ASM_OVLP_OPTIONS"};
    my $binPath = %$env{"BinPath"};
    my $threads = %$cfg{"THREADS"};
    my $volFile = "$volDir/reads_info.txt";


    my $jobMkVol = Job->new (
        name => "altr_mk_vol",
        ifiles => [$trimReads],
        ofiles => [$volFile],
        gfiles => [$volFile],
        mfiles => ["$workDir/asm_read_list.txt"],
        cmds => ["echo $trimReads > $workDir/asm_read_list.txt",
                 "$binPath/oc2mkdb $volDir $workDir/asm_read_list.txt"],
        msg => "making vol for aligning trimmed reads",
    );

    my $jobAlVol = Job->new (
        prefunc => sub($) {
            my ($job) = @_;
        
            my $count = getFileFirstItem($volFile, 0);
            
            for (my $i=0; $i<$count; $i=$i+1) {
                my $jobSub = Job->new(
                    name => "altr_al_vol_$i", 
                    ifiles => [$volFile],
                    ofiles => ["$volDir/pm_result_$i"],
                    gfiles => ["$volDir/pm_result_$i"],
                    mfiles => [],
                    cmds => ["$binPath/oc2asmpm $options -t $threads $volDir $i $volDir/pm_result_$i"],
                    msg => "aligning volumn $i for aligning trimmed reads",
                );
                push @{$job->pjobs}, $jobSub;
                push @{$job->ofiles}, "$volDir/pm_result_$i";
            }

        },

        name => "altr_al_vol",
        ifiles => [$volFile],
        ofiles => [], # prefunc
        mfiles => [],
        pjobs => [], # prefunc
        msg => "aligning vol for aligning trimmed reads",
    );

    my $jobCatVol = Job->new (
        prefunc => sub ($) {
            my ($job) = @_;
            my $subPmStr = join(" ", @{$jobAlVol->ofiles});
            $job->ifiles($jobAlVol->ofiles);
            $job->cmds(["cat $subPmStr | $binPath/pigz  -f -p $threads - -c > $pmm4"]);
        },

        name => "altr_cat_vol",
        ifiles => [], # prefunc
        ofiles => [$pmm4],
        gfiles => [$pmm4],
        mfiles => [],
        cmds => [], # prefunc
        msg => "catenating vol for aligning trimmed reads",
    );

    serialRunJobs($env, $cfg, Job->new(
        name => "altr_job",
        ifiles => [$trimReads],
        ofiles => [$pmm4],
        mfiles => [$volDir, "$volDir/pm_result_*"],
        jobs => [$jobMkVol, $jobAlVol, $jobCatVol],
        msg => "aligning trimmd reads",
    ));
}

sub runTrimAlignReads($$) {
    my ($env, $cfg) = @_;

    if (%$cfg{"TRIM_METHOD"} eq "accurate") {
        runTrimAccurate($env, $cfg);
        runAlignReads($env, $cfg);
    } elsif (%$cfg{"TRIM_METHOD"} eq "accurate0") {
        runTrimAccurate0($env, $cfg);
        runAlignReads($env, $cfg);
    }else {
        runTrimBasesFast($env, $cfg);
    
    }

}

sub runAssemble($$) {
    my ($env, $cfg) = @_;

    my $prjDir = %$env{"WorkPath"} . "/" .%$cfg{"PROJECT"};
    my $workDir = "$prjDir/4-fsa";
    mkdir $workDir;
 
    my $script = "$prjDir/scripts/assemble.sh";
    my $pmm4 = "$prjDir/3-assembly/pm.m4.gz";
    my $trimReads = "$prjDir/trimReads.fasta.gz";
    my $filterm4 = "$workDir/filter.m4";
    my $contigs = "$workDir/contigs.fasta";

    my $binPath = %$env{"FsaBinPath"}; 
    my $threads = %$cfg{"THREADS"};
    my $filterOptions = %$cfg{"FSA_OL_FILTER_OPTIONS"};
    if (%$cfg{"GENOME_SIZE"}) {
        $filterOptions = $filterOptions . " --genome_size=" . %$cfg{"GENOME_SIZE"};
    }
    my $assembleOptions = %$cfg{"FSA_ASSEMBLE_OPTIONS"};
    my $job = Job->new(
        name => "ass_job",
        ifiles => [$pmm4, $trimReads],
        ofiles => [$filterm4, $contigs],
        gfiles => [$filterm4, $contigs],
        mfiles => [],
        cmds => ["$binPath/fsa_ol_filter $pmm4 $filterm4 --thread_size=$threads --output_directory=$workDir $filterOptions\n", 
                 "$binPath/fsa_assemble $filterm4 --read_file=$trimReads --thread_size=$threads --output_directory=$workDir $assembleOptions"],
,
        msg => "assembling",
    );
    serialRunJobs($env, $cfg, $job);
}

sub runAlignContigs($$) {
    my ($env, $cfg) = @_;
    
    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};
    my $workDir = "$prjDir/5-align_contigs";
    mkdir($workDir);

    my $rawreadList = "$prjDir/1-consensus/raw_reads/raw_read_list.txt"; # The garbage reads have been removed
    my $contigs = "$prjDir/4-fsa/contigs.fasta";
    my $ctg2ctg = "$workDir/ctg2ctg.m4a";
    my $rawread2ctg = "$workDir/rawread2ctg.m4a";

    my $binPath = %$env{"BinPath"};
    my $threads = %$cfg{"THREADS"};
    my $volDir = "$workDir/pm_dir";
    my $volFile = "$volDir/reads_info.txt";

    # ctg2ctg short names
    my $ctg2ctgVolDir = "$workDir/temp";
    mkdir($ctg2ctgVolDir);

    my $jobCtg2ctg = Job->new(
        name => "r2c_ctg2ctg",
        ifiles => [$contigs],
        ofiles => ["$ctg2ctg.gz"],
        mfiles => [$ctg2ctgVolDir],
        gfiles => ["$ctg2ctg.gz"],
        cmds => ["$binPath/oc2SplitCtgs $contigs $ctg2ctgVolDir/ctg_reads.fasta",
                 "echo \"$ctg2ctgVolDir/ctg_reads.fasta\"> $ctg2ctgVolDir/ctg_read_list.txt",
                 "$binPath/oc2mkdb $ctg2ctgVolDir $ctg2ctgVolDir/ctg_read_list.txt",
                 "$binPath/oc2pm -t $threads -j 0 $ctg2ctgVolDir $ctg2ctgVolDir/tmp_candidates.txt",
                 "$binPath/oc2FixCanInfo $ctg2ctgVolDir/ctg_reads.fasta $ctg2ctgVolDir/tmp_candidates.txt $ctg2ctgVolDir/candidates.txt",
                 "$binPath/oc2ctgpm -t $threads $contigs $ctg2ctgVolDir/candidates.txt $ctg2ctg",
                 "$binPath/pigz  -f -p $threads $ctg2ctg"],
        msg => "algining contigs to contigs",
    );

    my $jobMkVol = Job->new (

        name => "r2c_mk_vol",
        ifiles => [$rawreadList],
        ofiles => [$volFile],
        gfiles => [$volFile],
        mfiles => ["$workDir/asm_read_list.txt"],
        cmds => ["$binPath/oc2mkdb $volDir $rawreadList"],
        msg => "making vol for aligning rawreads to contigs",
    );

    my $jobAlVol = Job->new (
        prefunc => sub($) {
            my ($job) = @_;
            my $count = getFileFirstItem($volFile, 0);
            
            for (my $i=0; $i<$count; $i=$i+1) {
                my $jobSub = Job->new(
                    name => "r2c_al_vol_$i", 
                    ifiles => [$volFile, $contigs],
                    ofiles => ["${rawread2ctg}_$i"],
                    gfiles => ["${rawread2ctg}_$i"],
                    mfiles => [],
                    cmds => ["$binPath/oc2rm_worker -t $threads $volDir $contigs ${rawread2ctg}_$i -mn $i $count"],
                    msg => "aligning volumn $i for aligning rawreads to contigs",
                );
                push @{$job->pjobs}, $jobSub;
                push @{$job->ofiles}, "${rawread2ctg}_$i";
            }
        },

        name => "r2c_al_vol",
        ifiles => [$volFile, $contigs],
        ofiles => ["$ctg2ctg.gz"], # prefunc
        mfiles => [],
        pjobs => [$jobCtg2ctg], # prefunc
        msg => "aligning vol for aligning rawreads to contigs",
    );

    my $jobCatVol = Job->new (
        prefunc => sub ($) {
            my ($job) = @_;
            my $subPmStr = join(" ", @{$jobAlVol->ofiles}[1..$#{$jobAlVol->ofiles}]);
            push @{$job->ifiles}, @{$jobAlVol->ofiles}[1..$#{$jobAlVol->ofiles}];
            $job->cmds(["cat $subPmStr | $binPath/pigz  -f -p $threads - -c > $rawread2ctg.gz"]);
        },

        name => "r2c_cat_vol",
        ifiles => [], # prefunc
        ofiles => ["$rawread2ctg.gz"],
        gfiles => ["$rawread2ctg.gz"],
        mfiles => [],   #prefunc
        cmds => [], # prefunc
        msg => "catenating vol for aligning rawreads to contigs",
    );

    my $job = Job->new(
        name => "r2c_job",
        ifiles => [$contigs, $rawreadList],
        ofiles => ["$ctg2ctg.gz", "$rawread2ctg.gz"],
        mfiles => [$ctg2ctgVolDir, "${rawread2ctg}_*", $volDir],
        jobs => [$jobMkVol, $jobAlVol, $jobCatVol],
        msg => "aligning rawreads and contigs to contigs",
    );

    serialRunJobs($env, $cfg, $job);
}

sub runBridgeContigs($$) {
    my ($env, $cfg) = @_;

    my $prjDir = %$env{"WorkPath"} . "/" .%$cfg{"PROJECT"};
    my $workDir = "$prjDir/6-bridge_contigs";
    mkdir($workDir);

    my $ctg2ctg = "$prjDir/5-align_contigs/ctg2ctg.m4a.gz";
    my $rawread2ctg = "$prjDir/5-align_contigs/rawread2ctg.m4a.gz";
    my $contigs = "$prjDir/4-fsa/contigs.fasta";
    my $rawreadList = "$prjDir/1-consensus/raw_reads/raw_read_list.txt"; # The garbage reads have been removed
    my $bridged = "$workDir/bridged_contigs.fasta";
    my $readinfo= "$prjDir/4-fsa/readinfos.gz";

    my $binPath = %$env{"BinPath"}; 
    my $threads = %$cfg{"THREADS"};
    my $options = %$cfg{"FSA_CTG_BRIDGE_OPTIONS"};

    my $job = Job->new(
        name => "br_job",
        ifiles => [$ctg2ctg, $rawread2ctg, $contigs, $rawreadList],
        ofiles => [$bridged],
        gfiles => [$bridged],
        mfiles => [],
        cmds => ["$binPath/fsa_ctg_bridge  $rawreadList $contigs $rawread2ctg $bridged --readinfo_fname=$readinfo --ctg2ctg_file=$ctg2ctg  --thread_size=$threads --output_directory=$workDir $options"],
        msg => "bridging contigs",
    );
    serialRunJobs($env, $cfg, $job);
}


sub runPolishContigs($$$$$$$) {
    my ($env, $cfg, $name, $contigs, $contigTiles, $reads, $output) = @_;

    my $prjDir = %$env{"WorkPath"} . "/" .%$cfg{"PROJECT"};
    my $workDir = "$prjDir/$output";
    mkdir($workDir);

    my $binPath = %$env{"BinPath"}; 
    my $thread = %$cfg{"THREADS"};

    my $volDirReads = "$workDir/pac_reads";
    my $volFileReads = "$volDirReads/reads_info.txt";
    my $volDirContigs = "$workDir/pac_contigs";
    my $volFileContigs = "$volDirContigs/reads_info.txt";


    my $rd2ctg = "$workDir/rd2ctg.m4";
    my $rd2ctgFiltered = "$workDir/rd2ctg_filtered.m4";
    my $polishedContigs = "$workDir/polished_contigs.fasta";

    my $readList = "$workDir/read_list.txt";
    my $contigList = "$workDir/contig_list.txt";

    my $jobMkVol = Job->new (

        name => "${name}_mk_vol",
        ifiles => [$reads, $contigs],
        ofiles => [$volFileReads, $volFileContigs],
        gfiles => [$volFileReads, $volFileContigs],
        mfiles => ["$workDir/asm_read_list.txt"],
        cmds => ["echo  \"$contigs\" > \"$contigList\" && $binPath/oc2mkdb $volDirContigs $contigList", 
                 "echo  \"$reads\" > \"$readList\" && $binPath/oc2mkdb $volDirReads $readList"],
        msg => "making vol for polishing contigs",
    );

    my $jobAlVol = Job->new (
        prefunc => sub($) {
            my ($job) = @_;
            my $count = getFileFirstItem($volFileReads, 0);
            
            for (my $i=0; $i<$count; $i=$i+1) {
                my $jobSub = Job->new(
                    name => "${name}_al_vol_$i", 
                    ifiles => [$volFileReads, $contigs],
                    ofiles => ["${rd2ctg}_$i"],
                    gfiles => ["${rd2ctg}_$i"],
                    mfiles => [],
                    cmds => ["$binPath/oc2rm_worker -t $thread -u 1 -z 20 $volDirReads $contigs ${rd2ctg}_$i -mn $i $count"],
                    msg => "aligning vol $i for polishing contigs",
                );
                push @{$job->pjobs}, $jobSub;
                push @{$job->ofiles}, "${rd2ctg}_$i";
            }
        },

        name => "${name}_al_vol",
        ifiles => [$volFileReads, $contigs],
        ofiles => [], # prefunc
        mfiles => [],
        pjobs => [], # prefunc
        msg => "aligning vol for polishing contigs",
    );

    my $jobCatVol = Job->new (
        prefunc => sub ($) {
            my ($job) = @_;
            my $subPmStr = join(" ", @{$jobAlVol->ofiles}[0..$#{$jobAlVol->ofiles}]);
            push @{$job->ifiles}, @{$jobAlVol->ofiles}[0..$#{$jobAlVol->ofiles}];
            $job->cmds(["cat $subPmStr > $rd2ctg"]);
        },

        name => "${name}_cat_vol",
        ifiles => [], # prefunc
        ofiles => [$rd2ctg],
        gfiles => [$rd2ctg],
        mfiles => [],   #prefunc
        cmds => [], # prefunc
        msg => "catenating vol for polishing contigs",
    );
    
    my $jobPolish = Job->new(
        name => "${name}_cns",
        ifiles => [$rd2ctg], # prefunc
        ofiles => [$polishedContigs],
        gfiles => [$polishedContigs],
        mfiles => [],   #prefunc
        cmds => ["$binPath/filter_m4 $volDirReads $contigTiles $rd2ctg $rd2ctgFiltered",
                 "$binPath/pm4 $volDirContigs $rd2ctgFiltered",
                 "$binPath/ctgcns $volDirReads $volDirContigs $thread $polishedContigs"], # prefunc
        msg => "polishing contigs",
    );

    my $job = Job->new(
        name => "${name}_job",
        ifiles => [$contigs, $contigTiles, $reads],
        ofiles => [$polishedContigs],
        mfiles => [],
        jobs => [$jobMkVol, $jobAlVol, $jobCatVol, $jobPolish],
        msg => "polishing contigs",
    );

    serialRunJobs($env, $cfg, $job);

}

sub statReadN50($$$$) {
    my ($env, $cfg, $seq, $msg) = @_;

    my $prjDir = %$env{"WorkPath"} . "/" .%$cfg{"PROJECT"};

    my $binPath = %$env{"BinPath"}; 

    plgdInfo("N50 of $msg: $seq");
    my $cmd = "$binPath/fsa_rd_tools n50 --ifname $prjDir/$seq";
    system($cmd);
}

my %cfg = ();
my %env = ();


sub cmdCorrect($) {
    my ($fname) = @_;

    %cfg = loadNecatConfig($fname);
    %env = loadNecatEnv(\%cfg);
    initializeNecatProject(\%cfg);

    my $prjDir = %env{"WorkPath"} . "/" .%cfg{"PROJECT"};
    my $isGz = %cfg{"COMPRESS"};
    my $finalCnsReads = $isGz ? "1-consensus/cns_final.fasta.gz" : "1-consensus/cns_final.fasta";

    runConsensus(\%env, \%cfg);

    statReadN50(\%env, \%cfg, $finalCnsReads, "corrected reads");

}

sub cmdAssemble($) {
    my ($fname) = @_;

    %cfg = loadNecatConfig($fname);
    %env = loadNecatEnv(\%cfg);
    initializeNecatProject(\%cfg);

    my $prjDir = %env{"WorkPath"} . "/" .%cfg{"PROJECT"};

    runConsensus(\%env, \%cfg);

    runTrimAlignReads(\%env, \%cfg);

    runAssemble(\%env, \%cfg);

    if (%cfg{"POLISH_CONTIGS"} == 1 or %cfg{"POLISH_CONTIGS"} eq "true") {
        runPolishContigs(\%env, \%cfg, "plctg0", "$prjDir/4-fsa/contigs.fasta", 
                     "$prjDir/4-fsa/contig_tiles", "$prjDir/trimReads.fasta.gz", "4-fsa");
        
        statReadN50(\%env, \%cfg, "4-fsa/polished_contigs.fasta", "polished contigs");
    } else {
    
        statReadN50(\%env, \%cfg, "4-fsa/contigs.fasta", "contigs");
    }
}

sub cmdBridge($) {
    my ($fname) = @_;

    %cfg = loadNecatConfig($fname);
    %env = loadNecatEnv(\%cfg);
    initializeNecatProject(\%cfg);

    my $prjDir = %env{"WorkPath"} . "/" .%cfg{"PROJECT"};

    runConsensus(\%env, \%cfg);

    runTrimAlignReads(\%env, \%cfg);

    runAssemble(\%env, \%cfg);

    runAlignContigs(\%env, \%cfg);
    runBridgeContigs(\%env, \%cfg);
    
    if (%cfg{"POLISH_CONTIGS"} == 1 or %cfg{"POLISH_CONTIGS"} eq "true") {
        runPolishContigs(\%env, \%cfg, "plctg1", "$prjDir/6-bridge_contigs/bridged_contigs.fasta", 
                     "$prjDir/4-fsa/contig_tiles", "$prjDir/trimReads.fasta.gz", "6-bridge_contigs");
    
        statReadN50(\%env, \%cfg, "6-bridge_contigs/polished_contigs.fasta", "polished contigs");
    } else {
    
        statReadN50(\%env, \%cfg, "6-bridge_contigs/bridged_contigs.fasta", "bridged contigs");
    }
}


sub cmdConfig($$) {
    my ($fname, $config) = @_;

    open(F, "> $fname") or die; 
    foreach my $item (@$config) {
        if ($item->[2] ne "") {
            print F "# $item->[2]\n";
        }
        print F "$item->[0]=$item->[1]\n";
    }
    close(F);
}


sub usage() {
    print "Usage: necat.pl correct|assemble|bridge|config cfg_fname\n".
          "    correct:     correct rawreads\n" .
          "    assemble:    generate contigs\n" .
          "    bridge:      bridge contigs\n" .
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
        } elsif ($cmd eq "bridge") {
            cmdBridge($cfgfname);
        } elsif ($cmd eq "config") {
            cmdConfig($cfgfname, \@defaultConfig);
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
    my $exitcode = $?;
    stopRunningScripts(\%env, \%cfg);
    exit($exitcode);
}
