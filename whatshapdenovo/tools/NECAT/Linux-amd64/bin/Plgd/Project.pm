package Plgd::Project;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(serialRunJobs parallelRunJobs loadConfig loadEnv initializeProject runScript runScripts detectGrid stopRunningScripts scriptEnv runSingleTask runPatternTask runMultiTask);

use strict;

use FindBin;
use lib $FindBin::RealBin;
use Cwd;
use File::Basename;

use Plgd::Utils;
use Plgd::Script;
use Plgd::GridPbs;
use Plgd::GridLsf;
use Plgd::GridSge;
use Plgd::GridSlurm;

use Class::Struct;

our %running = ();

my $WAITING_FILE_TIME = 60;

sub loadConfig($$) {
    my ($fname, $cfg) = @_;
    open(F, "<$fname") or die "cann't open file: $fname, $!";

    while(<F>) {
        $_ =~ s/^\s+|\s+$//g;
        if (substr($_, 0, 1) ne "#") {
            
            my @items  = split("=", $_, 2);
            $items[1] =~s/^\s*"|"\s*$//g;
            $cfg->{$items[0]} = trim($items[1]);
        }
    }
}


sub loadEnv($) {
    my ($cfg) = @_;
    my %env = {};
    $env{"WorkPath"} = getcwd();
    $env{"OntsaBinPath"} = $FindBin::RealBin;
    $env{"FsaBinPath"} = $FindBin::RealBin;

    if (%$cfg{"USE_GRID"}) {
        detectGrid(\%env);
    }

    $env{"running"} = ();
    return %env;
}

struct Job => {
    prefunc => '$',
    name => '$',
    ifiles => '@',
    ofiles => '@',
    gfiles => '@',
    mfiles => '@',
    cmds => '@',
    jobs => '@',
    pjobs => '@',
    funcs => '@',
    msg => '$',
};

sub serialRunJobs {
    my ($env, $cfg, @jobs) = @_;

    foreach my $job (@jobs) {
        runJob($env, $cfg, $job);
    }
}

sub parallelRunJobs {
    my ($env, $cfg, @jobs) = @_;
    
    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};

    # check which job should be run
    my @running = ();
    my @scripts = ();
    foreach my $job (@jobs) {
        
        if (scalar @{$job->funcs} > 0 || scalar @{$job->jobs} > 0) {
            plgdError("Only cmds can run parallel.");
        }

        my $script = "$prjDir/scripts/" . $job->name . ".sh";
        
        requireFiles(@{$job->ifiles});
        if (filesNewer($job->ifiles, $job->ofiles) or not isScriptSucc($script)) {
            unlink @{$job->ofiles};

            writeScript($script, scriptEnv($env), @{$job->cmds});
            push @scripts, $script;
            push @running, $job;
        } else {
            plgdInfo("Skip ". $job->msg . " for outputs are newer.") if ($job->msg);
        }
    }
    
    
    if (scalar @scripts > 0) {
        foreach my $job (@running) {
            plgdInfo("Parallelly start " . $job->msg . ".") if ($job->msg);
        }

        runScripts($env, $cfg, \@scripts);

        foreach my $job (@running) {

            waitRequiredFiles($WAITING_FILE_TIME, @{$job->ofiles});
            
            if (%$cfg{"CLEANUP"} == 1) {
                deleteFiles(@{$job->mfiles});
            }

            plgdInfo("End " .$job->msg. ".") if ($job->msg);
        }
    }

}

sub runJob ($$$) {
    my ($env, $cfg, $job) = @_;
    $job->prefunc->($job) if ($job->prefunc);
    
    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};
   
    my $script = "$prjDir/scripts/" . $job->name. ".sh";

    requireFiles(@{$job->ifiles});
    if (filesNewer($job->ifiles, $job->ofiles) or not isScriptSucc($script)) {
        deleteFiles(@{$job->gfiles}) if ($job->gfiles); 
        deleteFiles("$script.done"); 

        plgdInfo("Start " . $job->msg . ".") if ($job->msg);

        if (scalar @{$job->cmds} > 0) {
            writeScript($script, scriptEnv($env), @{$job->cmds});
            runScript($env, $cfg, $script);
        } elsif (scalar @{$job->funcs} > 0) {
            foreach my $f (@{$job->funcs}) {
                $f->($env, $cfg);
            }
            echoFile("$script.done", "0");
        } elsif (scalar @{$job->jobs} > 0) {
            foreach my $j (@{$job->jobs}) {
                runJob($env, $cfg, $j);
            }
            echoFile("$script.done", "0");
        } elsif (scalar @{$job->pjobs} > 0) {
            parallelRunJobs($env, $cfg, @{$job->pjobs});
            echoFile("$script.done", "0");
        } else {
            pldgWarn("It is an empty job");
            # die "never come here"
        }

        waitRequiredFiles($WAITING_FILE_TIME, @{$job->ofiles});
        if (%$cfg{"CLEANUP"} == 1) {
            deleteFiles(@{$job->mfiles});
        }

        plgdInfo("End " .$job->msg . ".") if ($job->msg);
    } else {
        plgdInfo("Skip ". $job->msg . " for outputs are newer.") if ($job->msg);
    
    }
}

sub scriptEnv($) {
    my ($env) = @_;

    my $ontsaBinPath = %$env{"OntsaBinPath"};
    my $fsaBinPath = %$env{"FsaBinPath"};

    return "export PATH=$ontsaBinPath:$fsaBinPath:\$PATH\n";
}


sub initializeProject($) {
    my ($configs) = @_;
    mkdir %$configs{"PROJECT"};
    mkdir %$configs{"PROJECT"} . "/scripts";
}



sub runScript($$$) {
    my ($env, $cfg, $script) = @_;
    _runScripts($env, $cfg, $script);
}


sub runScripts($$$) {
    my ($env, $cfg, $scripts) = @_;
    
    _runScripts($env, $cfg, @$scripts);
}

sub _runScripts {
    my ($env, $cfg, @scripts) = @_;
    
    if (%$cfg{"USE_GRID"} eq "true" and %$env{"GridEngine"} ) {
        _runScriptsGrid($env, $cfg, \@scripts);
    } else {
        #foreach my $script (@scripts) {
        #    runScriptLocal($script);
        #}
        _runScriptLocal($env, $cfg, \@scripts);
    }
}

sub _runScriptLocal($$$) {
    my ($env, $cfg, $scripts) = @_;

    my @undone = @$scripts;
    until ((scalar @undone) == 0) {
        my @failed = ();
        foreach my $script (@undone) {
            runScriptLocal($script);
            my $r = getScriptReturn($script);
            if ($r != 0) {
                plgdWarn("Failed to run script, $r, $script");
                push @failed, $script;

                $env->{"ScriptError"} += 1;
                if ($env->{"ScriptError"} > $cfg->{"MAX_SCRIPT_ERROR"}) {
                    plgdError("Reached to maximum number of script errors");
                }
            }
        }

        @undone = @failed;
    }
}

sub detectGrid($) {
    my ($env) =  @_;

    my $r = detectPbs();
    $r = detectLsf() if (not $r);
    $r = detectSge() if (not $r);
    $r = detectSlurm() if (not $r);
    $$env{"GridEngine"} = $r;
}


sub waitScriptsGrid($$$$) {
    my ($env, $cfg, $running, $part) = @_;

    my @scripts = keys %$running;
    
    my @finished = ();
    until (@finished ~~ @scripts) {
        @finished = ();
        foreach my $s (@scripts) {
            my $jobid = $running{$s};
            my $state = checkScriptGrid($env, $cfg, $s, $jobid);
            if ($state eq "" or $state eq "C") {
                if (waitScript($s, 60, 5, 1)) {
                    push @finished, $s;
                } else {
                    plgdWarn("Failed to get script result, id=$jobid, $s");
                    push @finished, $s;
                }
            } else {
                sleep(5);
            }
        }
        last if ($part and @finished > 0);        
    }
    return @finished;
}


sub submitScriptGrid($$$) {
    my ($env, $cfg, $script) = @_;
    if (%$env{"GridEngine"} eq "PBS") {
        return submitScriptPbs($script, %$cfg{"THREADS"}, %$cfg{"MEMORY"}, %$cfg{"GRID_OPTIONS"});
    } elsif (%$env{"GridEngine"} eq "SGE") {
        return submitScriptSge($script, %$cfg{"THREADS"}, %$cfg{"MEMORY"}, %$cfg{"GRID_OPTIONS"});
    } elsif (%$env{"GridEngine"} eq "LSF") {
        return submitScriptLsf($script, %$cfg{"THREADS"}, %$cfg{"MEMORY"}, %$cfg{"GRID_OPTIONS"});
    } elsif (%$env{"GridEngine"} eq "Slurm") {
        return submitScriptSlurm($script, %$cfg{"THREADS"}, %$cfg{"MEMORY"}, %$cfg{"GRID_OPTIONS"});
    } else {
        plgdError("Not support Grid ". %$env{"GridEngine"});
    }
}

sub checkScriptGrid($$$$) {
    my ($env, $cfg, $script, $jobid) = @_;
    
    my $state = ""; # "R" (running), "Q" (Queue), "C" (Complete) ""(Unknown)
    if (%$env{"GridEngine"} eq "PBS") {
        $state = checkScriptPbs($script, $jobid);
    } elsif (%$env{"GridEngine"} eq "SGE") {
        $state = checkScriptSge($script, $jobid);
    } elsif (%$env{"GridEngine"} eq "LSF") {
        $state = checkScriptLsf($script, $jobid); 
    } elsif (%$env{"GridEngine"} eq "Slurm") {
        $state = checkScriptSlurm($script, $jobid)
    } else {
        plgdError("Not support Grid ". %$env{"GridEngine"});
    }

    return $state
}

sub incAndCheckScriptError($$) {
    my ($env, $cfg) = @_;

    $env->{"ScriptError"} += 1;
    if ($env->{"ScriptError"} > $cfg->{"MAX_SCRIPT_ERROR"}) {
        plgdError("Reached to maximum number of script errors");
    }
}

sub _runScriptsGrid($$$) {
    my ($env, $cfg, $scripts) = @_;

    
    my $node = %$cfg{"GRID_NODE"};
    my @undone = @$scripts;

    until ((scalar @undone) == 0) {
        my @failed = ();

        foreach my $s (@$scripts) {
            plgdInfo("Run script $s");
            my $r = submitScriptGrid($env, $cfg, $s);
            if ($r) {
                $running{$s} = $r;

                if ($node > 0 and (keys %running) >= $node) {
                    my @finished = waitScriptsGrid($env, $cfg, \%running, 1);

                    foreach my $i (@finished) {
                        my $r = getScriptReturn($i);
                        if ($r != 0) {
                            plgdWarn("Failed to run script, $r, $i");
                            push @failed, $i;

                            incAndCheckScriptError($env, $cfg);
                        }
                        delete $running{$i};
                    }
                }

            } else {
                plgdWarn("Failed to submit script $s");
                push @failed, $s;
                incAndCheckScriptError($env, $cfg);
            }
        }

        my @finished = waitScriptsGrid($env, $cfg, \%running, 0);
        foreach my $i (@finished) {
            my $r = getScriptReturn($i);
            if ($r != 0) {
                plgdWarn("Failed to run script, $r, $i");
                push @failed, $i;

                incAndCheckScriptError($env, $cfg);
            }
            delete $running{$i};
        }

        @undone = @failed;
    }
    
    
}

sub runScriptsGrid($$$) {
    my ($env, $cfg, $scripts) = @_;

    
    my $node = %$cfg{"GRID_NODE"};

    foreach my $s (@$scripts) {
        plgdInfo("Run script $s");
        my $r = submitScriptGrid($env, $cfg, $s);
        plgdError("Failed to submit script $s") if (not $r);

        $running{$s} = $r;
        my $rsize = keys %running;
	    if ($node > 0 and (keys %running) >= $node) {
            my @finished = waitScriptsGrid($env, $cfg, \%running, 1);
	        foreach my $i (@finished) {
                delete $running{$i};
            }
            checkScripts(@finished);
        }
        
    }
    my @finished = waitScriptsGrid($env, $cfg, \%running, 0);
    foreach my $i (@finished) {
        delete $running{$i};
    }
    checkScripts(@finished);
    
    
}

sub stopScirptGrid($$$$) {
    
    my ($env, $cfg, $script, $jobid) = @_;
    if (%$env{"GridEngine"} == "PBS") {
        return stopScriptPbs($jobid);
    } elsif (%$env{"GridEngine"} == "SGE") {
        return stopScriptSge($jobid);
    } elsif (%$env{"GridEngine"} == "LFS") {
        return stopScriptLfs($jobid);
    } elsif (%$env{"GridEngine"} == "Slurm") {
        return stopScriptSlurm($jobid);
    } else {
        plgdError("Not support Grid ". %$env{"GridEngine"});
    }
}

sub stopRunningScripts($$) {
    my ($env, $cfg) = @_;

    foreach my $i (keys %running) {
        stopScirptGrid($env, $cfg, $i, $running{$i});
        delete $running{$i};
    }

}

1;
