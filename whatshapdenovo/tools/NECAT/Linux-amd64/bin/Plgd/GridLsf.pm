package Plgd::GridLsf;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(detectLsf submitScriptLsf stopScriptLsf checkScriptLsf);

use strict;

use File::Basename;
use Plgd::Utils;

sub detectLsf () {
    my $path = `which bsub 2> /dev/null`;
    $path = trim($path);

    if (not $path eq "") {
        plgdInfo("Found LSF, which is $path");
        return "LSF"
    } else {
        return undef;
    }
}

sub submitScriptLsf($$$$) {
    plgdWarn("TODO: The code for Lsf isn't tested");
    my ($script, $thread, $memory, $options) = @_;

    my $jobName = basename($script);

    my $cmd = "bsub ";
    $cmd = $cmd . " -J $jobName";                                           # name
    $cmd = $cmd . " -R span[hosts=1] -n $thread"  if ($thread > 0);         # thread
    $cmd = $cmd . " -M $memory" if ($memory > 0);                           # memory
    $cmd = $cmd . " -o $script.log";                                        # output
    $cmd = $cmd . " -e $script.log";                                        # output
    $cmd = $cmd . " $options";                                              # other options
    $cmd = $cmd . " $script";                                               # script
    plgdInfo("Sumbit command: $cmd");    
    my $result = `$cmd`;

    my @items = split(" ", $result);
    if (scalar @items >= 2) {
        return $items[1];
    } else {
        plgdInfo("Failed to sumbit command");
    }
}


sub stopScriptLsf($) {
    plgdWarn("TODO: The code for Lsf isn't tested");
    
    my ($job) = @_;
    my $cmd = "bkill $job";
    plgdInfo("Stop script: $cmd");
    `$cmd`;
}

sub checkScriptLsf($$) {
    plgdWarn("TODO: The code for Lsf isn't tested");
    my ($script, $jobid) = @_;
    my $state = "";
    open(F, "bjobs -l $jobid |");
    while (<F>) {
        my @items = split(" ", $_);
        if (scalar @items >= 3 and $items[0] eq $jobid) {
            if ($items[2] eq "RUN") {
                $state = "R"; 
            } elsif ($items[2] eq "PEND") {
                $state = "Q";
            } elsif ($items[2] eq "PROV") {
                $state = "Q";
            } elsif ($items[2] eq "PSUSP") {
                $state = "Q";
            } elsif ($items[2] eq "USUSP") {
                $state = "Q";
            } elsif ($items[2] eq "SSUSP") {
                $state = "Q";
            } elsif ($items[2] eq "DONE") {
                $state = "C";
            } elsif ($items[2] eq "EXIT") {
                $state = "C";
            } elsif ($items[2] eq "UNKWN") {
                $state = "";
            } elsif ($items[2] eq "WAIT") {
                $state = "Q";
            } elsif ($items[2] eq "ZOMBI") {
                $state = "";
            } else {
                $state = "";
            }
            last;
        }
        
    }
    close(F);
    return $state;
}



