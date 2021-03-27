package Plgd::GridSge;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(detectSge submitScriptSge stopScriptSge checkScriptSge);


use strict;

use File::Basename;
use Plgd::Utils;


sub detectSge () {
    if (defined($ENV{'SGE_ROOT'})) {
        plgdInfo("Found Sun Grid Engine, which is " . $ENV{'SGE_ROOT'});
        return "SGE";
    } else {
        return undef;
    }
}


sub submitScriptSge($$$$) {

    my ($script, $thread, $memory, $options) = @_;

    my $jobName = basename($script);

    my $cmd = "qsub -cwd";
    $cmd = $cmd . " -N $jobName";                           # name
    $cmd = $cmd . " -pe smp $thread" if ($thread > 0);      # thread
    $cmd = $cmd . " -l vf=$memory" if ($memory > 0);        # memory
    $cmd = $cmd . " -o $script.log -j yes";                 # output
    $cmd = $cmd . " $options";                              # other options
    $cmd = $cmd . " $script";                               # script
    plgdInfo("Sumbit command: $cmd");    
    my $result = `$cmd`;
    my @items = split(" ", $result);
    if (scalar @items >= 3) {
        return $items[2];
    } else {
        plgdInfo("Failed to sumbit command");
    }
}


sub stopScriptSge($) {
    my ($job) = @_;
    my $cmd = "qdel $job";
    plgdInfo("Stop script: $cmd");
    `$cmd`;
}

sub checkScriptSge($$) {
    my ($script, $jobid) = @_;
    my $state = "";
    open(F, "qstat |");
    while (<F>) {
        my @items = split(" ", $_);
        if (scalar @items >= 5 and $items[0] eq $jobid) {
            if (grep {$_ eq $items[4]} ("qw", "hqw", "hRwq")) {
                $state = "Q"; 
            } elsif (grep {$_ eq $items[4]} ("r", "t", "Rr", "Rt")) {
                $state = "R"; 
            } elsif (grep {$_ eq $items[4]} ("s", "ts", "S", "tS", "T", "tT", "Rs", "Rts", "RS", "RtS", "RT", "RtT")) {
                $state = "Q";
            } elsif (grep {$_ eq $items[4]} ("Eqw", "Ehqw", "EhRqw", "dr", "dt", "dRr", "dRt", "ds", "dS", "dT", "dRs", "dRS", "dRT")) {
                $state = "C";
            } else {
                $state = "";
            }
            last;
        }
    }
    close(F);
    return $state;
}
