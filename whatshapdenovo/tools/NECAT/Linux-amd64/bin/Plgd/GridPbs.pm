package Plgd::GridPbs;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(detectPbs submitScriptPbs stopScriptPbs checkScriptPbs);

use Cwd;
use File::Basename;

use Plgd::Utils;


our $isPro = "";
our $version = "";
our $VERSION = '1.00';

sub detectPbs() {
    my $path = `which pbsnodes 2> /dev/null`;
    $path = trim($path);

    if (not $path eq "") {

        open(F, "pbsnodes --version 2>&1 |");
        while (<F>) {
            if (m/pbs_version\s+=\s+(.*)/) {
                $isPro   =  1;
                $version = $1;
            }
            if (m/Version:\s+(.*)/) {
                $version = $1;
            }
        }
        close(F);
    
        if ($isPro == 0) {
            plgdInfo("Found PBS/Torque '$version', which is $path");
            return "PBS";
        } else {
            plgdInfo("Found PBS/Pro '$version', which is $path");
            return "PBS";
        }

    } else {
        return undef;
    } 
}


sub submitScriptPbs($$$$) {
    
    my ($script, $thread, $memory, $options) = @_;

    my $jobName = basename($script);

    my $cmd = "qsub -j oe";
    $cmd = $cmd . " -d `pwd`" if ($isPro == 0); 
    $cmd = $cmd . " -N $jobName";                               # name
    $cmd = $cmd . " -l nodes=1:ppn=$thread" if ($thread > 0);   # thread
    $cmd = $cmd . " -l mem=$memory" if ($memory > 0);           # memory
    $cmd = $cmd . " -o $script.log";                            # output
    $cmd = $cmd . " $options";                                  # other options
    $cmd = $cmd . " $script";                                   # script
    plgdInfo("Sumbit command: $cmd");    
    my $result = `$cmd`;

    if (not $result eq "") {
        return trim($result);
    } else {
        plgdInfo("Failed to sumbit command");
    }
}

sub stopScriptPbs($) {
    my ($job) = @_;
    my $cmd = "qdel $job";
    plgdInfo("Stop script: $cmd");
    `$cmd`;
}

sub checkScriptPbs($$) {
    my ($script, $jobid) = @_;
    my $state = "";
    open(F, "qstat |");
    while (<F>) {
        my @items = split(" ", $_);
        if (scalar @items >= 6 and $jobid =~ /$items[0]/) {
            $state = $items[4];
            break;
        }
        
    }
    close(F);
    return $state;
} 



