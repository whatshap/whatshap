package Plgd::Script;

require Exporter;

@ISA = qw(Exporter);
@EXPORT =qw(isScriptDone writeScript writeScripts runScriptLocal isScriptSucc isScriptPatternSucc getScriptReturn waitScript checkScripts);

use strict;
use Plgd::Utils;

sub isScriptDone($) {
    my ($script) = @_;
    return ((-e $script . ".done") and ((stat($script))[9] <= (stat($script.".done"))[9]));
}

# write script command to file
# 0: file path
# 1: command string
sub writeScript {
    my ($fname, $env, @cmds) = @_;
    
    plgdDebug("Write Script, $fname");
    #if (! -e $fname) {
    {
        open(F, "> $fname") or die;
        print F "#!/bin/bash\n\n";
        print F "$env";

        print F "retVal=0\n";

        my $wrapCmds = wrapCommands(@cmds);
        print F "$wrapCmds\n";

        print F "echo \$retVal > $fname.done\n";
        close(F);

        chmod(0755 & ~umask(), $fname);
    } 
}

# 
sub writeScripts($$$) {
    my ($pattern, $env, $cmds) = @_;
    
    plgdDebug("Write Scripts. The pattern is $pattern");

    my $size = scalar @$cmds;
    my @fnames = ();
    for (my $i = 0; $i < scalar @$cmds; $i = $i + 1) {
        $fnames[$i] = sprintf($pattern, $i);
        writeScript($fnames[$i], $env, @$cmds[$i]);
    }
    return @fnames;
}

sub isScriptSucc($) {
    my ($script) = @_;
    return ((-e $script . ".done") and ((stat($script))[9] <= (stat($script.".done"))[9]) and (getScriptReturn($script) == 0));
}

sub isScriptPatternSucc($$) {
    my ($pattern, $count) = @_;
    for (my $i = 0; $i<$count; $i = $i + 1) {
        my $script = sprintf($pattern, $i);
        if (not isScriptSucc($script)) {
            return 0;
        }
    }
    return 1;
}


# wait the scripts is over
# 1: script files
# 2: waiting time
# 3: interval time
# 4: be silent
sub waitScript($$$$) {
    my ($script, $waitTime, $interval, $silent) = @_;
 
    my $startTime = time();
 
    while (not isScriptDone($script)) {
        if ($waitTime > 0 and time() - $startTime > $waitTime) {
            return 0;
        }

        if (not $silent) {
            plgdInfo("Wait script fininshed $script");
        }
        sleep($interval);
    }
    return 1;
}

sub waitScripts($$$) {
    my ($scripts, $waitTime, $silent) = @_;

    my $sleepTime = 1;
    my $startTime = time();

    my $done = 0;    
    until ($done) {
        $done = 1;
        foreach my $s (@$scripts) {
            if (not isScriptDone($s)) {

                if (not $silent) {
                    plgdInfo("Wait script fininshed $s");
                }
                $done = 0;
                last;
            }
        }
        if ($waitTime > 0 and time() - $startTime > $waitTime) {
            return 0;
        }
        sleep($sleepTime);
        $sleepTime = $sleepTime*2 < 60 ? $sleepTime*2 : 60;
    }
    return 1;
}

sub waitCheckScript($$$) {
    my ($script, $waitTime, $silent) = @_;
   
    if (waitScript($script, $waitTime, 10, $silent)) {
        checkScripts($script);
    } else {
        plgdError("Failed to wait script $script");
    }
}

sub checkScripts {
    foreach my $script (@_) {
        my $retCode = getScriptReturn($script);
        if ($retCode != 0) {
            plgdError("Failed to run script, $retCode, $script");
        } 
    }
}


sub queryScript($) {
    foreach my $script (@_) {
        my $retCode = getScriptReturn($script);
        if ($retCode != 0) {
            plgdError("Failed to run script, $retCode, $script");
        } 
    }
}

sub runScriptLocal($) {
    my ($script) = @_;
    plgdInfo("Run script: $script 2>&1 |tee $script.log");
    my $r = system($script . " 2>&1 | tee $script.log");
    if ($r != 0) {
        `echo $r > $script.done`; # 
    }
    #checkScripts($script);
}


sub getScriptReturn($) {
    my ($script) = @_;
    my $retVal = 127;
    if (-e "$script.done") {
        open F, "< $script.done" or die;
        while(<F>){
            $retVal = 0 + $_; # Transfer string to number;
            last;
        }
    }
    return $retVal;
}

sub wrapCommands {
    my $str = "";
    foreach my $c (@_) {
        #$str = $str . 
        #       "if [ \$retVal -eq 0 ]; then\n" .
        #       "  $c\n" .
        #       "  temp_result=\$?\n" .
        #       "  if [ \$retVal -eq 0 ]; then\n".
        #       "    retVal=\$temp_result\n" .
        #       "  fi\n" .
        #       "fi\n";
        $str = $str . 
               "if [ \$retVal -eq 0 ]; then\n" .
               "  $c\n" .
               "  temp_result=(\${PIPESTATUS[*]})\n" .
               "  for i in \${temp_result[*]} \n" .
               "  do\n" .
               "    if [ \$retVal -eq 0 ]; then\n" .
               "      retVal=\$i\n" .
               "    else\n" .
               "      break\n" .
               "    fi\n" .
               "  done\n".
               "fi\n";
    }
    return $str;
}

