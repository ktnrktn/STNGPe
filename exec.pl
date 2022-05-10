#!/usr/bin/perl
$numproc = 8;

$control = 1; #if $control is set to 1, the normal condition is set. Otherwise, the abnormal condition is set.

if($control == 1){
    #DA modulation
    $DAmod_D1 = 1.5;
    $DAmod_D2 = 0.5;

    #DA depletion
    $DAdep_g2s = 1.0; #GPe -> STN
    $DAdep_g2g = 1.0; #GPe -> GPe
    $DAdep_c2s = 1.0; #Ctx -> STN
    $DAdep_c2p = 1.0; #Ctx -> PGPe
    $DAdep_c2a = 1.0; #Ctx -> AGPe
    $DAdep_s2p = 1.0; #STN -> PGPe
    $DAdep_s2a = 1.0; #Str -> AGPe
    $DAdep_fistr = 1.0; #iStr
    $DAdep_fdstr = 1.0; #dStr

    $frate_ctx = 1.5;
    $frate_istr = 0.35;
    $frate_dstr = 0.3;

    $DAdep_amp_ctx = 1.0;
    $DAdep_amp_str = 1.0;
    
}else{
    #DA modulation
    $DAmod_D1 = 1.0;
    $DAmod_D2 = 1.0;

    #DA depletion
    $DAdep_g2s = 1.5; #GPe -> STN
    $DAdep_g2g = 1.4; #GPe -> GPe
    $DAdep_c2s = 0.5; #Ctx -> STN
    $DAdep_c2p = 1.0; #Ctx -> PGPe
    $DAdep_c2a = 1.0; #Ctx -> AGPe
    $DAdep_s2p = 0.7; #STN -> PGPe
    $DAdep_s2a = 3.0; #Str -> AGPe
    $DAdep_fistr = 3.0; #iStr
    $DAdep_fdstr = 2.0; #dSTr

    $frate_ctx = 1.5;
    $frate_istr = 1.2;
    $frate_dstr = 0.8;

    $DAdep_amp_ctx = 0.7;
    $DAdep_amp_str = 1.0;
}

#scaling of STN-GPe connection
$scale = 1.0;

#number of connections of background inputs
$num_ctx2s = 100;
$num_ctx2p = 100;
$num_ctx2a = 100;
$num_istr = 100;
$num_dstr = 100;

#amplitudes of Ctx-STN and Str-GPe inputs
$amp_ctx    = $DAdep_amp_ctx*1.0;  # amplitude of Ctx rhythm
$amp_str    = $DAdep_amp_str*1.0;  # amplitude of Str rhythm

for($i=1; $i<=$#ARGV; $i++){

    if($ARGV[$i] eq '-a'){ # for amplitude
	$amp_ctx = $ARGV[$i+1];
	$amp_str = $ARGV[$i+2];
    }

    if($ARGV[$i] eq '-n'){ # for the number of background inputs
	$num_ctx2s = $num_ctx2s*$ARGV[$i+1];
	$num_ctx2p = $num_ctx2p*$ARGV[$i+1]*0.4;
	$num_ctx2a = $num_ctx2a*$ARGV[$i+1]*0.5;
	$num_istr = $num_istr*$ARGV[$i+2];
	$num_dstr = $num_dstr*$ARGV[$i+3];
    }
    
    if($ARGV[$i] eq '-s'){ # for striatal inputs
	$scale = $ARGV[$i+1];
	$frate_istr = $scale*$frate_istr;
	$frate_dstr = $scale*$frate_dstr;
    }
}

#simulation parameters
$time = 15000; # simulation time

$gA_stn2pgpe = $DAdep_s2p*5.4;     # from STN to PGPe;
$gA_stn2agpe = 1.4;                # from STN to PGPe;

$gA_ctx2stn = $DAdep_c2s*3.75;      # from CTX to STN
$gA_ctx2pgpe = $DAdep_c2s*3.75;     # from CTX to PGPe
$gA_ctx2agpe = $DAdep_c2s*3.75*1.5; # from CTX to AGPe

$gN_stn2pgpe = $DAdep_s2p*2.7;      # NMDA from STN to PGPe
$gN_stn2agpe = 0.9;                 # NMDA from STN to AGPe

$gG_pgpe2pgpe = $DAdep_g2g*4.38;    # from PGPe to PGPe
$gG_pgpe2agpe = $DAdep_g2g*6.67;    # from PGPe to AGPe
$gG_agpe2pgpe = $DAdep_g2g*1.33;    # from PGPe to PGPe
$gG_agpe2agpe = $DAdep_g2g*1.33;    # from PGPe to AGPe

$gG_gpe2stn = $DAdep_g2s*6.82;      # from PGPe to STN
$gG_str2pgpe = 3.33;                # from D2-STR to PGPe
$gG_str2agpe = $DAdep_s2a*0.33;     # from D1-STR to AGPe

$freq_ctx   = 1.0;                  # frequency of CTX rhythm
$freq_str   = 1.0;                  # frequency of STR rhythm
$phase_diff = 0.0;  # phase difference between CTX and STR

$scale_D1  = $DAmod_D1;             # scale for U of STN
$scale_D2  = $DAmod_D2;             # scale for U of GP

$seed_syn   = -1;   # random seed for synaptic parameters
$seed_init  = -10;  # random seed for initial states
$seed_run   = -100; # random seed for run

    
system "nohup mpirun -np $numproc ../stngpMPI $time $gA_stn2pgpe $gA_stn2agpe $gA_ctx2stn $gA_ctx2pgpe $gA_ctx2agpe $gN_stn2pgpe $gN_stn2agpe $gG_pgpe2pgpe $gG_pgpe2agpe $gG_agpe2pgpe $gG_agpe2agpe $gG_gpe2stn $gG_str2pgpe $gG_str2agpe $frate_ctx $frate_istr $frate_dstr $num_ctx2s $num_ctx2p $num_ctx2a $num_istr $num_dstr $amp_ctx $freq_ctx $amp_str $freq_str $phase_diff $scale_D1 $scale_D2 $seed_syn $seed_init $seed_run > nohup.log";

exit;



