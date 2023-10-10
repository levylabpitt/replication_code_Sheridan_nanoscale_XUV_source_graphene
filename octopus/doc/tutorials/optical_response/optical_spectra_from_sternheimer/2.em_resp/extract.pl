#!/bin/perl
use Math::Trig;
$c = 137.035999679;

$ls = `ls -d em_resp/freq*`;
@list = ($ls =~/em\_resp\/freq_([0-9.]*)/g);
@freqs = sort {$a <=> $b} @list;

print "# Energy (eV)    Re alpha         Im alpha    Cross-section (A^2)\n";
format =
@###.###      @####.#######   @####.########   @####.########
$energy, $Re_alpha    ,  $Im_alpha,         $cross_section
.

foreach(@freqs)
{
    if ($_ eq "0.0000")
    {
        $Im_alpha = 0;
    }
    else
    {
        $crossfile = `cat em_resp/freq_$_/cross_section`;
        @crossbits = split(' ', $crossfile);
        $energy = $crossbits[26];
        $cross_section = $crossbits[27]; # isotropic average
        $Im_alpha = $c * $cross_section / (4 * pi * $energy);
    }
 
    $alphafile = `cat em_resp/freq_$_/alpha`;
    @alphabits = split(' ', $alphafile);
    $Re_alpha = $alphabits[15];
 
    write;
}
