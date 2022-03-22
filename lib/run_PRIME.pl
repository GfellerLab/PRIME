
################
# Written by David Gfeller.

# The product is provided free of charge for academic users and, therefore, on an "as is" basis, without warranty of any kind.

# For any question or commercial use, please ask david.gfeller@unil.ch.

# Copyright (2022) David Gfeller.
################



use strict;
use Getopt::Long;


my ($alleles, $input, $dir, $output_file, $lib_dir, $MixMHCpred_dir);

GetOptions ("alleles=s" => \$alleles,    # allele
	    "input=s"   => \$input,      # input file
	    "dir=s"   => \$dir,      # current path
	    "output=s"   => \$output_file,      # output directory
	    "lib=s"   => \$lib_dir,      # lib dir
	    "mix=s"   => \$MixMHCpred_dir)      # MixMHCpred package
    or die("Error in command line arguments\n");


##########
# Check MixMHCpred version
##########

my $MixMHCpred_version=`$MixMHCpred_dir -h | head -1`;

$MixMHCpred_version =~ s/\r?\n$//;

if(index($MixMHCpred_version, "MixMHCpred2.1")>=0 || index($MixMHCpred_version, "MixMHCpred2.0.2")>=0){
    print "\n######\n";
    print "WARNING:\n";
    print "It appears that you are using an old version of MixMHCpred ($MixMHCpred_version):\n";
    print "$MixMHCpred_dir\n";
    print "Make sure you use version 2.2 or above.\n";
    print "#######\n\n";
}


##########
# Check allele name
##########


my @allele_input=split(",", $alleles);

my @allele_list=();
my $h;
my $nh=scalar @allele_input;
my $N=20;

foreach $h (@allele_input){

    my $th=$h;
    #For human alleles, simplify the nomenclature
    if(substr($th, 0, 4) eq "HLA-"){
	$th=substr($th, 4, length($th)-4); # Remove the "HLA-"
    }
    my $i=index($th, "*");
    if($i>0){
	$th=substr($th, 0, $i).substr($th, $i+1, length($th)); # Remove the "*"
    }
    if( substr($th, 3, 1) eq ":" && length($th)==6 ){
	$th=substr($th, 0, 3).substr($th, 4, length($th)-3); # Remove the ":"
    }
    if(substr($th, 0, 3) eq "H-2"){
	$th="H2".substr($th, 3,  length($th)-2); # Make it compatible with mouse alleles (future work)
    } 
    push @allele_list, $th;
}


###########
# Check input file
###########

my @peptide=();
my @seq=();
my $l;
my $p;
my $s;

open IN, "$input", or die "No input\n";
while($l=<IN>){
    $l =~ s/\r?\n$//;
    if(substr($l, 0, 1) ne ">" && $l ne ""){
	push @peptide, $l;
    }
}
close IN;

if(scalar @peptide == 0){
    die "Empty peptide file\n";
}

my $Lmin=8;
my $Lmax=14;

my %map=("A", 0, "C", 1, "D", 2, "E", 3, "F", 4, "G", 5, "H", 6, "I", 7, "K", 8, "L", 9, "M", 10, "N", 11, "P", 12, "Q", 13, "R", 14, "S", 15, "T", 16, "V", 17, "W", 18, "Y", 19);
my ($ct1, $ct2, $Lp);
my @lg=();
my @lg_pres=(); # Keep track whether peptides of this lenght are found in input
for(my $le=$Lmin; $le<=$Lmax; $le++){
    $lg_pres[$le]=0;
}

foreach $p (@peptide){
    @seq=split('', $p);
    $ct2=0;
    foreach $s (@seq){
	if(!exists $map{$s}){
	    die "Error in input sequence - unknown amino acids $s: $p\n";
	} 
	$ct2++;
    }
    if($ct2 < $Lmin || $ct2 > $Lmax){
	die "Incompatible peptide length: $p\t$ct2. \nOnly peptides of length $Lmin to $Lmax are supported\n";
    } 
}

#######################################
# Load information about the alleles
#######################################

my @allele_list_pred=();

my %pres=();
my %maph=();

open IN, "$lib_dir/alleles_mapping.txt", or die;
while($l=<IN>){
    $l =~ s/\r?\n$//;
    if($l ne ""){
	my @a=split(' ', $l);
	$maph{$a[0]}=$a[1];
	$pres{$a[0]}=1;
    }
}
close IN;

for(my $i=0; $i<$nh; $i++){
    if(!exists $pres{$allele_list[$i]}){
	die "Predictions cannot be run for $allele_list[$i]\n";
    } else{
	push @allele_list_pred,$maph{$allele_list[$i]};
	
    }
}


################################
# Run the binding predictions
################################

my $rd=int(rand(1000000));


#Run MixMHCpred

my $th=$allele_list_pred[0];
if(scalar @allele_list_pred > 1){
    for(my $i=1; $i<scalar @allele_list_pred; $i++){
	$th=$th.",".$allele_list_pred[$i];
    }
}
system("$MixMHCpred_dir -i $input -o $lib_dir/../tmp/MixMHCpred_$rd.txt -a $th");

###########
# Run PRIME
###########

my $nh=scalar @allele_list;
my $Np=scalar @peptide;
my $all_list=$allele_list[0];
my $all_list_pred=$allele_list_pred[0];

for(my $i=1; $i<$nh; $i++){
    $all_list=$all_list." ".$allele_list[$i];
    $all_list_pred=$all_list_pred." ".$allele_list_pred[$i];
}

#print "$lib_dir/PRIME.x $lib_dir $rd $nh $Np $all_list $all_list_pred $output_file $MixMHCpred_dir $input\n";
if(-e $output_file){
    print "Overriding existing output file: $output_file\n";
    system("rm $output_file");
}
system("$lib_dir/PRIME.x $lib_dir $rd $nh $Np $all_list $all_list_pred $output_file $MixMHCpred_dir $input");

if(! -e $output_file){
    print "PRIME failed...\n";
}
system("rm $lib_dir/../tmp/MixMHCpred_$rd.txt");
