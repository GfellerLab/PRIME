
################
# Written by David Gfeller.

# The product is provided free of charge for academic users and, therefore, on an "as is" basis, without warranty of any kind.

# For any question or commercial use, please ask david.gfeller@unil.ch.

# Copyright (2019) David Gfeller.
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
# Check allele name
##########

my @allele_input=split(",", $alleles);
my @allele_list=();
my $h;
my $nh=scalar @allele_input;

foreach $h (@allele_input){
    
    if(substr($h, 0, 4) eq "HLA-"){
	$h=substr($h, 4, length($h)-4);
    }
    if(substr($h, 1, 1) eq "*"){
	$h=substr($h, 0, 1).substr($h, 2, length($h)-2);
    }
    if(substr($h, 3, 1) eq ":"){
	$h=substr($h, 0, 3).substr($h, 4, length($h)-3);
    }
    push @allele_list, $h;
}

###########
# Check input file
###########

my @peptide=();
my @peptide_num=([]);
my @seq=();
my $l;
my $p;
my $s;

open IN, "$input", or die "No input\n";
while($l=<IN>){
    $l =~ s/\r?\n$//;
    chomp($l);
    if(substr($l, 0, 1) ne ">" && $l ne ""){
	push @peptide, $l;
    }
}
close IN;

if(scalar @peptide == 0){
    print "Empty peptide file\n";
    exit;
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
	    print "Error in input sequence - unknown amino acids $s: $p\n";
	    exit;
	} 
	$ct2++;
    }
    if($ct2 < $Lmin || $ct2 > $Lmax){
	print "Incompatible peptide length: $p\t$ct2. \nOnly peptides of length $Lmin to $Lmax are supported\n";
	exit;
    } 
}

#######################################
# Load information about the alleles
#######################################


my $t;
my $cond=0;

my %pres=();
my %maph=();

open IN, "$lib_dir/alleles_mapping.txt", or die;
while($l=<IN>){
    $l =~ s/\r?\n$//;
    my @a=split(' ', $l);
    $maph{$a[0]}=$a[1];
}
close IN;


my @a=split("\/", $MixMHCpred_dir);
my $path="";
for(my $i=0; $i<(scalar @a)-1; $i++){
    $path=$path."/".$a[$i];
}


open IN, "$lib_dir/alleles.txt", or die;
while($l=<IN>){
    $l =~ s/\r?\n$//;
    chomp($l);
    $pres{$l}=1;
}
close IN;
my @b;

my @allele_list_pred=();
for(my $i=0; $i<$nh; $i++){
    if(!exists $pres{$allele_list[$i]}){
	print "Predictions cannot be run for $allele_list[$i]\n";
	exit;
    } else{
	$cond=1;
	if(! exists $maph{$allele_list[$i]}){
	    push @allele_list_pred, $allele_list[$i];
	} else {
	    push @allele_list_pred,$maph{$allele_list[$i]};
	}
    }
}

#######################################
# Load the coefficients of the logistic regression
#######################################

my @coef=();
open IN, "$lib_dir/coef.txt", or die;
#open IN, "../../train/coef.txt", or die;
while($l=<IN>){

    chomp($l);
    @a=split(' ', $l);
    push @coef, $a[1];
    
}


################################
# Run the MixMHCpred predictions
################################

my $rd=int(rand(1000000));

while(-e "../tmp/MixMHCpred_$rd.txt"){
    print "fail\n"; 
}


#Run MixMHCpred
my $th=$allele_list_pred[0];
if(scalar @allele_list_pred > 1){
    for(my $i=1; $i<scalar @allele_list_pred; $i++){
	$th=$th.",".$allele_list_pred[$i];
    }
}
system("$MixMHCpred_dir -i $input -o $lib_dir/../tmp/MixMHCpred_$rd.txt -a $th");

###########################
# Run the PRIME predictions
###########################

#Load the thresholds for different Pvalues
my @thresh=([]);
my @pval=();
my $Npval;
for(my $i=0; $i<scalar @allele_list_pred; $i++){
    
    $Npval=0;
    open IN, "$lib_dir/PerRank/$allele_list_pred[$i].txt", or die;
    while($l=<IN>){
	chomp($l);
	@a=split(' ', $l);
	$pval[$Npval]=$a[0];
	$thresh[$i][$Npval]=$a[1];
	$Npval++;
    }
}

my $rank;
my $al;
my $al_pred;
my $N=20;
my $score=0;
my $Pval;
my @Pval_list=();
my @score_list=();
my @rank_list=();
my @fr=();
my @aa=();
my @pos=();
my $h;
open IN, "$lib_dir/../tmp/MixMHCpred_$rd.txt";
open OUT, ">$output_file";
printf OUT "####################\n";
printf OUT "# Output from PRIME (v1.0)\n";
printf OUT "# Alleles: $th\n";
printf OUT "#\n";
printf OUT "# PRIME is freely available for academic users.\n";
printf OUT "# Private companies should contact eauffarth\@licr.org at the Ludwig Institute for Cancer Research Ltd for commercial licenses.\n";
printf OUT "#\n";
printf OUT "# To cite PRIME, please refer to:\n";
printf OUT "# The immunogenicity scale of amino acids. xxx (2019)\n";
printf OUT "#\n";
printf OUT "####################\n";

printf OUT "Peptide\t%%Rank_bestAllele\tScore_bestAllele\t%%RankBinding_bestAllele\tBestAllele";
for(my $i=0; $i<scalar @allele_list; $i++){
    printf OUT "\t%%Rank_$allele_list[$i]\tScore_$allele_list[$i]\t%%RankBinding_$allele_list[$i]";
}
print OUT "\n";

for(my $i=0; $i<12; $i++){
    $l=<IN>;
}

while ($l=<IN>) {

    #if(substr($l, 0, 1) ne "#" && substr($l, 0, 7) ne "Peptide"){

    chomp($l);
    @a=split(' ', $l);
    
    my $t=0;
    @score_list=();
    @Pval_list=();
    @rank_list=();
    for (my $j=5; $j<scalar @a; $j=$j+2) {
	
	$rank=-log($a[$j]);
	#print "$a[$j]\t$rank\n";
	push @rank_list, $a[$j];
	$al_pred=$allele_list_pred[$t];
	    
	####
	# Compute the PRIME score with each allele and the corresponding %Rank
	####
	    
	for (my $i=0; $i<$N; $i++) {
	    $fr[$i]=0;
	}
	@aa=split('', $a[0]);
	@pos=find_pos($al_pred, length($a[0]));
	   
	for (my $i=0; $i<$N; $i++) {
	    $fr[$i]=0;
	}

	my $sc=(scalar @pos);
	foreach $p (@pos) {
	    $fr[$map{$aa[$p]}]=$fr[$map{$aa[$p]}]+1.0/$sc;
	}
	    
	#print "@fr\t$rank\n";
	$score=$coef[0];
	for (my $i=0; $i<$N; $i++) {
	    $score=$score+$fr[$i]*$coef[$i+1];
	}
	$score=$score+$rank*$coef[$N+1];
	$score=1/(1+exp(-$score));
	    
	push @score_list, $score;
	    
	#Compute the %Rank
	my $k=0;
	while ($score<$thresh[$t][$k] && $k<$Npval) {
	    $k++;
	}

	if ($k==$Npval) {
	    $Pval=100;
	} else {
	    $Pval=$pval[$k];
	}
	push @Pval_list, $Pval;
	    
	$t++;
    }
    #######################
    #Find the lowest Pvalue (this will determine which is the predicted allele)
    #This means that we take the PRIME score with each allele, and then take the best allele based on %Rank, unlike MixMHCpred where the best score is used to take the predicted allele, and %Rank are computed for the combination of allele
    #######################
    
    my $min_Pval=1000;
    my $min_pos=-1;
    for (my $i=0; $i<scalar @Pval_list; $i++) {
	if ($Pval_list[$i]<$min_Pval) {
	    $min_Pval=$Pval_list[$i];
	    $min_pos=$i;
		    
	}
    }
    
    printf OUT "$a[0]\t%.2f\t%.6f\t%.2f\t$allele_list[$min_pos]",  $Pval_list[$min_pos], $score_list[$min_pos], $rank_list[$min_pos];
    for (my $i=0; $i<scalar @Pval_list; $i++) {
	printf OUT "\t%.2f\t%.6f\t%.2f", $Pval_list[$i], $score_list[$i], $rank_list[$i];
    }
    print OUT "\n";
    #}
}
close IN;
close OUT;

#system("cp $lib_dir/../tmp/MixMHCpred_$rd.txt $output_dir/MixMHCpred.txt");
system("rm $lib_dir/../tmp/MixMHCpred_$rd.txt");

sub find_pos(){

    my $al_pred=$_[0];
    my @pos=();
    my $le=$_[1];
    if($al_pred eq "B0801" || $al_pred eq "B1401" || $al_pred eq "B1402" || $al_pred eq  "B3701" || $al_pred eq "A6802"){  # 6
	push @pos, 3;
	for(my $j=5; $j<$le-1; $j++){
	    push @pos, $j;
	}
    } elsif($al_pred eq "A0201" || $al_pred eq "A0206" || $al_pred eq "A0211" ){  # 2
	for(my $j=4; $j<$le-1; $j++){
	    push @pos, $j;
	}
    }elsif($al_pred eq "A0203" || $al_pred eq "A0204" || $al_pred eq "A0205" ){   # 3
	push @pos, 4;
	for(my $j=6; $j<$le-1; $j++){
	    push @pos, $j;
	}
    }elsif($al_pred eq "A2501" || $al_pred eq "A2601" ){  # 4
	push @pos, 3;
	push @pos, 4;
	for(my $j=6; $j<$le-1; $j++){
	    push @pos, $j;
	}
    }elsif($al_pred eq "A2902" || $al_pred eq "C1502"){  # 5
	for(my $j=3; $j<$le-3; $j++){
	    push @pos, $j;
	}
	push @pos, $le-2;
    } else{
	for(my $j=3; $j<$le-1; $j++){ # 1
	    push @pos, $j;
	}
    }
    
    return(@pos);
    
}
