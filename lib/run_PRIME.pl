
################
# Written by David Gfeller.

# The product is provided free of charge for academic users and, therefore, on an "as is" basis, without warranty of any kind.

# For any question or commercial use, please ask david.gfeller@unil.ch.

# Copyright (2019) David Gfeller.
################

use strict;
use Getopt::Long;


my ($alleles, $input, $dir, $output_file, $lib_dir, $MixMHCpred_dir, $PredBinding, $NetMHCpan_dir, $affinity_threshold);

GetOptions ("alleles=s" => \$alleles,    # allele
	    "input=s"   => \$input,      # input file
	    "dir=s"   => \$dir,      # current path
	    "output=s"   => \$output_file,      # output directory
	    "lib=s"   => \$lib_dir,      # lib dir
	    "mix=s"   => \$MixMHCpred_dir,      # MixMHCpred package
	    "net=s"   => \$NetMHCpan_dir,      # NetMHCpan package
	    "thr=s"   => \$affinity_threshold,      # Threshodl on affinity
	    "pred=s"   => \$PredBinding)      # Use NetMHCpan instead of MixMHCpred
    or die("Error in command line arguments\n");

##########
# Check allele name
##########

my @allele_input=split(",", $alleles);


my @allele_list=();
my $h;
my $nh=scalar @allele_input;

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
    push @allele_list, $th;
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
	    #print "Error in input sequence - unknown amino acids $s: $p\n";
	    die "Error in input sequence - unknown amino acids $s: $p\n";
	    #exit;
	} 
	$ct2++;
    }
    if($ct2 < $Lmin || $ct2 > $Lmax){
	#print "Incompatible peptide length: $p\t$ct2. \nOnly peptides of length $Lmin to $Lmax are supported\n";
	die "Incompatible peptide length: $p\t$ct2. \nOnly peptides of length $Lmin to $Lmax are supported\n";
	#exit;
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
    my @a=split(' ', $l);
    $maph{$a[0]}=$a[1];
    $pres{$a[0]}=1;
}
close IN;


open IN, "$lib_dir/alleles.txt", or die;
while($l=<IN>){
    $l =~ s/\r?\n$//;
    $maph{$l}=$l;
    $pres{$l}=1;
}
close IN;

if($PredBinding eq "MixMHCpred"){
    
    for(my $i=0; $i<$nh; $i++){
	if(!exists $pres{$allele_list[$i]}){
	    die "Predictions cannot be run for $allele_list[$i]\n";
	} else{
	    push @allele_list_pred,$maph{$allele_list[$i]};
	    
	}
    }
} elsif($PredBinding eq "NetMHCpan"){
    
    my %pres=();
    open IN, "$lib_dir/allelenames_NetMHCpan.txt", or die;
    while($l=<IN>){
	$l =~ s/\r?\n$//;
	my @a=split(' ', $l);
	$pres{$a[0]}=1;
    }
    close IN;
    
    foreach $h (@allele_list){
	my $th="";
	if(length($h)==5 && (substr($h, 0, 1) eq "A" || substr($h, 0, 1) eq "B" || substr($h, 0, 1) eq "C" || substr($h, 0, 1) eq "G" )){
	    $th="HLA-".substr($h, 0, 3).":".substr($h, 3, 2);
	} else {
	    $th=$h;
	}
	if(!exists $pres{$th}){
	    die "Predictions cannot be run for $th\n";
	}
	push @allele_list_pred, $h;	
    }
}


#######################################
# Load the coefficients of the logistic regression
#######################################

my @coef=();
open IN, "$lib_dir/coef.txt", or die;
while($l=<IN>){
    
    $l =~ s/\r?\n$//;
    my @a=split(' ', $l);
    push @coef, $a[1];
    
}


################################
# Run the binding predictions
################################

my $rd=int(rand(1000000));


#Run MixMHCpred

if($PredBinding eq "MixMHCpred"){

    my $th=$allele_list_pred[0];
    if(scalar @allele_list_pred > 1){
	for(my $i=1; $i<scalar @allele_list_pred; $i++){
	    $th=$th.",".$allele_list_pred[$i];
	}
    }
    system("$MixMHCpred_dir -i $input -o $lib_dir/../tmp/MixMHCpred_$rd.txt -a $th");
    
} elsif($PredBinding eq "NetMHCpan"){

    my $th="";
    for(my $i=0; $i<scalar @allele_list_pred; $i++){
	$h=$allele_list_pred[$i];
	if(length($h)==5 && (substr($h, 0, 1) eq "A" || substr($h, 0, 1) eq "B" || substr($h, 0, 1) eq "C" || substr($h, 0, 1) eq "E" || substr($h, 0, 1) eq "G" )){
	    $th="HLA-".substr($h, 0, 3).":".substr($h, 3, 2);
	} else {
	    $th=$h;
	}
	#print "$h\t$th\n";
	system("NetMHCpan -f $input -a $th -p -BA >  $lib_dir/../tmp/NetMHCpan_$h\_$rd.txt");
    }
}

###########################
# Run the PRIME predictions
###########################

#Load the thresholds for different Pvalues
my @thresh=([]);
my @pval=();
my $Npval;
for(my $i=0; $i<scalar @allele_list_pred; $i++){
    
    $Npval=0;
    if(-e "$lib_dir/PerRank/$allele_list_pred[$i].txt"){
	open IN, "$lib_dir/PerRank/$allele_list_pred[$i].txt", or die;
    } elsif (-e "$lib_dir/PerRank/$maph{$allele_list_pred[$i]}.txt") {
	open IN, "$lib_dir/PerRank/$maph{$allele_list_pred[$i]}.txt", or die; # This is the case if we use NetMHCpan with alleles not in pre-computed
    } else {
	open IN, "$lib_dir/PerRank/A0201.txt", or die; # This is the case when using NetMHCpan with alleles that cannot be mapped to any of the precomputed ones.
    } while($l=<IN>){
	$l =~ s/\r?\n$//;
	my @a=split(' ', $l);
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
my @fr=();
my @aa=();
my @pos=();
my $th=$allele_input[0];

for(my $i=1; $i<scalar(@allele_input); $i++){
    $th=$th.",".$allele_input[$i];
}

#print "$PredBinding\t@allele_input\t@allele_list\t@allele_list_pred\n";


open OUT, ">$output_file";
printf OUT "####################\n";
printf OUT "# Output from PRIME (v1.0)\n";
printf OUT "# Alleles: $th\n";
printf OUT "# Affinity predictions: $PredBinding";
if($PredBinding eq "MixMHCpred"){
    print OUT " ($MixMHCpred_dir)\n";
} else {
    print OUT " ($NetMHCpan_dir)\n";
}
printf OUT "# PRIME is freely available for academic users.\n";
printf OUT "# Private companies should contact eauffarth\@licr.org at the Ludwig Institute for Cancer Research Ltd for commercial licenses.\n";
printf OUT "#\n";
printf OUT "# To cite PRIME, please refer to:\n";
printf OUT "# Schmidt et al. Prediction of neo-epitopes’ immunogenicity reveals TCR recognition determinants and immunoediting mechanisms, BioRxiv (2020).\n";
printf OUT "#\n";
printf OUT "####################\n";

printf OUT "Peptide\t%%Rank_bestAllele\tScore_bestAllele\t%%RankBinding_bestAllele\tBestAllele";
for(my $i=0; $i<scalar @allele_input; $i++){
    printf OUT "\t%%Rank_$allele_input[$i]\tScore_$allele_input[$i]\t%%RankBinding_$allele_input[$i]";
}
print OUT "\n";

my @pep;
my @rank_list=([]);

my $ct;

if($PredBinding eq "MixMHCpred"){
    open IN, "$lib_dir/../tmp/MixMHCpred_$rd.txt", or die;
    for(my $i=0; $i<12; $i++){
	$l=<IN>;
    }
    $ct=0;
    while ($l=<IN>) {
	$l =~ s/\r?\n$//;
	my @a=split(' ', $l);
	push @pep, $a[0];
	for (my $j=5; $j<scalar @a; $j=$j+2) {
	    push @{$rank_list[$ct]}, $a[$j];
	   
	}
	$ct++;
    }
    close IN;
    
} elsif ($PredBinding eq "NetMHCpan"){
    
    for(my $i=0; $i<scalar @allele_list_pred; $i++){
	$ct=0;
	my $pos_rank;
	open IN, "$lib_dir/../tmp/NetMHCpan_$allele_list_pred[$i]_$rd.txt", or die;
	while($l=<IN>){
	    $l =~ s/\r?\n$//;
	    my @a=split(' ', $l);
	    if($a[0] eq "Pos"){
		for(my $j=0; $j<scalar (@a); $j++){
		    if($a[$j] eq "%Rank"){
			$pos_rank=$j;
		    }
		}
	    }
	    if($a[0]==1){
		push @pep, $a[2];
		$rank_list[$ct][$i] = $a[$pos_rank];
		$ct++;
	    }
	}
	close IN;
    }   
}

#Precompute the positions for each allele and each length
my @pos_all=([[]]);
for (my $j=0; $j<scalar @allele_list_pred; $j++) {
    for(my $n=$Lmin; $n<=$Lmax; $n++){
	@{$pos_all[$j][$n]}=find_pos($allele_list_pred[$j], $n);
    }
}

for(my $n=0; $n<$ct; $n++){
   
    @score_list=();
    @Pval_list=();
    my $lg=length($pep[$n]);
    for (my $j=0; $j<scalar @allele_list_pred; $j++) {
	
	if($rank_list[$n][$j] <= $affinity_threshold){
	    $al_pred=$allele_list_pred[$j];
	    
	    ####
	    # Compute the PRIME score with each allele and the corresponding %Rank
	    ####
	    
	    for (my $i=0; $i<$N; $i++) {
		$fr[$i]=0;
	    }
	    @aa=split('', $pep[$n]);
	    @pos=@{$pos_all[$j][$lg]};
	    
	    for (my $i=0; $i<$N; $i++) {
		$fr[$i]=0;
	    }
	    
	    my $sc=(scalar @pos);
	    foreach $p (@pos) {
		$fr[$map{$aa[$p]}]=$fr[$map{$aa[$p]}]+1.0/$sc;
	    }
	    
	        
	    $score=$coef[0];
	    for (my $i=0; $i<$N; $i++) {
		$score=$score+$fr[$i]*$coef[$i+1];
	    }
	    $score=$score-log($rank_list[$n][$j])*$coef[$N+1];
	    $score=1/(1+exp(-$score));
	} else {
	    $score=0;
	}
	    
	push @score_list, $score;
	    
	#Compute the %Rank
	my $k=0;
	while ($score<$thresh[$j][$k] && $k<$Npval) {
	    $k++;
	}

	if ($k==$Npval) {
	    $Pval=100;
	} else {
	    $Pval=$pval[$k];
	}
	push @Pval_list, $Pval;
	#print "$pep[0]\t$score\t$Pval\n";
	#die;
	
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
    
    printf OUT "$pep[$n]\t%.2f\t%.6f\t%.2f\t$allele_input[$min_pos]",  $Pval_list[$min_pos], $score_list[$min_pos], $rank_list[$n][$min_pos];
    for (my $i=0; $i<scalar @Pval_list; $i++) {
	printf OUT "\t%.2f\t%.6f\t%.2f", $Pval_list[$i], $score_list[$i], $rank_list[$n][$i];
    }
    print OUT "\n";

}

close IN;
close OUT;

if($PredBinding eq "MixMHCpred"){
   system("rm $lib_dir/../tmp/MixMHCpred_$rd.txt");
} elsif($PredBinding eq "NetMHCpan"){
   system("rm $lib_dir/../tmp/NetMHCpan_*_$rd.txt");
}

sub find_pos(){

    my $al_pred=$_[0];
    my @pos=();
    my $le=$_[1];
    if($al_pred eq "B0801" || $al_pred eq "B1401" || $al_pred eq "B1402" || $al_pred eq  "B3701" || $al_pred eq "A6802" || $al_pred eq "H2-Dd" || $al_pred eq "H-2-Dd" || $al_pred eq "H2-Db" || $al_pred eq "H-2-Db" || $al_pred eq "H2-Kb" || $al_pred eq "H-2-Kb"){  # 6
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
