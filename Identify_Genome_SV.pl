use strict;


###################################################
## perl Identify_Genome_SV.pl blat_MCScanX.result.txt
open IN,"@ARGV[0]";
open OUT,">@ARGV[0].gene.sv.txt";
open OUT1,">@ARGV[0].small.nocolinearity.scaffold";
open OUT2,">@ARGV[0].false.positive.txt";
open OUT3,">@ARGV[0].duplication";
open OUT4,">@ARGV[0].large.inversion";
my @line;
my %chr_number;
my %nn_id_info;
my %small_scaffold;
my %big_scaffold;
my %colinearity;
my %c_sort;
my %nocolinearity;
my %deletion;
my %inversion;
my %duplications;
my %translocations;
my %false_positive;

while(<IN>){
	chomp;
	@info=split(/\t/,$_);
	@line=(@line,$_);
	if(exists $chr_number{@info[10]}){$chrnumber=$chr_number{@info[10]};$chrnumber+=1;$chr_number{@info[10]}=$chrnumber;}else{$chr_number{@info[10]}=1;}	
	$nn_id_info{@info[1]}=$_;
}

## small scaffold - deletions
@key_chr_number=keys %chr_number;
@key_nn_id_info=keys %nn_id_info;
for($i=0;$i<@key_chr_number;$i++){if($chr_number{@key_chr_number[$i]} < 4){$small_scaffold{@key_chr_number[$i]}=1;print OUT1 "@key_chr_number[$i]\tsmall_scafold\n";}else{$big_scaffold{@key_chr_number[$i]}=1;}}
for($i=0;$i<@key_nn_id_info;$i++){
	@nn_small=split(/\t/,$nn_id_info{@key_nn_id_info[$i]});
	if(exists $small_scaffold{@nn_small[10]}){if(@nn_small[18] <= 0.4){print OUT "@nn_small[1]\tdeletion\tsmall_scaffold\tdeletion\n";$deletion{@nn_small[1]}=1;}}
				}

## w/ or w.o colinearity scaffold - deletions 
my @line_big;
for($i=0;$i<@line;$i++){@original_list=split(/\t/,@line[$i]);if(exists $small_scaffold{@original_list[10]}){}else{@line_big=(@line_big,@line[$i]);}}


@id=split(/\t/,@line_big[0]);
my @jc; my $CID=@id[10]; my %jcgroup; my $jcmax="NONE";my $mainsca;
for($i=0;$i<@line_big;$i++){
	@jclist=split(/\t/,@line_big[$i]);
	if(@jclist[10] eq $CID){@jc=(@jc,@line_big[$i]);}
	else{
		for($n=0;$n<@jc;$n++){@jcid=split(/\t/,@jc[$n]);
			if(exists $jcgroup{@jcid[0]}){$addjcid=$jcgroup{@jcid[0]};$addjcid+=1;$jcgroup{@jcid[0]}=$addjcid;}else{$jcgroup{@jcid[0]}=1;}	
			}
		@jcnum=keys %jcgroup;
		for($m=0;$m<@jcnum;$m++){$genenum=@jc;$jcp=$jcgroup{@jcnum[$m]}/$genenum;if($jcp > 0.5 && $jcp > $jcmax){$mainsca=@jcnum[$m];$jcmax=$jcp;}}
		if($jcmax > 0.5){$colineartiy{$CID}=$mainsca;print OUT1 "$CID\t$mainsca\tcolinearity\n";}else{$nocolinearity{$CID}=1;print OUT1 "$CID\tnocolinearity\n";}
		$jcmax="NONE";%jcgroup=();@jc=();$mainsca=();$CID=@jclist[10];@jc=(@jc,@line_big[$i]);
	}
}
## last one
		for($n=0;$n<@jc;$n++){@jcid=split(/\t/,@jc[$n]);
                        if(exists $jcgroup{@jcid[0]}){$addjcid=$jcgroup{@jcid[0]};$addjcid+=1;$jcgroup{@jcid[0]}=$addjcid;}else{$jcgroup{@jcid[0]}=1;}
                        }
                @jcnum=keys %jcgroup;
                for($m=0;$m<@jcnum;$m++){$genenum=@jc;$jcp=$jcgroup{@jcnum[$m]}/$genenum;if($jcp > 0.5 && $jcp > $jcmax){$mainsca=@jcnum[$m];$jcmax=$jcp;}}
                if($jcmax > 0.5){$colineartiy{$CID}=$mainsca;print OUT1 "$CID\t$mainsca\tcolinearity\n";}else{$nocolinearity{$CID}=1;print OUT1 "$CID\tnocolinearity\n";}
                $jcmax="NONE";%jcgroup=();@jc=();$mainsca=();$CID=@jclist[10];@jc=(@jc,@line_big[$i]);

my @line_big_col;my @line_big_nocol;
for($i=0;$i<@line_big;$i++){@j_big=split(/\t/,@line_big[$i]);if(exists $nocolinearity{@j_big[10]}){@line_big_nocol=(@line_big_nocol,@line_big[$i]);}else{@line_big_col=(@line_big_col,@line_big[$i]);}}

#nocolinearity scaffold deletions and translocation
for($i=0;$i<@line_big_nocol;$i++){@j_nocol=split(/\t/,@line_big_nocol[$i]);
	if(@j_nocol[18] <= 0.4){print OUT "@j_nocol[1]\tdeletion\tnocolinearity\tdeletion\n";$deletion{@j_nocol[1]}=1;}

	if(@j_nocol[23] eq "#N/A" && @j_nocol[24] eq "#N/A"){
		print OUT "@j_nocol[1]\ttranslocation\tnocolinearity\ttranslocation\n";$translocations{@j_nocol[1]}=1;
	}
	elsif(@j_nocol[23] eq "#N/A" && @j_nocol[24] ne "#N/A"){
		if(@j_nocol[0] eq @j_nocol[19]){
			if(@j_nocol[20] > @j_nocol[3] || @j_nocol[21] < @j_nocol[2]){
				print OUT "@j_nocol[1]\tintra-duplication\tnocolinearity\tduplication\n";
				if(exists $duplications{@j_nocol[14]}){}else{$duplications{@j_nocol[14]}=1;}

			}
			else{print OUT2 "@j_nocol[1]\tfalse_positive\tnocolinearity\tfalse_positive\n";$false_positive{@j_nocol[1]}=1;}
		}
		else{
			print OUT "@j_nocol[1]\tinter-duplication\tnocolinearity\tduplication\n";
			if(exists $duplications{@j_nocol[14]}){}else{$duplications{@j_nocol[14]}=1;}
			
		}
	}
}

#gene sort in colinearity scaffold
@big_col_id=split(/\t/,@line_big_col[0]);
my $s_same=0;$s_dif=0;my @j_sort_col;my $sid=@big_col_id[10];
for($i=0;$i<@line_big_col;$i++){
	@s_col=split(/\t/,@line_big_col[$i]);
	if(@s_col[10] eq $sid){@j_sort_col=(@j_sort_col,@line_big_col[$i]);}
	else{
		for($n=0;$n<@j_sort_col;$n++){
			@jscol=split(/\t/,@j_sort_col[$n]);
			if(@jscol[4] eq @jscol[13]){$s_same++;}else{$s_dif++}
		}
		if($s_same >= $s_dif){$c_sort{$sid}="same";}else{$c_sort{$sid}="dif";}
		@j_sort_col=();$s_same=0;$s_dif=0;$sid=@s_col[10];@j_sort_col=(@j_sort_col,@line_big_col[$i]);
	}
}
#last one
		for($n=0;$n<@j_sort_col;$n++){
                        @jscol=split(/\t/,@j_sort_col[$n]);
                        if(@jscol[4] eq @jscol[13]){$s_same++;}else{$s_dif++}
                }
                if($s_same >= $s_dif){$c_sort{$sid}="same";}else{$c_sort{$sid}="dif";}


## colinearity scaffold sv
## deletions AlignmentRatio <= 0.4
## inter-scaffold translocation
## gene dirction inversion
## intra-scaffold translocation (w/ inversion or w/o inversion)
## duplication
## false_positive
my $sline;
for($i=0;$i<@line_big_col;$i++){
	@j_col=split(/\t/,@line_big_col[$i]);
	if(@j_col[18] <= 0.4){print OUT "@j_col[1]\tdeletion\tcolinearity\tdeletion\n";$deletion{@j_col[1]}=1;}

	if(@j_col[0] ne $colineartiy{@j_col[10]}){
		print OUT "@j_col[1]\tinter-scaffold-translocation\tcolinearity\ttranslocation\n";$translocations{@j_col[1]}=1;
	}
	else{
		if(@j_col[4] eq @j_col[13]){$sline="same";}else{$sline="dif";}
		if($sline ne $c_sort{@j_col[10]}){
			print OUT "@j_col[1]\tinversion\tcolinearity\tinversion\n";$inversion{@j_col[1]}=1;}

		if(@j_col[23] eq "#N/A" && @j_col[24] eq "#N/A"){
			print OUT "@j_col[1]\tintra-scaffold-translocation\tcolinearity\ttranslocation\n";$translocations{@j_col[1]}=1;
		}
		elsif(@j_col[23] eq "#N/A" && @j_col[24] ne "#N/A"){
			if(@j_col[0] eq @j_col[19]){
				if(@j_col[20] > @j_col[3] || @j_col[21] < @j_col[2]){
					print OUT "@j_col[1]\tintra-duplication\tcolinearity\tduplication\n";
					if(exists $duplications{@j_col[14]}){}else{$duplications{@j_col[14]}=1;}
				}
				else{print OUT2 "@j_col[1]\tfalse_positive\tcolinearity\tfalse_positive\n";$false_positive{@j_col[1]}=1;}
			}
			else{
				print OUT "@j_col[1]\tinter-duplication\tnocolinearity\tduplication\n";
				if(exists $duplications{@j_col[14]}){}else{$duplications{@j_col[14]}=1;}
			}
		}
	}
	$sline=();
}

@key_duplications=keys %duplications;
for($i=0;$i<@key_duplications;$i++){
	print OUT3 "@key_duplications[$i]\n";
}

my $star;my $end;my @j_largei;
for($i=0;$i<@line_big_col;$i++){
	@j_for_largei=split(/\t/,@line_big_col[$i]);
	if(exists $inversion{@j_for_largei[1]}){@j_largei=(@j_largei,@line_big_col[$i]);}
	else{
		if(@j_largei > 1){
			@largei_star=split(/\t/,@j_largei[0]);$star=@largei_star[2];
			@largei_end=split(/\t/,@j_largei[-1]);$end=@largei_end[3];
			$largei_loci=$end-$star;
			if($largei_loci > 100000){for($n=0;$n<@j_largei;$n++){print OUT4 "@j_largei[$n]\n";}print OUT4 "\n";}
		}
		$star=();$end=();@j_largei=();
	}
}
