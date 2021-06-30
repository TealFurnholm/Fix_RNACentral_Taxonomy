#use warnings;

#Get RNAcentral data
qx{ wget http://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/id_mapping.tsv.gz};
qx{ zcat id_mapping.tsv.gz | awk -F'\t' '{print $1"\t"$4"\t"$5}' | sort -u > rnacentral_ids.txt};


#Fix and deduplicate RNAcentral sequences 
qx{wget http://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_species_specific_ids.fasta.gz};
qx{comics reformat --fixjunk iupacton overwrite=t in=rnacentral_species_specific_ids.fasta.gz out=rnacentral_iupacfix.fasta.gz};
$cmd = "grep -oiP \"^\\S+\"";
qx{zcat rnacentral_iupacfix.fasta.gz | $cmd > rnacentral_hedrfix.fasta};
qx{gzip rnacentral_hedrfix.fasta};
qx{comics dedupe sort=length absorbrc=f absorbmatch=f exact=f threads=20 absorbcontainment=f overwrite=t in=rnacentral_hedrfix.fasta.gz out=rnacentral_sort.fasta.gz};
qx{comics dedupe sort=length absorbcontainment=t mergenames=t threads=20 mergedelimiter=+ exact=f overwrite=t in=rnacentral_sort.fasta.gz out=rnacentral_dedup.fasta.gz};

print "INPUT TAXONOMY\n";
open(INTAX, "TAXONOMY_DB_2021.txt")||die "unable to open TAXONOMY_DB_2021.txt: $!\n";
while(<INTAX>){
	if($_ !~ /\w/){next;}
        $_ =~ s/[\r\n]//;
        $_=uc($_);
        @stuff = split("\t", $_, -1);
        $tid =  shift(@stuff);
	$lin = join(";",@stuff);
	$PHY{$tid}=$lin;
}

print "INPUT RNACENTRAL IDS\n";
open(INIDS, "rnacentral_ids.txt")||die "unable to open rnacentral_ids.txt: $!\n";
while(<INIDS>){
	if($_ !~ /\w/){next;}
        $_ =~ s/[\r\n]//;
        $_=uc($_);
        (my $id, my $tid, my $type) = split("\t", $_, -1);
	$ID_TIDS{$id}{$tid}=$type;
	if($id !~ /^U/){print "weird id $id\n";}
}


$/=">";
print "INPUT RNA SEQS\n";
$inrna = "rnacentral_dedup.fasta.gz";
open(INRNA, "gunzip -c $inrna |") or die "cannot open $inrna: $!";
open(OUTSEQ,">", "rnacentral_clean.fasta")||die "unable to open rnacentral_clean.fasta: $!\n";
while(<INRNA>){
	if($_ !~ /\w/){next;}
	$_=uc($_);
	@stuff=split("\n",$_);
	$header=shift(@stuff);
	$seq=join("",@stuff);
	$seq =~ s/[^A-Z]+//g;
	if(length($seq) < 100 || length($seq) > 10000){next;}
	%TYPES=(); %TIDS=();
	@IDS=split('\+',$header);
	$top='';
	foreach my $xid (@IDS){
		if($top eq ''){$top=$xid;}
		$xid=~/^([^\_]+)/;
		$id=$1; 
		foreach my $tid (keys %{$ID_TIDS{$id}}){
			if($PHY{$tid}=~/^(BACTERIA|ARCHAEA|MONA|EUKARYOTA)/){
				$type=$ID_TIDS{$id}{$tid};
				$TIDS{$tid}++;
				$TYPES{$type}++;
			}
		}
	}
	
	if(keys %TIDS < 1){
		foreach my $xid (@IDS){
                	$xid=~/^([^\_]+)/;
                	$id=$1;
			foreach my $tid (keys %{$ID_TIDS{$id}}){
				if($PHY{$tid}=~/\w/){
					$type=$ID_TIDS{$id}{$tid};
					$TIDS{$tid}=1; 
					$TYPES{$type}++;
				}
			}
		}
	}

	#MAKE LCA
	@PHYL=(); @TIDS=(); 
	foreach my $tid (sort{$TIDS{$b}<=>$TIDS{$a}} keys %TIDS){ 
		push(@PHYL,$PHY{$tid}); 
		push(@TIDS,$tid); 
		$pc=@PHYL; if($pc>=25){last;} 
	}
	$LCA=MakeLCA(@PHYL);
	$tids=join(";",@TIDS);

	#REDUCE TYPES
        @TYPES=(); $max=0;
        foreach my $type (sort{$TYPES{$b}<=>$TYPES{$a}} keys %TYPES){
		if($TYPES{$type}>$max){$max=$TYPES{$type};} 
		if($TYPES{$type}>$max*0.75){push(@TYPES,$type);} 
	}
	$types=join(";", @TYPES);

	print OUTSEQ ">".$top."|".$types."|".$tids."|".$LCA."\n$seq\n";
}

#zip up the file
qx{gzip rnacentral_clean.fasta};


sub MakeLCA{
        my @ARR=();
        %seen=();
        @array1 = @_;
        @array1 = grep { !$seen{$_}++ } @array1;

        #get the kingdoms, JIC lineage is NCA
        %LET=();
        foreach my $lin (@array1){
                if($lin !~ /^(BACTERIA|ARCHAEA|MONA|EUKARYOTA)/i){next;}
                $lin =~ /^(.)/;
                $LET{$1}=1;
        }
        @LET=();
        foreach my $let (sort(keys %LET)){push(@LET,$let);}
        $let=join("",@LET);

        $len1 = @array1;
        if($len1 == 1){$LCA = $array1[0]; }
        elsif($len1 > 1){
                $first = $array1[0];
                @levels = split(";", $first);
                for($i=0; $i<=$#levels; $i++){
                        $alevel=$levels[$i];
                        @matched = grep(/\Q$alevel\E/i, @array1);
                        $len2 = @matched;
                        if($len2 == $len1){push(@ARR, $alevel);}
                        else{last;}
                }
                $len3 = @ARR;
                if($len3 > 1){$LCA = join(";", @ARR);}
                elsif($len3==1){$LCA = $ARR[0];}
                else{$LCA = "NCA"; }
        }
        else{$LCA = "NCA"; }

        #add kingdoms to NCA
        if($LCA eq "NCA"){ $LCA.="-".$let; }
        return($LCA);
}
	
	

