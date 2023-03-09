#!/usr/bin/perlq  
 
###### 

use Getopt::Long;
use Pod::Usage;

my $refseq; # gene annotation - UCSC
my $output; # length of sequence read
my $minisoform; # minimum number of isoforms per gene

GetOptions('r=s'=>\$refseq,'o=s'=>\$output,'m=i'=>\$minisoform);

if((!($refseq))||(!($output))){
    pod2usage();
}

open(REF, $refseq);
open (RRR, ">$output");

######### load transcript annotation
my %knownGene = ();
while (<REF>) {
  chomp($_);
  my @transcript = split(/\t/);
  if($transcript[1] =~ /NR/){
      next;
  }
  my $name  = $transcript[12];
  my $tran =  $transcript[1];  
     $knownGene{$name}{$tran}{"chrom"} = $transcript[2];
     $knownGene{$name}{$tran}{"strand"} = $transcript[3];
     $knownGene{$name}{$tran}{"txStart"} = $transcript[4];
     $knownGene{$name}{$tran}{"txEnd"} = $transcript[5];
     $knownGene{$name}{$tran}{"exonCount"} = $transcript[8];
     $knownGene{$name}{$tran}{"exonStarts"} = $transcript[9];
     $knownGene{$name}{$tran}{"exonEnds"} = $transcript[10];
}


########## load isoform annotation
my %isoGene = ();
my %isoStart = ();
my %isoEnd = ();

foreach my $name  (keys %knownGene )
{
  foreach my $tran (keys %{$knownGene{$name}})
  {

       my $i_start = $knownGene{$name}{$tran}{"txStart"};
       my $i_end = $knownGene{$name}{$tran}{"txEnd"};
          $isoGene{$name} = $isoGene{$name}.$tran.",";  
       if($isoStart{$name} == NULL )
         {$isoStart{$name} =  $i_start;}
       else
       {
           if($isoStart{$name} > $i_start){ $isoStart{$name} = $i_start; }
       }
       if($isoEnd{$name} == NULL )
         {$isoEnd{$name} =  $i_end;}
       else
       {
          if($isoEnd{$name} < $i_end){ $isoEnd{$name} = $i_end;}
       }
  }
}


################## process isoform information
foreach my $ID  (keys %isoGene )
{
   
    my @genename = split(/,/, $isoGene{$ID});
    my $size = @genename;
    my $i_chrom = $knownGene{$ID}{$genename[0]}{"chrom"};
    my $i_strand = $knownGene{$ID}{$genename[0]}{"strand"};
    my $i_start = $isoStart{$ID};
    my $i_end = $isoEnd{$ID};
    my %ISO_Index = ();

    if($size >= $minisoform)  ### you can specify the number of the isoform per gene here
    {
       print RRR "$ID\t$i_chrom\t$i_strand\t$i_start\t$i_end\t";
       for(my $j=0; $j<= $#genename; $j++)
       {
          my $name = $genename[$j];
          print RRR "$name,";
          my @start = split(/,/, $knownGene{$ID}{$name}{"exonStarts"});
          my @end = split(/,/, $knownGene{$ID}{$name}{"exonEnds"}); 
          for (my $ijk=0; $ijk<= $#start; $ijk++)
          {
               my $sss = $start[$ijk];
               my $eee = $end[$ijk];
               for (my $abc=$sss; $abc<=$eee; $abc++)
               {$ISO_Index{$abc}{$j} = 1;}
          } # ijk
       } # j 

       print RRR "\n";   

       my %NEW_EXON =();
       my $CCC =0; 
       my $pre_POS = $i_start-10;
       my @pre_Index = ();
       for(my $j=0; $j <= $#genename; $j++)
       {$pre_Index{$j}=0;}
     
       for my $ijk (sort {$a<=>$b} keys %ISO_Index)
       {
           my $tot =0;
           my @cur_Index=();
           for(my $j=0; $j<=$#genename; $j++)
           {
               my $name = $genename[$j]; 
               my $value = exists $ISO_Index{$ijk}{$j} ? $ISO_Index{$ijk}{$j} : 0;
               $cur_Index[$j] = $value; 
               if($cur_Index[$j] != $pre_Index[$j]) 
               {$tot = $tot +1;} 
           }
           my $move = $ijk - $pre_POS;
           if($move != 1)
           {
               $NEW_EXON{$CCC}{"start"} =  $ijk;
               $NEW_EXON{$CCC}{"Index"} =  [@cur_Index];
               if($CCC > 0)
               {
                  $NEW_EXON{$CCC-1}{"end"} =  $pre_POS;
               }
               @pre_Index = @cur_Index;
               $CCC = $CCC+1;
           } 
           else{        
               if($tot >0)
               {
                  $NEW_EXON{$CCC}{"start"} =  $ijk;
                  $NEW_EXON{$CCC-1}{"end"} =  $ijk-1;
                  $NEW_EXON{$CCC}{"Index"} =  [@cur_Index];
                  @pre_Index = @cur_Index;
                  $CCC = $CCC+1;
               }    
           }
           $pre_POS = $ijk;
       } # ijk
       $NEW_EXON{$CCC-1}{"end"} =  $i_end;

       ### print data structure 
       for my $CCC (sort {$a<=>$b} keys %NEW_EXON) 
       {
           print RRR "$ID\t$i_chrom\t$i_strand\t";
           my $sss = $NEW_EXON{$CCC}{"start"};
           print RRR  "$sss\t";
           my $eee =  $NEW_EXON{$CCC}{"end"};
           print RRR  "$eee\t";
           my $Read_C = 0;
           #for (my $abc=$sss; $abc<=$eee; $abc++)
           #{
           #   my $value = exists $genome{$i_chrom}{$abc} ? $genome{$i_chrom}{$abc} : 0;  
           #  $Read_C =$Read_C+$value;
           #}
           #print RRR "$Read_C\t";
           my @index_get = $NEW_EXON{$CCC}{"Index"}; 
           for my $INDEX ( 0 .. $#{$NEW_EXON{$CCC}{"Index"}}){
           my $index_get = $NEW_EXON{$CCC}{"Index"}[$INDEX];  
           print RRR  "$index_get,"; }
           print RRR  "\n";
	   #print "haha\n";
       } ### for CCC 

    } # if size
 } #i  : Cluster ID -1 




close(SAM);
close(REF);
close(RRR);

=head1 SYNOPSIS

    -r ---RefSeqAnnotation file

    -o ---The file name that you want to save the results
