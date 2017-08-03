#!usr/bin/perl
use strict;
use warnings;
use Bio::Graphics;
use Bio::SeqFeature::Generic;
use Bio::Graphics::Glyph::segments;
use Bio::Graphics::Glyph::heat_map;
use Bio::Graphics::Glyph::generic;

my $file= shift @ARGV;	
my $filename =shift @ARGV;
my $exon=shift @ARGV;
my $hmmer=shift @ARGV;
my $one_file=shift @ARGV;
my $cluster=shift @ARGV;
my @files;
my %high;
my %low;
my $difference;
my $highest=0;
if ($one_file)
  {
  @files=split(",",$filename)
  }
else {push(@files,$filename);}
for (@files)
  {
  my $high=0;
  my $low=100000000;
  my $atm=`cut -f 3 $file$_.gff3 | sort -g | tail -n 1`;
  if ($atm>$high){$high=$atm;}
  $atm=`cut -f 2 $file$_.gff3 | sort -g | tail -n 1`;
  if ($atm>$high){$high=$atm;}
  $atm=`cut -f 2 $file$_.gff3 | sort -g | head -n 1`;
  if ($atm<$low){$low=$atm-1;}
  $atm=`cut -f 3 $file$_.gff3 | sort -g | head -n 1`;
  if ($atm<$low){$low=$atm-1;}
  $difference=$high-$low;
  $high{$_}=$difference;
  $low{$_}=$low;
  if ($difference>$highest){$highest=$difference;}
  }
$difference=int($highest/10000+1)*10000;
my $panel=Bio::Graphics::Panel->new(-length => $difference,
				      -width => 900,
				      -pad_left=>100,
				      -pad_right=>50,
				      -key_style => "left");
				      
my $full_length=Bio::SeqFeature::Generic->new(-start=>0,-end=>$difference);
$panel->add_track($full_length,
		  -glyph=>'arrow',
		  -tick=>2,
		  -fgcolor=>'black',
		  -double=>1);
		  
my @track;
for (my $i=0;$i<scalar(@files);$i++)
  {
  my $name=$files[$i];
  if (!$exon)
    {
    my $scaf=`cut -f 1 $file$files[$i].gff3| head -n 1`;
    my $short=substr $files[$i],0,4;
    chomp $scaf;
    $scaf="$short\_$scaf";
    $track[$i]=$panel->add_track(
			  -glyph =>'heat_map',
			  -label => 1,
			  -strand_arrow=>1,
			  -key=>"$scaf",
			  -min_score=>0,
			  -max_score=>50,
			  -start_color=>"green",
			  -end_color=>"red"
			  );
    open my $IN, "<$file$files[$i].gff3";
    $name=substr $name,5;
    while (my $line=<$IN>)
	{
	next if $line=~/^\#/;
	chomp $line;
	my ($scaffold,$start,$end,$dir,$ID,$type,$motifs,$clusters)=split ("\t",$line);
	my $ID_name="$ID";
	if ($cluster){$ID_name="$ID_name|$clusters";}
	if ($ID eq $name){$ID_name="ORIGINIAL_$ID_name";}
	if ($hmmer and $ID ne $motifs){$ID_name="$ID_name|$motifs";}
	my $score;
	if ($type eq "mRNA"){$score=0;}
	elsif ($type eq "tRNA" or $type eq "ncRNA" or $type eq "miRNA" or $type eq "rRNA" or $type eq "snRNA"){$score=30;}
	elsif ($type eq "enhancer" or $type eq "silencer"){$score=40;}
	elsif ($type eq "ori" or $type eq "origin_of_replication"){$score=45;}
	else {$score=50;}
	my $feature = Bio::SeqFeature::Generic->new(-seq_id=>"$ID_name",
						    -start=>$start-$low{$files[$i]},
						    -end=>$end-$low{$files[$i]},
						    -strand=>$dir,
						    -score=>$score
						    );
	$track[$i]->add_feature($feature);
      }
    }
  else
    {
    my $scaf=`cut -f 1 $file$files[$i].gff3| head -n 1`;
    my $short=substr $files[$i],0,4;
    chomp $scaf;
    $scaf="$short\_$scaf";
    $track[$i]=$panel->add_track(-glyph =>'segments',
			    -label => 1,
			    -strand_arrow=>1,
			    -key=>"$scaf"
			    );
    open my $IN, "<$file$files[$i]\_exon.gff3";
    $name=substr $name,5;
    my %exon_hash;
    while (my $line=<$IN>)
      {
      next if $line=~/^\#/;
      chomp $line;
      my ($scaffold,$start,$end,$dir,$ID,$type,$motifs,$clusters)=split ("\t",$line);
      next if ($exon_hash{$ID});
      $exon_hash{$ID}=1;
      my $exons_start=`grep $ID $file$files[$i]\_exon.gff3 | cut -f 2`;
      my @exons_start=split ("\n",$exons_start);
      my $exons_end=`grep $ID $file$files[$i]\_exon.gff3 | cut -f 3`;
      my @exons_end=split ("\n",$exons_end);
      my @exons;
      #my $ID_name="$ID|$type";
      my $ID_name="$ID";
      if ($cluster){$ID_name="$ID_name|$clusters";}
      if ($ID eq $name){$ID_name="ORIGINIAL_$ID_name";}
      if ($hmmer and $ID ne $motifs){$ID_name="$ID_name|$motifs";}
      for (my $j=0;$j<scalar(@exons_end);$j++)
	{
	$exons[$j]= Bio::Graphics::Feature->new(-start=>$exons_start[$j]-$low{$files[$i]},-stop=>$exons_end[$j]-$low{$files[$i]},-type=>'exon');
	}
      my $feature  = Bio::Graphics::Feature->new(-name=>$ID_name,-segments=>[@exons],-type=>'gene');
      $track[$i]->add_feature($feature);
      }
    }
  }

print $panel->png;