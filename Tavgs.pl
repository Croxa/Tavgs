#1 usr/bin/perl
use strict;
use warnings;
use Getopt::Long;	#This module is used for the options in the shell


##########################################################################################################
# Here are the global variables used in the programm
my $gff_file;
my $protein;
my $protein_file;
my $length;
my $number_of_genes='';
my $number_length;
my $file_protein;
my $out_file;
my $protein_data;
my $compare_gff;
my $compare_protein;
my $evalue;
my $delete="";
my @delete;
my $exon="0";
my @blast;
my $help;
my $hmmer="";
my @hmmer;
my $hmmer_prot;
my $pfam;
my $one_file="0";
my $all_files="";
my $name;
my $hmmer_e;
my @blast_file;
my @fas;
my $db_name;
my $cluster="0";
my @cluster;
my $cluster_prot;
my $first_arg=$ARGV[0];
my $only_mRNA;
########################################################################################################
# This part checks for the given options and gives every variable the expected value
GetOptions ('delete' => \$delete,'adapt=s{3}' => \@hmmer, 'exon' => \$exon,'number:i' => \$number_of_genes,'length:i' => \$length,'file=s' => \$protein_file,'protein=s' => \$protein,'help|?' => \$help,'out=s' => \$out_file,'blast_protein=s{5}' => \@blast,'compact' => \$one_file,'blast_file=s{4}' => \@blast_file,'make_cluster_db:s' => \$db_name,'cluster=s{2}' => \@cluster,'genes' =>\$only_mRNA);
if ($help)								#will only print the help menu and then stop everything
  {
  help($help);
  die "\n\n";
  }
if ($db_name){make_cluster_db();die"\n\n";}
if (@cluster)
  {
  $cluster=$cluster[0];
  $cluster_prot=$cluster[1];
  }
if (!@blast_file){$gff_file=$first_arg;if (!$gff_file){die "Missing options or parameters. Use -h option for help!\n";}}		#The first parameter has to be a gff-file unless the option blast_file is active
else
  {
  if (!($first_arg=~/^-/))
    {$gff_file=$first_arg;}
  else{$gff_file=$blast_file[1];}
  }
if (!$number_of_genes and !$length){$length=25000;}										#default length is 25000 bp if nothing else is given
  #some options should not be used together because they won't work or are useless together
if($number_of_genes and $length){die "You cannot use -l AND -n!\n";}
if (@blast and $protein_file){die "You cannot use -f AND -blast_p!\n";}
if (@blast_file and $protein_file){die "You cannot use -f AND -blast_f!\n";}
if ($number_of_genes eq "0" or $number_of_genes){$number_length="-n";if($number_of_genes==0){$number_of_genes=10;}}
elsif ($length){$number_length="-l";if($length==1){$length=25000;}}
if (!$protein and !$protein_file){$protein=$ARGV[1];}
if ($protein){$file_protein='-p';}
elsif ($protein_file){$file_protein='-f';}
if (!$number_of_genes and !$length){$length=25000;}
if (@hmmer) 
  {
  $hmmer="_hmmer";
  $hmmer_prot=$hmmer[0];
  $pfam=$hmmer[1];
  $hmmer_e=$hmmer[2];
  if ($pfam!~/.hmm/){die "Please give a .hmm file for hmmscan!\n"}
  }
if (@blast)
  {
  $protein=$blast[0];
  $protein_data=$blast[1];
  $compare_gff=$blast[2];
  $compare_protein=$blast[3];
  $evalue=$blast[4]
  }
if (@blast_file)
  {
  $protein_data=$blast_file[0];
  $compare_gff=$blast_file[1];
  $compare_protein=$blast_file[2];
  $evalue=$blast_file[3];
  my $fas=`grep \">\" $blast_file[0]`;
  @fas=split("\n",$fas);
  }
if (!$out_file)						#If no path for the files is given, one will be arranged
  {
  if (!$out_file and !$help and !$protein_file and !@blast_file){$out_file="$protein/";}
  elsif (!$out_file and !$help and $protein_file and !@blast_file){$out_file="list_of_proteins/";}
  elsif (!$out_file and !$help and !$protein_file and @blast_file){$out_file="protein_list/";}
  if ($out_file)
    {
    if($gff_file)
      {
      $gff_file=~/([\w.-]*)$/;
      my $short=substr $1,0,4;
      $out_file="$short\_$out_file";
      }
    }
  }
mkdir ("$out_file");


######################################################################################################


if (!-e $gff_file){die "gff_file does not exist!\n"}
if (@blast)
  {
  if (!-e $protein_data){die "$protein_data does not exist!\n";}
  if (!-e $compare_gff){die "$compare_gff does not exist!\n";}
  if (!-e $compare_protein){die "$compare_protein does not exist!\n";}
  blast();
  }
elsif (@fas)
  {
  for my $file_prot(@fas)
    {
    $file_prot=substr $file_prot,1;
    $protein=$file_prot;
    blast();
    }
  }
else
  {
  make_file();
  }
if ($one_file and $all_files)
  {
  make_display($out_file,$all_files);
  }
if ($delete)
  {unlink @delete;}





##########################################################################################################
sub make_file			# The main sub for making the files which willbe used by Tavgs_display.pl. Different options will enable different files here.
  {
  my @proteins;
  if ($file_protein eq "-p")
    {
    push @proteins, $protein;
    }
  elsif($file_protein eq "-f")
    {
    open my $IN, "<$protein_file" or die "$protein_file does not exist\n";
    while (my $line=<$IN>)
      {
      chomp $line;
      push @proteins, $line;
      }
    close $IN;
    }
  else {die "-f or -p required\n";}
  $gff_file=~/([\w.-]*)$/;
  my $data=$1;
  $name=substr $data,0,4;
  for my $list (@proteins)
    {
    my $constant_name=$list;
    my $gff_line=`grep '$list\[;, ]' $gff_file | awk '{if (\$3==\"mRNA\" || \$3==\"tRNA\" || \$3==\"ncRNA\" || \$3==\"miRNA\" || \$3==\"rRNA\" || \$3==\"snRNA\" || \$3==\"enhancer\" || \$3==\"silencer\" || \$3==\"ori\" || \$3==\"origin_of_replication\")print \$0}'`;
    if (!$gff_line){$gff_line=`grep '$list\$' $gff_file | awk '{if (\$3==\"mRNA\" || \$3==\"tRNA\" || \$3==\"ncRNA\" || \$3==\"miRNA\" || \$3==\"rRNA\" || \$3==\"snRNA\" || \$3==\"enhancer\" || \$3==\"silencer\" || \$3==\"ori\" || \$3==\"origin_of_replication\")print \$0}'`;}	# The given element can be of of these features: mRNA, tRNA, ncRNA, miRNA, rRNA, snRNA, enhancer, silencer, ori, origin_of_replication, has to be in the thirf column of the gff-file

    if (!$gff_line)
      {
      $list=~/([\dA-Za-z.]*)$/;
      $list=$1;
      if ($list)
	{
	$gff_line=`grep '$list\[;, ]' $gff_file | awk '{if (\$3==\"mRNA\")print \$0}'`;
	if (!$gff_line){`grep '$list\$' $gff_file | awk '{if (\$3==\"mRNA\")print \$0}'`;}
	}
      }
    if (!$gff_line){print "$list does not exist!";}
    next if (!$gff_line);
    chomp $gff_line;
    my @gff_line=split("\t",$gff_line);
    my $start=$gff_line[3];
    my $end=$gff_line[4];
    my $scaffold=$gff_line[0];
    my $nearby_genes=`grep $scaffold $gff_file | awk '{if (\$3==\"mRNA\" || \$3==\"tRNA\" || \$3==\"ncRNA\" || \$3==\"miRNA\" || \$3==\"rRNA\" || \$3==\"snRNA\" || \$3==\"enhancer\" || \$3==\"silencer\" || \$3==\"ori\" || \$3==\"origin_of_replication\")print \$0}'`;      # this line greps every element on the scaffold, the later options decide which will be printed into the file
    my @nearby_genes=split("\n",$nearby_genes);
    $list="$name\_$list";
    open my $OUT, ">$out_file$list.gff3" or die "$out_file does not exist!\n";
    push (@delete,"$out_file$list.gff3");
    if ($number_length eq "-l")		#This is by default enabled and only disabled when the option number_of_genes is enabled. This prints out every element which is in the given range (default 25000)
      {
      for (@nearby_genes)
	{
	my @column=split("\t",$_);
	if ($column[3]>$start-$length and $column[4]<$end+$length)
	  {
	  $column[8]=~/ID=([\w.]*)[;,]/;
	  my $ID=$1;
	  if (!$only_mRNA){print $OUT join("\t",$column[0],$column[3],$column[4],$column[6],$ID,$column[2]),"\n";}
	  elsif($only_mRNA and $column[2] eq "mRNA"){print $OUT join("\t",$column[0],$column[3],$column[4],$column[6],$ID,$column[2]),"\n";}
	  }
	}
      }
    elsif ($number_length eq "-n") 			#This option will count to the X'th gene upstream and downstream and will print every element in between
      {
      my @all_genes;
      $all_genes[0]=$gff_line;
      my %start_all_genes;
      $start_all_genes{$gff_line}=$start;  
      my $ori_gene;
      my $start_num;
      my $end_num;
      my @frontgenes;
      my @backgenes;
      for (@nearby_genes)
	{
	next if ($_ eq $gff_line);
	my @column=split("\t",$_);
	for (my $i=0;$i<scalar(@all_genes);$i++)
	  {
	  if ($start_all_genes{$all_genes[$i]}>$column[3])
	    {
	    splice @all_genes,$i,0,$_;
	    $start_all_genes{$_}=$column[3];
	    $i=scalar(@all_genes);
	    }
	  elsif ($start_all_genes{$all_genes[scalar(@all_genes)-$i-1]}<$column[3])
	    {

	    if (scalar(@all_genes)-$i+1>scalar(@all_genes))
	      {
	      push @all_genes,$_;
	      }
	    else 
	      {
	      splice @all_genes,scalar(@all_genes)-$i,0,$_;
	      }
	    $start_all_genes{$_}=$column[3];
	    $i=scalar(@all_genes);
	    }	  
	  }
	} 
      for (my $i=0;$i<scalar(@all_genes);$i++)
	{
	my @columns=split("\t",$all_genes[$i]);
	if ($columns[2] eq "mRNA")
	  {
	  if (!$ori_gene){ push (@frontgenes,$i);}
	  else {push (@backgenes,$i);}
	  }
	if ($all_genes[$i] eq $gff_line)
	  {
	  $ori_gene=$i;
	  }
	}
      
      if (scalar(@frontgenes)-$number_of_genes<1)
	{$start_num=0;}
      else{$start_num=$frontgenes[scalar(@frontgenes)-$number_of_genes-1];}
      if ($number_of_genes>=scalar(@backgenes))
	{$end_num=scalar(@all_genes)-1;}
      else{$end_num=$backgenes[$number_of_genes-1];}
      for (my $i=$start_num;$i<=$end_num;$i++)
	{
	my @column=split("\t",$all_genes[$i]);
	$column[8]=~/ID=([\w.]*)[;,]/;
	my $ID=$1;
	if (!$only_mRNA){print $OUT join("\t",$column[0],$column[3],$column[4],$column[6],$ID,$column[2]),"\n";}
	elsif($only_mRNA and $column[2] eq "mRNA"){print $OUT join("\t",$column[0],$column[3],$column[4],$column[6],$ID,$column[2]),"\n";}
	}
      }
    close $OUT;
    if ($hmmer)			#If hmmer is active this will compare every protein to the hmm databank and will print the found domains to the file. Uses hmmscan in shell
      {
      open my $IN, "<$out_file$list.gff3" or die "Couldn't open $out_file$list.gff3\n";
      unlink "$out_file$list$hmmer.gff3";
      open $OUT, ">>$out_file$list$hmmer.gff3" or die "Couldn't open $out_file$list$hmmer.gff3\n" ;
      while (my $line=<$IN>)
	{
	chomp $line;
	my @column=split("\t",$line);
	my $protein_atm=$column[4];
	my $variable=2;
	my $sequence="";
	open my $IN_prot,"<$hmmer_prot";
	while (my $lines=<$IN_prot>)
	  {
	  if ($lines=~/($protein_atm)$/)
	    {
	    $variable=0;
	    }
	  if ($lines=~/^>/)
	    {
	    $variable++;
	    }
	  next if ($variable>1);
	  next if ($lines=~/^>/);
	  chomp $lines;
	  $sequence=join("",$sequence,$lines);
	  }
	close $IN_prot;
	if ($sequence)
	  {
	  open my $OUT_prot, ">$out_file$name$protein_atm.fas" or die "Couldn't open $out_file$name$protein_atm.fas\n";
	  print $OUT_prot ">$protein_atm\n$sequence\n";
	  #if (!$sequence){die "$protein_atm could not be found in the $hmmer_prot (hmmscan)\n";}
	  close $OUT_prot;
	  system "hmmscan --tblout $out_file$name$hmmer\_pfam_result_$protein_atm.txt $pfam $out_file$name$protein_atm.fas > out.log";
	  unlink "out.log";
	  push (@delete,"$out_file$name$hmmer\_pfam_result_$protein_atm.txt","$out_file$name$protein_atm.fas");
	  open my $IN_pfam, "<$out_file$name$hmmer\_pfam_result_$protein_atm.txt";
	  $variable=0;
	  my @motifs;
	  my $counter=0;
	  while (my $lines=<$IN_pfam>)
	    {
	    next if ($lines=~/^#/);
	    my @columns=split /\s+/,$lines;
	    if ($columns[4]<$hmmer_e){push (@motifs,$columns[0]);}
	    #push (@motifs,$columns[0]);
	    }
	  if (@motifs){print $OUT $line,"\t",join("|",@motifs),"\n";}
	  else {print $OUT $line,"\t",$column[4],"\n";}
	  }
	else 
	  {
	  print $OUT $line,"\t",$column[4],"\n";
	  }
	}
      close $IN;
      close $OUT;
      rename "$out_file$list$hmmer.gff3","$out_file$list.gff3";  
      }
      
    if ($cluster)		#Uses a cluster databank to give the elements new names. Uses blast to find the best hit and therefor the best cluster
      {
      my %cluster_hash;
      open my $IN, "<$cluster";
      my $blast_comp;
      while (my $line=<$IN>)
	{
	chomp $line;
	if ($line=~/^#/){$blast_comp=substr $line,1;chomp $blast_comp}
	next if ($line=~/^#/);
	my @column=split("\t",$line);
	$cluster_hash{$column[0]}=$column[1];
	}
      close $IN;
      open $IN, "<$out_file$list.gff3";
      open my $OUT, ">$out_file$list\_cluster.gff3";
      push(@delete,"$out_file$list\_cluster.gff3");

      while (my $line=<$IN>)
	{
	chomp $line;
	my @column=split("\t",$line);
	my $prot_atm=$column[4];
	if ($column[5] eq "mRNA")
	  {
	  open my $IN_prot,"<$cluster_prot";
	  my $variable=2;
	  my $sequence="";
	  while (my $lines=<$IN_prot>)
	    {
	    chomp $lines;
	    if ($lines=~/($prot_atm)$/)
	      {
	      $variable=0;
	      }
	    if ($lines=~/^>/)
	      {
	      $variable++;
	      }
	    next if ($variable>1);
	    next if ($lines=~/^>/);
	    chomp $lines;
	    $sequence=join("",$sequence,$lines);
	    }
	  open my $OUT_prot, ">$out_file$prot_atm.fas";
	  print $OUT_prot ">$prot_atm\n$sequence\n";
	  close $OUT_prot;
	  push (@delete, "$out_file$prot_atm.fas");
	  close $IN_prot;
	  }
	my $cluster_atm=$column[4];
	if ($column[5] eq "mRNA")
	  {
	  system "blastp -query $out_file$prot_atm.fas -db cluster_db/blast_db/$blast_comp\_db -outfmt 6 -out $out_file$prot_atm\_clusterblast.txt -evalue 1e-10";
	  push (@delete,"$out_file$prot_atm\_clusterblast.txt");
	  $cluster_atm=`head -n 1 $out_file$prot_atm\_clusterblast.txt | cut -f 2`;
	  chomp $cluster_atm;
	  }
	if (!$cluster_hash{$cluster_atm}){$cluster_atm="No_cluster";}
	else {$cluster_atm=$cluster_hash{$cluster_atm};}
	if ($hmmer){print $OUT "$line\t$cluster_atm\n";}
	else {print $OUT "$line\t$column[4]\t$cluster_atm\n";}
	}
      close $IN;
      close $OUT;
      rename ("$out_file$list\_cluster.gff3","$out_file$list.gff3");
      }
    if ($exon)			#Will print out the exons instead of just the length of the element, if no exons are found, it will just use the regular element
      {
      open my $IN_exon, "<$out_file$list.gff3";
      open my $OUT_exon, ">$out_file$list\_exon.gff3";
      push (@delete,"$out_file$list\_exon.gff3");
      while (my $line=<$IN_exon>)
	{
	chomp $line;
	my @column=split("\t",$line);
	my $exons=`grep '$column[4]\[\\;\\,\ \]' \"$gff_file\" | awk '{if (\$3==\"exon\")print \$0}'`;
	chomp $exons;
	if (!$exons){$exons=`grep $column[4]\$ \"$gff_file\" | awk '{if (\$3==\"exon\")print \$0}'`;}
	if (!$exons){print $OUT_exon "$line\n";}
	my @exons=split("\n",$exons);
	if (@exons)
	  {
	  for (@exons)
	    {
	    my @exon_column=split("\t",$_);
	    my $print;
	    if (!$column[6])
	      {
	      $print=join ("\t",$column[0],$exon_column[3],$exon_column[4],$column[3],$column[4],$column[5],$column[4]);
	      }
	    else
	      {
	      $print=join ("\t",$column[0],$exon_column[3],$exon_column[4],$column[3],$column[4],$column[5],$column[6]);
	      }
	    if ($cluster)
	      {$print="$print\t$column[7]";}
	    print $OUT_exon $print,"\n";
	    }
	  }
	}
      close $IN_exon;
      close $OUT_exon;
      }
    if (!$one_file){make_display($out_file,$list)}	#if compact is not active this will produce a png-file
    else 						#if compact is active this will store the filenames and gives them in one big step to Tavgs_display.pl
      {
      if ($all_files){$all_files="$all_files,$list";}
      else {$all_files=$list;}
      }
    }
  }

  
#######################################################################################################################################


sub make_display
  {
  my $directory=shift;
  my $file=shift;
  if (!$hmmer){$hmmer="0";}
  my $time=get_time();
  system "perl Tavgs_display.pl $directory \"$file\" $exon $hmmer $one_file $cluster > $directory$time\_$name\_display.png";	#Tavgs_display.pl will produce the png-file with all the given parameters
  }

  
#######################################################################################################################################


sub blast			#blasts a protein and will produce a png-file for every found protein
  {
  open my $IN, "<$protein_data";
  my $variable=2;
  my $sequence="";
  while (my $line=<$IN>)
    {
    if ($line=~/($protein)$/)
      {
      $variable=0;
      }
    if ($line=~/^>/)
      {
      $variable++;
      }
    next if ($variable>1);
    next if ($line=~/^>/);
    chomp $line;
    $sequence=join("",$sequence,$line);
    }
  close $IN;
  open my $OUT, ">$out_file$protein.fas";
  print $OUT ">$protein\n$sequence\n";
  close $OUT;
  if (!$sequence){die "$protein could not be found in the $protein_data\n";}
  $file_protein="-p";
  if (@hmmer)
    {
    if(@blast){$hmmer_prot=$blast[1];}
    elsif(@blast_file){$hmmer_prot=$blast_file[0];}
    }
  if (@cluster)
    {
    if(@blast){$cluster_prot=$blast[1];}
    elsif(@blast_file){$cluster_prot=$blast_file[0];}
    }
  if (!@blast_file or !($first_arg=~/^-/)){make_file();}
  if ($hmmer){$hmmer_prot=$hmmer[0];}
  if ($cluster){$cluster_prot=$cluster[1];}
  mkdir "db/";
  $compare_protein=~/([\w.-]*)$/;
  my $db_name=$1;
  system "makeblastdb -in $compare_protein -out db/$db_name\_db -dbtype prot" unless -e "db/$db_name\_db.phr";
  system "blastp -query $out_file$protein.fas -out $out_file$protein\_blast.txt -db db/$db_name\_db -outfmt 6 -evalue $evalue";
  push (@delete,"$out_file$protein\_blast.txt","$out_file$protein.fas");
  open $IN, "<$out_file$protein\_blast.txt"; 
  $gff_file=$compare_gff;
  while (my $line=<$IN>)
    {
    my @column=split("\t",$line);
    $protein=$column[1];
    make_file();
    }
  close $IN;
  }

  
##############################################################################################################
  
  
sub get_time		#Will get a timestamp to use for the file_names, to not overwrite the files 
  {
  my @time=localtime(time);
  my $timestamp=$time[5]+1900;
  $time[4]=$time[4];
  $timestamp="$timestamp\-$time[4]\-$time[3]_$time[2]:$time[1]:$time[0]";
  return $timestamp;
  }
  
  
################################################################################################################


sub make_cluster_db		#Will make a file with clusters, the fasta file gets blasted against itself and every hit is put into the same cluster
  {
  mkdir "cluster_db";
  mkdir "cluster_db/blast_db";
  $db_name=~/([\w.-]*)$/;
  my $out_name=$1;
  system "makeblastdb -in $db_name -out cluster_db/blast_db/$out_name\_db -dbtype prot" unless -e "cluster_db/blast_db/$out_name\_db.phr";
  #system "blastp -query $db_name -out cluster_db/out_$out_name -db cluster_db/blast_db/$out_name\_db -outfmt 6 -evalue 1e-40 &";
  my $while=1;
  my $time=0;
  my $one=`grep \">\" $db_name | wc -l`;chomp $one;
  my $prozent=0;
  my $counter=0;
  sleep 5;
  open my $fh, "<cluster_db/out_$out_name";
  while($while)
    {
    if((stat($fh))[9]+120<time){$while=0;}
    my $mom=`cut -f 1 cluster_db/out_$out_name | sort | uniq | wc -l`;chomp $mom;
    $prozent=(int(($mom/$one)*100))/100;
    if ($prozent>=$counter*0.01)
      {print $counter,"\%\t";
      $counter++;}
    }
  close $fh;
  $out_name="out_$out_name";
  open my $IN, "<cluster_db/$out_name";
  open my $OUT, ">cluster_db/$out_name\_cluster_db.txt";
  print $OUT "\#",substr $out_name,4;
  my $cluster_nr=0;
  while (my $line=<$IN>)
    {
    chomp $line;
    my @column=split("\t",$line);
    my $control=`grep -P '$column[0]\\t' cluster_db/$out_name\_cluster_db.txt`;
    if (!$control)
      {
      $cluster_nr++;
      my $column1=`awk '{if (\$1=="$column[0]") print \$2}' cluster_db/$out_name`;
      my @column1=split ("\n",$column1);
      for (@column1)
	{	
	my $column2=`awk '{if (\$1=="$_") print \$2}' cluster_db/$out_name`;
	my @column2=split("\n",$column2);
	for (@column2)
	  {
	  my $control2=`grep -P '$_\\t' cluster_db/$out_name\_cluster_db.txt`;
	  if (!$control2)
	    {
	    print $OUT "\n$_\tcluster_$cluster_nr"
	    }
	  }    
	}
      }
    }
  $out_name=substr $out_name,4;
  unlink "cluster_db/out_$out_name";
  }
  
  
################################################################################################################


sub help		# prints out the help menu
  {
    print "This program uses the Bio::Graphics perl module for visualisation of genome features. Make sure you have it installed (guide for installing: http://bioperl.org/howtos/BioGraphics_HOWTO.html or https://nsaunders.wordpress.com/2006/06/25/an-introduction-to-biographics/). Furthermore this program uses \"gff3\" and \"fasta/fas\" files. Please make sure, that your the last column in your gff files have \"ID=...\" and that your fasta files are all the same with your mRNA-ID at the end of the description line. If you want to use the \"-p\" option, make sure that you have command line blastp installed. If you want to use the \"-a\" option, make sure that you have command line hmmscan installed (hmmer packet). Genes are displayed as green, RNAs except mRNAs are yellow, enhancer/silencer are yellow/orange and origins of replications are red. \n\nDefault:\ngff3-file protein-ID 			-	Produces a png file of the protein and every other protein within 25000 bp upstream and downstream, using the gff3-file\n\n";
    print "Options: gff_file -option parameter(s)\n\n";
    print "-o directory/ 				-	--out: The files will be saved in this directory, use\"/\" at the end, the files will be named automatically. Will make a directory with the protein-ID as default. Works with every other option.\n\n";
    print "-d 					-	--delete: If this is active, every file that is created except the png file will be deleted. Works with every other option.\n\n";
    print "-e 					-	--exon: If this is active, the png file will show the gene seperated by exons but it will not show different colors for different features. Works with every other option.\n\n";
    print "-a protein-file PFAM-files e-value	-	--adapt: If this is active, the png file will show the gene names by motifs found in the protein. Needs \"hmmscan\" to work. Works with every other option.\n\n";
    print "-g 					-	--genes: If this is active, only the regions which code for mRNAs are shown in the png file. Works with every other option\n\n";
    print "-p protein-ID				-	--protein: The programm will use this protein to create to png file. This is used in the default settings. Does not work with -f, -blast_f or -blast_p.\n\n";
    print "-f protein_file				-	--file: Reads a file which has a protein-ID in each line and produces a png file for each protein. Does not work with -p, -blast_f or -blast_p.\n\n";
    print "-n number_of_genes			-	--number of genes: The png file will show the next X number of genes upstream and downstream of the protein on the scaffold. If no parameter is given, X=10. Does not work with -l.\n\n";
    print "-l length				-	--length: The png file show the genes within X bp upstream and downstream of the protein on the scaffold. This is used in the default settings. If no parameter is given, length=25000. Does not work with -n.\n\n";
    print "-blast_p protein-ID protein-file comparing-gff3 comparing-protein-file evalue	-	--blast_protein: This option will blast the protein from a protein file against a \"comparing-protein-file\" (fasta-format) and will create png files for every hit. The higher the evalue, the less and better are the hits. Does not work with -p, -f or -blast_f. \n\n";
    print "-blast_f protein_file comparing-gff3 comparing-protein-file evalue		-	--blast_file: This option will blast all proteins from a protein file against a \"comparing-protein-file\" (fasta-format) and will create png files for every hit. The higher the evalue, the less and better are the hits. Does not work with -p or -f, -blast_p. \n\n";
    print "-co 					-	--compact: If this is active, multiple png files will be transformed into one png file with every region. Only has an effect with -b or -f.\n\n";
    print "-h/? 					-	--help:shows the manual with all options \n\n";
    print "-m fasta-file				-	--make_cluster_db: Make a databank for the -cl option. The fasta file will blast with itself to set comparing clusters. This can take a while.\n\n";
    print "-cl cluster_db fasta_file		-	--cluster: Replace the names of the features with clusters out of a cluster_db. The fasta file has to be with the comparing proteins.\n\n"
  }