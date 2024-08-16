#!/usr/bin/perl

use strict;
use Getopt::Long;
use Cwd 'abs_path';

my $help;
my $version;
my $primer;
my $inbam;
my $outbam;
my $type;
my $split;
my $thread = 10;
my $max_mismatch = 100;

GetOptions(
  'help' => \$help,
  'version' => \$version,
  'primer=s' => \$primer,
  'in=s' => \$inbam,
  'out=s' => \$outbam,
  'type=s' => \$type,
  'split=i' => \$split,
  'thread=i' => \$thread,
  'max_mis=i' => \$max_mismatch
);

if ($help) {
  print "
Basic usage:

  AtomSeqTrimmer -in <in.bam> -out <out.bam> -primer <primer position BED> -type <pe|se>

Required:
  -in\t\tInput BAM/SAM file after aligning.
  -out\t\tOutput coordinate sorted BAM file after clipping primers.
  -primer\tBED file containing position of primers.
  -type\t\tPair-end or single-end read.[pe|se].

Optional:
  -help\t\tPrint this message.
  -version\tPrint the version message.
  -thread\tNumber of threads to launch. Default = 10
  -split\tWhether split the input BAM into two BAM (plus and minus strand). 1/0 means yes/no. Default=0
  -max_mis\tMax mismatch in a read. Default=100

  The BED file of primers should contain at least first six colum:
  1. Chormosome ID  (consistent with reference genome);
  2. Start of primer (0-based);
  3. End of primer (1-based);
  4. Primer ID;
  5. Anything;
  6. Targeted strand (+/-).

";
  exit;
}

if ($version) {
  print "AtomPrimerTrimmer v1.0.2\n";
  exit;
}

my %para = (0=>'-in', 1=>'-out', 2=>'-primer', 3=>'-type');
my @para = ($inbam, $outbam, $primer, $type);
my @empty_para = ();
for (my $i = 0; $i < @para; $i++) {
  unless ($para[$i]) {
    push @empty_para, $para{$i};
  }
}
if (@empty_para && !$help && !$version) {
  die "Error: Can not get the parameter: @empty_para\nUse \"-h\" to see usage document\n";
}

my %type_option = ('pe' => 1, 'se' => 1);
unless ($type_option{$type}) {
  die "Error: \'$type\' is not in '-type' option list\n\n";
}

my %split_option = ('1' => 1, '0' => 1);
if ($split && !$split_option{$split}) {
  die "Error: \'$split\' is not in '-split' option list\n\n";
}

if ($type eq 'se' && $split eq '1') {
  warn "Warining: '-split' must be '0' when '-type' is 'se'\n\n";
}

$inbam = abs_path($inbam);
$outbam = abs_path($outbam);
$primer = abs_path($primer);
(my $outdir = $outbam) =~ s/(.+)\/.+/$1/;

my %chr2pos = ();
open I,"$primer";
while(<I>){
  chomp;
  my @section = split /\t/,$_;
  $section[1] = $section[1] + 1;
  if ($chr2pos{$section[0]}{$section[5]}) {
    $chr2pos{$section[0]}{$section[5]} = "$chr2pos{$section[0]}{$section[5]}\,$section[1]\-$section[2]";
  } else {
    $chr2pos{$section[0]}{$section[5]} = "$section[1]\-$section[2]";
  }
}
close I;

if ($type eq 'pe') { &pair_end; }
else { &single_end; }



sub single_end {

system "samtools view -H $inbam > $outdir/tmp.sam";
open O,">>$outdir/tmp.sam";

open I,"samtools view -F 256 $inbam |";
LINE: while(<I>) {
  next LINE if $_ =~ /^@/;
  chomp;
  my $rawline = $_;
  my @section = split /\t/, $rawline;
  my $cigar = $section[5];
  my $flag = $section[1];
  my $chr = $section[2];
  my $read_start = $section[3];
  my $read_end = $read_start - 1 + &compute_ref_length($cigar); 
  (my $strand, my $readtype) = (&find_strand_read($flag))[0,1];
  if ($readtype==1 || $readtype==2) {
    die "Error: Input contains pair-end read\n\n";
  }

  my $read_seq = $section[9];
  my $read_qua = $section[10];
  my $head_S_count = 0;
  if ($strand eq "+" && $cigar =~ /^(\d+)S/) { $head_S_count = $1;}
  elsif ($strand eq "-" && $cigar =~ /(\d+)S$/) { $head_S_count = $1;}
  next LINE if $head_S_count > 2;
    
  my $target_start = 0;
  my $target_end = 0;
  my @pos = split /\,/, $chr2pos{$chr}{$strand};
  PRIMER_CYCLE: foreach my $pos (@pos) {
    (my $primer_start, my $primer_end) = split /\-/, $pos;
    (my $diff_head, my $diff_tail) = &compare_read_and_primer($strand,$read_start,$read_end,$primer_start,$primer_end);
    if ($diff_head <= 4 && $diff_head >= -1 && $diff_tail > 0) {
      $target_start = $primer_start;
      $target_end = $primer_end;
      last PRIMER_CYCLE;
    }
  }

  if ($target_start == 0 || $target_end == 0) {
    print O "$rawline\n";
    next LINE;
  }

  my $MD;
  for (my $i=11; $i<@section; $i++) {
    if ($section[$i] =~ /MD:Z:/) { $MD = $section[$i]; }
  }
  my $mismatch_count = &compute_mutation_count($cigar, $MD);
  next LINE if $mismatch_count > $max_mismatch;

  my $clip_length = &compute_clip_length(1, $strand, $read_start, $read_end, $target_start, $target_end);
  my $newcigar = &clip_cigar($cigar,1,$strand,$clip_length);
  next LINE if $newcigar eq 'NA';

  my $newread_length = &compute_read_length($newcigar);
  (my $newread_seq, my $newread_qua) = &clip_seq_qua($read_seq,$read_qua,$newread_length,1,$strand);
  my $index_md = -1;
  my $index_nm = -1;
  for (my $i=11; $i<@section; $i++) {
    if ($section[$i] =~ /MD:Z:/) { $index_md = $i; }
    elsif ($section[$i] =~ /NM:i:/) { $index_nm = $i; }
  }
  if ($index_md > 10) { $section[$index_md] = &clip_MD($section[$index_md],1,$strand,$clip_length); }
  if ($index_md > 10 && $index_nm > 10) {
    $section[$index_nm] = &compute_edit_distance($newcigar, $section[$index_md]);
    $section[$index_nm] = "NM:i:" . $section[$index_nm];
  }

  if ($strand eq "+") {
    $section[3] = &edit_pos($section[3],$clip_length);
  }

  $section[5] = $newcigar;
  $section[9] = $newread_seq;
  $section[10] = $newread_qua;
  my $newline = join "\t", @section;
  print O "$newline\n";
}
close I;
close O;
system "samtools view -@ $thread -bS $outdir/tmp.sam | samtools sort -O bam -@ $thread > $outbam";
system "rm $outdir/tmp.sam";

}



sub pair_end {

my %readid2count = ();
open I,"samtools view $inbam | ";
while(<I>){
  next LINE if $_ =~ /^@/;
  my @section = split /\t/,$_;
  my $readtype = &find_strand_read($section[1]);
  if ($readtype == 0) {
    die "Error: Input contains single-end read\n\n";
  }
  $readid2count{$section[0]}++;
}
close I;

open I,"samtools sort -n -@ $thread $inbam | samtools view | ";
open O,">$outdir/tmp";
while(<I>){
  my @section = split /\t/,$_;
  if ($readid2count{$section[0]} == 2) {
    print O "$_";
  }
}
close I;
close O;

if ($split) {
  system "samtools view -H $inbam > $outdir/tmp_plus.sam";
  system "samtools view -H $inbam > $outdir/tmp_minus.sam";
  open O1,">>$outdir/tmp_plus.sam";
  open O2,">>$outdir/tmp_minus.sam";
} else {
  system "samtools view -H $inbam > $outdir/tmp.sam";
  open O,">>$outdir/tmp.sam";
}

open I,"$outdir/tmp";
LINE: while(my $line1=<I>){
  my $line2=<I>;
  my @section1 = split /\t/, $line1;
  my @section2 = split /\t/, $line2;
  (my $line1_strand, my $line1_readtype) = &find_strand_read($section1[1]);
  (my $line2_strand, my $line2_readtype) = &find_strand_read($section2[1]);
  next LINE if $line1_strand eq $line2_strand || $line1_readtype == $line2_readtype;
  
  my @first;
  my @second;
  if ($line1_readtype == 1) {
    @first = @section1;
    @second = @section2;  
  }
  else {
    @first = @section2;
    @second = @section1;
  }

  my $cigar1 = $first[5];
  my $flag1 = $first[1];
  my $chr1 = $first[2];
  my $read1_start = $first[3];
  my $read1_end = $read1_start - 1 + &compute_ref_length($cigar1);
  my $strand1 = (&find_strand_read($flag1))[0];
  my $read1_seq = $first[9];
  my $read1_qua = $first[10];
  my $cigar2 = $second[5];
  my $flag2 = $second[1];
  my $read2_start = $second[3];
  my $read2_end = $read2_start - 1 + &compute_ref_length($cigar2);
  my $strand2;
  if ($strand1 eq "+") { $strand2 = "-"; }
  else { $strand2 = "+"; }
  my $read2_seq = $second[9];
  my $read2_qua = $second[10];

  my $head_S_count = 0;
  if ($strand1 eq "+" && $cigar1 =~ /^(\d+)S/) { $head_S_count = $1;}
  elsif ($strand1 eq "-" && $cigar1 =~ /(\d+)S$/) { $head_S_count = $1;}
  next LINE if $head_S_count > 2;

  my $target_start = 0;
  my $target_end = 0;
  my @pos = split /\,/, $chr2pos{$chr1}{$strand1};
  PRIMER_CYCLE: foreach my $pos (@pos) {
    (my $primer_start, my $primer_end) = split /\-/, $pos;
    (my $diff_head, my $diff_tail) = &compare_read_and_primer($strand1,$read1_start,$read1_end,$primer_start,$primer_end);
    if ($diff_head <= 4 && $diff_head >= -1 && $diff_tail > 0) {
      $target_start = $primer_start;
      $target_end = $primer_end;
      last PRIMER_CYCLE;
    }
  }
  next LINE if $target_start==0 || $target_end==0;

  my $MD1; my $MD2;
  for (my $i=11; $i<@first; $i++) {
    if ($first[$i] =~ /MD:Z:/) { $MD1 = $first[$i]; }
  }
  for (my $i=11; $i<@second; $i++) {
    if ($second[$i] =~ /MD:Z:/) { $MD2 = $second[$i]; }
  }
  my $mismatch_count1 = &compute_mutation_count($cigar1, $MD1);
  my $mismatch_count2 = &compute_mutation_count($cigar2, $MD2);
  next LINE if $mismatch_count1 > $max_mismatch || $mismatch_count2 > $max_mismatch;

  my $clip_length1 = &compute_clip_length(1, $strand1, $read1_start, $read1_end, $target_start, $target_end);
  my $clip_length2 = &compute_clip_length(2, $strand2, $read2_start, $read2_end, $target_start, $target_end);
  my $newcigar1 = &clip_cigar($cigar1,1,$strand1,$clip_length1);
  my $newcigar2 = &clip_cigar($cigar2,2,$strand2,$clip_length2);
  next LINE if $newcigar1 eq 'NA' || $newcigar2 eq 'NA';

  my $newread1_length = &compute_read_length($newcigar1);
  my $newread2_length = &compute_read_length($newcigar2);
  (my $newread1_seq, my $newread1_qua) = &clip_seq_qua($read1_seq,$read1_qua,$newread1_length,1,$strand1);
  (my $newread2_seq, my $newread2_qua) = &clip_seq_qua($read2_seq,$read2_qua,$newread2_length,2,$strand2);

  my $index1_mc = -1;
  my $index1_md = -1;
  my $index1_nm = -1;
  my $index2_mc = -1;
  my $index2_md = -1;
  my $index2_nm = -1;
  for (my $i=11; $i<@first; $i++) {
    if ($first[$i] =~ /MC:Z:/) { $index1_mc = $i; }
    elsif ($first[$i] =~ /MD:Z:/) { $index1_md = $i; }
    elsif ($first[$i] =~ /NM:i:/) { $index1_nm = $i; }
  }
  for (my $i=11; $i<@second; $i++) {
    if ($second[$i] =~ /MC:Z:/) { $index2_mc = $i; }
    elsif ($second[$i] =~ /MD:Z:/) { $index2_md = $i; }
    elsif ($second[$i] =~ /NM:i:/) { $index2_nm = $i; }
  }
  if ($index1_mc > 10) { $first[$index1_mc] = "MC:Z:$newcigar2"; }
  if ($index2_mc > 10) { $second[$index2_mc] = "MC:Z:$newcigar1"; }
  if ($index1_md > 10) { $first[$index1_md] = &clip_MD($first[$index1_md],1,$strand1,$clip_length1); }
  if ($index2_md > 10) { $second[$index2_md] = &clip_MD($second[$index2_md],2,$strand2,$clip_length2); }
  if ($index1_nm > 10 && $index1_md > 10) {
    $first[$index1_nm] = &compute_edit_distance($newcigar1, $first[$index1_md]);
    $first[$index1_nm] = "NM:i:$first[$index1_nm]";
  }
  if ($index2_nm > 10 && $index2_md > 10) {
    $second[$index2_nm] = &compute_edit_distance($newcigar2, $second[$index2_md]);
    $second[$index2_nm] = "NM:i:$second[$index2_nm]";
  }

  if ($strand1 eq "+") {
    $first[3] = &edit_pos($first[3],$clip_length1);
    $second[7] = $first[3];
    $second[3] = &edit_pos($second[3],$clip_length2);
    $first[7] = $second[3];
  }

  my $newstart1 = $first[3];
  my $newend1 = $newstart1 - 1 + &compute_ref_length($newcigar1);
  my $newstart2= $second[3];
  my $newend2 = $newstart2 - 1 + &compute_ref_length($newcigar2);
  my @newrange = ($newstart1, $newend1, $newstart2, $newend2);
  @newrange = sort {$a <=> $b} @newrange;
  my $fragment_length = $newrange[-1] - $newrange[0] + 1;
  if ($strand1 eq "+") {
    $first[8] = $fragment_length;
    $second[8] = $fragment_length * (-1);
  } else {
    $first[8] = $fragment_length * (-1);
    $second[8] = $fragment_length;
  }

  $first[5] = $newcigar1;
  $first[9] = $newread1_seq;
  $first[10] = $newread1_qua;
  $second[5] = $newcigar2;
  $second[9] = $newread2_seq;
  $second[10] = $newread2_qua;
  my $newline1 = join "\t", @first;
  my $newline2 = join "\t", @second;
  if ($split) {
    if ($strand1 eq '+') {
      print O1 "$newline1";
      print O1 "$newline2";
    } else {
      print O2 "$newline1";
      print O2 "$newline2";
    }
  } else {
    print O "$newline1";
    print O "$newline2";
  }
}
close I;

if ($split) {
  close O1;
  close O2;
  (my $outbam1 = $outbam) =~ s/\.bam/_plus\.bam/;
  (my $outbam2 = $outbam) =~ s/\.bam/_minus\.bam/;
  system "samtools view -@ $thread -bS $outdir/tmp_plus.sam | samtools sort -O bam -@ $thread > $outbam1";
  system "samtools view -@ $thread -bS $outdir/tmp_minus.sam | samtools sort -O bam -@ $thread > $outbam2";
  system "rm $outdir/tmp $outdir/tmp_plus.sam $outdir/tmp_minus.sam";
} else {
  close O;

  system "samtools view -@ $thread -bS $outdir/tmp.sam | samtools sort -O bam -@ $thread > $outbam";
  system "rm $outdir/tmp $outdir/tmp.sam";
}

}



sub find_strand_read {
  my $flag = shift @_;
  my @coefficient = ();
  while ($flag > 0) {
    my $remain = $flag % 2;
    push @coefficient, $remain;
    $flag = ($flag - $remain) / 2;
  }
  my $strand;
  my $read = 0;
  if ($coefficient[4] == 1) { $strand = '-'; }
  else { $strand = '+'; }
  if ($coefficient[6] == 1) { $read = 1; }
  elsif ($coefficient[7] == 1) { $read = 2; }
  return ($strand, $read);
}



sub compare_read_and_primer {
  (my $strand, my $read_start, my $read_end, my $primer_start, my $primer_end) = @_;
  my $diff_head;
  my $diff_tail;
  if ($strand eq '+') {
    $diff_head = $read_start - $primer_start;
    $diff_tail = $read_end - $primer_end;
  } else {
    $diff_head = $primer_end - $read_end;
    $diff_tail = $primer_start - $read_start;
  }
  return ($diff_head, $diff_tail);
}



sub compute_ref_length {
  my $cigar = shift @_;
  my %count = ('M'=>0, 'D'=>0);
  foreach my $type (keys %count) {
    if ($cigar =~ /$type/) {
      my @section = split /$type/, $cigar;
      foreach my $section (@section) {
	if ($section =~ /(\d+)$/) {
	  $count{$type} = $count{$type} + $1;
	}
      }
    }
  }
  my $length = $count{'M'} + $count{'D'};
  return $length;
}



sub compute_read_length {
  my $cigar = shift @_;
  my %count = ('M'=>0, 'I'=>0, 'S'=>0);
  foreach my $type (keys %count) {
    if ($cigar =~ /$type/) {
      my @section = split /$type/, $cigar;
      foreach my $section (@section) {
        if ($section =~ /(\d+)$/) {

          $count{$type} = $count{$type} + $1;
        }
      }
    }
  }
  my $length = $count{'M'} + $count{'I'} + $count{'S'};
  return $length;
}



sub compute_edit_distance {
  (my $cigar, my $MD) = @_;
  my $NM = 0;
  if ($cigar =~ /I/) {
    my @subcigar = split /I/, $cigar;
    foreach my $subcigar (@subcigar) {
      if ($subcigar =~ /(\d+)$/) {
	$NM = $NM + $1;
      }
    }
  }
  $MD =~ s/MD:Z://;
  while ($MD =~ /[A-Z]/) {
    $MD =~ s/^\d+\^?([A-Z]+)//;
    $NM = $NM + length($1);
  }
  return $NM;
}



sub compute_clip_length {
  (my $read_type, my $strand, my $read_start, my $read_end, my $target_start, my $target_end) = @_;
  my $clip_length;
  if (($read_type == 1 && $strand eq "+") || ($read_type == 2 && $strand eq "-")) {
    $clip_length = $target_end - $read_start + 1;
  }
  elsif (($read_type == 1 && $strand eq "-") || ($read_type == 2 && $strand eq "+")) {
    $clip_length = $read_end - $target_start + 1;
  }
}



sub clip_cigar {
  (my $rawcigar, my $read_type, my $strand, my $clip_length) = @_;
  my $newcigar = "";
  if (($read_type == 1 && $strand eq "+") || ($read_type == 2 && $strand eq "-")) {
    if ($clip_length > 0) {
      my $count = 0;
      CIGAR_CYCLE: while ($rawcigar) {
        $rawcigar =~ s/^(\d+)([A-Z])//;
        if ($2 eq 'M' || $2 eq 'D') {
          $count = $count + $1;
          if ($count > $clip_length) {
            my $over = $count - $clip_length;
            $newcigar = $over . $2 . $rawcigar;
            last CIGAR_CYCLE;
          }
        }
      }
    } else {
      $newcigar = $rawcigar;
    }
  }

  elsif (($read_type == 1 && $strand eq "-") || ($read_type == 2 && $strand eq "+")) {
    if ($clip_length > 0) {
      my $count = 0;
      CIGAR_CYCLE: while ($rawcigar) {
        $rawcigar =~ s/(\d+)([A-Z])$//;
        if ($2 eq 'M' || $2 eq 'D') {
          $count = $count + $1;
          if ($count > $clip_length) {
            my $over = $count - $clip_length;
            $newcigar = $rawcigar .$over . $2;
            last CIGAR_CYCLE;
          }
        }
      }
    } else {
      $newcigar = $rawcigar;
    }
  }
  unless ($newcigar =~ /\w/) { $newcigar = 'NA';}
  return $newcigar;
}



sub clip_MD {
  (my $rawMD, my $read_type, my $strand, my $clip_length) = @_;
  $rawMD =~ s/MD:Z://;
  my $newMD;
  if (($read_type == 1 && $strand eq "+") || ($read_type == 2 && $strand eq "-")) {
    if ($clip_length > 0) {
      my $count = 0;
      MD_CYCLE: while ($rawMD) {
        if ($rawMD =~ /^\d+/) {
          $rawMD =~ s/^(\d+)//;
          $count = $count + $1;
          if ($count >= $clip_length) {
            my $over = $count - $clip_length;
            $newMD = $over . $rawMD;
            last MD_CYCLE;
          }
        }
        elsif ($rawMD =~ /^\^?[A-Z]+/) {
          $rawMD =~ s/^(\^?)([A-Z]+)//;
          $count = $count + length($2);
          if ($count > $clip_length) {
            my $over = $count - $clip_length;
            my $remain_word = substr($2, -$over);
            $newMD = 0 . $1 . $remain_word . $rawMD;
            last MD_CYCLE;
          }
          elsif ($count == $clip_length) {
            $newMD = $rawMD;
            last MD_CYCLE;
          }
        }
      }
    } else {
      $newMD = $rawMD;
    }
  }

  elsif (($read_type == 1 && $strand eq "-") || ($read_type == 2 && $strand eq "+")) {
    if ($clip_length > 0) {
      my $count = 0;
      MD_CYCLE: while ($rawMD) {
        if ($rawMD =~ /\d+$/) {
          $rawMD =~ s/(\d+)$//;
          $count = $count + $1;
          if ($count >= $clip_length) {
            my $over = $count - $clip_length;
            $newMD = $rawMD . $over;
            last MD_CYCLE;
          }
        }
        elsif ($rawMD =~ /\^?[A-Z]+$/) {
          $rawMD =~ s/(\^?)([A-Z]+)$//;
          $count = $count + length($2);
          if ($count > $clip_length) {
            my $over = $count - $clip_length;
            my $remain_word = substr($2, -$over);
            $newMD = $rawMD . $1 . $remain_word . 0;
            last MD_CYCLE;
          }
          elsif ($count == $clip_length) {
            $newMD = $rawMD;
            last MD_CYCLE;
          }
        }
      }
    } else {
      $newMD = $rawMD;
    }
  }
  unless ($newMD =~ /\w/) { $newMD = 'NA';}
  $newMD = "MD:Z:$newMD";
  return $newMD;
}



sub clip_seq_qua {
  (my $seq, my $qua, my $length, my $read_type, my $strand) = @_;
  my $newseq;
  my $newqua;
  if (($read_type == 1 && $strand eq "+") || ($read_type == 2 && $strand eq "-")) {
    $newseq = substr($seq, -$length);
    $newqua = substr($qua, -$length);
  }
  elsif (($read_type == 1 && $strand eq "-") || ($read_type == 2 && $strand eq "+")) {
    $newseq = substr($seq, 0, $length);
    $newqua = substr($qua, 0, $length);
  }
  return ($newseq, $newqua);
}



sub edit_pos {
  (my $rawpos, my $clip_length) = @_;
  my $newpos;
  if ($clip_length < 0) { $clip_length = 0;}
  $newpos = $rawpos + $clip_length;
  return $newpos;
}



sub compute_mutation_count {
  (my $cigar, my $MD) = @_;
  $MD =~ s/MD:Z://;
  $cigar = $cigar . "T";
  my $count_I = 0; 
  my $count_D = 0;
  my $count_M = 0;
  if ($cigar =~ /I/) { $count_I = (split /I/, $cigar) - 1; }
  if ($cigar =~ /D/) { $count_D = (split /D/, $cigar) - 1; }
  $MD =~ s/\^([A-Z]+)//g;
  my @MD = split //, $MD;
  foreach my $base (@MD) {
    if ($base =~ /[A-Z]/) { $count_M++; }
  }
  my $mismatch_count = $count_I + $count_D + $count_M;
  return $mismatch_count;
}

  

