char *PerlScriptName[]={"rec_sum.pl","count.pl","p\
rocess_list.pl","make_license.pl","CCsed.script","\
msa2bootstrap.pl","t_coffee_dpa","t_coffee_dpa2","\
tc_generic_method.pl","rnapdb2protpdb.pl","generic\
_method.tc_method","clustalw_method.tc_method","ex\
tract_from_pdb","install.pl","clean_cache.pl","nat\
ure_protocol.pl","mocca","dalilite.pl","wublast.pl\
","blastpgp.pl","ncbiblast_lwp.pl","wublast_lwp.pl\
","RNAplfold2tclib.pl","fasta_seq2RNAplfold_templa\
tefile.pl","fasta_seq2hmmtop_fasta.pl","fasta_seq2\
consan_aln.pl","clustalw_aln2fasta_aln.pl","seq2na\
me_seq.pl","seq2intersection.pl","msf_aln2fasta_al\
n.pl","msa.pl","upp.pl","clustalo.pl","dca.pl","bl\
ast_aln2fasta_aln.pl","blast_xml2fasta_aln.pl","fa\
sta_aln2fasta_aln_unique_name.pl","newick2name_lis\
t.pl","excel2fasta.pl","nameseq2fasta.pl","any_fil\
e2unix_file.pl","EndList"};char *PerlScriptFile[]=\
{"use File::Copy;\nuse Env qw(HOST);\nuse Env qw(H\
OME);\nuse Env qw(USER);\n$x_field=0;\n$y_field=1;\
\n$y_field_set=1;\n$nyf=1;\n\n$interval=0;\n$file=\
\"stdin\";\n\n$print_avg=1;\n$print_sd=0;\n$print_\
sum=0;\n$print_n=0;\nforeach $value ( @ARGV)\n    \
{\n	if ($value ne $ARGV[$np]) \n	    {\n	    ;\n	 \
   }\n	elsif($value eq \"-s\")\n	     {\n	       $\
step=$ARGV[++$np];\n	       $np++;\n	     }\n	elsi\
f($value eq \"-print_all\")\n	    {\n	    $print_s\
d=$print_avg=$print_n=$print_sum=1;\n	    $np++;\n\
	    }\n	elsif($value eq \"-print_sum\")\n	    {\n\
	    $print_sum=1;\n	    $print_avg=0;\n	    $np++\
;\n	    }\n	elsif($value eq \"-print_n\")\n	    {\\
n	    $print_n=1;\n	    $print_avg=0;\n	    $np++;\
\n	    }\n	elsif($value eq \"-print_avg\")\n	    {\
\n	    $print_avg=1;\n	    $print_avg=0;\n	    $np\
++;\n	    }\n	elsif($value eq \"-sd\")\n	    {\n	 \
   $print_sd=1;\n	    $print_avg=0;\n	    $np++;\n\
	    }\n	elsif($value eq \"-h\")\n	    {\n	    $he\
ader=1;\n	    $np++;\n	    }\n	elsif ($value eq \"\
-i\")\n	    {\n	    $interval= $ARGV[++$np];\n	   \
 $np++;\n    	    }\n	elsif ($value eq \"-r2\")\n	\
    {\n	      $r2=1;\n	      \n	      $np  = $ARGV\
[++$np];\n	      $nsim= $ARGV[++$np];\n	      $np+\
+;\n    	    }\n	elsif ($value eq \"-r\")\n	    {\\
n	    $min= $ARGV[++$np];\n	    $max= $ARGV[++$np]\
;\n	    $np++;\n    	    }\n	\n	elsif ($value eq \\
"-x\")\n	    {\n	    $x_field= $ARGV[++$np]-1;\n	 \
   $np++;\n    	    }\n	elsif ($value eq \"-y\")\n\
	    {\n	    $nyf=0;  \n	    while ($ARGV[$np+1] &\
& !($ARGV[$np+1]=~/\\-/))\n	      {\n		$y_field[$n\
yf++]=$ARGV[++$np]-1;\n		$y_field_set=1;\n	      }\
\n\n	    $np++;\n    	    }\n	elsif ($value eq \"-\
file\")\n	    {\n	    $file= $ARGV[++$np];\n	    $\
file_set=1;\n	    $np++;\n    	    }       \n	elsi\
f ( $value eq \"h\" ||  $value eq \"-h\" || $value\
 eq \"-H\" || $value eq \"-help\" || $value eq \"h\
elp\")\n	  {\n	    print STDOUT \"data_analyse: An\
alyse and discretization of data\\n\";\n	    print\
 STDOUT \"       -file:    <file containing the da\
ta to analyze>,.<def=STDIN>\\n\";\n	    print STDO\
UT \"       -x: <field containing the X>,.........\
......<Def=0>\\n\";\n	    print STDOUT \"       -y\
: <field containing the Y>,...............<Def=1>\\
\n\";\n	    print STDOUT \"       -i:<Interval siz\
e on the X>,...............<Def=0>\\n\";\n	    pri\
nt STDOUT \"       -i:<0:only one interval>\\n\";\\
n	    print STDOUT \"       -r:<Range of the X>\\n\
\";\n	    print STDOUT \"       -s:<Step on the  X\
, 0 means non sliding bins>\\n\";\n	    print STDO\
UT \"       -sd: print standard deviation on the Y\
\";\n	    print STDOUT \"       -h  : print column\
 header \\n\";\n	    exit (0);\n	  }\n	elsif ($val\
ue=~/-/)\n	  {\n	    print \"$value is not a valid\
 FLAG[FATAL]\\n\";\n	    exit (0);\n	   } \n	elsif\
 ($list eq \"\") \n	    {\n	    $file=$ARGV[$np];\\
n	    $np++;\n	    }\n	\n	\n      }\n\n\n\n\n\nif \
($file eq \"stdin\")\n	{\n	$remove_file=1;\n	$file\
=\"tmp$$\";\n	open (F, \">$file\");\n	while (<STDI\
N>)\n		{\n		print F $_;\n		}\n	close (F);\n	 \n	;}\
\n\n\n\nif ($interval && $step)\n  {\n    my $nl;\\
n    open(F,$file);\n    while (<F>)\n      {\n	$l\
ine=$_;\n	\n	if (!/\\S/){next;}\n	@list=($line=~/(\
\\S+)/g);\n	$val{$nl}{x}=$list[$x_field];\n	$val{$\
nl}{y}=$list[$y_field[0]];\n	$nl++\n      }\n    c\
lose (F);\n    \n    for (my $a=$min; $a<($max+$in\
terval); $a+=$step)\n      {\n	my ($avgx, $avgy, $\
cn);\n	\n	my $rmin=$a-$interval;\n	my $rmax=$a;\n	\
$cn=0;\n	for (my $b=0; $b<$nl; $b++)\n	  {\n	    m\
y $x=$val{$b}{x};\n	    my $y=$val{$b}{y};\n	    i\
f ($x<=$rmax && $x>=$rmin)\n	      {\n		$avgx+=$x;\
\n		$avgy+=$y;\n		$cn++;\n		$tcn++;\n		$val{$b}{us\
ed}=1;\n	      }\n	  }\n	if ($cn)\n	  {\n	    $avg\
x/=$cn;\n	    $avgy/=$cn;\n	  }\n	printf \"%.3f %.\
3f %.3f\\n\", $avgx, $avgy, $avgx-$avgy;\n      }\\
n    for (my $a=0; $a<$nl; $a++)\n      {\n	if ( !\
$val{$a}{used})\n	  {\n	    print \"---$val{$a}{x}\
; $val{$a}{y}\\n\";\n	  }\n      }\n  }\nelse\n  {\
\n    if ($interval && $max)\n      {\n	$interval_\
size=($max-$min)/$interval;\n      }\n    elsif ($\
interval)\n      {\n	open(F,$file);  \n	my $set_ma\
x=0;\n	my $set_min=0;\n	while (<F>)\n	  {\n	    my\
 $v=$_;\n	    chomp($v);\n	    print \"--$v--\";\n\
	    \n	    if ($v<$min ||!$set_min){$set_min=1;$m\
in=$v;}\n	    if ($v>$max ||!$set_max){$set_max=1;\
$max=$v;}\n	  }\n	close (F);\n	print \"$min $max u\
uuu\";\n	$interval_size=($max-$min)/$interval;\n  \
    }\n    open(F,$file);  \n    while (<F>)\n    \
  {\n	$line=$_;\n	if (!/\\S/){next;}\n	@list=($lin\
e=~/(\\S+)/g);\n	\n	if ($interval==0){$bin=0;}\n	e\
lse{$bin=int (($list[$x_field]-$min)/($interval_si\
ze));}\n	\n	\n	if ($bin && $bin==$interval){$bin--\
;}\n	for ( $a=0; $a<$nyf; $a++)\n	  {\n	    $sum{$\
a}{$bin}+=$list[$y_field[$a]];\n	    $sum2{$a}{$bi\
n}+=$list[$y_field[$a]]*$list[$y_field[$a]];\n	   \
 $n{$a}{$bin}++;\n	  }\n      }\n    \n    if (!$i\
nterval){$interval=1;}\n    for ( $a=0; $a<$interv\
al; $a++)\n      {\n	printf ( \"%4d %4d \", $inter\
val_size*$a, $interval_size*($a+1));\n	for ( $b=0;\
 $b<$nyf; $b++)	\n	  {\n	    $i=$interval*$a;\n	  \
  if ( $n{$b}{$a}==0)\n	      {\n		$avg=0;\n		$sd=\
0;\n	      }\n	    else\n	      {\n		$avg=$sum{$b}\
{$a}/$n{$b}{$a};\n		$sd=sqrt($sum2{$b}{$a}*$n{$b}{\
$a}-$sum{$b}{$a}*$sum{$b}{$a})/($n{$b}{$a}*$n{$b}{\
$a});\n	      }\n	    if ($print_n) {printf \"%15.\
4f \", $n{$b}{$a};}\n	    if ($print_sum){printf \\
"%15.4f \", $sum{$b}{$a};}\n	    if ($print_avg){p\
rintf \"%15.4f \", $avg}\n	    if ($print_sd) {pri\
ntf \"%15.4f \", $sd;}\n	  }\n	printf (\"\\n\");\n\
      }\n  }\n\nif ( $remove_file){unlink $file;}\\
n","use File::Copy;\nuse Env qw(HOST);\nuse Env qw\
(HOME);\nuse Env qw(USER);\n\nforeach $v (@ARGV){$\
cl.=$v;}\n\n\nif ( $cl=~/-k(\\d+)/){$k=$1;}\nelse \
{$k=1;}\nif ( $cl=~/-w(\\d+)/){$w=$1;}\nelse {$w=-\
1;}\nif ( $cl=~/-p(\\d+)/){$p=$1;}\nelse {$p=-1;}\\
n\nwhile (<STDIN>)\n  {\n    @l=($_=~/(\\S+)/g);\n\
    $v=$l[$k-1];\n    if ( !$h{$v}){@ll=($v, @ll);\
}\n    \n    if ( $w==-1)\n      {$h{$v}++;}\n    \
else\n      {$h{$v}+=$l[$w-1];}\n\n    if ($p!=-1)\
{$print{$v}=$l[$p-1];}\n\n  }\nforeach $v (@ll)\n \
 {\n    print \"$v $print{$v} $h{$v}\\n\";\n  }\n"\
,"\nuse Env qw(HOST);\nuse Env qw(HOME);\nuse Env \
qw(USER);\n$random_tag=int (rand 10000)+1;\n$uniqu\
e_prefix=\"$$.$HOST.$random_tag\";\n$queue=\"disti\
llery.and.mid\";\n$monitor=0;\n$stderr_file=\"/dev\
/null\";\n$stdio_file=\"/dev/null\";\n$log_file=\"\
/dev/null\";\n$pause_time=0;\n$max_sub_jobs=60;\n$\
min_sub_jobs=30;\n$output_all=0;\n$var='\\$';\n\nf\
oreach $value ( @ARGV)\n    {\n	if ($value ne $ARG\
V[$np]) \n	    {\n	    ;\n	    }\n	elsif ($value e\
q \"-max_sub_jobs\")\n	    {\n	    $max_sub_jobs= \
$ARGV[++$np];\n	    $np++;\n    	    }	\n	elsif ($\
value eq \"-min_sub_jobs\" )\n	    {\n	    $min_su\
b_jobs= $ARGV[++$np];\n	    $np++;\n    	    }\n	e\
lsif ($value eq \"-para\")\n	    {\n	    $para=1;\\
n	    $monitor=1;\n	    $np++;\n    	    }\n	elsif\
 ($value eq \"-monitor\") \n	    {\n	    $monitor=\
1;\n	    $np++;\n	    }\n	elsif ($value eq \"-no_m\
onitor\") \n	    {\n	    $monitor=0;\n	    $np++;\\
n	    }\n	elsif ($value eq \"-queue\")\n	    {\n	 \
   $queue=$ARGV[++$np];\n	    $np++;\n	    }	\n	el\
sif ($value eq \"-stderr_file\")\n	    {\n	    $st\
derr_file=$ARGV[++$np];\n	    $np++;\n	    }\n	els\
if ($value eq \"-stdio_file\")\n	    {\n	    $stdi\
o_file=$ARGV[++$np];\n	    $np++;\n	    }\n	elsif \
($value eq \"-output_all\")\n	    {\n	    $output_\
all=1;\n	    $np++;\n	    }\n	elsif ($value eq \"-\
pause\") \n	    {\n	    $pause_time=$ARGV[++$np];\\
n	    $np++;\n	    }\n	elsif ($value eq \"-log\")\\
n	      {\n	       $log=1;\n	       \n	       if (\
$ARGV[$np+1]=~/\\-\\S+/) \n	          {\n		  $log_\
file=\"stderr\";\n	          }\n	       else \n	  \
        {\n		  $log_file=$ARGV[++$np]; \n		  ++$np\
;\n		 \n	          }\n	      }\n	elsif ( $value eq\
 \"-com\")\n	    {\n		\n		if (!$ARGV[$np+1]=~/^\\'\
/) { $com=$ARGV[++$np];}\n		else {$com=$ARGV[++$np\
];}\n\n	     $np++;\n	    }\n	elsif ( $value eq \"\
-check\")\n	  {\n	    \n	    if (!$ARGV[$np+1]=~/^\
\\'/) { $check=$ARGV[++$np];}\n	    else {$check=$\
ARGV[++$np];}\n	    $np++;\n	  }\n	elsif ($com eq \
\"\") \n	    {\n	    $com_set=1;\n	    $com=$ARGV[\
$np];\n	    \n	    $np++;\n	    }\n	elsif ($list e\
q \"\") \n	    {\n	    $list_set=1;\n	    $list=$A\
RGV[$np];\n	    $np++;\n	    }\n	elsif ( $var_set \
eq \"\")\n	    {\n	    $var_set=1;\n	    $var=$ARG\
V[$np];\n	    $np++;\n	    }\n	}\n\n\n\n\nif ( $co\
m eq \"\"){print \"You Need to Provide a Command [\
FATAL]\\n\";\n	      die;\n	     }\n\n\n\nif ($lis\
t_set==0) \n    {\n    $x= int (rand 100000)+1;\n \
   $tmp_file_name=\"tmp_file_$x\";\n    open ( TMP\
, \">$tmp_file_name\");\n    while (<STDIN>)\n    \
  {\n	print TMP $_;\n      }\n    close (TMP);\n  \
  open (F, $tmp_file_name);\n    }\nelse \n    {\n\
    open (F, $list);\n    }\n\nif ($para==0) \n   \
 {\n\n     @tc_list= <F>;\n     close (F); \n     \
\n     foreach $val(@tc_list) \n	    {\n	      \n	\
      \n	      \n	      $loc_com=$com;\n	      if \
($check){$loc_check=$check;}\n	      \n	      @i_v\
al=($val=~/([^\\s]+)/g);\n	      \n	      if ( $#i\
_val==0)\n		{\n		  if ($check){$loc_check=~s/$var/\
$i_val[0]/g;}\n		  $loc_com=~s/$var/$i_val[0]/g;\n\
		}\n	      else\n		{\n		  for ($n=1; $n<=$#i_val+\
1;$n++ )\n		    {\n		      \n		      $sub=\"$var$n\
\";\n		      \n		      $loc_com=~s/$sub/$i_val[$n-\
1]/g;\n		      if ($check){$loc_check=~s/$var/$i_v\
al[0]/g;}\n		    }\n		}\n	      if ( $check && -e \
$loc_check)\n		{\n		  print STDERR \"skipping $loc\
_com...\\n\";\n		  }\n	      else\n		{\n		  system\
 \"$loc_com\";\n		}\n	    }\n    exit;\n    }\n\ne\
lsif ($para==1) \n    {\n    print STDERR \"do par\
allel execution of: \\\"$com $list\\\"\\n\";\n    \
\n    if ($log==1) \n	{\n	if ($log_file eq \"stdou\
t\" || $log_file eq \"stderr\" ) \n		{\n		$log_fil\
e=\"\";\n	        }\n\n        else \n		{\n		syste\
m \"echo LOG FILE> $log_file\";\n		\n	        }\n	\
}\n    else	\n	{\n	open ( OUT, \">/dev/null\");\n	\
}\n	\n    \n    $id=0;\n    $n_sub=0;\n    while (\
$val=<F>) \n	    {	    	    \n	    $job_log[$id]=\\
"$HOME/tmp/$unique_prefix.$id.log_file\";\n	    \n\
	    $job=$unique_prefix.\"_$id\";\n	    open (JOB\
, \">$job\");\n	    \n	    $loc_com=$com;\n	    ch\
op $val;\n\n	    $loc_com=~s/\\$/$val/g;\n	 \n	   \
 print JOB \"#!/bin/csh\\n\";\n	    print JOB \"#\\
\$ -cwd\\n\";\n	    print JOB \"#\\$ -N $unique_pr\
efix\\n\";\n	    if ($queue && !($queue eq \" \"))\
 {print JOB \"#\\$ -l $queue\\n\";}\n	    print JO\
B \"#\\n\";	    \n            print JOB \"$loc_com\
\\n\";\n	    print JOB \"echo FINISHED  >> $job_lo\
g[$id]\\n\";\n	    print JOB \"pwd\\n\";\n	    \n	\
    close (JOB);\n	    if ( $output_all==1)\n		{\n\
		system \"qsub $job >  $unique_prefix\";		\n	    \
    }\n	    else\n		{system \"qsub $job -e $stderr\
_file -o $stdio_file >$unique_prefix\";	        \n\
	        } \n\n\n\n	    print STDERR \"$id: $outpu\
t_all\\n\";\n	    $n_sub++;\n	    if ( $max_sub_jo\
bs && $n_sub==$max_sub_jobs) \n		{\n		$n_sub=monit\
or_process($min_sub_jobs,@job_log); 		 \n		\n	    \
    }	\n	   \n            unlink $unique_prefix;\n\
	    sleep $pause_time;\n	    $id++;\n	    }\n\n  \
  close (OUT);\n    close (F);\n\n    print STDERR\
 \"Your $id Jobs Have Been Submited (NAME=$unique_\
prefix)\\n\";\n    monitor_process (0, @job_log);\\
n    foreach $file(@job_log) {if (-e $file) {unlin\
k($file);}}\n    \n    }\n\nsub monitor_process ( \
@job_list)\n    {\n    my (@job_list)=@_;\n    my \
$min_sub_jobs=shift (@job_list);\n    my $n_sub_jo\
bs;\n    my $finished;\n    my $n=0;\n\n    $n_sub\
_jobs=-1;\n    $finished=0;\n    print STDERR \"\\\
nMonitor Batch: [$min_sub_jobs]\";\n       \n    w\
hile (!$finished && (($n_sub_jobs>$min_sub_jobs)||\
 $n_sub_jobs==-1) ) \n	{\n	$finished=1;\n	$n_sub_j\
obs=0;\n	$n=0;\n	foreach $file (@job_list)\n	     \
   {\n	\n		if (-e $file){;}\n		else \n		    {\n		 \
   $finished=0; $n_sub_jobs++;\n	            }\n	 \
       }\n	system \"sleep 1\";\n        }\n    \n \
   return $n_sub_jobs;\n    }\n    \n    \nif ($tm\
p_file_name){unlink($tmp_file_name);}\n","\n\nfore\
ach ($np=0; $np<=$#ARGV; $np++)\n    {\n    $value\
=$ARGV[$np];\n\n    if ($value eq \"-file\")\n    \
  {\n      $file= $ARGV[++$np];\n      }\n    elsi\
f ($value eq \"-type\")\n      {\n        $type= $\
ARGV[++$np];\n      }\n    elsif ($value eq \"-ins\
titute\")\n      {\n        $institute= $ARGV[++$n\
p];\n      }\n    elsif ($value eq \"-author\")\n \
     {\n        $author= $ARGV[++$np];\n      }\n \
   elsif ($value eq \"-date\")\n      {\n        $\
date= $ARGV[++$np];\n      }\n     elsif ($value e\
q \"-program\")\n      {\n        $program= $ARGV[\
++$np];\n      }\n    elsif ($value eq \"-email\")\
\n      {\n        $email= $ARGV[++$np];\n      }\\
n    else\n      {\n	print \"$value is an unkown a\
rgument[FATAL]\\n\";\n	exit (1);\n      }\n  }\n\n\
\n\nopen F, $file || die;\nprint $INSTITUTE;\nif (\
 $type eq \"c\"){print \"/************************\
******COPYRIGHT NOTICE****************************\
***/\\n\";}\nif ( $type eq \"perl\"){print \"#####\
#########################COPYRIGHT NOTICE#########\
#####################/\\n\";}\nif ( $type eq \"txt\
\"){print \"-------------------------------COPYRIG\
HT NOTICE------------------------------/\\n\";}\n\\
n\nwhile (<F>)\n  {\n  s/\\$INSTITUTE/$institute/g\
;\n  s/\\$AUTHOR/$author/g;\n  s/\\$DATE/$date/g;\\
n  s/\\$PROGRAM/$program/g;  \n  s/\\$EMAIL/$email\
/g;  \n  if ( $type eq \"txt\"){print $_;}\n  elsi\
f ($type eq \"c\"){chop $_; print \"\\/*$_*\\/\\n\\
";}\n  elsif ($type eq \"perl\"){print \"\\#$_\";}\
\n}\nclose (F);\nif ( $type eq \"c\"){print \"/***\
***************************COPYRIGHT NOTICE*******\
************************/\\n\";}\nif ( $type eq \"\
perl\"){print \"##############################COPY\
RIGHT NOTICE##############################/\\n\";}\
\nif ( $type eq \"txt\"){print \"-----------------\
--------------COPYRIGHT NOTICE--------------------\
----------/\\n\";}\n\n","\nwhile (<>)	\n	{\n	s/\\=\
cc/123456789/g;\n	s/\\bcc/\\$\\(CC\\)/g;\n	s/12345\
6789/\\=cc/g;\n	print $_;\n	}\n\n","$version=\"1.0\
0\";\n$rseed= int(rand(100000))+1;\n\n\nif ( $#ARG\
V==-1)\n  {\n    print \"msa2bootstrap -i <input_f\
ile> -input <seq|msa|matrix|tree> -n <N-Boostrap> \
-o <outtree> -tmode <nj|upgma|parsimony|ml> -dmode\
 <kimura> -alignpg <t_coffee | muscle | clustalw> \
-rtree <file> -stype <prot|cdna|dna> -recompute -s\
ystem <cygwin|unix>\";\n    print \"\\n\\t-i: inpu\
t file, can be sequneces, msa, matrix, trees, type\
 is specified via -input\";\n    print \"\\n\\t-in\
put: Type of input data\";\n    print \"\\n\\t\\tm\
sa: msa in fasta format\";\n    print \"\\n\\t\\ts\
eq: compute an msa with -alignpg\";\n    print \"\\
\n\\t\\tmatrix: phylipp distance matrix fed direct\
ly to method -tmode [caveat: tmode=nj or upgma]\";\
\n    print \"\\n\\t\\ttree: list of newick trees \
directly fed to consence in order to generate a bo\
otstraped tree\";\n    \n    print \"\\n\\t-n: num\
ber of bootstrap replicates\";\n    print \"\\n\\t\
-o: name of the output tree. Files are not overwri\
tten. Use -recompute to overwrite existing file\";\
\n    print \"\\n\\t-tmode: tree mode: nj|upgma|pa\
rsimony|ml\";\n    print \"\\n\\t-dmode: distance \
mode\";\n    print \"\\n\\t-alignpg: program for a\
ligning sequences (t_coffee=default)\";\n    print\
 \"\\n\\t-rtree: replicate tree file (default: no \
file)\";\n    print \"\\n\\t-rmsa: replicate msa f\
ile (default: no file)\";\n    print \"\\n\\t-rmat\
: replicate matrix file (default: no file)\";\n   \
 print \"\\n\\t-stype: sequence type: protein, dna\
 or cdna\";\n    print \"\\n\\t-recompute: force f\
iles to be overwritten\";\n    print \"\\n\\t-syst\
em: cygwin|unix\";\n      \n\n    \n    &my_exit (\
EXIT_FAILURE);\n  }\nforeach $arg (@ARGV){$command\
.=\"$arg \";}\n\nprint \"CLINE: $command\\n\";\n$t\
hreshold=100;\n$trim_msa=0;\n$stype=\"prot\";\npri\
nt \"msa2bootstrap \";\n\n$system=\"cygwin\";\nif(\
($command=~/\\-system (\\S+)/))\n  {\n    $system=\
$1;\n    if ( $system eq \"cygwin\")\n      {\n	$e\
xec_extension=\".exe\";\n      }\n    elsif ( $sys\
tem eq \"unix\")\n      {\n	$exec_extension=\"\";\\
n	print \"system=Unix\";die;\n      }\n    else\n \
     {\n	print \"msa2boostrap: -system=$system is \
an unknown mode [FATAL]\\n\"; die;\n      }\n    \\
n    print \"-system $system \";\n  }\nif(($comman\
d=~/\\-stype (\\S+)/))\n  {\n    $stype=$1;\n  }\n\
print \"-stype=$stype \";\n\n\n\nif(($command=~/\\\
-i (\\S+)/))\n  {\n    $msa=$1;\n    print \"-i $m\
sa \";\n  }\n\nif(($command=~/\\-rtree (\\S+)/))\n\
  {\n    $rtree=$1;\n    print \"-rtree=$rtree \";\
\n  }\n\nif(($command=~/\\-rmsa (\\S+)/))\n  {\n  \
  $rmsa=$1;\n  }\nif(($command=~/\\-rmat (\\S+)/))\
\n  {\n    $rmat=$1;\n  }\n$input=\"seq\";\nif(($c\
ommand=~/\\-input (\\S+)/))\n  {\n    $input=$1;\n\
  }\nprint \"-input=$input \";\n\n$dmode=\"kimura\\
";\nif(($command=~/\\-dmode (\\S+)/))\n  {\n    $d\
mode=$1;\n  }\nprint \"-dmode=$dmode \";\n$alignpg\
=\"muscle\";\nif(($command=~/\\-alignpg (\\S+)/))\\
n  {\n    $alignpg=$1;\n  }\nprint \"-alignpg=$dmo\
de \";\n\n$tmode=\"nj\";\nif(($command=~/\\-tmode \
(\\S+)/))\n  {\n    $tmode=$1;\n  }\nprint \"-tmod\
e=$tmode \";\n$recompute=0;\nif(($command=~/\\-rec\
ompute/))\n  {\n    $recompute=1;\n    print \"-re\
compute \";\n  }\n\n$out=$msa;\n$out=~s/\\..*//;\n\
$out.=\".bph\";\nif(($command=~/\\-o (\\S+)/))\n  \
{\n    $out=$1;\n    \n  }\nprint \"-out=$out \";\\
nif (-e $out && !$recompute)\n  {\n    print \"\\n\
No Computation Required $out already exists\\n\";\\
n    &my_exit (EXIT_SUCCESS);\n    \n  }\n\n$n=100\
;\nif(($command=~/\\-n (\\d+)/))\n  {\n    $n=$1;\\
n  }\nprint \"-n=$n \";\n$seed=3;\nif(($command=~/\
\\-s (\\d+)/))\n  {\n    $seed=$1;\n  }\nprint \"-\
s=$seed \";\n\nif(($command=~/\\-run_name (\\d+)/)\
)\n  {\n    $suffix=$1;\n  }\nelse\n  {\n    $msa=\
~/([^.]+)/;\n    $suffix=$1;\n  }\nprint \"-run_na\
me=$suffix\\n\";\n\n\nif ( $input eq \"seq\")\n  {\
\n    $seq=$msa;\n    $msa=\"$suffix.prot_msa\";\n\
    \n    if ($stype eq \"cdna\")\n      {\n	$cdna\
_seq=$seq;\n	$clean_cdna_seq=&vtmpnam();\n	$seq=&v\
tmpnam();\n	`t_coffee -other_pg seq_reformat -in $\
cdna_seq -action +clean_cdna >$clean_cdna_seq`;\n	\
`t_coffee -other_pg seq_reformat -in $clean_cdna_s\
eq -action +translate >$seq`;\n	\n      }\n\n    i\
f (!-e $msa || $recompute)\n      {\n	print \"\\n#\
####   Compute an MSA With $alignpg\\n\";\n	\n	if \
( $alignpg eq \"t_coffee\")\n	  {`$alignpg $seq -o\
utfile=$msa >/dev/null 2>/dev/null`;}\n	elsif ( $a\
lignpg eq \"muscle\")\n	  {\n	    `$alignpg -in $s\
eq > $msa 2>/dev/null`;\n	  }\n	elsif ( $alignpg e\
q \"clustalw\")\n	  {\n	    `$alignpg -infile=$seq\
 -outfile=$msa -quicktree >/dev/null 2>/dev/null`;\
\n	  }\n	elsif ( $align eq \"mafft\")\n	  {\n	    \
`$alignpg $seq > $msa >/dev/null 2>/dev/null`;\n	 \
 }\n	else\n	  {\n	    `$alignpg -in=$seq -outfile=\
$msa`;\n	  }\n      }\n    if (!-e $msa)\n      {\\
n	print \"\\nError: $alignpg Could Not produce the\
 MSA $msa [FATAL]\\n\";\n      }\n\n    if ($stype\
 eq \"cdna\")\n      {\n	$msa2=\"$suffix.cdna_msa\\
";\n	`t_coffee -other_pg seq_reformat -in $clean_c\
dna_seq -in2 $msa -action +thread_dna_on_prot_aln \
-output fasta_aln  >$msa2`;\n	$msa=$msa2;\n      }\
\n    \n    $input=\"msa\";\n  }\n\n\n\n$seqboot_o\
=&vtmpnam();\n$seqboot_c=&vtmpnam();\n\n$protdist_\
o=&vtmpnam();\n$protdist_c=&vtmpnam();\nif ( $inpu\
t eq \"msa\")\n  {\n    if ($tmode eq \"nj\" || $t\
mode eq \"upgma\"){$input=\"matrix\";}\n    \n    \
$lmsa= &vtmpnam ();\n    `t_coffee -other_pg seq_r\
eformat -in $msa -output phylip_aln > $lmsa`;\n   \
 \n    if ( -e \"outfile\"){unlink (\"outfile\");}\
\n    # run seqboot\n  \n    if ( $n>1)\n      {\n\
	print \"Run SeqBoot .....\";\n	open (F, \">$seqbo\
ot_c\");\n	print F \"$lmsa\\nR\\n$n\\nY\\n$seed\\n\
\";\n	close (F);\n	`seqboot$exec_extension  < $seq\
boot_c`;\n	if ( -e \"outfile\"){ print \"[OK]\\n\"\
;}\n	else { print \"[FAILED]\\n\";&my_exit (EXIT_F\
AILURE);}\n	`mv outfile $seqboot_o`;\n      }\n   \
 else\n      {\n	`cp $lmsa $seqboot_o`;\n      }\n\
\n    if ($rmsa){`cp $seqboot_o $rmsa`;}\n    \n  \
  if ($tmode eq \"nj\" || $tmode eq \"upgma\")\n  \
    {\n	if ( $stype eq \"prot\")\n	  {\n	    # run\
 protdist\n	    print \"Run Protdist [dmode=$dmode\
]\";\n	    if ($dmode eq \"kimura\")\n	      {\n		\
$dmode=\"P\\nP\\nP\";\n	      }\n	    else\n	     \
 {\n		print \"\\n$dmode is an unknown mode for Pro\
tdist [FATAL:msa2bootstrap.pl]\\n\";\n		&my_exit (\
EXIT_FAILURE);\n	      }\n	    open (F, \">$protdi\
st_c\");\n	    if ($n>1){print F \"$seqboot_o\\n$d\
mode\\nM\\nD\\n$n\\nY\\n\";}\n	    else {printf F \
\"$seqboot_o\\n$dmode\\nY\\n\";}\n	    close (F);\\
n	    `protdist$exec_extension  < $protdist_c`;\n	\
    if ( -e \"outfile\"){ print \"[OK]\\n\";}\n	  \
  else { print \"[FAILED]\\n\";&my_exit (EXIT_FAIL\
URE);}\n	    `mv outfile $protdist_o`;\n	 \n	  }\n\
	elsif ( $stype eq \"cdna\" || $stype eq \"dna\")\\
n	  {\n	    print \"Run dnadist [dmode=default\";\\
n	    open (F, \">$protdist_c\");\n	    if ($n>1){\
print F \"$seqboot_o\\nM\\nD\\n$n\\nY\\n\";}\n	   \
 else {printf F \"$seqboot_o\\nY\\n\";}\n	    clos\
e (F);\n	    `protdist$exec_extension  < $protdist\
_c`;\n	    if ( -e \"outfile\"){ print \"[OK]\\n\"\
;}\n	    else { print \"[FAILED]\\n\";&my_exit (EX\
IT_FAILURE);}\n	    `mv outfile $protdist_o`;\n	  \
}\n      }\n  }\nelsif ( $input eq \"matrix\")\n  \
{\n    $protdist_o=&vtmpnam();\n    print \"MSA: $\
msa\\n\";\n    `cp $msa $protdist_o`;\n    $n=1;\n\
  }\n\n\n\n\n\n$nb_o=&vtmpnam();\n$nb_c=&vtmpnam()\
;\nif ($input eq \"matrix\" && $tmode ne \"parsimo\
ny\" && $tmode ne \"ml\")\n  {\n    print \"Run ne\
ighbor [tmode=$tmode]\";\n\n    if ($tmode eq \"nj\
\")\n      {\n	$tmode=\"\\nN\\nN\";\n      }\n    \
elsif ( $tmode eq \"upgma\")\n      {\n	$tmode = \\
"\\nN\";\n      }\n    else\n      {\n	print \"\\n\
 ERROR: $tmode is an unknown tree computation mode\
\\n\";\n	&my_exit (EXIT_FAILURE);\n      }\n\n    \
open (F, \">$nb_c\");\n    if ($n>1){print F \"$pr\
otdist_o$tmode\\nM\\n$n\\n$seed\\nY\\n\";}\n    el\
se {print F \"$protdist_o$tmode\\nY\\n\";}\n    cl\
ose (F);\n\n    `neighbor$exec_extension  < $nb_c`\
;\n    if ( -e \"outtree\"){ print \"[Neighbor OK]\
\\n\";}\n    else { print \"[FAILED]\\n\";&my_exit\
 (EXIT_FAILURE);}\n    `mv outtree $nb_o`;\n    un\
link (\"outfile\");\n  }\nelsif ($input eq \"msa\"\
 && $tmode eq \"parsimony\")\n  {\n    if ( -e \"o\
utfile\"){unlink (\"outfile\");}\n    if ( -e \"ou\
ttree\"){unlink (\"outtree\");}\n    \n    if ($st\
ype eq \"prot\")\n      {\n	print \"Run protpars [\
tmode=$tmode]\";\n	open (F, \">$nb_c\");\n	if ($n>\
1){print F \"$seqboot_o\\nM\\nD\\n$n\\n$seed\\n10\\
\nY\\n\";}\n	else {print F \"$seqboot_o\\nY\\n\";}\
\n	close (F);\n	`protpars$exec_extension  < $nb_c`\
;\n      }\n    elsif ( $stype eq \"dna\" || $styp\
e eq \"cdna\")\n      {\n	print \"Run dnapars [tmo\
de=$tmode]\";\n	open (F, \">$nb_c\");\n	if ($n>1){\
print F \"$seqboot_o\\nM\\nD\\n$n\\n$seed\\n10\\nY\
\\n\";}\n	else {print F \"$seqboot_o\\nY\\n\";}\n	\
close (F);\n	`dnapars$exec_extension  < $nb_c`;\n \
     }\n    if ( -e \"outtree\"){ print \"[OK]\\n\\
";}\n    else { print \"[FAILED]\\n\";&my_exit (EX\
IT_FAILURE);}\n    `mv outtree $nb_o`;\n   unlink \
(\"outfile\");\n  }\nelsif ($input eq \"msa\" && $\
tmode eq \"ml\")\n  {\n    if ( -e \"outfile\"){un\
link (\"outfile\");}\n    if ( -e \"outtree\"){unl\
ink (\"outtree\");}\n    \n    if ($stype eq \"pro\
t\")\n      {\n	print \"Error: ML impossible with \
Protein Sequences [ERROR]\";\n	&my_exit (EXIT_FAIL\
URE);\n      }\n    elsif ( $stype eq \"dna\" || $\
stype eq \"cdna\")\n      {\n	print \"Run dnaml [t\
mode=$tmode]\";\n	open (F, \">$nb_c\");\n	if ($n>1\
){print F \"$seqboot_o\\nM\\nD\\n$n\\n$seed\\n10\\\
nY\\n\";}\n	else {print F \"$seqboot_o\\nY\\n\";}\\
n	close (F);\n	`dnaml$exec_extension  < $nb_c`;\n \
     }\n    if ( -e \"outtree\"){ print \"[OK]\\n\\
";}\n    else { print \"[FAILED]\\n\";&my_exit (EX\
IT_FAILURE);}\n    `mv outtree $nb_o`;\n   unlink \
(\"outfile\");\n  }\n\n\nelse\n  {\n    `cp $msa $\
nb_o`;\n    $n=2;\n  }\n\nif ($rmsa && -e $seqboot\
_o){print \"\\nOutput List of $n Replicate MSA: $r\
msa\\n\";`cp $seqboot_o $rmsa`;}\nif ($rmat && -e \
$protdist_o){print \"\\nOutput List of $n Replicat\
e MATRICES: $rmat\\n\";`cp $protdist_o $rmat`;}\ni\
f ($rtree && -e $nb_o){print \"\\nOutput List of $\
n Replicate TREES: $rtree\\n\";`cp $nb_o $rtree`;}\
\n\n\n\n$con_o=&vtmpnam();\n$con_c=&vtmpnam();\nif\
 ($n >1)\n  {\n    print \"Run Consense.....\";\n \
   open (F, \">$con_c\");\n    print F \"$nb_o\\nY\
\\n\";\n    close (F);\n    `consense$exec_extensi\
on  < $con_c`;\n    if ( -s \"outtree\"  > 0) { pr\
int \"[OK]\\n\";}\n    else { print \"[FAILED]\\n\\
";&my_exit (EXIT_FAILURE);}\n    `mv outtree $con_\
o`;\n    unlink (\"outfile\");\n  }\nelse\n  {\n  \
  `cp $nb_o $con_o`;\n  }\n\n\n`cp $con_o $out`;\n\
if ( !-e $out)\n  {\n    print \"Tree Computation \
failed [FAILED]\\n\";\n    &my_exit (EXIT_FAILURE)\
;\n  }\nelsif ($n>1)\n  {\n    print \"\\nOutput B\
ootstrapped Tree: $out\\n\";\n    $avg=`t_coffee -\
other_pg seq_reformat -in $out -action +avg_bootst\
rap`;\n    $avg=~s/\\n//g;\n    print \"$avg\\n\";\
\n  }\nelse\n  {\n    print \"\\nOutput Tree: $out\
\\n\";\n  }\n\nopen (F, \"$out\");\nwhile (<F>)\n \
 {\n    \n    $tree.=$_;\n  }\nclose (F);\n$tree=~\
s/\\n//g;\nprint \"BPH: $tree\\n\";\n\n\n&my_exit \
(EXIT_SUCCESS);\n\nsub my_exit \n  {\n    my $m=@_\
[0];\n    &clean_vtmpnam();\n    exit ($m);\n  }\n\
sub vtmpnam \n  {\n    my $file;\n\n\n    $ntmp++;\
\n    $file=\"tmp4msa2bootstrap.$rseed.$$.$ntmp\";\
\n    \n    push (@tmpfile, $file);\n    return $f\
ile;\n  }\nsub clean_vtmpnam \n  {\n    my $t;\n  \
  foreach $t (@tmpfile)\n      {\n	if ( -e $t){unl\
ink ($t)};\n      }\n  }\n","use Env;\n$seq_reform\
at=\"t_coffee -other_pg seq_reformat \";\n$Version\
Tag=\"1.00\";\n$step=1;\n$unset=\"\";\n$scoreT1=$s\
coreT2=$nseqT=$dp_limit=$unset;\n@tl=();\nchomp($t\
c_version=`t_coffee -version`);$tc_version=~s/PROG\
RAM: //;\n\n\nprint STDERR \"\\n******************\
***********************************************\";\
\nprint STDERR \"\\n*           HIGH LEVEL PROGRAM\
: T-COFFEE_DPA Version $VersionTag\";\nprint STDER\
R \"\\n*           LOW  LEVEL PROGRAM: $tc_version\
 \";\nprint STDERR \"\\n**************************\
***************************************\";\n\nif (\
!@ARGV)\n  {\n    print \"t_coffee_dpa accepts eve\
ry t_coffee_flag.\\nType t_coffee to obtain a list\
\\n\";\n    print \"Requires $TC_VERSION\\n\";\n  \
  print \"Requires \";\n    print \"t_coffee_dpa s\
pecific flags:\\n\";\n    print \"\\t-dpa_master_a\
ln....................Master alignment: provided O\
R computed\\n\";\n    print \"\\t-dpa_master_aln..\
..................By default, Computed with t_coff\
ee -very_fast\\n\";\n    print \"\\t-dpa_master_al\
n=<file>.............Use file, (must be an aln in \
Fasta or ClustalW\\n\";\n    print \"\\t-dpa_maste\
r_aln=<program>..........Compute aln with pg -in s\
eq -out aln`\\n\";\n    print \"\\t-dpa_maxnseq...\
....................Maximum number of sequences in\
 subgroups\\n\";\n    print \"\\t-dpa_min_score1..\
..................Minimum Id for two sequences to \
be grouped in ref_aln\\n\";\n    print \"\\t-dpa_m\
in_score2....................Minimum Id within a s\
ubgroup\\n\";\n    print \"\\t-dpa_debug..........\
...............Keep Tmp File (for debug purpose)\\\
n\\n\";\n    \n    exit (0);\n  }\nforeach $arg (@\
ARGV)\n  {\n    $arg_list.=\" $arg\";\n  }\n$arg_l\
ist=~s/[=,;]/ /g;\n\n\n($seq0, $arg_list)=&extract\
_val_from_arg_list(\"^\",$arg_list, \"SPLICE\",\"u\
nset\");\n($seq1, $arg_list)=&extract_val_from_arg\
_list(\"-seq\",$arg_list, \"SPLICE\",\"unset\");\n\
($seq2, $arg_list)=&extract_val_from_arg_list(\"-i\
n\",$arg_list, \"KEEP\",\"unset\");\n($seq3, $arg_\
list)=&extract_val_from_arg_list(\"-infile\",$arg_\
list, \"SPLICE\",\"unset\");\n($prf,  $arg_list)=&\
extract_val_from_arg_list(\"-profile\",$arg_list, \
\"SPLICE\",\"unset\");\n\n$gl{'Seq'}=$seq=&vtmpnam\
();#file containing all the sequences\n\n   #1-rem\
ove sequences from -in\nif ( $arg_list =~/\\-in\\b\
/)\n  {\n    my $save, $name;\n    while($arg_list\
=~/\\-in\\b[^-]+(\\bS[\\w.]+)/)\n      {\n	$name=$\
1;$name=~s/^.//;\n	if ( !-e $name){$save.=\" S$nam\
e \";}\n\n	$arg_list=~s/S$name/ /;\n      }\n    $\
arg_list=~s/\\-in\\b/\\-in $save /;\n  }\n   #2-pr\
epare \n\nif (!($arg_list=~/\\-outorder/))\n  {\n \
   \n    $output_cl .=\" -outorder=$seq\";\n  }\n@\
output_flag=(\"-output\",\"-outfile\", \"-run_name\
\", \"-outorder\"); \nforeach $v1 (@output_flag)\n\
  {\n    ($v2, $arg_list)=&extract_val_from_arg_li\
st($v1,$arg_list, \"SPLICE\",\"unset\");\n    if (\
$v2 ne \"\")\n      {\n\n	if ($v1 eq \"-run_name\"\
){$run_name=$v2;$output_cl .=\" $v1 $v2 \";}\n	els\
if ( $v1 eq \"-outorder\")\n	  {\n	    if ( $v2 eq\
 \"input\"){$v2=$seq;}\n	    $outorder=$v2;$output\
_cl .=\" $v1 $v2 \";\n	  }\n	else\n	  {\n	    $out\
put_cl .=\" $v1 $v2 \";\n	  }\n      }\n }\n\n\n($\
dpa_master_aln, $arg_list)  =&extract_val_from_arg\
_list(\"-dpa_master_aln\",$arg_list, \"SPLICE\", \\
"t_coffee\");\n$dpa_master_aln=~s/\\s//g;\n($nseqT\
, $arg_list)           =&extract_val_from_arg_list\
(\"-dpa_maxnseq\",$arg_list, \"SPLICE\", 30);\n($s\
coreT1, $arg_list)         =&extract_val_from_arg_\
list(\"-dpa_min_score1\",$arg_list, \"SPLICE\", 80\
);\n($scoreT2, $arg_list)         =&extract_val_fr\
om_arg_list(\"-dpa_min_score2\"    ,$arg_list, \"S\
PLICE\", 30);\n($dpa_limit, $arg_list)       =&ext\
ract_val_from_arg_list(\"-dpa_limit\"        ,$arg\
_list, \"SPLICE\", 0);\n($dpa_delta_id, $arg_list)\
    =&extract_val_from_arg_list(\"-dpa_delta_id\" \
       ,$arg_list, \"SPLICE\", 1);\n($dpa_debug, $\
arg_list)       =&extract_val_from_arg_list(\"-dpa\
_debug\"           ,$arg_list, \"SPLICE\", 0);\n\n\
\n$in_seq=$seq0.\" \".$seq1.\" \".$seq2.\" \".$seq\
3;\n$in_prf=(($prf ne $unset)?\"$prf \":\"\");\n&e\
xit_dpa (($in_seq eq \"\" && $in_prf eq \"\")?1:0,\
 \"ERROR: You did not Provide any sequences. Use t\
he -seq flag [FATAL: t_coffee_dpa]\\n\", EXIT_FAIL\
URE);\n\n\nprint STDERR \"\\nSTART DPA COMPUTATION\
\";\n\n\n\nif ($in_seq=~/\\S+/)\n  {\n    \n    pr\
int STDERR \"\\n Step $step: Gather all the sequen\
ces into the tmp file: [$seq]\";$step++;	\n    &my\
_system (\"t_coffee $in_seq -convert -quiet -outpu\
t fasta_seq -outfile=$seq -maxnseq 0\");\n  }\n\ni\
f ( !-e $seq){$seq=\"\";}\n\nif ($in_prf=~/\\S+/)\\
n  {\n    $seq_in_type=\"profile\"; \n    $seq.= $\
in_prf; \n  }\nif ($seq eq \"\"){ &exit_dpa (1, \"\
\\nERROR: No Sequence FOund. Provide Sequences wit\
h the -seq flag [FATAL: t_coffee_dpa]\", EXIT_FAIL\
URE);}\n\n \n\nif ( $run_name)\n  {\n    $suffix=$\
run_name;\n  }\nelsif ($in_seq=~/\\b(S[\\w.]+\\b)/\
)\n  {\n    my $suffix1, $sufffix2;\n    $suffix1=\
$suffix2=$1;\n    $suffix2=~s/^S//;\n    if ( -e $\
suffix1){$suffix=$suffix1;}\n    elsif ( -e $suffi\
x2){$suffix=$suffix2;}\n    else\n      {\n	$suffi\
x=&vtmpnam();	\n      }\n    $suffix=~s/\\.\\w+//;\
\n  }\n\nelse\n  {\n    $suffix=&vtmpnam();\n  }\n\
\n\nif (!$run_name){$output_cl.=\" -run_name $suff\
ix \";}\n\n\n$gl{'Tree'}=&seq2dpa_tree ($seq, \"$s\
uffix.dpadnd\");\n\nprint STDERR \"\\n Step $step:\
 Prepare guide tree: $gl{'Tree'}\";$step++;\n\npri\
nt STDERR \"\\n Step $step: Identify and Align Clo\
sely Related Groups\";$step++;\n%gl=&make_one_pass\
 (0, $scoreT1,\"Align\",%gl);\n\nprint STDERR \"\\\
n Step $step: Make Multiple Group Alignment\";$ste\
p++;\nwhile (!%gl ||$gl{'Ng'}>$nseqT)\n  {\n    %g\
l=&make_one_pass ($nseqT, $scoreT2,\"t_coffee\",%g\
l);\n    if ( $gl{'Newgroups'}==0){$scoreT2--;}   \
 \n  }\nprint STDERR \"\\n Step $step: Make The Fi\
nal Alignment\";$step++;\n\n\n$arg_list .=$output_\
cl;\n\n\n%gl=&tree2group (0,0, %gl);\n$gl{$gl{'0'}\
{'File'}}{'Output'}=\"\";\n$a=0;\n&align_groups (\\
"t_coffee\",'0', $arg_list, \" \", %gl);\n\n\n\nif\
 ( !$dpa_keep_tmpfile){&clean_tmp_file (@tl);}\n\n\
\n\nsub seq2dpa_tree \n  {\n    my $seq=@_[0];\n  \
  my $newtree=@_[1];\n    my $aln=&vtmpnam ();\n\n\
    &my_system (\"t_coffee -special_mode quickaln \
-in $seq -outfile $aln -quiet\");\n    &my_system \
(\"$seq_reformat -in $aln -action +aln2tree +tree2\
dpatree -output newick >$newtree\");\n    return $\
newtree;\n  }	\nsub seq2dpa_tree_old \n  {\n    my\
 $aln=@_[0];\n    my $newtree=@_[1];\n    \n    \n\
    &my_system(\"$seq_reformat -in $aln -action +s\
eq2dpatree -output newick > $newtree\");\n    retu\
rn $newtree;\n  }\nsub aln2dpa_tree \n  {\n    my \
$aln=@_[0];\n    my $newtree=&vtmpnam();\n    \n  \
  &my_system(\"$seq_reformat -in $aln -action +aln\
2tree +tree2dpatree -output newick > $newtree\");\\
n    return $newtree;\n  }\nsub group_file2ngroups\
\n  {\n    my $file=@_[0];\n    my $n;\n    \n    \
open ( F, $file);\n    while (<F>)\n      {\n	$n+=\
/\\>/;\n      }\n    close (F);\n    return $n;\n \
 }\n\nsub make_one_pass\n  {\n    my ($N, $ID,$pg,\
 %gl)=@_;\n    my $a;\n\n    %gl=&tree2group ($N,$\
ID,%gl);\n    if (!$gl{'Newgroups'}){return %gl;}\\
n    else\n      {\n	for ( $a=0; $a< $ng; $a++)\n	\
  {\n	    if ($gl{$gl{$a}{'File'}}{'Ng'}>1){&displ\
ay_group($a, %gl);}\n	    &align_groups ($pg, $a, \
$arg_list, \" -quiet=quiet \", %gl);\n	  }\n	retur\
n %gl;\n      }\n  }\n\nsub tree2group \n  {\n    \
my ($N, $ID, %gl)=@_;\n    my $prefix=&vtmpnam();\\
n    my $group_file=&vtmpnam();\n    my $file;\n  \
  my $oldtree=&vtmpnam();\n    my $n;\n    my $tre\
e;\n\n\n    if ( $gl{'Ng'}==1){return %gl;}\n    $\
tree=$gl{'Tree'}; \n    \n    #1 extract the group\
s\n    &my_system (\"$seq_reformat -in $tree -acti\
on +tree2group $N $ID $prefix > $group_file\");\n \
   $n=group_file2ngroups($group_file);\n    \n    \
\n    $gl{'Newgroups'}=1;\n    if ( $n==$gl{'Ng'})\
\n      {\n	$gl{'Newgroups'}=0;\n	return %gl;\n   \
   }\n    $gl{'Iteration'}++;\n    $gl{'MaxNseq'}=\
$N;$gl{'MinID'}=$ID;\n    $gl{'GroupFile'}=$group_\
file;$gl{'Ng'}=$ng=0;\n    #2 Process the group li\
st into the hash\n    open (F, $group_file);\n    \
while (<F>)\n      {\n	$gl{'File'}.=$_;\n	if (/\\>\
/)\n	  {\n	    $line=$_;\n	    $line=~s/\\>//;\n	 \
   @list=($line=~/(\\S+)/g);\n	    $file=$gl{$ng}{\
'File'}=shift @list;\n	    $gl{$file}{'Output'}=$f\
ile;\n	    \n	    $gl{$file}{'Ng'}=$#list+1;\n	   \
 if ($gl{$file}{'Ng'}>1){ $gl{$file}{'Tlist'}=$gl{\
$file}{'Alist'}=\"(\";}\n	    foreach $l (@list)\n\
	      {\n	\n		$gl{$file}{'List'}.=\" $l \";\n		\n\
		if (!$gl{$l}{'Tlist'})\n		  {\n		    $gl{$l}{'Tl\
ist'}=\"$l\";\n		    $gl{$l}{'Alist'}=\"$l\";\n		 \
   $gl{$l}{'Nseq'}=1;\n		    $gl{$l}{'Ng'}=1;\n		 \
 }\n		$gl{$file}{'Tlist'}.=\"$gl{$l}{'Tlist'},\";\\
n		$gl{$file}{'Alist'}.=\"$gl{$l}{'Tlist'}|\";\n		\
$gl{$file}{'Nseq'}+=$gl{$l}{'Nseq'};\n	      }\n	 \
   \n\n	    chop($gl{$file}{'Tlist'});chop($gl{$fi\
le}{'Alist'});\n	    if ($gl{$file}{'Ng'}>1){$gl{$\
file}{'Tlist'}.=\")\"; $gl{$file}{'Alist'}.=\");\"\
;}\n	    $ng++;\n	  }	\n      }\n    $gl{'Ng'}=$ng\
;\n    close (F);\n    \n    #3 Update the old tre\
e with the new groups\n    $gl{'Tree'}=&vtmpnam();\
\n    &my_system (\"$seq_reformat -in $tree -actio\
n +collapse_tree $group_file -output newick > $gl{\
'Tree'}\");\n    \n    return %gl;\n  }\n\nsub dis\
play_group \n  {\n    my ($g,%gl)=@_;\n    my $f;\\
n    \n    if ( $g==-1)\n      {\n	print STDERR \"\
\\nIteration $gl{'Iteration'} [MaxN=$gl{'MaxNseq'}\
][MinID=$gl{'MinID'}]\";\n      }\n    else\n     \
 {\n\n	$f=$gl{$g}{'File'};\n	$action=($gl{$f}{'Ng'\
}==1 || $gl{'Iteration'}==1)?\"KEEP  \":\"ALIGN \"\
;\n        print STDERR \"\\n\\t[$action][MaxN=$gl\
{'MaxNseq'}][MinID=$gl{'MinID'}][File $f][Nseq=$gl\
{$f}{'Nseq'}][Ngroups=$gl{$f}{'Ng'}][$gl{$f}{'Alis\
t'}]\";\n      }\n  }\n      \n\n\nsub align_group\
s\n  {\n    my ($pg, $g, $arg, $extra_arg,%gl)=@_;\
\n    my $f;\n    my $Output,$Outflag;\n    \n    \
\n    $f=$gl{$g}{'File'};\n    $Output=($gl{$f}{'O\
utput'});\n    \n    if ( $pg eq \"Align\")\n     \
 {\n	if ( !-e $f)\n	  {\n	    $command=\"$seq_refo\
rmat -in $gl{'Seq'}  -action +extract_aln $gl{'Gro\
upFile'}\";\n	    if ($gl{$f}{'Ng'}>1)\n	      {\n\
		&my_system ($command);\n		$command=\"t_coffee -s\
pecial_mode quick_aln  S$f -outfile=$Output -quiet\
\";\n	      }\n	  }\n	else \n	  {$command=\"\";}\n\
      }\n    elsif ( -e $f)\n      {	\n	$Outflag=(\
$Output)?\"-outfile=$Output\":\"\";\n	$command=\"$\
pg -infile $f $Outflag -quiet stdout $arg $extra_a\
rg -maxnseq 0 -convert -quiet stdout\";\n      }\n\
    elsif ( $gl{$f}{'Ng'}==1)\n      {\n	$action=(\
$dpa_debug)?\"cp\":\"mv\";\n	$command=\"$action $g\
l{$f}{'List'} $Output\";\n      }\n    else\n     \
 {\n	$Outflag=($Output)?\"-outfile=$Output\":\"\";\
\n	$command=\"$pg -profile $gl{$f}{'List'} $Outfla\
g $arg $extra_arg -maxnseq 0\";\n      }\n    \n  \
  &my_system ($command);\n    return $outfile;\n  \
}\n    \nsub my_system \n  {\n    my $command=@_[0\
];\n    my $force=@_[1];\n    my $status;\n\n    i\
f ( $dpa_debug) {print STDERR \"\\nCOMMAND: $comma\
nd\";}\n    $status=system ($command);\n\n    if (\
!$force)\n       {\n	 &exit_dpa (($status==1), \"F\
ailed in Command:\\n$command\\n[FATAL: t_coffee_dp\
a]\\n\", EXIT_FAILURE);\n       }\n    \n    retur\
n $status;\n  }\n\nsub vtmpnam\n  {\n    my $prefi\
x=@_[0];\n    my $tmp_file_name;\n\n    $tmp_prefi\
x=($prefix)?$prefix:\"dpa_tmp_file_$$\";\n   \n   \
 $tmp_count++;\n    $tmp_file_name=\"$tmp_prefix\"\
.\"$tmp_count\";\n    $tl[$#tl+1]=$tmp_file_name;\\
n    return $tmp_file_name;\n  }\n\nsub clean_tmp_\
file\n  {\n\n    my $list;\n    my $file;\n    \n \
   if ($dpa_debug){return;}\n    $list=vtmpnam();\\
n    `ls -1 | grep $tmp_prefix>$list`;\n    \n    \
open (F,$list);\n    while ( <F>)\n      {\n	$file\
=$_;\n	chop $file;\n	if ( -e $file){unlink $file;}\
\n      }\n    close (F);\n    unlink $list;\n  }\\
n\n\nsub exit_dpa\n  {\n  my $condition=@_[0];\n  \
my $error_msg=@_[1];\n  my $exit_value=@_[2];\n  i\
f ( $condition)\n    {\n      print \"$error_msg\\\
n\";\n      exit ($exit_value);\n    }\n  else\n  \
  {\n      return;\n    }\n  \n}\nsub extract_val_\
from_arg_list\n  {\n    my $arg=@_[0];\n    my $ar\
g_list=@_[1];\n    my $keep_flag=@_[2];\n    my $d\
efault_value=@_[3];\n    my $val=\"\";\n    \n    \
#protect\n    $arg_list=~s/\\s-/ \\@/g;\n    $arg=\
~s/-/\\@/g;\n    \n    #search\n    if ($arg eq \"\
^\")\n      {\n	$arg_list=~/^([^@]*)/;\n	$val=$1;\\
n      }\n    else\n      {$arg_list=~/$arg ([^@]*\
)/;$val=$1;}\n    \n    #remove trailing spaces\n \
   $val=~s/\\s*$//;\n    \n    #remove the parsed \
sequence if needed\n    if (($val ne \"\") && $kee\
p_flag ne \"KEEP\")\n      {\n	if ( $arg eq \"^\")\
{$arg_list=~s/$val/ /;}\n	else {$arg_list=~s/($arg\
 [^@]*)/ /;}\n      }\n	\n    #unprotect\n    $arg\
_list=~s/\\@/-/g;\n    $arg=~s/\\@/-/g;\n    \n   \
 if (($val eq \"\") && $default_value ne \"unset\"\
){$val=$default_value;}\n    \n    return $val, $a\
rg_list;\n  }\n$program=\"T-COFFEE (dev_shivashaml\
oo@20190813_14:24)\";\n\n","\n$DEBUG=1;\n$dpa_nseq\
=10;\n$dpa_sim=0;\nif (!@ARGV)\n  {\n    `t_coffee\
`;\n    exit (0);\n  }\nforeach $arg (@ARGV)\n  {\\
n    $arg_list.=\" $arg\";\n  }\n$max_nseq=10;\n($\
seq0, $arg_list)=&extract_val_from_arg_list(\"^\",\
$arg_list);\n($seq1, $arg_list)=&extract_val_from_\
arg_list(\"-seq\",$arg_list);\n($seq2, $arg_list)=\
&extract_val_from_arg_list(\"-in\",$arg_list, \"KE\
EP\");\n($seq3, $arg_list)=&extract_val_from_arg_l\
ist(\"-infile\",$arg_list);\n$in_seq=$seq0.\" \".$\
seq1.\" \".$seq2.\" \".$seq3;\n\n$seq=vtmpnam();\n\
`t_coffee $in_seq -convert -output fasta_seq -outf\
ile=$seq`;\n\n\n($dpa_nseq, $arg_list)=&extract_va\
l_from_arg_list(\"-dpa_nseq\",$arg_list);\n($maste\
r_aln, $arg_list)=&extract_val_from_arg_list(\"-ma\
ster_aln\",$arg_list);\n($sim_matrix, $arg_list)=&\
extract_val_from_arg_list(\"-sim_matrix\",$arg_lis\
t);\n($core_seq, $arg_list)=&extract_val_from_arg_\
list(\"-core_seq\",$arg_list);\n($dpa_sim, $arg_li\
st)=&extract_val_from_arg_list(\"-dpa_sim\",$arg_l\
ist);\n($run_name, $arg_list)=&extract_val_from_ar\
g_list(\"-run_name\",$arg_list);\n($output, $arg_l\
ist)=&extract_val_from_arg_list(\"-output\",$arg_l\
ist);\n\n\n\nif (!$sim_mat && !$master_aln)#Comput\
e the fast alignment\n  {\n    $ref_aln=vtmpnam();\
\n    `t_coffee -seq=$seq -very_fast -outfile=$ref\
_aln -quiet`;\n    \n  }\n\nif (!$sim_mat)\n  {\n \
   $sim_mat=vtmpnam();\n    `seq_reformat -in $ref\
_aln -output sim > $sim_mat`;\n  }\n\nif ( !$core_\
seq)\n  {\n    $core_seq=vtmpnam();\n    `seq_refo\
rmat -in $ref_aln -action +trimTC N$max_nseq -outp\
ut fasta_seq > $core_seq`;\n  }\n@core_name=`seq_r\
eformat -in $core_seq -output name `; \n\n@tot_nam\
e=`seq_reformat -in $seq -output name `;\n\nforeac\
h $s (@core_name){$s=~s/\\s//g;$hcore{$s}=1;}\nfor\
each $s (@tot_name){$s=~s/\\s//g;}\nprint STDERR \\
"T-Coffee_dpa:\\n\";\nprint STDERR \"\\tTOTAL  SEQ\
: @tot_name\\n\";\nprint STDERR \"\\tCHOSEN SEQ: @\
core_name\\n\";\n\n\n\nopen (F, $sim_mat);\nwhile \
( <F>)\n  {\n    @l=($_=~/(\\b[\\S]+\\b)/g);\n    \
if (($l[0] eq \"TOP\" || $l[0] eq \"BOT\"))\n     \
 {\n	$s1=$l[1];$s2=$l[2];$v=$l[3];\n	if ($hcore{$s\
1} && !$hcore{$s2})\n	  {\n	    if (!$hseq{$s2}{\"\
sim\"} || $v>$hseq{$s2}{\"sim\"})\n	      {\n		$hs\
eq{$s2}{\"sim\"}=$v;$hseq{$s2}{\"seq\"}=$s1;\n	   \
   }\n	  }\n      }\n  }\nclose (F);\nforeach $s (\
@tot_name)\n  {\n\n    if ( !$hseq{$s}{\"seq\"}){;\
}\n    else\n      {\n	$s2=$hseq{$s}{\"seq\"};\n	$\
v=$hseq{$s}{\"sim\"};\n		\n	if ($v>$dpa_sim)\n	  {\
\n	    $hseq{$s}{'used'}=1;\n	    $seq_list{$s2}{$\
seq_list{$s2}{'nseq'}++}=$s;\n	  }\n      }\n  }\n\
foreach $s (@core_name){$seq_list{$s}{$seq_list{$s\
}{'nseq'}++}=$s;$hseq{$s}{'used'}=1;}\nforeach $s \
(@tot_name){if (!$hseq{$s}{'used'}){$seq_list{'unu\
sed'}{$seq_list{'unused'}{'nseq'}++}=$s;}}\n\n\n$n\
=0;\nforeach $s (@core_name)\n  {\n    $ng++;\n   \
 $n=$seq_list{$s}{'nseq'};\n    for (@g_list=(), $\
a=0; $a<$n; $a++){@g_list=(@g_list,$seq_list{$s}{$\
a});}\n\n    $g_seq=vtmpnam();\n    $g_aln=vtmpnam\
();\n    \n    print STDERR \"Group $ng: $#g_list \
Seq: @g_list: \";\n    \n    \n    `seq_reformat -\
in $seq -action +lower +keep_name +extract_seq  @g\
_list -output fasta_seq > $g_seq`;\n    \n    \n  \
  if ( $#g_list==0)\n      {\n	print STDERR \"[No \
aln]\\n\";\n	$g_aln=$g_seq;\n      }\n    elsif ($\
#g_list<$max_nseq) \n      {\n	print STDERR \"[t_c\
offee]\\n\";\n	`t_coffee $g_seq -outfile=$g_aln -q\
uiet $arg_list`;\n      }\n    else\n      {\n	pri\
nt STDERR \"[t_coffee_dpa]\\n\";\n	`t_coffee_dpa2 \
$g_seq -outfile=$g_aln $arg_list -sim_matrix $sim_\
matrix -dpa_nseq $dpa_nseq`;\n      }\n    @profil\
e_list=(@profile_list, $g_aln);\n  }\n\n\nprint \"\
UNUSED $seq_list{'unused'}{'nseq'}\";\n\nif ($seq_\
list{'unused'}{'nseq'})\n    {\n      $prf=vtmpnam\
();\n      \n      `t_coffee -profile @profile_lis\
t $arg_list -outfile=$prf -quiet`;\n      $n=$seq_\
list{\"unused\"}{'nseq'};\n      $new_seq=vtmpnam(\
);\n      $new_prf=vtmpnam();\n      for ($a=0; $a\
<$n-1; $a++)\n	{\n	  $s=$seq_list{\"unused\"}{$a};\
\n	  print STDERR \"\\nADD Sequence $s\";\n	  \n	 \
 `seq_reformat -in $seq -action +lower +keep_name \
+extract_seq $s  -output fasta_seq > $new_seq`;\n	\
  `t_coffee -profile $prf $new_seq $arg_list -outf\
ile=$new_prf`;\n	  `cp $new_prf $prf`;\n	}\n      \
$s=$seq_list{\"unused\"}{$a};\n      `seq_reformat\
 -in $seq -action +lower +keep_name +extract_seq $\
s  -output fasta_seq > $new_seq`;\n      @profile_\
list=($prf, $new_seq);\n    }\n    \n      \nif ($\
run_name){$arg_list.=\" -run_name $run_name\";}\ne\
lse \n  {\n    $in_seq=~/([\\w-]+)/;\n    $arg_lis\
t.=\" -run_name $1\";\n  }\nif ( $output){$arg_lis\
t.=\" -output $output \";}\n\n`t_coffee -profile @\
profile_list $arg_list`;\n\n\n&clean (@tmp_file_li\
st);\n\n\nsub vtmpnam\n  {\n    my $tmp_file_name;\
\n    $tmp_name_counter++;\n    $tmp_file_name=\"t\
mp_file_$tmp_name_counter\\_Pid$$\";\n    $tmp_fil\
e_list[$ntmp_file++]=$tmp_file_name;\n    return $\
tmp_file_name;\n  }\nsub clean\n  {\n  my @fl=@_;\\
n  my $file;\n  return;\n\n  foreach $file ( @fl)\\
n    {\n      if ( -e $file){unlink($file);}\n    \
}\n}\nsub extract_val_from_arg_list\n  {\n    my $\
arg=@_[0];\n    my $arg_list=@_[1];\n    my $keep_\
flag=@_[2];\n    #protect\n    $arg_list=~s/\\s-/ \
\\@/g;\n    $arg=~s/-/\\@/g;\n    \n    #search\n \
   if ($arg eq \"^\")\n      {\n	$arg_list=~/^([^@\
]*)/;\n	$val=$1;\n      }\n    else\n      {$arg_l\
ist=~/$arg ([^@]*)/;$val=$1;}\n    \n    #remove t\
he parsed sequence if needed\n    if ($val && $kee\
p_flag ne \"KEEP\")\n      {\n	if ( $arg eq \"^\")\
{$arg_list=~s/$val/ /;}\n	else {$arg_list=~s/($arg\
 [^@]*)/ /;}\n      }\n	\n    #unprotect\n    $arg\
_list=~s/\\@/-/g;\n    $arg=~s/\\@/-/g;\n    \n   \
 return $val, $arg_list;\n  }\n\n","use Env;\nuse \
FileHandle;\nuse Cwd;\nuse File::Path;\nuse Sys::H\
ostname;\n\n\nour $PIDCHILD;\nour $ERROR_DONE;\nou\
r @TMPFILE_LIST;\nour $EXIT_FAILURE=1;\nour $EXIT_\
SUCCESS=0;\n\nour $REFDIR=getcwd;\nour $EXIT_SUCCE\
SS=0;\nour $EXIT_FAILURE=1;\n\nour $PROGRAM=\"tc_g\
eneric_method.pl\";\nour $CL=$PROGRAM;\n\nour $CLE\
AN_EXIT_STARTED;\nour $debug_lock=$ENV{\"DEBUG_LOC\
K\"};\nour $debug_generic_method=$ENV{\"DEBUG_GENE\
RIC_METHOD\"};\nour $LOCKDIR=$ENV{\"LOCKDIR_4_TCOF\
FEE\"};\nif (!$LOCKDIR){$LOCKDIR=getcwd();}\nour $\
ERRORDIR=$ENV{\"ERRORDIR_4_TCOFFEE\"};\nour $ERROR\
FILE=$ENV{\"ERRORFILE_4_TCOFFEE\"};\n&set_lock ($$\
);\nif (isshellpid(getppid())){lock4tc(getppid(), \
\"LLOCK\", \"LSET\", \"$$\\n\");}\nour %RECODE;\no\
ur $RECODE_N;\n\n\n\n\nour $BLAST_MAX_NRUNS=2;\nou\
r $COMMAND;\nour $PIDCHILD;\n\n$REF_EMAIL=\"\";\n$\
tmp_dir=\"\";\n$init_dir=\"\";\n\n\n$test=0;\nif (\
$test==1)\n  {\n    $SERVER=\"NCBI\";\n    $query=\
$ARGV[0];\n    $hitf=$ARGV[1];\n    %s=read_fasta_\
seq($query);\n    @sl=keys(%s);\n    &blast_xml2pr\
ofile (\"xx\", $s{$sl[0]}{seq},$maxid,$minid,$minc\
ov, $hitf);\n    myexit ($EXIT_FAILURE);\n  }\n\nf\
oreach $v(@ARGV){$cl.=\"$v \";}\n$COMMAND=$cl;\n($\
mode)=&my_get_opt ( $cl, \"-mode=\",1,0);\n\n($A)=\
(&my_get_opt ( $cl, \"-name1=\",0,0));\n($B)=(&my_\
get_opt ( $cl, \"-name2=\",0,0));\n($TMPDIR)=(&my_\
get_opt ( $cl, \"-tmpdir=\",0,0));\n($CACHE)=(&my_\
get_opt ( $cl, \"-cache=\",0,0));\n($SERVER)=((&my\
_get_opt ( $cl, \"-server=\",0,0)));\n($EMAIL)=((&\
my_get_opt ( $cl, \"-email=\",0,0)));\n\nif (!$A){\
$A=\"A\";}\nif (!$B){$B=\"B\";}\n\n\nif (!$TMPDIR)\
\n  {\n    $HOME=$ENV{HOME};\n    if ($ENV{TMP_4_T\
COFFEE}){$TMPDIR=$ENV{TMP_4_TCOFFEE};}\n    else{$\
TMPDIR=\"$HOME/.t_coffee/tmp/\";}\n  }\nif ( ! -d \
$TMPDIR)\n  {\n    mkdir $TMPDIR;\n  }\nif ( ! -d \
$TMPDIR)\n  {\n    print \"ERROR: Could not create\
 temporary dir: $TMPDIR\\n\";\n    myexit ($EXIT_F\
AILURE);\n  }\n\n$EMAIL=~s/XEMAILX/\\@/g;\nif (!$E\
MAIL)\n  {\n    if ($ENV{EMAIL_4_TCOFFEE}){$EMAIL=\
$ENV{EMAIL_4_TCOFFEE};}\n    elsif ($ENV{EMAIL}){$\
EMAIL=$ENV{EMAIL};}\n    else {$EMAIL=$REF_EMAIL;}\
\n  }\n\n($maxid,$minid,$mincov,$trim)=(&my_get_op\
t ( $cl, \"-maxid=\",0,0, \"-minid=\",0,0,\"-minco\
v=\",0,0, \"-trim=\",0,0));\nif (!$cl=~/\\-maxid\\\
=/){$maxid=95;}\nif (!$cl=~/\\-minid\\=/){$minid=3\
5;}\nif (!$cl=~/\\-mincov\\=/){$mincov=80;}\nif (!\
$cl=~/\\-trim\\=/){$trim;}\n\n\n\n\nif ($mode eq \\
"seq_msa\")\n  {\n    &seq2msa($mode,&my_get_opt (\
 $cl, \"-infile=\",1,1, \"-method=\",1,2, \"-param\
=\",0,0,\"-outfile=\",1,0, \"-database=\",0,0));\n\
  }\nelsif ( $mode eq \"tblastx_msa\")\n  {\n    &\
seq2tblastx_lib ($mode,&my_get_opt ( $cl, \"-infil\
e=\",1,1, \"-outfile=\",1,0));\n  }\nelsif ( $mode\
 eq \"tblastpx_msa\")\n  {\n    &seq2tblastpx_lib \
($mode,&my_get_opt ( $cl, \"-infile=\",1,1, \"-out\
file=\",1,0));\n  }\nelsif ( $mode eq \"thread_pai\
r\")\n  {\n    &seq2thread_pair($mode,&my_get_opt \
( $cl, \"-infile=\",1,1, \"-pdbfile1=\",1,1, \"-me\
thod=\",1,2,\"-param=\",0,0, \"-outfile=\",1,0, ))\
;\n  }\nelsif ( $mode eq \"pdbid_pair\")\n  {\n   \
 &seq2pdbid_pair($mode,&my_get_opt ( $cl, \"-pdbfi\
le1=\",1,0, \"-pdbfile2=\",1,0, \"-method=\",1,2,\\
"-param=\",0,0, \"-outfile=\",1,0, ));\n  }\nelsif\
 ( $mode eq \"pdb_pair\")\n  {\n    &seq2pdb_pair(\
$mode,&my_get_opt ( $cl, \"-pdbfile1=\",1,1, \"-pd\
bfile2=\",1,1, \"-method=\",1,2,\"-param=\",0,0, \\
"-outfile=\",1,0, ));\n  }\nelsif ( $mode eq \"rna\
pdb_pair\")\n{\n    &seq2rnapdb_pair($mode,&my_get\
_opt ( $cl, \"-pdbfile1=\",1,1, \"-pdbfile2=\",1,1\
, \"-method=\",1,2,\"-param=\",0,0, \"-outfile=\",\
1,0, ));\n}\nelsif ( $mode eq \"profile_pair\")\n \
 {\n     &seq2profile_pair($mode,&my_get_opt ( $cl\
, \"-profile1=\",1,1, \"-profile2=\",1,1, \"-metho\
d=\",1,2,\"-param=\",0,0, \"-outfile=\",1,0 ));\n \
 }\nelsif ($mode eq \"pdb_template_test\")\n  {\n \
   &blast2pdb_template_test ($mode,&my_get_opt ( $\
cl, \"-infile=\",1,1));\n\n  }\nelsif ($mode eq \"\
psi_template_test\")\n  {\n    &psiblast2profile_t\
emplate_test ($mode,&my_get_opt ( $cl, \"-seq=\",1\
,1,\"-blast=\",1,1));\n\n  }\n\nelsif ( $mode eq \\
"pdb_template\")\n  {\n    &blast2pdb_template ($m\
ode,&my_get_opt ( $cl, \"-infile=\",1,1, \"-databa\
se=\",1,0, \"-method=\",1,0, \"-outfile=\",1,0,\"-\
pdb_type=\",1,0));\n  }\n\nelsif ( $mode eq \"prof\
ile_template\")\n  {\n\n    &seq2profile_template \
($mode,&my_get_opt ( $cl, \"-infile=\",1,1, \"-dat\
abase=\",1,0, \"-method=\",1,0, \"-outfile=\",1,0)\
);\n  }\nelsif ( $mode eq \"psiprofile_template\")\
\n  {\n    &seq2profile_template ($mode,&my_get_op\
t ( $cl, \"-infile=\",1,1, \"-database=\",1,0, \"-\
method=\",1,0, \"-outfile=\",1,0));\n  }\nelsif ( \
$mode eq \"RNA_template\")\n  {\n    &seq2RNA_temp\
late ($mode,&my_get_opt ( $cl, \"-infile=\",1,1,\"\
-pdbfile=\",1,1,\"-outfile=\",1,0));\n  }\nelsif (\
 $mode eq \"tm_template\")\n  {\n    &seq2tm_templ\
ate ($mode,&my_get_opt ( $cl, \"-infile=\",1,1,\"-\
arch=\",1,1,\"-psv=\",1,1, \"-outfile=\",1,0));\n \
 }\nelsif ( $mode eq \"psitm_template\")\n  {\n   \
 &seq2tm_template ($mode,&my_get_opt ( $cl, \"-inf\
ile=\",1,1, \"-arch=\",1,1,\"-psv=\",1,1, \"-outfi\
le=\",1,0,\"-database=\",1,0));\n  }\nelsif ( $mod\
e eq \"ssp_template\")\n  {\n    &seq2ssp_template\
 ($mode,&my_get_opt ( $cl, \"-infile=\",1,1,\"-seq\
=\",1,1,\"-obs=\",1,1, \"-outfile=\",1,0));\n  }\n\
elsif ( $mode eq \"psissp_template\")\n  {\n    &s\
eq2ssp_template ($mode,&my_get_opt ( $cl, \"-infil\
e=\",1,1,\"-seq=\",1,1,\"-obs=\",1,1, \"-outfile=\\
",1,0));\n  }\n\n\n\nelse\n  {\n    myexit(flush_e\
rror( \"$mode is an unknown mode of tc_generic_met\
hod.pl\"));\n  }\nmyexit ($EXIT_SUCCESS);\n\n\nsub\
 seq2ssp_template\n  {\n  my ($mode, $infile,$gor_\
seq,$gor_obs,$outfile)=@_;\n  my %s, %h;\n  my $re\
sult;\n  my (@profiles);\n  &set_temporary_dir (\"\
set\",$infile,\"seq.pep\");\n  %s=read_fasta_seq (\
\"seq.pep\");\n\n\n  open (R, \">result.aln\");\n\\
n  #print stdout \"\\n\";\n  foreach $seq (keys(%s\
))\n    {\n\n      open (F, \">seqfile\");\n      \
$s{$seq}{seq}=uc$s{$seq}{seq};\n      print (F \">\
$s{$seq}{name}\\n$s{$seq}{seq}\\n\");\n      close\
 (F);\n      $lib_name=\"$s{$seq}{name}.ssp\";\n  \
    $lib_name=&clean_file_name ($lib_name);\n\n   \
   if ($mode eq \"ssp_template\"){&seq2gor_predict\
ion ($s{$seq}{name},$s{$seq}{seq}, \"seqfile\", $l\
ib_name,$gor_seq, $gor_obs);}\n      elsif ($mode \
eq \"psissp_template\")\n	{\n	  &seq2msa_gor_predi\
ction ($s{$seq}{name},$s{$seq}{seq},\"seqfile\", $\
lib_name,$gor_seq, $gor_obs);\n	}\n\n      if ( !-\
e $lib_name)\n	{\n	  myexit(flush_error(\"GORIV fa\
iled to compute the secondary structure of $s{$seq\
}{name}\"));\n	  myexit ($EXIT_FAILURE);\n	}\n    \
  else\n	{\n	  print stdout \"!\\tProcess: >$s{$se\
q}{name} _E_ $lib_name \\n\";\n	  print R \">$s{$s\
eq}{name} _E_ $lib_name\\n\";\n	}\n      unshift (\
@profiles, $lib_name);\n    }\n  close (R);\n  &se\
t_temporary_dir (\"unset\",$mode, $method,\"result\
.aln\",$outfile, @profiles);\n}\n\nsub seq2tm_temp\
late\n  {\n  my ($mode,$infile,$arch,$psv,$outfile\
,$db)=@_;\n  my %s, %h;\n  my $result;\n  my (@pro\
files);\n  &set_temporary_dir (\"set\",$infile,\"s\
eq.pep\");\n  %s=read_fasta_seq (\"seq.pep\");\n\n\
\n  open (R, \">result.aln\");\n\n  #print stdout \
\"\\n\";\n  foreach $seq (keys(%s))\n    {\n      \
open (F, \">seqfile\");\n      print (F \">$s{$seq\
}{name}\\n$s{$seq}{seq}\\n\");\n      close (F);\n\
      $lib_name=\"$s{$seq}{name}.tmp\";\n      $li\
b_name=&clean_file_name ($lib_name);\n\n      if (\
$mode eq \"tm_template\")\n	{\n	  &safe_system (\"\
t_coffee -other_pg fasta_seq2hmmtop_fasta.pl -in=s\
eqfile -out=$lib_name -arch=$arch -psv=$psv\");\n	\
}\n      elsif ( $mode eq \"psitm_template\")\n	{\\
n	  &seq2msa_tm_prediction ($s{$seq}{name},$s{$seq\
}{seq}, $db, \"seqfile\", $lib_name,$arch, $psv);\\
n	}\n      if ( !-e $lib_name)\n	{\n	  myexit(flus\
h_error(\"hmmtop failed to compute the secondary s\
tructure of $s{$seq}{name}\"));\n	  myexit ($EXIT_\
FAILURE);\n	}\n      else\n	{\n	  print stdout \"!\
\\tProcess: >$s{$seq}{name} _T_ $lib_name\\n\";\n	\
  print R \">$s{$seq}{name} _T_ $lib_name\\n\";\n	\
}\n      unshift (@profiles, $lib_name);\n    }\n \
 close (R);\n  &set_temporary_dir (\"unset\",$mode\
, $method,\"result.aln\",$outfile, @profiles);\n}\\
n\n\n\nsub seq2RNA_template\n  {\n    \n    my ($m\
ode, $infile, $pdbfile, $outfile)=@_;\n    my %s, \
%h ;\n    my $result;\n    my (@profiles);\n    my\
 ($seq_mode, $pdb_mode, $pwd);\n    \n    #use $se\
q_mode to estimate the template of sequences WITHO\
UT a PDB\n    #use $pdb_mode to estimate the templ\
ate of sequences WITH    a PDB\n\n    $seq_mode=$E\
NV{\"SEQ2TEMPLATE4_F_\"};\n    $pdb_mode=$ENV{\"PD\
B2TEMPLATE4_F_\"};\n    \n    if (!$pdb_mode){$pdb\
_mode=\"find_pair-p\";}\n    if (!$seq_mode){$seq_\
mode=\"RNAplfold\";}\n    \n    my $cwd = cwd();\n\
    &set_temporary_dir (\"set\",$infile,\"seq.pep\\
");\n    %s=read_fasta_seq (\"seq.pep\");\n    %pd\
b_template_h = &read_template_file($pdbfile);\n   \
 my $pdb_chain;\n    \n       \n    open (R, \">re\
sult.aln\");\n    #print stdout \"\\n\";\n    fore\
ach $seq (keys(%s))\n      {\n	\n	open (F, \">seqf\
ile\");\n	print (F \">$s{$seq}{name}\\n$s{$seq}{se\
q}\\n\");\n	close (F);\n	$pdb_chain = $pdb_templat\
e_h{$seq};\n	$lib_name=\"$s{$seq}{name}.rfold\";\n\
	$lib_name=&clean_file_name ($lib_name);\n	if ($pd\
b_template_h{$seq} eq \"\")\n	  {\n	    if    ($se\
q_mode eq \"RNAplfold\"){RNAplfold2lib (\"seqfile\\
", \"$lib_name\");}\n	    elsif ($seq_mode eq \"no\
\"){$lib_name=0;}\n	    else\n	      {\n		myexit(a\
dd_error (EXIT_FAILURE,$$,$$,getppid(), \"seq2RNA_\
template failure::method $seq_mode not available f\
or sequences without PDB structures\"));\n	      }\
\n	  }\n	elsif ($pdb_template_h{$seq} ne \"\")\n	 \
 {\n	    my $pdbf;\n	    if (-e \"$cwd/$pdb_chain\\
"){$pdbf=\"$cwd/$pdb_chain\";}\n	    else {$pdbf=\\
"$CACHE$pdb_chain\";}\n	    \n\n	    if($pdb_mode \
eq \"x3dna-ssr\")\n	      {\n		x3dnassr2lib (\"seq\
file\", \"$pdbf\", \"$lib_name\");\n	      }\n	   \
 elsif ($pdb_mode eq \"find_pair-p\")\n	      {\n	\
	x3dna_find_pair2lib (\"seqfile\", \"$pdbf\", \"$l\
ib_name\", \"find_pair -p\");\n	      }\n	    elsi\
f ($pdb_mode eq \"find_pair\")\n	      {\n		x3dna_\
find_pair2lib (\"seqfile\", \"$pdbf\", \"$lib_name\
\", \"find_pair\");\n	      }\n	    elsif ($pdb_mo\
de eq \"RNAplfold\")\n	      {\n		RNAplfold2lib (\\
"seqfile\", \"$lib_name\");\n	      }\n	    elsif \
($pdb_mode eq \"no\"){$lib_name=0;}\n	    else\n	 \
     {\n		myexit(add_error (EXIT_FAILURE,$$,$$,get\
ppid(), \"seq2RNA_template failure::Could not find\
 method $pdb_mode\"));\n	      }\n	  }\n	if ($lib_\
name)\n	  {\n	    print stdout \"!\\tProcess: >$s{\
$seq}{name} _F_ $lib_name\\n\";\n	    print R \">$\
s{$seq}{name} _F_ $lib_name\\n\";\n	    unshift (@\
profiles, $lib_name);\n	  }\n      }\n    close (R\
);\n    &set_temporary_dir (\"unset\",$mode, $meth\
od,\"result.aln\",$outfile, @profiles);\n  }\n\n\n\
\nsub psiblast2profile_template_test\n  {\n  my ($\
mode, $seq,$blast)=@_;\n  my %s, %h, ;\n  my ($res\
ult,$psiblast_output,$profile_name,@profiles);\n  \
my $trim=0;\n  my $maxid=100;\n  my $minid=0;\n  m\
y $mincov=0;\n  my $maxcov=100;\n\n  %s=read_fasta\
_seq ($seq);\n  open (R, \">result.aln\");\n\n  #p\
rint stdout \"\\n\";\n  foreach $seq (keys(%s))\n \
   {\n\n      open (F, \">seqfile\");\n      print\
 (F \">$A\\n$s{$seq}{seq}\\n\");\n      close (F);\
\n      $psiblast_output=$blast;\n      if ( -e $p\
siblast_output)\n	{\n	  %profile=blast_xml2profile\
($s{$seq}{name}, $s{$seq}{seq},$maxid, $minid,$min\
cov,$psiblast_output);\n\n\n\n	  $profile_name=\"$\
s{$seq}{name}.prf\";\n	  $profile_name=&clean_file\
_name ($profile_name);\n	  unshift (@profiles, $pr\
ofile_name);\n	  output_profile ($profile_name, \\\
%profile, $trim);\n	  print stdout \"!\\tProcess: \
>$s{$seq}{name} _R_ $profile_name [$profile{n} Seq\
.] [$SERVER/blast/$db][$CACHE_STATUS]\\n\";\n	  pr\
int R \">$s{$seq}{name} _R_ $profile_name\\n\";\n	\
}\n    }\n  close (R);\n\n  die;\n}\nsub seq2profi\
le_template\n    {\n      my ($mode, $infile, $db,\
 $method, $outfile)=@_;\n      if    ($method eq \\
"psiblast\"){return psiblast2profile_template ($mo\
de, $infile, $db, $method, $outfile);}\n      elsi\
f ($method eq \"blastp\")   {return psiblast2profi\
le_template ($mode, $infile, $db, $method, $outfil\
e);}\n      elsif ($method eq \"hh\")      {return\
 hh2profile_template ($mode, $infile, $db, $method\
, $outfile);}\n    }\n\nsub psiblast2profile_templ\
ate\n  {\n  my ($mode, $infile, $db, $method, $out\
file)=@_;\n  my %s, %h, ;\n  my ($result,$psiblast\
_output,$profile_name,@profiles);\n  &set_temporar\
y_dir (\"set\",$infile,\"seq.pep\");\n  %s=read_fa\
sta_seq (\"seq.pep\");\n  open (R, \">result.aln\"\
);\n\n  #print stdout \"\\n\";\n  foreach $seq (ke\
ys(%s))\n    {\n      open (F, \">seqfile\");\n   \
   print (F \">$A\\n$s{$seq}{seq}\\n\");\n      cl\
ose (F);\n      $psiblast_output=&run_blast ($s{$s\
eq}{name},$method, $db, \"seqfile\",\"outfile\");\\
n\n      if ( -e $psiblast_output)\n	{\n	  my %pro\
file=blast_xml2profile($s{$seq}{name}, $s{$seq}{se\
q},$maxid, $minid,$mincov,$psiblast_output);\n	  u\
nlink ($psiblast_output);\n	  \n	  $profile_name=\\
"$s{$seq}{name}.prf\";\n	  $profile_name=&clean_fi\
le_name ($profile_name);\n	  unshift (@profiles, $\
profile_name);\n	  output_profile ($profile_name, \
\\%profile, $trim);\n	  \n	  print stdout \"!\\tPr\
ocess: >$s{$seq}{name} _R_ $profile_name [$profile\
{n} Seq.] [$SERVER/blast/$db][$CACHE_STATUS]\\n\";\
\n	  print R \">$s{$seq}{name} _R_ $profile_name\\\
n\";\n	  \n	  \n	}\n      \n    }\n  close (R);\n \
 \n  \n\n  &set_temporary_dir (\"unset\",$mode, $m\
ethod,\"result.aln\",$outfile, @profiles);\n}\n\ns\
ub hh2profile_template\n  {\n\n  #for each sequenc\
e, build a profile, in FASTA, with ungapped querry\
 on top  \n  my ($mode, $infile, $db, $method, $ou\
tfile)=@_;\n  my %s, %h, ;\n  my ($result,$psiblas\
t_output,$profile_name,@profiles);\n  &set_tempora\
ry_dir (\"set\",$infile,\"seq.pep\");\n  %s=read_f\
asta_seq (\"seq.pep\");\n  open (R, \">result.aln\\
");\n  \n  my $hh=$ENV{\"HHSEARCH_4_TCOFFEE\"};\n \
 if (!$hh)\n    {\n      print \"ERROR: HHSEARCH_4\
_TCOFFEE is not set\\n\";\n      myexit ($EXIT_FAI\
LURE);\n    }\n  \n  #print stdout \"\\n\";\n  for\
each $seq (keys(%s))\n    {\n      my ($profile_na\
me, $nseq);\n      open (F, \">seqfile\");\n      \
print (F \">$A\\n$s{$seq}{seq}\\n\");\n      close\
 (F);\n      \n      #This function should input a\
 querry and a database and return as output a fast\
a MSA with quesry on top\n      $profile_name=\"$s\
{$seq}{name}.prf\";\n      $profile_name=&clean_fi\
le_name ($profile_name);\n      unshift (@profiles\
, $profile_name);\n      \n      \n      safe_syst\
em  (\"$hh -name=$s{$seq}{name} -method=search -db\
=$db -seq=seqfile -outfile=$profile_name\");\n    \
  if (-e $profile_name){$nseq=fasta2nseq($profile_\
name);}\n      \n      print stdout \"!\\tProcess:\
 >$s{$seq}{name} _R_ $profile_name [$nseq Seq.] [$\
method/$db][$CACHE_STATUS]\\n\";\n      print R \"\
>$s{$seq}{name} _R_ $profile_name\\n\";\n    }\n  \
close (R);\n  &set_temporary_dir (\"unset\",$mode,\
 $method,\"result.aln\",$outfile, @profiles);\n}\n\
\nsub blast2pdb_template_test\n    {\n      my ($m\
ode,$infile)=@_;\n      my ($maxid,$minid,$mincov)\
;\n      $maxid=100;\n      $minid=0;\n      $minc\
ov=0;\n\n      print \"$infile\\n\";\n\n      %p=b\
last_xml2profile($s{$seq}{name}, $s{$seq}{seq},$ma\
xid, $minid,$mincov,$infile);\n      $c=1;\n      \
print stdout \"!\\tProcess: >$s{$seq}{name} [$SERV\
ER/blast/$db][$CACHE_STATUS]\\n\";\n      while (!\
$found && $c<$p{n})\n	{\n	  $pdbid=&id2pdbid($p{$c\
}{identifyer});\n	  if ( length ($pdbid)>5){$pdbid\
=id2pdbid($p{$c}{definition});}\n\n	  if ( length \
($pdbid)>5)\n	    {\n	      myexit(add_error (EXIT\
_FAILURE,$$,$$,getppid(), \"BLAST_FAILURE::Could N\
ot Parse PDBID ($p{$c}{identifyer},$p{$c}{definiti\
on})\"));\n	    }\n\n\n	  if (!&pdb_is_released($p\
dbid))\n	    {\n	      print stdout \"\\t\\t**$pdb\
id [WARNIG: PDB NOT RELEASED or WITHDRAWN]\\n\";\n\
	      $c++;\n	    }\n	  elsif (!&pdb_has_right_ty\
pe ($pdbid,$type))\n	    {\n	      my $ptype=&pdb2\
type ($pdbid);\n	      my $etype=&type2etype($type\
);\n\n	      print stdout \"\\t\\t**$pdbid [$ptype\
 cannot be used (expected: $etype)]\\n\";\n	      \
$c++;\n	    }\n	  else\n	    {\n	      $found=1;\n\
	    }\n	}\n\n      if ($found)\n	{\n	  print stdo\
ut \"\\t\\t >$s{$seq}{name} _P_ $pdbid\\n\";\n	}\n\
      else\n	{\n	  print stdout \"\\t\\t >$s{$seq}\
{name} No Template Selected\\n\";\n	}\n      die;\\
n    }\nsub blast2pdb_template\n  {\n  my ($mode, \
$infile, $db, $method, $outfile,$type)=@_;\n  my %\
s, %h, ;\n  my ($result,$blast_output);\n  &set_te\
mporary_dir (\"set\",$infile,\"seq.pep\");\n  %s=r\
ead_fasta_seq (\"seq.pep\");\n  open (R, \">result\
.aln\");\n\n\n  #print stdout \"\\n\";\n  foreach \
$seq (keys(%s))\n    {\n      my $c;\n      my $fo\
und;\n\n      open (F, \">seqfile\");\n      print\
 (F \">$A\\n$s{$seq}{seq}\\n\");\n      close (F);\
\n\n      $blast_output=&run_blast ($s{$seq}{name}\
,$method, $db, \"seqfile\",\"outfile\");\n\n      \
%p=blast_xml2profile($s{$seq}{name}, $s{$seq}{seq}\
,$maxid, $minid,$mincov,$blast_output);\n      unl\
ink ($blast_output);\n\n      $c=1;\n      print s\
tdout \"!\\tProcess: >$s{$seq}{name} [$SERVER/blas\
t/$db][$CACHE_STATUS]\\n\";\n      while (!$found \
&& $c<$p{n})\n	{\n	  $pdbid=&id2pdbid($p{$c}{ident\
ifyer});\n	  if ( length ($pdbid)>5){$pdbid=id2pdb\
id($p{$c}{definition});}\n\n	  if ( length ($pdbid\
)>5)\n	    {\n	      myexit(add_error (EXIT_FAILUR\
E,$$,$$,getppid(), \"BLAST_FAILURE::Could Not Pars\
e PDBID ($p{$c}{identifyer},$p{$c}{definition})\")\
);\n	    }\n\n\n	  if (!&pdb_is_released($pdbid))\\
n	    {\n	      print stdout \"\\t\\t**$pdbid [PDB\
 NOT RELEASED or WITHDRAWN]\\n\";\n	      $c++;\n	\
    }\n	  elsif (!&pdb_has_right_type ($pdbid,$typ\
e))\n	    {\n	      my $ptype=&pdb2type ($pdbid);\\
n	      my $etype=&type2etype($type);\n\n	      pr\
int stdout \"\\t\\t**$pdbid [$ptype cannot be used\
 (expected: $etype)]\\n\";\n	      $c++;\n	    }\n\
	  else\n	    {\n	      $found=1;\n	    }\n	}\n\n \
     if ($found)\n	{\n	  print R \">$s{$seq}{name}\
 _P_ $pdbid\\n\";\n	  print stdout \"\\t\\t >$s{$s\
eq}{name} _P_ $pdbid\\n\";\n	}\n      else\n	{\n	 \
 print R \">$s{$seq}{name}\\n\";\n	  print stdout \
\"\\t\\t >$s{$seq}{name} No Template Selected\\n\"\
;\n	}\n    }\n  close (R);\n  &set_temporary_dir (\
\"unset\",$mode, $method,\"result.aln\",$outfile);\
\n}\nsub type2etype\n  {\n    my $type=shift;\n   \
 my $etype;\n\n    if ( $type=~/n/){$etype.=\"NMR \
\";}\n    if ( $type=~/d/){$etype.=\"diffraction \\
";}\n    if ( $type=~/m/){$etype.=\"model \";}\n  \
  return $etype;\n  }\nsub pdb2type\n  {\n     my \
$pdb=shift;\n     my $f=vtmpnam();\n\n     my $val\
ue= &safe_system (\"t_coffee -other_pg extract_fro\
m_pdb -model_type $pdb > $f\");\n     my $r=&file2\
string ($f);\n     chomp($r);\n     return $r;\n  \
 }\nsub pdb_has_right_type\n  {\n    my $pdb=shift\
;\n    my $type=shift;\n\n    my $f=vtmpnam();\n\n\
    my $value= &safe_system (\"t_coffee -other_pg \
extract_from_pdb -model_type $pdb > $f\");\n    my\
 $r=&file2string ($f);\n    chomp($r);\n\n\n    if\
 ( $r eq \"NMR\" && $type=~/n/){return 1;}\n    el\
sif ( $r eq \"diffraction\" && $type=~/d/){return \
1;}\n    elsif ( $r eq \"model\" && $type=~/m/){re\
turn 1;}\n    else {return 0;}\n  }\nsub pdb_is_re\
leased\n  {\n    my $pdb=shift;\n    my $f=vtmpnam\
();\n\n    $value= &safe_system (\"t_coffee -other\
_pg extract_from_pdb -is_released_pdb_name $pdb > \
$f\");\n    my $r=&file2string ($f);\n    chomp($r\
);\n    return $r;\n  }\nsub blast_msa\n  {\n    m\
y ($blast,$infile,$db,$outfile)=@_;\n    my ($a, %\
s1, %s, %qs, %qs1);\n    my $seqfile;\n    my $SEQ\
=new FileHandle;\n    my $seqfile=\"seqfile\";\n  \
  my @txt;\n\n\n    %s1=&read_fasta_seq ($db);\n  \
  %s=&fasta_hash2index_hash(%s1);\n    %qs1=&read_\
fasta_seq ($infile);\n    %qs=&fasta_hash2index_ha\
sh(%qs1);\n\n\n    #&safe_system (\"formatdb -i $d\
b\");\n    if ($blast eq \"blastp\"){&safe_system \
 (\"blastall -i $infile -d $db -m7 -p blastp -o io\
\");}\n    elsif ($blast eq \"blastn\"){&safe_syst\
em  (\"blastn -query $infile -db $db -outfmt 5 -wo\
rd_size 4 -out io\");}\n\n    &set_blast_type (\"i\
o\");\n\n\n    my %FB=&xml2tag_list (\"io\", \"Ite\
ration\");\n    open (F, \">$outfile\");\n    prin\
t F \"! TC_LIB_FORMAT_01\\n\";\n    print F \"$s{n\
}\\n\";\n    for ( my $a=0; $a<$s{n}; $a++)\n     \
 {\n	print F \"$s{$a}{name} $s{$a}{len} $s{$a}{seq\
}\\n\";\n      }\n\n\n    for ( my $a=0; $a<$FB{n}\
; $a++)\n      {\n	my %p=blast_xml2profile ($qs{$a\
}{name}, $qs{$a}{seq},100, 0, 0, $FB{$a}{body});\n\
	my $query=$p{0}{name};\n	my $i= $s1{$query}{order\
}+1;\n	for (my $b=1; $b<$p{n}; $b++)\n	  {\n	    m\
y $l=length ($p{$b}{Qseq});\n	    my $hit=$p{$b}{d\
efinition};\n	    my $Qstart=$p{$b}{Qstart};\n	   \
 my $Hstart=$p{$b}{Hstart};\n	    my $identity=$p{\
$b}{identity};\n	    my @lrQ=split (//,$p{$b}{Qseq\
});\n	    my @lrH=split (//,$p{$b}{Hseq});\n\n	   \
 my $j= $s1{$hit}{order}+1;\n	    #if ( $j==$i){ne\
xt;}\n	    printf F \"# %d %d\\n\", $i, $j;\n	    \
#  print  F \"\\n$p{$b}{Qseq} ($Qstart)\\n$p{$b}{H\
seq} ($Hstart)\";\n	    for ($c=0; $c<$l; $c++)\n	\
      {\n		my $rQ=$lrQ[$c];\n		my $rH=$lrH[$c];\n	\
	my $n=0;\n\n		if ($rQ ne \"-\"){$n++, $Qstart++;}\
\n		if ($rH ne \"-\"){$n++; $Hstart++;}\n\n		if ( \
$n==2)\n		  {\n		    printf F \"\\t%d %d %d\\n\", \
$Qstart-1, $Hstart-1,$identity;\n		  }\n	      }\n\
	  }\n      }\n    print F \"! SEQ_1_TO_N\\n\";\n \
   close (F);\n    return $output;\n  }\n\nsub bla\
st_msa_old\n  {\n    my ($infile,$outfile)=@_;\n  \
  my ($a, %seq);\n    %s1=&read_fasta_seq ($infile\
);\n    foreach $s (keys (%s1))\n      {\n	$i=$s1{\
$s}{order};\n	$s{$i}{name}=$s;\n	$s{$i}{seq}=$s1{$\
s}{seq};\n	$s{$i}{len}=length( $s{$i}{seq});\n	$s{\
n}++;\n      }\n    &safe_system (\"formatdb -i $i\
nfile\");\n    &safe_system (\"blastall -i $infile\
 -d $infile -m7 -o io\");\n    &set_blast_type (\"\
io\");\n\n    %FB=&xml2tag_list (\"io\", \"Iterati\
on\");\n\n    open (F, \">$outfile\");\n    print \
F \"! TC_LIB_FORMAT_01\\n\";\n    print F \"$s{n}\\
\n\";\n    for ( $a=0; $a<$s{n}; $a++)\n      {\n	\
print F \"$s{$a}{name} $s{$a}{len} $s{$a}{seq}\\n\\
";\n      }\n    for ( $a=0; $a<$FB{n}; $a++)\n   \
   {\n	%p=blast_xml2profile ($s{$a}{name}, $s{$a}{\
seq},100, 0, 0, $FB{$a}{body});\n	for ($b=1; $b<$p\
{n}; $b++)\n	  {\n	    my $l=length ($p{$b}{Qseq})\
;\n	    my $hit=$p{$b}{definition};\n	    my $Qsta\
rt=$p{$b}{Qstart};\n	    my $Hstart=$p{$b}{Hstart}\
;\n	    my $identity=$p{$b}{identity};\n	    my @l\
rQ=split (//,$p{$b}{Qseq});\n	    my @lrH=split (/\
/,$p{$b}{Hseq});\n	    my $i= $s1{$s{$a}{name}}{or\
der}+1;\n	    my $j= $s1{$hit}{order}+1;\n	    #if\
 ( $j==$i){next;}\n	    printf F \"# %d %d\\n\", $\
i, $j;\n	    #  print  F \"\\n$p{$b}{Qseq} ($Qstar\
t)\\n$p{$b}{Hseq} ($Hstart)\";\n	    for ($c=0; $c\
<$l; $c++)\n	      {\n		my $rQ=$lrQ[$c];\n		my $rH\
=$lrH[$c];\n		my $n=0;\n\n		if ($rQ ne \"-\"){$n++\
, $Qstart++;}\n		if ($rH ne \"-\"){$n++; $Hstart++\
;}\n\n		if ( $n==2)\n		  {\n		    printf F \"\\t%d\
 %d %d\\n\", $Qstart-1, $Hstart-1,$identity;\n		  \
}\n	      }\n	  }\n      }\n    print F \"! SEQ_1_\
TO_N\\n\";\n    close (F);\n    return $output;\n\\
n  }\n\nsub seq2msa\n  {\n    my ($mode, $infile, \
$method, $param, $outfile,$database)=@_;\n    &set\
_temporary_dir (\"set\",$infile,\"seq.pep\", $data\
base, \"db.pep\");\n    $param.=\" >/dev/null 2>&1\
 \";\n\n\n    #make sure test.pep is in FASTA\n   \
 &safe_system (\"t_coffee -other_pg seq_reformat -\
in seq.pep -output fasta_seq > x\");\n    `mv x se\
q.pep`;\n\n    if ( $method eq \"blastp\")\n      \
{\n	&blast_msa (\"blastp\",\"seq.pep\",$database,\\
"result.aln\");\n      }\n    elsif ( $method eq \\
"blastn\")\n      {\n	&blast_msa (\"blastn\",\"seq\
.pep\",$database,\"result.aln\");\n      }\n\n    \
elsif ( $method eq \"muscle\")\n      {\n	`muscle \
-in seq.pep -out result.aln $param`;\n      }\n   \
 elsif ( $method eq \"probcons\")\n      {\n	`prob\
cons seq.pep >result.aln 2>/dev/null`;\n      }\n \
   elsif ( $method eq \"mafft\")\n      {\n	`mafft\
 --quiet --localpair --maxiterate 1000 seq.pep> re\
sult.aln  2>/dev/null`\n      }\n    elsif ( $meth\
od=~/prank/)\n      {\n	`$method -d=seq.pep -o=res\
ult.aln -quiet 2>/dev/null`;\n	`mv result.aln.1.fa\
s result.aln`;\n      }\n    elsif ($method eq \"c\
lustalo\")\n      {\n	`clustalo -i seq.pep > resul\
t.aln`;\n      }\n    else\n      {\n	`$method -in\
file=seq.pep -outfile=result.aln`;\n      }\n\n   \
 &set_temporary_dir (\"unset\",$mode, $method,\"re\
sult.aln\",$outfile);\n    myexit ($EXIT_SUCCESS);\
\n  }\n\nsub seq2thread_pair\n  {\n    my ($mode, \
$infile, $pdbfile1, $method, $param, $outfile)=@_;\
\n    &set_temporary_dir (\"set\",$infile,\"seq.pe\
p\",$pdbfile1,\"struc.pdb\");\n    if ($method eq \
\"fugueali\")\n      {\n	#Env Variable that need t\
o be defined for Fugue\n	if (!$ENV{FUGUE_LIB_LIST}\
){$ENV{FUGUE_LIB_LIST}=\"DUMMY\";}\n	if (!$ENV{HOM\
STRAD_PATH})  {$ENV{HOMSTRAD_PATH}=\"DUMMY\";}\n	i\
f (!$ENV{HOMS_PATH}){$ENV{HOMS_PATH}=\"DUMMY\";}\n\
\n	`joy struc.pdb >x 2>x`;\n	&check_file(\"struc.t\
em\", \"Joy failed [FATAL:$PROGRAM/$method]\");\n	\
`melody -t struc.tem >x 2>x`;\n	&check_file(\"stru\
c.tem\", \"Melody failed [FATAL:$PROGRAM/$method]\\
");\n	`fugueali -seq seq.pep -prf struc.fug -print\
 > tmp_result.aln`;\n\n	&check_file(\"tmp_result.a\
ln\", \"Fugue failed [FATAL:$PROGRAM/$method]\");\\
n	&safe_system (\"t_coffee -other_pg seq_reformat \
-in tmp_result.aln -output fasta_aln >result.aln\"\
);\n      }\n    elsif ( $method eq \"t_coffee\")\\
n      {\n	&safe_system (\"t_coffee -in Pstruc.pdb\
 Sseq.pep Mslow_pair -outfile result.aln -quiet\")\
;\n      }\n    else\n      {\n	&safe_system (\"$m\
ethod -infile=seq.pep -pdbfile1=struc.pdb -outfile\
=result.aln $param>x 2>x\");\n      }\n    &set_te\
mporary_dir (\"unset\",$mode,$method,\"result.aln\\
",$outfile);\n    myexit ($EXIT_SUCCESS);\n  }\nsu\
b seq2pdbid_pair\n  {\n    my ($mode, $pdbfile1, $\
pdbfile2, $method, $param, $outfile)=@_;\n    my (\
$name);\n\n\n    &set_temporary_dir (\"set\");\n  \
  $name=$pdbfile1.\" \".$pdbfile2;\n\n    if (    \
&cache_file(\"GET\",\"\",\"$name\",\"$method\",\"d\
ali\",$outfile,\"EBI\"))\n      {return $outfile;}\
\n    else\n      {\n	if ($method eq \"daliweb\")\\
n	  {\n	    $pdbfile1=~/(....)(.)/;\n	    $id1=$1;\
 $c1=$2;\n\n	    $pdbfile2=~/(....)(.)/;\n	    $id\
2=$1; $c2=$2;\n\n	    $command=\"t_coffee -other_p\
g dalilite.pl --pdb1 $id1 --chainid1 $c1 --pdb2 $i\
d2 --chainid2 $c2 --email=$EMAIL  >dali_stderr 2>d\
ali_stderr\";\n	    $dali=`$command`;\n\n	    open\
 (F, \"dali_stderr\");\n	    while (<F>)\n	      {\
\n		if ( /JobId: dalilite-(\\S+)/)\n		{\n		  $jobi\
d=$1;\n		}\n	      }\n	    close (F);\n	    unlink\
 (\"dali_stderr\");\n\n	    $output1=\"dalilite-$j\
obid.txt\";\n	    if ( -e $output1)\n	      {\n		u\
nlink ($output1);\n		&url2file (\"http://www.ebi.a\
c.uk/Tools/es/cgi-bin/jobresults.cgi/dalilite/dali\
lite-$jobid/aln.html\", \"output2\");\n\n		if ( -e\
 \"output2\")\n		  {\n		    my ($seq1, $seq2);\n		\
    $seq1=$seq2=\"\";\n\n		    open (F, \"output2\\
");\n		    while (<F>)\n		      {\n			$l=$_;\n			i\
f ( $l=~/Query\\s+(\\S+)/)\n			  {\n			    $seq1.=\
$1;\n			  }\n			elsif ( $l=~/Sbjct\\s+(\\S+)/)\n		\
	  {\n			    $seq2.=$1;\n			  }\n		      }\n		    \
close (F);\n		    unlink (\"output2\");\n		    if \
($seq1 ne \"\" && $seq2 ne \"\")\n		      {\n			$o\
utput3=\">$A\\n$seq1\\n>$B\\n$seq2\\n\";\n			$outp\
ut3=~s/\\./-/g;\n			open (F, \">result.aln\");\n		\
	print F \"$output3\";\n			close (F);\n		      }\n\
		  }\n	      }\n	  }\n      }\n    &cache_file(\"\
SET\",\"\",\"$name\",\"$method\",\"dali\",\"result\
.aln\",\"EBI\");\n    &set_temporary_dir (\"unset\\
",$mode, $method, \"result.aln\",$outfile);\n    m\
yexit ($EXIT_SUCCESS);\n  }\nsub seq2pdb_pair\n  {\
\n    my ($mode, $pdbfile1, $pdbfile2, $method, $p\
aram, $outfile)=@_;\n\n    &set_temporary_dir (\"s\
et\",$pdbfile1,\"pdb1.pdb\",$pdbfile2,\"pdb2.pdb\"\
);\n    if ($method eq \"t_coffee\")\n      {\n	&s\
afe_system (\"t_coffee -in Ppdb1.pdb Ppdb2.pdb -qu\
iet -outfile=result.aln\");\n      }\n    elsif ( \
$method eq \"DaliLite\")\n      {\n	if ( &safe_sys\
tem (\"DaliLite -pairwise pdb1.pdb pdb2.pdb >tmp1\\
")==$EXIT_SUCCESS)\n	  {\n	     my ($seq1, $seq2);\
\n	     $seq1=$seq2=\"\";\n\n	     open (F, \"tmp1\
\");\n	     while (<F>)\n	       {\n		 $l=$_;\n		 \
if ( $l=~/Query\\s+(\\S+)/)\n		   {\n		     $seq1.\
=$1;\n		   }\n		 elsif ( $l=~/Sbjct\\s+(\\S+)/)\n	\
	   {\n		     $seq2.=$1;\n		   }\n	       }\n	    \
 close (F);\n	     unlink (\"tmp1\");\n	     if ($\
seq1 ne \"\" && $seq2 ne \"\")\n	       {\n		 my $\
output3=\">$A\\n$seq1\\n>$B\\n$seq2\\n\";\n		 $out\
put3=~s/\\./-/g;\n		 open (F, \">result.aln\");\n	\
	 print F \"$output3\";\n		 close (F);\n	       }\\
n	   }\n	else\n	  {\n	    print \"ERROR: DalLite f\
ailed to align the considered structures[tc_generi\
c_method.pl]\\n\";\n	  }\n      }\n    elsif ( $me\
thod eq \"TMalign\")\n      {\n	if ( &safe_system \
(\"TMalign pdb1.pdb pdb2.pdb >tmp1\")==$EXIT_SUCCE\
SS)\n	  {\n	    `tail -4 tmp1 > tmp2`;\n\n	    ope\
n (F, \"tmp2\");\n	    while (<F>)\n	      {\n		un\
shift(@l, $_);\n	      }\n	    close (F);\n	    op\
en (F, \">result.aln\");\n	    $l[3]=~s/[^a-zA-Z0-\
9-]/\\-/g;\n	    $l[1]=~s/[^a-zA-Z0-9-]/\\-/g;\n	 \
   print F \">$A\\n$l[3]\\n>$B\\n$l[1]\\n\";\n	   \
 close (F);\n	  }\n	else\n	  {\n	    print \"ERROR\
: TMalign failed to align the considered structure\
s[tc_generic_method.pl]\\n\";\n	    `rm result.aln\
 >/dev/null 2>/dev/null`;\n	  }\n      }\n    elsi\
f ( $method eq \"mustang\")\n      {\n	if ( &safe_\
system (\"mustang -i pdb1.pdb pdb2.pdb -F fasta >/\
dev/null 2>/dev/null\")==$EXIT_SUCCESS)\n	  {\n	  \
  `mv results.afasta result.aln`;\n	  }\n	else\n	 \
 {\n	    print \"ERROR: mustang failed to align th\
e considered structures[tc_generic_method.pl]\\n\"\
;\n	    `rm result.aln >/dev/null 2>/dev/null`;\n	\
  }\n      }\n    else\n      {\n	if ( &safe_syste\
m (\"$method -pdbfile1=pdb1.pep -pdbfile2=pdb2.pdb\
 -outfile=result.aln $param>x 2>x\")==$EXIT_SUCCES\
S)\n	  {\n	    `mv results.afasta result.aln`;\n	 \
 }\n	else\n	  {\n	    print \"ERROR: $method faile\
d to align the considered structures[tc_generic_me\
thod.pl]\\n\";\n	    `rm result.aln >/dev/null 2>/\
dev/null`;\n	  }\n      }\n    &set_temporary_dir \
(\"unset\",$mode, $method, \"result.aln\",$outfile\
);\n    myexit ($EXIT_SUCCESS);\n  }\n\nsub seq2rn\
apdb_pair\n  {\n    my ($mode, $pdbfile1, $pdbfile\
2, $method, $param, $outfile)=@_;\n    \n    if ($\
method eq \"runsara.py\")\n      {\n	my $path=$ENV\
{PATH};\n	\n	if ($ENV{X3DNA_4_SARA}){$ENV{PATH}=\"\
$ENV{X3DNA_4_SARA}:$path\";}\n	\n	open(TMP,\"<$pdb\
file1\");\n	my $count = 0;\n	my $line;\n	while (<T\
MP>)\n	  {\n	    $line = $_;\n	    if ($count ==1)\
\n	      {\n		last;\n	      }\n	    $count += 1;\n\
	  }\n	\n	\n	$chain1 = substr($line,length($line)-\
3,1);\n	\n	close TMP;\n	open(TMP,\"<$pdbfile2\");\\
n	my $count = 0;\n	while (<TMP>)\n	  {\n	    $line\
 = $_;\n	    if ($count ==1)\n	      {\n		last;\n	\
      }\n	    $count += 1;\n	  }\n	$chain2 = subst\
r($line,length($line)-3,1);\n	close TMP;\n	\n	$tmp\
_file=&vtmpnam();\n	\n	safe_system(\"runsara.py $p\
dbfile1 $chain1 $pdbfile2 $chain2 -s -o $tmp_file \
--limitation 5000 > /dev/null 2> /dev/null\");\n	i\
f ($ENV{X3DNA_4_SARA}){$ENV{PATH}=$path;}\n	\n	ope\
n(TMP,\"<$tmp_file\") or die \"cannot open the sar\
a tmp file:$!\\n\";\n	open(OUT,\">$outfile\") or d\
ie \"cannot open the $outfile file:$!\\n\";\n	\n	m\
y $switch = 0;\n	my $seqNum = 0;\n	foreach my $lin\
e (<TMP>)\n	  {\n	    next unless ($line=~/SARAALI\
/);\n	    if ($line=~/>/)\n	      {\n		$switch =0;\
\n		print OUT \">seq$seqNum\\n\";\n		$seqNum++;\n	\
      }\n	    if ($switch < 2){\n	      $switch++;\
\n	      next;\n	    }\n	    \n	    if ($line =~/R\
EMARK\\s+SARAALI\\s+([^\\*]+)\\*/)\n	      {\n		my\
 $string = $1;\n		print OUT \"$string\\n\";\n	    \
  }\n	  }\n	close TMP;\n	close OUT;\n	unlink($tmp_\
file);\n      }\n  }\nsub seq2profile_pair\n  {\n \
   my ($mode, $profile1, $profile2, $method, $para\
m, $outfile)=@_;\n    \n    \n    if ($method eq \\
"clustalw\")\n      {\n	`clustalw -profile1=$profi\
le1 -profile2=$profile2 -outfile=$outfile`;\n     \
 }\n    elsif ( $method eq \"clustalo\")\n      {\\
n	\n	`clustalo --p1 $profile1 --p2 $profile2 -o $o\
utfile --force`;\n      }\n    elsif ( $method eq \
\"hhalign\")\n      {\n	hhalign ( $profile1,$profi\
le2,$outfile,$param);\n      }\n    else\n      {\\
n	`$method -profile1=$profile1 -profile2=$profile2\
 -outfile=$outfile $param> /dev/null 2>/dev/null`;\
\n      }\n    myexit ($EXIT_SUCCESS);\n  }\n\nsub\
 pg_is_installed\n  {\n    my @ml=@_;\n    my ($r,\
 $p, $m);\n    my $supported=0;\n\n    my $p=shift\
 (@ml);\n    if ($p=~/::/)\n      {\n	if (safe_sys\
tem (\"perl -M$p -e 1\")==$EXIT_SUCCESS){return 1;\
}\n	else {return 0;}\n      }\n    else\n      {\n\
	$r=`which $p 2>/dev/null`;\n	if ($r eq \"\"){$r=0\
;}\n	else {$r=1;}\n\n	if ($r==0 && is_blast_packag\
e ($p)){return pg_is_installed (\"legacy_blast.pl\\
");}\n	else {return $r;}\n      }\n  }\n\nsub is_b\
last_package\n  {\n    my $p=shift;\n    if ( $p=~\
/blastp/){return 1;}\n    elsif ($p=~/blastall/){r\
eturn 1;}\n    elsif ($p=~/blastn/){return 1;}\n  \
  elsif ($p=~/blastx/){return 1;}\n    elsif ($p=~\
/formatdb/){return 1;}\n    else {return 0;}\n  }\\
n\nsub check_internet_connection\n  {\n    my $int\
ernet;\n    my $tmp;\n    &check_configuration ( \\
"wget\");\n\n    $tmp=&vtmpnam ();\n\n    if     (\
&pg_is_installed    (\"wget\")){`wget www.google.c\
om -O$tmp >/dev/null 2>/dev/null`;}\n    elsif  (&\
pg_is_installed    (\"curl\")){`curl www.google.co\
m -o$tmp >/dev/null 2>/dev/null`;}\n\n    if ( !-e\
 $tmp || -s $tmp < 10){$internet=0;}\n    else {$i\
nternet=1;}\n    if (-e $tmp){unlink $tmp;}\n\n   \
 return $internet;\n  }\nsub check_pg_is_installed\
\n  {\n    my @ml=@_;\n    my $r=&pg_is_installed \
(@ml);\n    if (!$r && $p=~/::/)\n      {\n	print \
STDERR \"\\nYou Must Install the perl package $p o\
n your system.\\nRUN:\\n\\tsudo perl -MCPAN -e 'in\
stall $pg'\\n\";\n      }\n    elsif (!$r)\n      \
{\n	myexit(flush_error(\"\\nProgram $p Supported b\
ut Not Installed on your system\"));\n      }\n   \
 else\n      {\n	return 1;\n      }\n  }\nsub set_\
temporary_dir\n  {\n    my @list=@_;\n    my $dir_\
mode, $a, $mode, $method;\n\n    $dir_mode=shift (\
@list);\n\n\n    if ( $dir_mode eq \"set\")\n     \
 {\n	$initial_dir=cwd();\n	if ( !$tmp_dir)\n	  {\n\
	    $rand=rand (100000);\n	    $tmp_dir=\"$TMPDIR\
/tmp4tcoffee_profile_pair_dir_$$\\_P_$rand\";\n	  \
}\n	if ( !-d $tmp_dir)\n	  {\n	    push (@TMPDIR_L\
IST, $tmp_dir);\n	    `mkdir $tmp_dir`;\n	  }\n\n	\
for ( $a=0; $a<=$#list; $a+=2)\n	      {\n		if (-e\
 $list[$a]){ `cp $list[$a] $tmp_dir/$list[$a+1]`;}\
\n	      }\n	chdir $tmp_dir;\n      }\n    elsif (\
 $dir_mode eq \"unset\")\n      {\n	$mode=shift (@\
list);\n	$method=shift (@list);\n\n	if (!-e $list[\
0])\n	  {\n	   myexit(flush_error(\"Program $metho\
d failed to produce $list[1]\" ));\n	    myexit ($\
EXIT_FAILURE);\n	  }\n	else\n	  {\n	    chdir $ini\
tial_dir;\n	    # `t_coffee -other_pg seq_reformat\
 -in $tmp_dir/$list[0] -output fasta_aln -out $tmp\
_dir/result2.aln`;\n	    `cp $tmp_dir/$list[0] $tm\
p_dir/result2.aln`;\n	    if ( $list[1] eq \"stdou\
t\")\n	      {\n		open (F, \"$tmp_dir/result2.aln\\
");\n		while (<F>){print $_;}close(F);\n	      }\n\
	    else\n	      {\n		`mv $tmp_dir/result2.aln $l\
ist[1]`;\n	      }\n	    shift (@list); shift (@li\
st);\n	    foreach $f (@list)\n	      {\n		if (-e \
(\"$tmp_dir/$f\")){`mv $tmp_dir/$f .`;}\n	      }\\
n	  }\n      }\n  }\n\n\n\n\nsub my_get_opt\n  {\n\
    my @list=@_;\n    my ($cl, $a, $argv, @argl);\\
n\n    \n    @argl=();\n    $cl=shift @list;\n    \
for ( my $a=0; $a<=$#list; $a+=3)\n      {\n	my $o\
ption=$list[$a];\n	my $optional=$list[$a+1];\n	my \
$status=$list[$a+2];\n	my $argv=\"\";\n	if ($cl=~/\
$option(\\S+)/){$argv=$1;}\n	@argl=(@argl,$argv);\\
n\n\n	#$optional:0=>optional\n	#$optional:1=>must \
be set\n	#$status: 0=>no requirement\n	#$status: 1\
=>must be an existing file\n	#$status: 2=>must be \
an installed package\n	\n\n	if ($optional==0){;}\n\
	elsif ( $optional==1 && $argv eq \"\")\n	  {\n	  \
  myexit(flush_error( \"ERROR: Option $option must\
 be set\"));\n	    myexit ($EXIT_FAILURE);\n	  }\n\
	if ($status==0){;}\n	elsif ($status ==1 && $argv \
ne \"\" && !-e $argv)\n	  {\n	    myexit(flush_err\
or( \"File [$argv] must exist\"));\n	    myexit ($\
EXIT_FAILURE);\n	  }\n	elsif ( $status==2 && $argv\
 ne \"\" && &check_pg_is_installed ($argv)==0)\n	 \
 {\n	    myexit(flush_error( \" $argv is not insta\
lled\"));\n	    myexit ($EXIT_FAILURE);\n	  }\n   \
   }\n    return @argl;\n    }\n\nsub check_file\n\
  {\n    my ($file, $msg)=@_;\n\n    if ( !-e $fil\
e)\n      {\n	myexit(flush_error(\"$msg\"));\n    \
  }\n    }\nsub hhalign\n  {\n    my ($aln1, $aln2\
, $outfile, $param)=@_;\n    my $hh=$ENV{\"HHALIGN\
_4_TCOFFEE\"};\n    \n    \n    if ($hh)\n      {\\
n	\n	#external_hhalign\n	# set via HHALIGN_4_TCOFF\
EE\n	#<pg> -profile1 <fasta_prf with seq1 top> -pr\
ofile2 <fasta profile with seq2 top> -outfile < fa\
sta alignmentof seq1 and 2 | tc_lib of seq 1 and 2\
>\n	\n	safe_system (\"$hh -method=align -profile1=\
$aln1 -profile2=$aln2 -outfile=$outfile\");\n     \
 }\n    else\n      {\n	&local_hhalign ($aln1, $al\
n2, $outfile, $param);\n      }\n  }\n\n    \n    \
\nsub local_hhalign\n  {\n    my ($aln1, $aln2, $o\
utfile, $param)=@_;\n    my $h1, $h2;\n\n    $h{0}\
{index}=0;\n    $h{1}{index}=1;\n\n    $h{0}{aln}=\
$aln1;\n    $h{1}{aln}=$aln2;\n\n\n\n    %{$h{0}}=\
aln2psi_profile (%{$h{0}});\n    %{$h{1}}=aln2psi_\
profile (%{$h{1}});\n\n    $param=~s/#S/ /g;\n    \
$param=~s/#M/\\-/g;\n    $param=~s/#E/\\=/g;\n\n\n\
\n    $command=\"hhalign -i $h{0}{a3m} -t $h{1}{a3\
m} -tc $outfile.tmp -rank 1 -mapt 0 $param\";\n   \
 `$command`;\n\n  #  `hhalign -i $h{0}{a3m} -t $h{\
1}{a3m} -tc $outfile.tmp -rank 1 -mapt 0 -gapf 0.8\
 -gapg 0.8`;\n\n\n    # To run global use the foll\
owing\n\n    open (I, \"$outfile.tmp\");\n    open\
 (O, \">$outfile\");\n    $h{0}{cons}=s/\\./x/g;\n\
    $h{1}{cons}=s/\\./x/g;\n\n    print O \"! TC_L\
IB_FORMAT_01\\n2\\n$h{0}{name} $h{0}{len} $h{0}{se\
q}\\n$h{1}{name} $h{1}{len} $h{1}{seq}\\n#1 2\\n\"\
;\n\n    while (<I>)\n      {\n	if (/(\\d+)\\s+(\\\
d+)\\s+(\\d+)/)\n	  {\n	    print O \"\\t$h{0}{$1}\
\\t$h{1}{$2}\\t$3\\n\";\n	  }\n      }\n    print \
O \"! SEQ_1_TO_N\\n\";\n\n    close (O);\n    clos\
e (I);\n  }\n\nsub aln2psi_profile\n  {\n    my (%\
h)=@_;\n    my ($aln,$i,$hv, $a, @c, $n);\n\n\n   \
 $i=$h{index};\n    $aln=$h{aln};\n\n    `cp $aln \
$$.hhh_aln`;\n    $command=\"t_coffee -other_pg se\
q_reformat -in $aln -output hasch\";\n    $hv=`$co\
mmand`;chomp ($hv);\n\n    $h{a2m}=\"$tmp/$hv.tmp4\
hhpred.a2m\";\n    $h{a3m}=\"$tmp/$hv.tmp4hhpred.a\
3m\";\n    if ( -e $h{a3m}){;}\n    else\n      {\\
n	$x=`which hhconsensus`;\n	`hhconsensus  -M 50 -i\
 $h{aln} -oa2m $h{a2m}`;\n	if (!-e $h{a2m})\n	  {\\
n	    print STDERR \"Program tc_generic_method.pl \
FAILED to run:\\n\\thhconsensus  -M 50 -i $h{aln} \
-oa2m $h{a2m}\";\n	    myexit ($EXIT_FAILURE);\n	 \
 }\n\n	`hhconsensus  -M 50 -i $h{aln} -oa3m $h{a3m\
}`;\n	if (!-e $h{a3m})\n	  {\n	    print STDERR \"\
Program tc_generic_method.pl FAILED to run:\\n\\th\
hconsensus  -M 50 -i $h{aln} -oa3m $h{a3m}\";\n	  \
  myexit ($EXIT_FAILURE);\n	  }\n       `buildali.\
pl $h{a3m} -n 1`;\n      }\n\n\n    $h{a2m_seq}=`h\
ead -n 2 $h{a2m} | grep -v \">\"`;chomp ($h{a2m_se\
q});\n    $h{a3m_seq}=`head -n 2 $h{a3m} | grep -v\
 \">\"`;chomp ($h{a3m_seq});\n    $h{cons}=$h{a2m_\
seq};\n    $h{seq}=`head -n 2 $h{aln} | grep -v \"\
>\"`;chomp ($h{seq});\n\n\n\n    @c=split (//, $h{\
cons});\n    $h{len}=$#c+1;\n    for ($n=0,$a=0, $\
b=0; $a<$h{len};$a++)\n      {\n	if ( $c[$a]=~/[A-\
Z]/)\n	  {\n	    $h{++$n}=++$b;\n\n	  }\n	elsif ( \
$c[$a]=~/[a-z\\.]/)\n	  {\n	    ++$b;\n	  }\n     \
 }\n\n    $name=`head -n 2 $h{aln} | grep \">\"`;\\
n    $name=~/\\>(\\S+)/;\n    $h{name}=$1;\n\n    \
`cp $h{a2m} $i.a2m`;\n    `cp $h{a3m} $i.a3m`;\n  \
  `cp $h{aln} $i.hh_aln`;\n\n    return %h;\n  }\n\
sub read_fasta_seq_index\n  {\n    my $f=@_[0];\n \
   my %hseq;\n    my (@seq, @com, @name);\n    my \
($a, $s,$nseq);\n\n    open (F, $f);\n    while (<\
F>)\n      {\n	$s.=$_;\n      }\n    close (F);\n\\
n\n    @name=($s=~/>(\\S*).*\\n[^>]*/g);\n\n    @s\
eq =($s=~/>.*.*\\n([^>]*)/g);\n    @com =($s=~/>\\\
S*(.*)\\n([^>]*)/g);\n\n\n    $nseq=$#name+1;\n\n \
   for ($a=0; $a<$nseq; $a++)\n      {\n	my $s;\n	\
my $n=$name[$a];\n	$hseq{$a}{name}=$n;\n	$seq[$a]=\
~s/[^A-Za-z]//g;\n	$hseq{$a}{order}=$a;\n	$hseq{$a\
}{seq}=$seq[$a];\n	$hseq{$a}{com}=$com[$a];\n\n   \
   }\n    return %hseq;\n  }\nsub read_fasta_seq\n\
  {\n    my $f=@_[0];\n    my %hseq;\n    my (@seq\
, @com, @name);\n    my ($a, $s,$nseq);\n\n    ope\
n (F, $f);\n    while (<F>)\n      {\n	$s.=$_;\n  \
    }\n    close (F);\n\n\n    @name=($s=~/>(\\S*)\
.*\\n[^>]*/g);\n\n    @seq =($s=~/>.*.*\\n([^>]*)/\
g);\n    @com =($s=~/>\\S*(.*)\\n([^>]*)/g);\n\n\n\
    $nseq=$#name+1;\n\n    for ($a=0; $a<$nseq; $a\
++)\n      {\n	my $s;\n	my $n=$name[$a];\n	$hseq{$\
n}{name}=$n;\n	$seq[$a]=~s/[^A-Za-z]//g;\n	$hseq{$\
n}{order}=$a;\n	$hseq{$n}{seq}=$seq[$a];\n	$hseq{$\
n}{com}=$com[$a];\n\n      }\n    return %hseq;\n \
 }\n\n\nsub read_fasta_aln\n  {\n    my $f=@_[0];\\
n    my %hseq;\n    my (@seq, @com, @name);\n    m\
y ($a, $s,$nseq);\n\n    open (F, $f);\n    while \
(<F>)\n      {\n	$s.=$_;\n      }\n    close (F);\\
n\n\n    @name=($s=~/>(\\S*).*\\n[^>]*/g);\n\n    \
@seq =($s=~/>.*.*\\n([^>]*)/g);\n    @com =($s=~/>\
\\S*(.*)\\n([^>]*)/g);\n\n\n    $nseq=$#name+1;\n\\
n    for ($a=0; $a<$nseq; $a++)\n      {\n	my $s;\\
n	my $n=$name[$a];\n	$hseq{$n}{name}=$n;\n	$seq[$a\
]=~s/[^A-Za-z-.()[\\]]//g;\n	$hseq{$n}{order}=$a;\\
n	$hseq{$n}{seq}=$seq[$a];\n	$hseq{$n}{com}=$com[$\
a];\n\n      }\n    return %hseq;\n  }\n\nsub reco\
de_name2\n{\n	my ($in)=shift;\n	my $mode=shift;\n\\
n	my %seq;\n	my $new_name;\n\n	if (! -e $in){retur\
n;}\n\n	#needed by ClustalOmega to avoid very long\
 names\n	open (INFILE, \"+<$in\");\n\n	my $line;\n\
\n	if ($mode eq \"code\")\n	{\n		chomp($line = <IN\
FILE>);\n		my $line_length = length($line);\n		$ne\
w_name=++$RECODE_N;\n		$new_name=\">$new_name\";\n\
		my $new_length = length($new_name);\n		$RECODE {\
$new_name}=$line;\n		for ($count = $new_length; $c\
ount < $line_length; $count++)\n		{\n			$new_name \
.= \" \";\n		}\n		$new_name=\"$new_name\\n\";\n		s\
eek INFILE, 0, 0\n			or die \"could not seek: $!\"\
;\n		print INFILE \"$new_name\";\n	}\n	else\n	{\n	\
	my $n_found = 0;\n		my $file_pos=0;\n		$file_pos=\
tell INFILE;\n		while (<INFILE>)\n		{\n			$line=$_\
;\n			$line =~ s/\\s*$//;\n\n			$old_name= $RECODE\
{$line};\n			if ($old_name ne \"\")\n			{\n				see\
k INFILE, $file_pos, 0\n					or die \"could not se\
ek: $!\";\n				print INFILE \"$old_name\\n\";\n			\
	$file_pos++;\n				if ($file_pos == 2)\n				{\n			\
		print \"stop\\n\";\n					break;\n				}\n			}\n		\
	$file_pos=tell INFILE;\n		}\n\n	}\n\n\n	close INF\
ILE;\n}\n\n\nsub recode_name\n{\n	my ($in)=shift;\\
n	my $mode=shift;\n	my $f=new FileHandle;\n	my %se\
q;\n	my $new_name;\n\n	if (! -e $in){return;}\n\n	\
#needed by ClustalOmega to avoid very long names\n\
	%seq=read_fasta_aln ($in);\n\n	open ($f, \">$in\"\
);\n	foreach my $s (keys(%seq))\n	{\n		if ($mode e\
q \"code\")\n		{\n			$new_name=++$RECODE_N;\n			$R\
ECODE {$new_name}=$seq{$s}{name};\n		}\n		else\n		\
{\n			$new_name=$RECODE{$seq{$s}{name}};\n		}\n		p\
rint $f \">$new_name\\n$seq{$s}{seq}\\n\";\n	}\n	c\
lose $f;\n}\n\nsub fasta_hash2index_hash\n  {\n   \
 my %s1=@_;\n    my %s;\n    foreach my $s (keys (\
%s1))\n      {\n	my $i=$s1{$s}{order};\n	$s{$i}{na\
me}=$s;\n	$s{$i}{seq}=$s1{$s}{seq};\n	$s{$i}{len}=\
length( $s{$i}{seq});\n	$s{n}++;\n      }\n    ret\
urn %s;\n  }\nsub file_contains\n  {\n    my ($fil\
e, $tag, $max)=(@_);\n    my ($n);\n    $n=0;\n\n \
   if ( !-e $file && ($file =~/$tag/)) {return 1;}\
\n    elsif ( !-e $file){return 0;}\n    else\n   \
   {\n	open (FC, \"$file\");\n	while ( <FC>)\n	  {\
\n	    if ( ($_=~/$tag/))\n	      {\n		close (FC);\
\n		return 1;\n	      }\n	    elsif ($max && $n>$m\
ax)\n	      {\n		close (FC);\n		return 0;\n	      \
}\n	    $n++;\n	  }\n      }\n    close (FC);\n   \
 return 0;\n  }\n\n\nsub file2string\n  {\n    my \
$f=@_[0];\n    my $string, $l;\n    open (F,\"$f\"\
);\n    while (<F>)\n      {\n\n	$l=$_;\n	#chomp (\
$l);\n	$string.=$l;\n      }\n    close (F);\n    \
$string=~s/\\r\\n//g;\n    $string=~s/\\n//g;\n   \
 return $string;\n  }\n\n\nsub tag2value\n  {\n\n \
   my $tag=(@_[0]);\n    my $word=(@_[1]);\n    my\
 $return;\n\n    $tag=~/$word=\"([^\"]+)\"/;\n    \
$return=$1;\n    return $return;\n  }\n\nsub hit_t\
ag2pdbid\n  {\n    my $tag=(@_[0]);\n    my $pdbid\
;\n\n    $tag=~/id=\"(\\S+)\"/;\n    $pdbid=$1;\n \
   $pdbid=~s/_//;\n    return $pdbid;\n  }\nsub id\
2pdbid\n  {\n    my $in=@_[0];\n    my $id;\n\n   \
 $in=~/(\\S+)/;\n    $id=$in;\n    $id=~s/PDB/pdb/\
g;\n\n    if ($id =~/pdb(.*)/){$id=$1;}\n    elsif\
 ( $id=~/(\\S+)\\s+mol:protein/){$id=$1;}\n    $id\
=~s/[:|��_]//g;\n    return $id;\n  }\nsub set\
_blast_type\n  {\n    my $file =@_[0];\n    if (&f\
ile_contains ($file,\"EBIApplicationResult\",100))\
{$BLAST_TYPE=\"EBI\";}\n    elsif (&file_contains \
($file,\"NCBI_BlastOutput\",100)) {$BLAST_TYPE=\"N\
CBI\";}\n    else\n      {\n	$BLAST_TYPE=\"\";\n  \
    }\n    return $BLAST_TYPE;\n  }\nsub is_valid_\
blast_xml\n    {\n      my $file=shift;\n      my \
$line;\n\n\n      if ( !-e $file) {return 0;}\n   \
   $line=&file2tail ($file,100);\n\n      if ( $li\
ne=~/<\\/EBIApplicationResult/ || $line=~/<\\/NCBI\
_BlastOutput/ || $line=~/<\\/BlastOutput/ ){return\
 1;}\n      return 0;\n    }\nsub file2blast_flavo\
r\n      {\n	my $file=shift;\n	if (&file_contains \
($file,\"EBIApplicationResult\",100)){return \"EBI\
\";}\n	elsif (&file_contains ($file,\"NCBI_BlastOu\
tput\",100)){return \"NCBI\";}\n	else {return \"UN\
KNOWN\";}\n      }\nsub blast_xml2profile\n  {\n  \
  my ($name,$seq,$maxid, $minid, $mincov, $file)=(\
@_);\n    my (%p, $a, $string, $n);\n\n\n\n    if \
($BLAST_TYPE eq \"EBI\" || &file_contains ($file,\\
"EBIApplicationResult\",100)){%p=ebi_blast_xml2pro\
file(@_);}\n    elsif ($BLAST_TYPE eq \"NCBI\" || \
&file_contains ($file,\"NCBI_BlastOutput\",100)){%\
p=ncbi_blast_xml2profile(@_);}\n    else\n      {\\
n	myexit(add_error ( $$,$$,getppid(), \"BLAST_FAIL\
URE::unkown XML\",$CL));\n      }\n    for ($a=0; \
$a<$p{n}; $a++)\n      {\n	my $name=$p{$a}{name};\\
n	$p{$name}{seq}=$p{$a}{seq};\n	$p{$name}{index}=$\
a;\n      }\n    return %p;\n  }\nsub ncbi_tblastx\
_xml2lib_file\n  {\n    my  ($outlib,$string)=(@_)\
;\n    my ($L,$l, $a,$b,$c,$d,$i,$nhits,@identifye\
rL);\n    my (%ITERATION);\n\n    open (F, \">>$ou\
tlib\");\n\n    $seq=~s/[^a-zA-Z]//g;\n    $L=leng\
th ($seq);\n\n    %ITERATION=xml2tag_list ($string\
, \"Iteration\");\n    for ($i=0; $i<$ITERATION{n}\
;$i++)\n      {\n	my ($qindex, $qlen, %hit, $strin\
g);\n	$string=$ITERATION{$i}{body};\n\n	$qindex=xm\
ltag2value($string,\"Iteration_iter-num\");\n	$qle\
n  =xmltag2value($string,\"Iteration_query-len\");\
\n	%hit=&xml2tag_list  ($string, \"Hit\");\n\n	for\
 ($a=0; $a<$hit{n}; $a++)\n	  {\n	    my ($string)\
;\n	    $string=$hit{$a}{body};\n\n	    $hindex=xm\
ltag2value($string,\"Hit_accession\")+1;\n	    if \
($hindex<=$qindex){next;}\n	    else  {print F  \"\
# $qindex $hindex\\n\";}\n\n\n	    $hlen=xmltag2va\
lue  ($string,\"Hit_len\");\n	    %HSP=&xml2tag_li\
st  ($string, \"Hsp\");\n\n	    for ($b=0; $b<$HSP\
{n}; $b++)\n	      {\n		my ($string, $qs,$qe,$qf,$\
hs,$he,$hf,$s, $d, $e);\n		$string=$HSP{$b}{body};\
\n\n		$qs=xmltag2value  ($string,\"Hsp_query-from\\
");\n		$qe=xmltag2value  ($string,\"Hsp_query-to\"\
);\n		$qf=xmltag2value  ($string,\"Hsp_query-frame\
\");\n\n		$hs=xmltag2value  ($string,\"Hsp_hit-fro\
m\");\n		$he=xmltag2value  ($string,\"Hsp_hit-to\"\
);\n		$hf=xmltag2value  ($string,\"Hsp_hit-frame\"\
);\n\n		$s=xmltag2value  ($string,\"Hsp_identity\"\
);\n		$l=xmltag2value  ($string,\"Hsp_align-len\")\
;\n		$s=int(($s*100)/$l);\n\n		if ($qf>0)\n		  {$r\
qs=$qs; $rqe=$qe;}\n		else\n		  {\n		    $rqe=($ql\
en-$qs)+1;\n		    $rqs=($qlen-$qe)+1;\n		  }\n\n		\
if ($hf>0)\n		  {$rhs=$hs; $rhe=$he;}\n		else\n		 \
 {\n		    $rhe=($hlen-$hs)+1;\n		    $rhs=($hlen-$\
he)+1;\n		  }\n		for ($d=0,$e=$rqs; $e<$rqe; $e++,\
$d++)\n		  {\n		    my ($r1,$r2);\n		    $r1=$e;\n\
		    $r2=$rhs+$d;\n		    print F \" $r1 $r2 $s 0\\
\n\";\n		  }\n	      }\n	  }\n      }\n    print F\
 \"! SEQ_1_TO_N\\n\";\n\n    close (F);\n    retur\
n %lib;\n  }\n\nsub ncbi_tblastpx_xml2lib_file\n  \
{\n    my  ($outlib,$string,%s)=(@_);\n    my ($L,\
$l, $a,$b,$c,$d,$i,$nhits,@identifyerL);\n    my (\
%ITERATION,%hdes, %qdes);\n\n    open (F, \">>$out\
lib\");\n\n    $seq=~s/[^a-zA-Z]//g;\n    $L=lengt\
h ($seq);\n\n    %ITERATION=xml2tag_list ($string,\
 \"Iteration\");\n    for ($i=0; $i<$ITERATION{n};\
$i++)\n      {\n	my ($qindex, $qlen, %hit, $string\
);\n	$string=$ITERATION{$i}{body};\n\n	$qdef=xmlta\
g2value($string,\"Iteration_query-def\");\n	%qdes=\
&tblastpx_name2description($qdef,%s);\n	$qlen  =xm\
ltag2value($string,\"Iteration_query-len\");\n	%hi\
t=&xml2tag_list  ($string, \"Hit\");\n\n	for ($a=0\
; $a<$hit{n}; $a++)\n	  {\n	    my ($string);\n	  \
  $string=$hit{$a}{body};\n	    $hdef=xmltag2value\
($string,\"Hit_def\");\n	    %hdes=&tblastpx_name2\
description($hdef,%s);\n	    if ($hdes{index}<=$qd\
es{index}){next;}\n	    else  {print F  \"# $qdes{\
index} $hdes{index}\\n\";}\n\n\n	    $hlen=xmltag2\
value  ($string,\"Hit_len\");\n	    %HSP=&xml2tag_\
list  ($string, \"Hsp\");\n\n	    for ($b=0; $b<$H\
SP{n}; $b++)\n	      {\n		my ($string, $l,$qs,$qe,\
$qf,$hs,$he,$hf,$s, $d, $e, @s1, @s2);\n		$string=\
$HSP{$b}{body};\n\n		$qs=xmltag2value  ($string,\"\
Hsp_query-from\");\n		$qe=xmltag2value  ($string,\\
"Hsp_query-to\");\n		$qf=$qdes{frame};\n		$qseq=xm\
ltag2value  ($string,\"Hsp_qseq\");\n\n		$hs=xmlta\
g2value  ($string,\"Hsp_hit-from\");\n		$he=xmltag\
2value  ($string,\"Hsp_hit-to\");\n		$hf=$hdes{fra\
me};\n		$hseq=xmltag2value  ($string,\"Hsp_hseq\")\
;\n\n		$s=xmltag2value  ($string,\"Hsp_identity\")\
;\n		$l=xmltag2value  ($string,\"Hsp_align-len\");\
\n		$s=int(($s*100)/$l);\n		@s1=tblastpx_hsp2coord\
inates($qseq,$qs,$qe,%qdes);\n		@s2=tblastpx_hsp2c\
oordinates($hseq,$hs,$he,%hdes);\n\n\n		for ($f=0;\
 $f<=$#s1; $f++)\n		  {\n		    if ($s1[$f]==-1 || \
$s2[$f]==-1){next;}\n		    else\n		      {\n			pri\
nt F \" $s1[$f] $s2[$f] $s 0\\n\";\n		      }\n		 \
 }\n	      }\n	  }\n      }\n    print F \"! SEQ_1\
_TO_N\\n\";\n\n    close (F);\n    return %lib;\n \
 }\nsub tblastpx_hsp2coordinates\n  {\n    my ($se\
q, $s, $e, %des)=@_;\n    my @list;\n    my @sa;\n\
    my @gap=(-1,-1,-1);\n\n    $s=$des{start}+3*($\
s-1);\n\n    if ($des{strand} eq \"d\"){;}\n    el\
se {$s=($des{length}-$s)+1;}\n\n    foreach $c (sp\
lit (//,$seq))\n      {\n	if ( $c eq '-'){push (@l\
ist,@gap);}\n	elsif ($des{strand} eq \"d\")\n	  {\\
n	    push(@list,$s++,$s++,$s++);\n	  }\n	else\n	 \
 {\n	    push(@list, $s--,$s--,$s--);\n	  }\n     \
 }\n    return @list;\n  }\n\nsub tblastpx_name2de\
scription\n  {\n    my ($name, %s)=@_;\n    my @at\
=split(\"__\", $name);\n    my %des;\n\n    $des{n\
ame}=$at[0];\n    $des{strand}=$at[1];\n\n    $des\
{start}=$at[2];\n    $des{end}=$at[3];\n    $des{l\
ength}=$at[4];\n    $des{index}=$s{$at[0]}{order}+\
1;\n    return %des;\n  }\nsub ncbi_blast_xml2prof\
ile\n  {\n    my ($name,$seq,$maxid, $minid, $minc\
ov, $string)=(@_);\n    my ($L,$l, $a,$b,$c,$d,$nh\
its,@identifyerL);\n\n\n    $seq=~s/[^a-zA-Z]//g;\\
n    $L=length ($seq);\n\n    #This is causing the\
 NCBI parser to fail when Iteration_query-def is m\
issing\n    #%query=&xml2tag_list ($string, \"Iter\
ation_query-def\");\n    #$name=$query{0}{body};\n\
\n    %hit=&xml2tag_list ($string, \"Hit\");\n\n\n\
    for ($nhits=0,$a=0; $a<$hit{n}; $a++)\n      {\
\n	my ($ldb,$id, $identity, $expectation, $start, \
$end, $coverage, $r);\n	my (%ID,%DE,%HSP);\n\n	$ld\
b=\"\";\n\n	%ID=&xml2tag_list ($hit{$a}{body}, \"H\
it_id\");\n	$identifyer=$ID{0}{body};\n\n	%DE=&xml\
2tag_list ($hit{$a}{body}, \"Hit_def\");\n	$defini\
tion=$DE{0}{body};\n\n	%HSP=&xml2tag_list ($hit{$a\
}{body}, \"Hsp\");\n	for ($b=0; $b<$HSP{n}; $b++)\\
n	  {\n	    my (%START,%END,%E,%I,%Q,%M);\n\n\n	  \
  %START=&xml2tag_list ($HSP{$b}{body}, \"Hsp_quer\
y-from\");\n	    %HSTART=&xml2tag_list ($HSP{$b}{b\
ody}, \"Hsp_hit-from\");\n\n	    %LEN=  &xml2tag_l\
ist ($HSP{$b}{body}, \"Hsp_align-len\");\n	    %EN\
D=  &xml2tag_list ($HSP{$b}{body}, \"Hsp_query-to\\
");\n	    %HEND=  &xml2tag_list ($HSP{$b}{body}, \\
"Hsp_hit-to\");\n	    %E=&xml2tag_list     ($HSP{$\
b}{body}, \"Hsp_evalue\");\n	    %I=&xml2tag_list \
    ($HSP{$b}{body}, \"Hsp_identity\");\n	    %Q=&\
xml2tag_list     ($HSP{$b}{body}, \"Hsp_qseq\");\n\
	    %M=&xml2tag_list     ($HSP{$b}{body}, \"Hsp_h\
seq\");\n\n	    for ($e=0; $e<$Q{n}; $e++)\n\n	   \
   {\n		$qs=$Q{$e}{body};\n		$ms=$M{$e}{body};\n\n\
		$expectation=$E{$e}{body};\n		$identity=($LEN{$e\
}{body}==0)?0:$I{$e}{body}/$LEN{$e}{body}*100;\n		\
$start=$START{$e}{body};\n		$end=$END{$e}{body};\n\
		$Hstart=$HSTART{$e}{body};\n		$Hend=$HEND{$e}{bo\
dy};\n\n		$coverage=($L)?(($end-$start)*100)/$L:0;\
\n\n		if ($identity>$maxid || $identity<$minid || \
$coverage<$mincov){next;}\n		@lr1=(split (//,$qs))\
;\n		@lr2=(split (//,$ms));\n		$l=$#lr1+1;\n		for \
($c=0;$c<$L;$c++){$p[$nhits][$c]=\"-\";}\n		for ($\
d=0,$c=0; $c<$l; $c++)\n		  {\n		    $r=$lr1[$c];\\
n		    if ( $r=~/[A-Za-z]/)\n		      {\n\n			$p[$n\
hits][$d + $start-1]=$lr2[$c];\n			$d++;\n		      \
}\n		  }\n		$Qseq[$nhits]=$qs;\n		$Hseq[$nhits]=$m\
s;\n		$QstartL[$nhits]=$start;\n		$HstartL[$nhits]\
=$Hstart;\n		$identityL[$nhits]=$identity;\n		$end\
L[$nhits]=$end;\n		$definitionL[$nhits]=$definitio\
n;\n		$identifyerL[$nhits]=$identifyer;\n		$commen\
t[$nhits]=\"$ldb|$identifyer [Eval=$expectation][i\
d=$identity%][start=$Hstart end=$Hend]\";\n		$nhit\
s++;\n	      }\n	  }\n      }\n\n\n    $profile{n}\
=0;\n    $profile{$profile{n}}{name}=$name;\n    $\
profile{$profile{n}}{seq}=$seq;\n    $profile {n}+\
+;\n\n    for ($a=0; $a<$nhits; $a++)\n      {\n	$\
n=$a+1;\n\n	$profile{$n}{name}=\"$name\\_$a\";\n	$\
profile{$n}{seq}=\"\";\n	$profile{$n}{Qseq}=$Qseq[\
$a];\n	$profile{$n}{Hseq}=$Hseq[$a];\n	$profile{$n\
}{Qstart}=$QstartL[$a];\n	$profile{$n}{Hstart}=$Hs\
tartL[$a];\n	$profile{$n}{identity}=$identityL[$a]\
;\n	$profile{$n}{definition}=$definitionL[$a];\n	$\
profile{$n}{identifyer}=$identifyerL[$a];\n	$profi\
le{$n}{comment}=$comment[$a];\n\n	for ($b=0; $b<$L\
; $b++)\n	  {\n	    if ($p[$a][$b])\n	      {\n		$\
profile{$n}{seq}.=$p[$a][$b];\n	      }\n	    else\
\n	      {\n		$profile{$n}{seq}.=\"-\";\n	      }\\
n	  }\n      }\n\n    $profile{n}=$nhits+1;\n    r\
eturn %profile;\n  }\nsub ebi_blast_xml2profile\n \
 {\n    my ($name,$seq,$maxid, $minid, $mincov, $s\
tring)=(@_);\n    my ($L,$l, $a,$b,$c,$d,$nhits,@i\
dentifyerL,$identifyer);\n\n\n\n    $seq=~s/[^a-zA\
-Z]//g;\n    $L=length ($seq);\n    %hit=&xml2tag_\
list ($string, \"hit\");\n\n    for ($nhits=0,$a=0\
; $a<$hit{n}; $a++)\n      {\n	my ($ldb,$id, $iden\
tity, $expectation, $start, $end, $coverage, $r);\\
n	my (%Q,%M,%E,%I);\n\n	$ldb=&tag2value ($hit{$a}{\
open}, \"database\");\n	$identifyer=&tag2value ($h\
it{$a}{open}, \"id\");\n\n	$description=&tag2value\
 ($hit{$a}{open}, \"description\");\n\n	%Q=&xml2ta\
g_list ($hit{$a}{body}, \"querySeq\");\n	%M=&xml2t\
ag_list ($hit{$a}{body}, \"matchSeq\");\n	%E=&xml2\
tag_list ($hit{$a}{body}, \"expectation\");\n	%I=&\
xml2tag_list ($hit{$a}{body}, \"identity\");\n\n\n\
	for ($b=0; $b<$Q{n}; $b++)\n	  {\n\n	    $qs=$Q{$\
b}{body};\n	    $ms=$M{$b}{body};\n\n	    $expecta\
tion=$E{$b}{body};\n	    $identity=$I{$b}{body};\n\
\n\n	    $start=&tag2value ($Q{$b}{open}, \"start\\
");\n	    $end=&tag2value ($Q{$b}{open}, \"end\");\
\n	    $startM=&tag2value ($M{$b}{open}, \"start\"\
);\n	    $endM=&tag2value ($M{$b}{open}, \"end\");\
\n	    $coverage=(($end-$start)*100)/$L;\n\n	   # \
print \"$id: ID: $identity COV: $coverage [$start \
$end]\\n\";\n\n	    if ($identity>$maxid || $ident\
ity<$minid || $coverage<$mincov){next;}\n	    # pr\
int \"KEEP\\n\";\n\n\n	    @lr1=(split (//,$qs));\\
n	    @lr2=(split (//,$ms));\n	    $l=$#lr1+1;\n	 \
   for ($c=0;$c<$L;$c++){$p[$nhits][$c]=\"-\";}\n	\
    for ($d=0,$c=0; $c<$l; $c++)\n	      {\n		$r=$\
lr1[$c];\n		if ( $r=~/[A-Za-z]/)\n		  {\n\n		    $\
p[$nhits][$d + $start-1]=$lr2[$c];\n		    $d++;\n	\
	  }\n	      }\n\n	    $Qseq[$nhits]=$qs;\n	    $H\
seq[$nhits]=$ms;\n	    $QstartL[$nhits]=$start;\n	\
    $HstartL[$nhits]=$Hstart;\n	    $identityL[$nh\
its]=$identity;\n	    $endL[$nhits]=$end;\n	    $d\
efinitionL[$nhits]=$definition;\n	    $identifyerL\
[$nhits]=$identifyer;\n	    $comment[$nhits]=\"$ld\
b|$identifyer [Eval=$expectation][id=$identity%][s\
tart=$startM end=$endM]\";\n	    $nhits++;\n	  }\n\
      }\n\n    $profile{n}=0;\n    $profile{$profi\
le{n}}{name}=$name;\n    $profile{$profile{n}}{seq\
}=$seq;\n    $profile {n}++;\n\n    for ($a=0; $a<\
$nhits; $a++)\n      {\n	$n=$a+1;\n	$profile{$n}{n\
ame}=\"$name\\_$a\";\n	$profile{$n}{seq}=\"\";\n	$\
profile{$n}{Qseq}=$Qseq[$a];\n	$profile{$n}{Hseq}=\
$Hseq[$a];\n	$profile{$n}{Qstart}=$QstartL[$a];\n	\
$profile{$n}{Hstart}=$HstartL[$a];\n	$profile{$n}{\
identity}=$identityL[$a];\n	$profile{$n}{definitio\
n}=$definitionL[$a];\n	$profile{$n}{identifyer}=$i\
dentifyerL[$a];\n	$profile{$n}{comment}=$comment[$\
a];\n\n	for ($b=0; $b<$L; $b++)\n	  {\n	    if ($p\
[$a][$b])\n	      {\n		$profile{$n}{seq}.=$p[$a][$\
b];\n	      }\n	    else\n	      {\n		$profile{$n}\
{seq}.=\"-\";\n	      }\n	  }\n      }\n    $profi\
le{n}=$nhits+1;\n\n    return %profile;\n  }\nsub \
output_profile\n  {\n    my ($outfile,$profileR, $\
trim)=(@_);\n    my ($a);\n    my %profile=%$profi\
leR;\n    my $P= new FileHandle;\n    my $tmp=vtmp\
nam();\n\n    open ($P, \">$tmp\");\n    for ($a=0\
; $a<$profile{n}; $a++)\n      {\n	print $P \">$pr\
ofile{$a}{name} $profile{$a}{comment}\\n$profile{$\
a}{seq}\\n\";\n      }\n    close ($P);\n\n    if \
( $trim)\n      {\n	\n	if ($trim>0)\n	  {\n	    &s\
afe_system (\"t_coffee -other_pg seq_reformat -in \
$tmp -action +trim _aln_n$trim\\_K1 -output fasta_\
aln -out $outfile\");\n	  }\n	else\n	  {\n	    &sa\
fe_system (\"t_coffee -other_pg seq_reformat -in $\
tmp -action +trim _aln_%%$trim\\_K1 -output fasta_\
aln -out $outfile\");\n	  }\n      }\n    else\n  \
    {\n	&safe_system (\"mv $tmp $outfile\");\n    \
  }\n    return;\n  }\nsub blast_xml2hit_list\n  {\
\n    my $string=(@_[0]);\n    return &xml2tag_lis\
t ($string, \"hit\");\n  }\nsub xmltag2value\n  {\\
n    my ($string_in, $tag)=@_;\n    my %TAG;\n    \
%TAG=xml2tag_list ($string_in, $tag);\n    return \
$TAG{0}{body};\n  }\n\nsub xml2tag_list\n  {\n    \
my ($string_in,$tag)=@_;\n    my $tag_in, $tag_out\
;\n    my %tag;\n\n    if (-e $string_in)\n      {\
\n	$string=&file2string ($string_in);\n      }\n  \
  else\n      {\n	$string=$string_in;\n      }\n  \
  $tag_in1=\"<$tag \";\n    $tag_in2=\"<$tag>\";\n\
    $tag_out=\"/$tag>\";\n    $string=~s/>/>##1/g;\
\n    $string=~s/</##2</g;\n    $string=~s/##1/<#/\
g;\n    $string=~s/##2/#>/g;\n    @l=($string=~/(\\
\<[^>]+\\>)/g);\n    $tag{n}=0;\n    $in=0;$n=-1;\\
n\n\n\n    foreach $t (@l)\n      {\n\n	$t=~s/<#//\
;\n	$t=~s/#>//;\n\n	if ( $t=~/$tag_in1/ || $t=~/$t\
ag_in2/)\n	  {\n\n	    $in=1;\n	    $tag{$tag{n}}{\
open}=$t;\n	    $n++;\n\n	  }\n	elsif ($t=~/$tag_o\
ut/)\n	  {\n\n\n	    $tag{$tag{n}}{close}=$t;\n	  \
  $tag{n}++;\n	    $in=0;\n	  }\n	elsif ($in)\n	  \
{\n\n	    $tag{$tag{n}}{body}.=$t;\n	  }\n      }\\
n\n    return %tag;\n  }\n\n\nsub seq2gor_predicti\
on\n  {\n    my ($name, $seq,$infile, $outfile, $g\
or_seq, $gor_obs)=(@_);\n    my ($l);\n\n    `gorI\
V -prd $infile -seq $gor_seq -obs $gor_obs > gor_t\
mp`;\n    open (GR, \">$outfile\");\n    open (OG,\
 \"gor_tmp\");\n\n    while (<OG>)\n      {\n\n	$l\
=$_;\n	if ($l=~/\\>/){print GR \"$l\";}\n	elsif ( \
$l=~/Predicted Sec. Struct./)\n	  {\n	    $l=~s/Pr\
edicted Sec. Struct\\.//;\n	    print GR \"$l\";\n\
	  }\n      }\n    close (GR);\n    close (OG);\n \
   return;\n  }\nsub seq2msa_tm_prediction\n  {\n \
   my ($name, $seq, $db, $infile, $outfile, $arch,\
 $psv)=(@_);\n    my (%p,%gseq,%R, $blast_output, \
%s, $l);\n    my $R2=new FileHandle;\n    my $meth\
od=\"psitm\";\n    my $SERVER=\"EBI\";\n\n    $bla\
st_output=&run_blast ($name,\"blastp\", $db, $infi\
le, \"outfile\");\n\n    if (&cache_file(\"GET\",$\
infile,$name,$method,$db,$outfile,$SERVER))\n     \
 {\n	print \"\\tPSITM: USE Cache\\n\";\n	return $o\
utfile;\n      }\n    else\n      {\n	$CACHE_STATU\
S=\"COMPUTE CACHE\";\n	%p=blast_xml2profile($name,\
$seq,$maxid, $minid,$mincov,$blast_output);\n\n\n	\
open (F, \">tm_input\");\n	for (my $a=0; $a<$p{n};\
 $a++)\n	  {\n	    my $s;\n\n	    $s=$p{$a}{seq};\\
n	    $s=uc($s);\n	    print F \">$p{$a}{name}\\n$\
s\\n\";\n	    #print stdout \">$p{$a}{name}\\n$s\\\
n\";\n	  }\n	close (F);\n	print \"\\tPSITM: kept  \
$p{n} Homologues for Sequence $p{0}{name}\\n\";\n	\
&safe_system (\"t_coffee -other_pg fasta_seq2hmmto\
p_fasta.pl -in=tm_input -out=$outfile -output=cons\
 -cov=70 -trim=95 -arch=$arch -psv=$psv\");\n	unli\
nk (\"tm_input\");\n	&cache_file(\"SET\",$infile,$\
name,$method,$db,$outfile,$SERVER);\n	return;\n   \
   }\n  }\n\n\nsub seq2msa_gor_prediction\n  {\n  \
  my ($name, $seq,$infile, $outfile, $gor_seq, $go\
r_obs)=(@_);\n    my (%p,%gseq,%R, $blast_output, \
%s, $l);\n    my $R2=new FileHandle;\n    my $db=\\
"uniprot\";\n    my $method=\"psigor\";\n    my $S\
ERVER=\"EBI\";\n\n    $blast_output=&run_blast ($n\
ame,\"blastp\", \"uniprot\", $infile, \"outfile\")\
;\n\n    if (&cache_file(\"GET\",$infile,$name,$me\
thod,$db,$outfile,$SERVER))\n      {\n	print \"\\t\
PSIGOR: USE Cache\\n\";\n	return $outfile;\n      \
}\n    else\n      {\n	$CACHE_STATUS=\"COMPUTE CAC\
HE\";\n	%p=blast_xml2profile($name,$seq,$maxid, $m\
inid,$mincov,$blast_output);\n\n\n	open (F, \">gor\
_input\");\n	for (my $a=0; $a<$p{n}; $a++)\n	  {\n\
	    my $s;\n\n	    $s=$p{$a}{seq};\n	    $s=uc($s\
);\n	    print F \">$p{$a}{name}\\n$s\\n\";\n	    \
#print stdout \">$p{$a}{name}\\n$s\\n\";\n	  }\n	c\
lose (F);\n	print \"\\tGORTM: kept  $p{n} Homologu\
es for Sequence $p{0}{name}\\n\";\n	&safe_system (\
\"t_coffee -other_pg fasta_seq2hmmtop_fasta.pl -in\
=gor_input -out=$outfile -output=cons -cov=70 -tri\
m=95 -gor_seq=$gor_seq -gor_obs=$gor_obs -mode=gor\
\");\n	unlink (\"tm_input\");\n	&cache_file(\"SET\\
",$infile,$name,$method,$db,$outfile,$SERVER);\n	r\
eturn;\n      }\n  }\n\n\n\nsub run_blast\n  {\n  \
  my ($name, $method, $db, $infile, $outfile, $run\
)=(@_);\n    if (!$run){$run=1;}\n    my $error_lo\
g=vtmpnam();\n    my $cl_db;\n    \n    if (&cache\
_file(\"GET\",$infile,$name,$method,$db,$outfile,$\
SERVER) && is_valid_blast_xml ($outfile))\n      {\
return $outfile;}\n    else\n      {\n	$CACHE_STAT\
US=\"COMPUTE CACHE\";\n	if ( $SERVER eq \"EBI_SOAP\
\")\n	  {\n	    &check_configuration (\"EMAIL\",\"\
SOAP::Light\",\"INTERNET\");\n\n	    $cl_method=$m\
ethod;\n	    if ($cl_method =~/wu/)\n	      {\n		$\
cl_method=~s/wu//;\n		if ( $cl_method eq \"psiblas\
t\")\n		  {\n		    add_warning($$,$$,\"PSI BLAST c\
annot be used with the wuBLAST Client. Use server=\
EBI Or server=LOCAL. blastp will be used instead\"\
);\n		    $cl_method=\"blastp\";\n		  }\n\n		$comm\
and=\"t_coffee -other_pg wublast.pl --email $EMAIL\
 $infile -D $db -p $cl_method --outfile $outfile -\
o xml>/dev/null 2>$error_log\";\n		&safe_system ( \
$command);\n		if (-e \"$outfile.xml\") {`mv $outfi\
le.xml $outfile`;}\n	      }\n	    else\n	      {\\
n		if ($cl_method eq \"psiblast\"){$cl_method =\"b\
lastp -j5\";}\n\n		$command=\"t_coffee -other_pg b\
lastpgp.pl --email $EMAIL $infile -d $db --outfile\
 $outfile -p $cl_method --mode PSI-Blast>/dev/null\
 2>$error_log\";\n		&safe_system ( $command);\n\n	\
	if (-e \"$outfile.xml\") {`mv $outfile.xml $outfi\
le`;}\n	      }\n	  }\n	elsif ($SERVER eq \"EBI_RE\
ST\" || $SERVER eq \"EBI\")\n	  {\n	    $cl_method\
=$method;\n	    &check_configuration(\"EMAIL\",\"X\
ML::Simple\", \"INTERNET\");\n	    if ($db eq \"un\
iprot\"){$db1=\"uniprotkb\";}\n	    else {$db1=$db\
;}\n\n	    \n	    if ($cl_method =~/wu/)\n	      {\
\n		$cl_method=~s/wu//;\n		if ( $cl_method eq \"ps\
iblast\"){$cl_method=\"blastp\";}\n\n		$command=\"\
t_coffee -other_pg wublast_lwp.pl --email $EMAIL -\
D $db1 -p $cl_method --outfile $outfile --align 5 \
--stype protein $infile>/dev/null 2>error_log\";\n\
	      }\n	    else\n	      {\n		if ( $cl_method =\
~/psiblast/){$cl_method =\"blastp -j5\";}\n		$comm\
and=\"t_coffee -other_pg ncbiblast_lwp.pl --email \
$EMAIL -D $db1 -p $cl_method --outfile $outfile --\
align 5 --stype protein $infile>/dev/null 2>$error\
_log\";\n		#DEBUG\n		#$command=\"t_coffee -other_p\
g ncbiblast_lwp.pl --email $EMAIL -D $db1 -p $cl_m\
ethod --outfile $outfile --align 5 --stype protein\
 $infile\";\n		\n		my $maxrun=5;#number of crashes\
 accepetd\n		my $nrun;\n		my $keep_going=1;\n		whi\
le ($keep_going)\n		  {\n		    \n		    #print \"--\
--- $command [$nrun]\\n\";\n		    $nrun++;\n		    \
$keep_going=0;\n		    &safe_system ( $command,5);\\
n		    \n		    my $success=0;\n		    $success =$su\
ccess || -e \"$outfile.out.xml\";\n		    $success \
=$success || -e \"$outfile.xml.xml\";\n		    $succ\
ess =$success || -e \"$outfile.out..xml\";\n		    \
$success =$success || -e \"$outfile.xml..xml\";\n	\
	    \n		    if (!$success && ($nrun<$maxrun || -e\
 \"$outfile.out.txt\"))\n		      {\n			$keep_going\
=1;\n			add_warning($$,$$,\"[ncbiblast_lwp.pl] [$c\
ommand] failed to produce xml output -- will ne tr\
ied again [$nrun]\");\n		      }\n		  }\n		\n		if \
(-e \"$outfile.out.xml\") {`mv $outfile.out.xml $o\
utfile`;}\n		elsif (-e \"$outfile.xml.xml\"){`mv $\
outfile.xml.xml $outfile`;}\n		elsif (-e \"$outfil\
e.out..xml\") {`mv $outfile.out..xml $outfile`;}\n\
		elsif (-e \"$outfile.xml..xml\"){`mv $outfile.xm\
l..xml $outfile`;}\n		else\n		  {\n		    add_warni\
ng($$,$$,\"[ncbiblast_lwp.pl] [$command] failed to\
 produce xml output\");\n		  }\n	      }\n	  }\n	e\
lsif ($SERVER eq \"NCBI\")\n	  {\n	    &check_conf\
iguration (\"INTERNET\");\n	    if ($db eq \"unipr\
ot\"){$cl_db=\"swissprot\";}\n	    else {$cl_db=$d\
b;}\n\n	    if ( $method eq \"psiblast\")\n	      \
{\n		add_warning($$,$$,\"PSI BLAST cannot be used \
with the NCBI BLAST Client. Use server=EBI Or serv\
er=LOCAL. blastp will be used instead\");\n		$cl_m\
ethod=\"blastp\";\n	      }\n	    else\n	      {\n\
		$cl_method=$method;\n	      }\n	      \n	    &ch\
eck_configuration ($cl_method);  \n	    $command=\\
"$cl_method -db $cl_db -query $infile -out $outfil\
e -outfmt 5 -remote\";\n	    &safe_system ($comman\
d);\n	  }\n	elsif ($SERVER =~/CLIENT_(.*)/)\n	  {\\
n	    my $client=$1;\n	    $command=\"$client -p $\
method -d $db -i $infile -o $outfile -m 7\";\n	   \
 &safe_system ($command);\n	  }\n	elsif ( $SERVER \
eq \"LOCAL_blastall\")\n	  {\n	    &check_configur\
ation (\"blastall\");\n	    if ($method eq \"blast\
p\")\n	      {\n		$command=\"blastall -d $db -i $i\
nfile -o $outfile -m7 -p blastp\";\n	      }\n	   \
 &safe_system ($command);\n	  }\n	elsif ( $SERVER \
eq \"LOCAL\")\n	  {\n	    if ($ENV{\"BLAST_DB_DIR\\
"}) {\n	    	$x=$ENV{\"BLAST_DB_DIR\"};\n			$cl_db\
=\"$x/$db\";\n	    }\n	    else{\n			$cl_db=$db;\n\
	    }\n\n        ##\n		## BLAST+ provide differen\
t binaries names and CLI options\n		## Use the 'le\
gacy_blast.pl' to keep compatibility with old blas\
t commands\n		##\n		$path=`which legacy_blast.pl 2\
>/dev/null`;  \n		$path=`dirname $path`; \n		chomp\
($path);\n	    if ($method eq \"blastp\"){\n			&ch\
eck_configuration(\"legacy_blast.pl\");\n			$comma\
nd=\"legacy_blast.pl blastpgp --path $path -d $cl_\
db -i $infile -o $outfile -m7 -j1\";\n	    }\n	   \
 elsif ($method eq \"psiblast\")\n	      {\n		&che\
ck_configuration(\"legacy_blast.pl\");\n		$command\
=\"legacy_blast.pl blastpgp --path $path -d $cl_db\
 -i $infile -o $outfile -m7 -j5\";\n	      }\n	   \
 elsif ($method eq \"blastn\")\n	      {\n		&check\
_configuration(\"legacy_blast.pl\");\n		$command=\\
"legacy_blast.pl blastall --path $path -p blastn -\
d $cl_db -i $infile -o $outfile -m7 -W6\";\n	     \
 }\n	    &safe_system ($command);\n	  }\n	else\n	 \
 {\n\n	    myexit(add_error (EXIT_FAILURE,$$,$$,ge\
tppid(), \"BLAST_FAILURE::UnknownServer\",$CL));\n\
	  }\n\n\n	#Check that everything went well\n\n	if\
 ( !-e $outfile || !&is_valid_blast_xml($outfile))\
\n	  {\n\n	    if ( -e $outfile)\n	      {\n		add_\
warning ($$,$$,\"Corrupted Blast Output (Run $run)\
\");\n		unlink($outfile);\n	      }\n	    if ( -e \
$error_log)\n	      {\n\n		my $error_msg=file2stri\
ng ($error_log);\n\n		if ( $error_msg =~/enter a v\
alid email/)\n		  {\n		    myexit(add_error (EXIT_\
FAILURE,$$,$$,getppid(), \"BLAST_FAILURE::Invalid_\
or_rejected_email::$EMAIL\", \"$command\"));\n		  \
}\n	      }\n	    if ( $run==$BLAST_MAX_NRUNS)\n	 \
     {\n\n		myexit(add_error (EXIT_FAILURE,$$,$$,g\
etppid(), \"BLAST_FAILURE::UnknownReason\", \"$com\
mand\"));\n	      }\n	    else\n	      {\n		my $ou\
t;\n		if ($SERVER eq \"NCBI\") {$SERVER=\"EBI\"; }\
\n		elsif ($SERVER eq \"EBI\"){$SERVER=\"NCBI\";}\\
n		add_warning ($$,$$,\"Blast for $name failed (Ru\
n: $run out of $BLAST_MAX_NRUNS. Use $SERVER)\");\\
n		$out=&run_blast ($name, $method, $db,$infile, $\
outfile, $run+1);\n		if ($SERVER eq \"NCBI\") {$SE\
RVER=\"EBI\"; }\n		elsif ($SERVER eq \"EBI\"){$SER\
VER=\"NCBI\";}\n		return $out;\n	      }\n	  }\n\n\
	&cache_file(\"SET\",$infile,$name,$method,$db,$ou\
tfile,$SERVER);\n	#system (\"cp $outfile ~/Dropbox\
/tmp/cedric.out\");\n	#die;\n	return $outfile;\n  \
    }\n  }\n\nsub cache_file\n  {\n    my ($cache_\
mode,$infile,$name,$method,$db, $outfile,$server)=\
(@_);\n    my $cache_file;\n    #Protect names so \
that they can be turned into legal filenames\n    \
$name=&clean_file_name ($name);\n\n    if ($db=~/\\
\//)\n      {\n	$db=~/([^\\/]+)$/;\n	$db=$1;\n    \
  }\n    $cache_file_sh=\"$name.$method.$db.$serve\
r.tmp\";\n    $cache_file=\"$CACHE/$name.$method.$\
db.$server.tmp\";\n\n    if ($infile ne \"\")\n   \
   {\n	$cache_file_infile_sh=\"$name.$method.$db.$\
server.infile.tmp\";\n	$cache_file_infile=\"$CACHE\
/$name.$method.$db.$server.infile.tmp\";\n      }\\
n\n    if ($cache_mode eq \"GET\")\n      {\n	if (\
$CACHE eq \"\" || $CACHE eq \"no\" || $CACHE eq \"\
ignore\"  || $CACHE eq \"local\" || $CACHE eq \"up\
date\"){return 0;}\n	elsif ( !-d $CACHE)\n	  {\n	 \
   print STDERR \"ERROR: Cache Dir: $CACHE Does no\
t Exist\";\n	    return 0;\n	  }\n	else\n	  {\n	  \
  if ( -e $cache_file && &fasta_file1_eq_fasta_fil\
e2($infile,$cache_file_infile)==1)\n	      {\n		`c\
p $cache_file $outfile`;\n		$CACHE_STATUS=\"READ C\
ACHE\";\n		return 1;\n	      }\n	  }\n      }\n   \
 elsif ($cache_mode eq \"SET\")\n      {\n	if ($CA\
CHE eq \"\" || $CACHE eq \"no\" || $CACHE eq \"ign\
ore\"  || $CACHE eq \"local\" || $CACHE eq \"updat\
e\"){return 0;}\n	elsif ( !-d $CACHE)\n	  {\n	    \
print STDERR \"ERROR: Cache Dir: $CACHE Does not E\
xist\";\n	    return 0;\n	  }\n	elsif (-e $outfile\
)\n	  {\n	    `cp $outfile $cache_file`;\n	    if \
($cache_file_infile ne \"\"){ `cp $infile $cache_f\
ile_infile`;}\n\n	    #functions for updating the \
cache\n	    #`t_coffee -other_pg clean_cache.pl -f\
ile $cache_file_sh -dir $CACHE`;\n	    #`t_coffee \
-other_pg clean_cache.pl -file $cache_file_infile_\
sh -dir $CACHE`;\n	    return 1;\n	  }\n      }\n \
   $CACHE_STATUS=\"COMPUTE CACHE\";\n    return 0;\
\n  }\nsub file1_eq_file2\n  {\n    my ($f1, $f2)=\
@_;\n    if ( $f1 eq \"\"){return 1;}\n    elsif (\
 $f2 eq \"\"){return 1;}\n    elsif ( !-e $f1){ret\
urn 0;}\n    elsif ( !-e $f2){return 0;}\n    elsi\
f ($f1 eq \"\" || $f2 eq \"\" || `diff $f1 $f2` eq\
 \"\"){return 1;}\n\n    return 0;\n  }\nsub clean\
_file_name\n  {\n    my $name=@_[0];\n\n    $name=\
~s/[^A-Za-z1-9.-]/_/g;\n    return $name;\n  }\nsu\
b url2file\n  {\n    my ($address, $out)=(@_);\n\n\
    if (&pg_is_installed (\"wget\"))\n	{\n	  retur\
n &safe_system (\"wget $address -O$out >/dev/null \
2>/dev/null\");\n	}\n    elsif (&pg_is_installed (\
\"curl\"))\n      {\n	return &safe_system (\"curl \
$address -o$out >/dev/null 2>/dev/null\");\n      \
}\n    else\n      {\n	myexit(flus_error(\"neither\
 curl nor wget are installed. Imnpossible to fectc\
h remote file\"));\n	exit ($EXIT_FAILURE);\n      \
}\n  }\nsub fasta_file1_eq_fasta_file2\n  {\n    m\
y ($f1, $f2)=@_;\n    my (%s1, %s2);\n    my @name\
s;\n    %s1=read_fasta_seq ($f1);\n    %s2=read_fa\
sta_seq ($f2);\n\n    @names=(keys (%s1));\n\n    \
foreach $n (keys(%s1))\n      {\n	if ($s1{$n}{seq}\
 ne $s2{$n}{seq}){return 0;}\n      }\n\n    forea\
ch $n (keys(%s2))\n      {\n	if ($s1{$n}{seq} ne $\
s2{$n}{seq}){return 0;}\n      }\n    return 1;\n \
 }\n\n\n\nsub read_template_file\n  {\n    my $pdb\
_templates = @_[0];\n    my $tmp=new FileHandle;\n\
    open ($tmp, \"<$pdb_templates\");\n    my %tem\
p_h;\n    my ($skip,$seq, $temp);\n\n    #supports\
 both a simple [seq pdb] format and the regular fa\
sta like template format\n    while (<$tmp>)\n    \
  {\n	\n	$line = $_;\n	if ($line=~/\\>(\\S+)\\s_._\
\\s(\\S+)/){$temp_h{$1}= $2;}\n	elsif ($line =~/(\\
\S+)\\s(\\S+)/){$temp_h{$1}= $2;}\n      }\n    cl\
ose($tmp);\n    return %temp_h;\n  }\n\n\n\n\n\n\n\
sub seq2tblastx_lib\n  {\n    my ($mode, $infile, \
$outfile)=@_;\n    my (%s, $method,$nseq);\n\n    \
$method=$mode;\n    &set_temporary_dir (\"set\",$i\
nfile,\"infile\");\n    %s=read_fasta_seq(\"infile\
\");\n\n\n    foreach $seq (keys(%s))\n      {\n	$\
slist[$s{$seq}{order}]=$s{$seq}{seq};\n	$sname[$s{\
$seq}{order}]=$s{$seq}{name};\n	$slen[$s{$seq}{ord\
er}]=length ($s{$seq}{seq});\n      }\n    $nseq=$\
#sname+1;\n    open (F, \">outfile\");\n    print \
F \"! TC_LIB_FORMAT_01\\n\";\n    print F \"$nseq\\
\n\";\n    for ($a=0; $a<$nseq;$a++)\n      {\n	pr\
int F \"$sname[$a] $slen[$a]  $slist[$a]\\n\"\n   \
   }\n    close (F);\n    &safe_system (\"formatdb\
 -i infile -p F\");\n    &safe_system (\"blastall \
-p tblastx -i infile -d infile -m 7 -S1>blast.outp\
ut\");\n\n    ncbi_tblastx_xml2lib_file (\"outfile\
\", file2string (\"blast.output\"));\n    &set_tem\
porary_dir (\"unset\",$mode, $method, \"outfile\",\
$outfile);\n    myexit ($EXIT_SUCCESS);\n    }\nsu\
b seq2tblastpx_lib\n  {\n    my ($mode, $infile, $\
outfile)=@_;\n    my (%s, $method,$nseq);\n    $me\
thod=$mode;\n    &set_temporary_dir (\"set\",$infi\
le,\"infile\");\n    %s=read_fasta_seq(\"infile\")\
;\n\n    foreach $seq (keys(%s))\n      {\n	$slist\
[$s{$seq}{order}]=$s{$seq}{seq};\n	$sname[$s{$seq}\
{order}]=$s{$seq}{name};\n	$slen[$s{$seq}{order}]=\
length ($s{$seq}{seq});\n      }\n    $nseq=$#snam\
e+1;\n    open (F, \">outfile\");\n    print F \"!\
 TC_LIB_FORMAT_01\\n\";\n    print F \"$nseq\\n\";\
\n    for ($a=0; $a<$nseq;$a++)\n      {\n	print F\
 \"$sname[$a] $slen[$a]  $slist[$a]\\n\"\n      }\\
n    close (F);\n    &safe_system(\"t_coffee -othe\
r_pg seq_reformat -in infile -output tblastx_db1 >\
 tblastxdb\");\n    &safe_system (\"formatdb -i tb\
lastxdb -p T\");\n    &safe_system (\"blastall -p \
blastp -i tblastxdb -d tblastxdb -m7 >blast.output\
\");\n    ncbi_tblastpx_xml2lib_file (\"outfile\",\
 file2string (\"blast.output\"), %s);\n    &set_te\
mporary_dir (\"unset\",$mode, $method, \"outfile\"\
,$outfile);\n    myexit ($EXIT_SUCCESS);\n    }\n\\
nsub x3dna_find_pair2lib\n      {\n      my ($seq,\
 $pdb, $lib, $pg)=@_;\n      my $outfile1=\"dssr-2\
ndstrs.dbn\";\n      my $outfile2=\"simple.output\\
";\n      my $f= new FileHandle;\n      my ($rnaSS\
,$pdbSS);\n      my $command;\n      my %s_pdb;\n \
     my %s_seq;\n      \n      #$pg: \"find_pair\"\
 OR \"find_pair -p\"\n      \n      if (!pg_is_ins\
talled (\"find_pair\"))\n	{\n	  add_warning ($$,$$\
, \"x3dna/find_pairs could not be used to extract \
RNA secondary structures. Secondary structures wil\
l be extracted by x3dna-ssr instead\");\n	  return\
 x3dnassr2lib ($seq, $pdb, $lib);\n	}\n      \n   \
   #get PDB sequence\n      safe_system (\"t_coffe\
e -other_pg extract_from_pdb $pdb -seq >$outfile1\\
");\n      \n      #get find_pair contacts\n      \
$command=\"$pg $pdb simple.output > /dev/null 2>/d\
ev/null\";\n      safe_system ($command);\n\n     \
 if (($command=~/find_pair -p/)){$outfile2=\"allpa\
irs.ana\";}\n      else {$outfile2=\"simple.output\
\";}\n      \n      if ( !-e $outfile2)\n	{\n	  my\
exit(flush_error(\"x3dna failed to compute the sec\
ondary structure file $outfile2 for $pdb\"));\n	  \
myexit ($EXIT_FAILURE);\n	}\n      \n\n      #Hand\
le situations when the pdb sequence differs from t\
he RNA sequence\n      #my @out=file2array($outfil\
e1);\n      %s_pdb=read_fasta_seq_index ($outfile1\
);\n      %s_seq=read_fasta_seq_index ($seq);\n   \
   my $rnaS=uc($s_seq{0}{seq});\n      my $pdbS=uc\
($s_pdb{0}{seq});\n      \n      my $vienna;\n    \
  my @lu;\n    \n      if ($rnaS ne $pdbS)\n	{\n	 \
 \n	  my ($rna,$pdb);\n	  $rnaSS=$rnaS;\n	  $pdbSS\
=$pdbS;\n	  $rnaSS=~s/T/U/g;\n	  $pdbSS=~s/T/U/g;\\
n	  ($rnaSS,$pdbSS)=nw ($rnaS, $pdbS);\n	  \n	  my\
 @rnaA =split (//,$rnaSS);\n	  my @pdbA=split (//,\
$pdbSS);\n	  my $l=@rnaA;\n	  \n	  #print \"\\n---\
 $s_seq{0}{name} $rnaSS\\n--- $s_seq{0}{name} $pdb\
SS\\n\\n\";\n	  \n	  for (my $b=0,my $a=0; $a<$l; \
$a++)\n	    {\n	      if   ($rnaA[$a] ne '-' && $p\
dbA[$a] ne '-'){$lu[++$pdb]=++$rna;}\n	      elsif\
($rnaA[$a] eq '-'){$lu[++$pdb]=-1;}\n	      elsif(\
$pdbA[$a] eq '-'){++$rna;}\n	    }\n	}\n      else\
\n	{\n	  for (my $a=0; $a<=length ($rnaS); $a++)\n\
	    {\n	      $lu[$a]=$a;\n	    }\n	}\n      my $\
l=length ($rnaS);\n      open ($f, \">$lib\");\n  \
    print $f \"! TC_LIB_FORMAT_01\\n\";\n      pri\
nt $f \"1\\n\";\n      print $f \"$s_seq{0}{name} \
$l $rnaS\\n\";\n      print $f \"!CMT: [SOURCE] >$\
s_seq{0}{name} 3D contact library Generated by $pg\
 (x3dna)\\n\";\n      print $f \"#1 1\\n\";\n     \
 \n      my $ne;\n      my @array=file2array($outf\
ile2);\n      for (my $a=0; $a<5; $a++){shift (@ar\
ray);}\n      while (!($array[0]=~/####/))\n	{\n	 \
 my $line= shift (@array);\n	  my @l=($line=~/(\\d\
+)/g);\n	  \n	 \n	  my $f1=$lu[$l[0]];\n	  my $s1=\
$lu[$l[1]];\n\n	  #print \"\\n$line\\n$l[0] --> $f\
1\\n$l[1] --> $s1\\n\\n\"; \n	  \n	  if (!$f1 || !\
$s1)\n	    {\n	      print \"\\n---- $rnaSS\\n----\
 $pdbSS\\n$line\\n[$l[0] --- $l[1]]<---->[$f1 --- \
$s1]\\n\";\n	      myexit(flush_error(\"Error whil\
e parsing s3dna::find_pair output\"));\n	    }\n	 \
 elsif ($f1==-1 || $s1==-1){;}\n	  elsif ($f1<$s1)\
{print $f \"$f1 $s1 100\\n\";}\n	  else {print $f \
\"$s1 $f1 100\\n\";$ne++;}\n	}\n      print $f \"!\
 SEQ_1_TO_N\\n\";\n      close ($f);\n      return\
;\n    }\nsub RNAplfold2lib\n  {\n    my ($seq, $l\
ib)=@_;\n    my $f= new FileHandle;\n    \n    &sa\
fe_system (\"t_coffee -other_pg RNAplfold2tclib.pl\
 -in=$seq -out=$lib\");\n    \n    if ( !-e $lib)\\
n	{\n	 myexit(flush_error(\"RNAplfold failed to co\
mpute the secondary structure of $s{$seq}{name}\")\
);\n	 myexit ($EXIT_FAILURE);\n       }\n    open \
($f, \">>$lib\");\n    print $f \"!CMT: [SOURCE] 2\
D contact library Generated by RNAPlfold (Vienna P\
ackage)\\n\";\n    close $f;\n    return;\n  }\nsu\
b x3dnassr2lib\n    {\n      my ($seq, $pdb, $lib)\
=@_;\n      my $outfile=\"dssr-2ndstrs.dbn\";\n   \
   my $f= new FileHandle;\n      \n\n      if (!pg\
_is_installed (\"x3dna-ssr\"))\n	{\n	  add_warning\
 ($$,$$, \"x3dna-ssr could not be used to extract \
RNA secondary structures. Secondary structures wil\
l be predicted ab-initio instead with RNAPlfold\")\
;\n	  return RNAplfold2lib ($seq,$lib);\n	}\n     \
 \n      safe_system (\"x3dna-ssr -i=$pdb >/dev/nu\
ll 2>/dev/null\");\n      if ( !-e $outfile)\n	{\n\
	  myexit(flush_error(\"x3dna-ssr failed to comput\
e the secondary structure file \"));\n	  myexit ($\
EXIT_FAILURE);\n	}\n\n      #Handle situations whe\
n the pdb sequence differs from the RNA sequence\n\
      @out=file2array($outfile);\n      my %s=read\
_fasta_seq ($seq);\n      my @names=keys (%s);\n  \
    my $rnaS=uc($s{$names[0]}{seq});\n      my $pd\
bS=uc($out[1]);\n      my $vienna;\n      \n      \
#x3dna returns non legitimate nucleotides\n       \
$pdbS=~s/[^AGCTU]//g;\n      \n      if ($rnaS ne \
$pdbS)\n	{\n	  my ($rna,$pdb);\n	  my $rnaSS=$rnaS\
;\n	  my $pdbSS=$pdbS;\n	  $rnaSS=~s/T/U/g;\n	  $p\
dbSS=~s/T/U/g;\n	  ($rnaSS,$pdbSS)=nw ($rnaSS, $pd\
bSS);\n	  my @rnaA =split (//,$rnaSS);\n	  my @pdb\
A=split (//,$pdbSS);\n	  my @SS=split (//, $out[2]\
);\n	  \n	  my $l=@rnaA;\n	  for (my $b=0,my $a=0;\
 $a<$l; $a++)\n	    {\n	      if   ($rnaA[$a] ne '\
-' && $pdbA[$a] ne '-'){$vienna.=$SS[$b++];}\n	   \
   elsif($rnaA[$a] eq '-'){$b++;}\n	      elsif($p\
dbA[$a] eq '-'){$vienna.='.';}\n	    }\n	}\n      \
else\n	{\n	  $vienna=$out[2];\n	}\n    \n\n      o\
pen ($f, \">seq\");\n      print $f \">$names[0]\\\
n$rnaS\\n\";\n      close $f;\n      \n      open \
($f, \">str\");\n      print $f \">$names[0]\\n$vi\
enna\\n\";\n      close $f;\n      \n      safe_sy\
stem (\"t_coffee -other_pg seq_reformat -in seq -i\
n2 str -output vienna2tc_lib >$lib\");\n      if (\
 !-e $lib)\n	    {\n	      myexit(flush_error(\"se\
q_reformat failed to convert your secondary struct\
ure\"));\n	      myexit ($EXIT_FAILURE);\n	    }\n\
      \n      open ($f, \">>$lib\");\n      print \
$f \"!CMT: [SOURCE] >$names[0] 2D Contact library \
generated by x3dna-ssr\\n\";\n      #print $f \"! \
Vienna_Format: >$names[0]\\n\";\n      #print $f \\
"! Vienna_Format: $vienna\\n\";\n      \n      clo\
se $f;\n      return;\n    }\n\n\nsub file2head\n \
     {\n	my $file = shift;\n	my $size = shift;\n	m\
y $f= new FileHandle;\n	my $line;\n	open ($f,$file\
);\n	read ($f,$line, $size);\n	close ($f);\n	retur\
n $line;\n      }\nsub file2tail\n      {\n	my $fi\
le = shift;\n	my $size = shift;\n	my $f= new FileH\
andle;\n	my $line;\n\n	open ($f,$file);\n	seek ($f\
,$size*-1, 2);\n	read ($f,$line, $size);\n	close (\
$f);\n	return $line;\n      }\n\n\nsub vtmpnam\n  \
    {\n	my $r=rand(100000);\n	my $f=\"file.$r.$$\"\
;\n	while (-e $f)\n	  {\n	    $f=vtmpnam();\n	  }\\
n	push (@TMPFILE_LIST, $f);\n	return $f;\n      }\\
n\nsub myexit\n  {\n    my $code=@_[0];\n    if ($\
CLEAN_EXIT_STARTED==1){return;}\n    else {$CLEAN_\
EXIT_STARTED=1;}\n    ### ONLY BARE EXIT\n    exit\
 ($code);\n  }\nsub set_error_lock\n    {\n      m\
y $name = shift;\n      my $pid=$$;\n\n\n      &lo\
ck4tc ($$,\"LERROR\", \"LSET\", \"$$ -- ERROR: $na\
me $PROGRAM\\n\");\n      return;\n    }\nsub set_\
lock\n  {\n    my $pid=shift;\n    my $msg= shift;\
\n    my $p=getppid();\n    &lock4tc ($pid,\"LLOCK\
\",\"LRESET\",\"$p$msg\\n\");\n  }\nsub unset_lock\
\n   {\n\n    my $pid=shift;\n    &lock4tc ($pid,\\
"LLOCK\",\"LRELEASE\",\"\");\n  }\nsub shift_lock\\
n  {\n    my $from=shift;\n    my $to=shift;\n    \
my $from_type=shift;\n    my $to_type=shift;\n    \
my $action=shift;\n    my $msg;\n\n    if (!&lock4\
tc($from, $from_type, \"LCHECK\", \"\")){return 0;\
}\n    $msg=&lock4tc ($from, $from_type, \"LREAD\"\
, \"\");\n    &lock4tc ($from, $from_type,\"LRELEA\
SE\", $msg);\n    &lock4tc ($to, $to_type, $action\
, $msg);\n    return;\n  }\nsub isshellpid\n  {\n \
   my $p=shift;\n    if (!lock4tc ($p, \"LLOCK\", \
\"LCHECK\")){return 0;}\n    else\n      {\n	my $c\
=lock4tc($p, \"LLOCK\", \"LREAD\");\n	if ( $c=~/-S\
HELL-/){return 1;}\n      }\n    return 0;\n  }\ns\
ub isrootpid\n  {\n    if(lock4tc (getppid(), \"LL\
OCK\", \"LCHECK\")){return 0;}\n    else {return 1\
;}\n  }\nsub lock4tc\n	{\n	  my ($pid,$type,$actio\
n,$value)=@_;\n	  my $fname;\n	  my $host=hostname\
;\n\n	  if ($type eq \"LLOCK\"){$fname=\"$LOCKDIR/\
.$pid.$host.lock4tcoffee\";}\n	  elsif ( $type eq \
\"LERROR\"){ $fname=\"$LOCKDIR/.$pid.$host.error4t\
coffee\";}\n	  elsif ( $type eq \"LWARNING\"){ $fn\
ame=\"$LOCKDIR/.$pid.$host.warning4tcoffee\";}\n\n\
	  if ($debug_lock)\n	    {\n	      print STDERR \\
"\\n\\t---lock4tc(tcg): $action => $fname =>$value\
 (RD: $LOCKDIR)\\n\";\n	    }\n\n	  if    ($action\
 eq \"LCHECK\") {return -e $fname;}\n	  elsif ($ac\
tion eq \"LREAD\"){return file2string($fname);}\n	\
  elsif ($action eq \"LSET\") {return string2file \
($value, $fname, \">>\");}\n	  elsif ($action eq \\
"LRESET\") {return string2file ($value, $fname, \"\
>\");}\n	  elsif ($action eq \"LRELEASE\")\n	    {\
\n	      if ( $debug_lock)\n		{\n		  my $g=new Fil\
eHandle;\n		  open ($g, \">>$fname\");\n		  print \
$g \"\\nDestroyed by $$\\n\";\n		  close ($g);\n		\
  safe_system (\"mv $fname $fname.old\");\n		}\n	 \
     else\n		{\n		  unlink ($fname);\n		}\n	    }\\
n	  return \"\";\n	}\n\nsub file2string\n	{\n	  my\
 $file=@_[0];\n	  my $f=new FileHandle;\n	  my $r;\
\n	  open ($f, \"$file\");\n	  while (<$f>){$r.=$_\
;}\n	  close ($f);\n	  return $r;\n	}\nsub file2ar\
ray\n	{\n	  my $file=@_[0];\n	  my $f=new FileHand\
le;\n	  my @r;\n	  open ($f, \"$file\");\n	  while\
 (<$f>){@r=(@r,$_);}\n	  close ($f);\n	  return @r\
;\n	}\nsub string2file\n    {\n    my ($s,$file,$m\
ode)=@_;\n    my $f=new FileHandle;\n\n    open ($\
f, \"$mode$file\");\n    print $f  \"$s\";\n    cl\
ose ($f);\n  }\n\nBEGIN\n    {\n      srand;\n\n  \
    $SIG{'SIGUP'}='signal_cleanup';\n      $SIG{'S\
IGINT'}='signal_cleanup';\n      $SIG{'SIGQUIT'}='\
signal_cleanup';\n      $SIG{'SIGILL'}='signal_cle\
anup';\n      $SIG{'SIGTRAP'}='signal_cleanup';\n \
     $SIG{'SIGABRT'}='signal_cleanup';\n      $SIG\
{'SIGEMT'}='signal_cleanup';\n      $SIG{'SIGFPE'}\
='signal_cleanup';\n\n      $SIG{'SIGKILL'}='signa\
l_cleanup';\n      $SIG{'SIGPIPE'}='signal_cleanup\
';\n      $SIG{'SIGSTOP'}='signal_cleanup';\n     \
 $SIG{'SIGTTIN'}='signal_cleanup';\n      $SIG{'SI\
GXFSZ'}='signal_cleanup';\n      $SIG{'SIGINFO'}='\
signal_cleanup';\n\n      $SIG{'SIGBUS'}='signal_c\
leanup';\n      $SIG{'SIGALRM'}='signal_cleanup';\\
n      $SIG{'SIGTSTP'}='signal_cleanup';\n      $S\
IG{'SIGTTOU'}='signal_cleanup';\n      $SIG{'SIGVT\
ALRM'}='signal_cleanup';\n      $SIG{'SIGUSR1'}='s\
ignal_cleanup';\n\n\n      $SIG{'SIGSEGV'}='signal\
_cleanup';\n      $SIG{'SIGTERM'}='signal_cleanup'\
;\n      $SIG{'SIGCONT'}='signal_cleanup';\n      \
$SIG{'SIGIO'}='signal_cleanup';\n      $SIG{'SIGPR\
OF'}='signal_cleanup';\n      $SIG{'SIGUSR2'}='sig\
nal_cleanup';\n\n      $SIG{'SIGSYS'}='signal_clea\
nup';\n      $SIG{'SIGURG'}='signal_cleanup';\n   \
   $SIG{'SIGCHLD'}='signal_cleanup';\n      $SIG{'\
SIGXCPU'}='signal_cleanup';\n      $SIG{'SIGWINCH'\
}='signal_cleanup';\n\n      $SIG{'INT'}='signal_c\
leanup';\n      $SIG{'TERM'}='signal_cleanup';\n  \
    $SIG{'KILL'}='signal_cleanup';\n      $SIG{'QU\
IT'}='signal_cleanup';\n\n      our $debug_lock=$E\
NV{\"DEBUG_LOCK\"};\n\n\n\n\n      foreach my $a (\
@ARGV){$CL.=\" $a\";}\n      if ( $debug_lock ){pr\
int STDERR \"\\n\\n\\n********** START PG: $PROGRA\
M *************\\n\";}\n      if ( $debug_lock ){p\
rint STDERR \"\\n\\n\\n**********(tcg) LOCKDIR: $L\
OCKDIR $$ *************\\n\";}\n      if ( $debug_\
lock ){print STDERR \"\\n --- $$ -- $CL\\n\";}\n\n\
\n\n\n    }\nsub flush_error\n  {\n    my $msg=shi\
ft;\n    return add_error ($EXIT_FAILURE,$$, $$,ge\
tppid(), $msg, $CL);\n  }\nsub add_error\n  {\n   \
 my $code=shift;\n    my $rpid=shift;\n    my $pid\
=shift;\n    my $ppid=shift;\n    my $type=shift;\\
n    my $com=shift;\n\n    $ERROR_DONE=1;\n    loc\
k4tc ($rpid, \"LERROR\",\"LSET\",\"$pid -- ERROR: \
$type\\n\");\n    lock4tc ($$, \"LERROR\",\"LSET\"\
, \"$pid -- COM: $com\\n\");\n    lock4tc ($$, \"L\
ERROR\",\"LSET\", \"$pid -- STACK: $ppid -> $pid\\\
n\");\n\n    return $code;\n  }\nsub add_warning\n\
  {\n    my $rpid=shift;\n    my $pid =shift;\n   \
 my $command=shift;\n    my $msg=\"$$ -- WARNING: \
$command\\n\";\n    print STDERR \"$msg\";\n    lo\
ck4tc ($$, \"LWARNING\", \"LSET\", $msg);\n  }\n\n\
sub signal_cleanup\n  {\n    print dtderr \"\\n***\
* $$ (tcg) was killed\\n\";\n    &cleanup;\n    ex\
it ($EXIT_FAILURE);\n  }\nsub clean_dir\n  {\n    \
my $dir=@_[0];\n    if ( !-d $dir){return ;}\n    \
elsif (!($dir=~/tmp/)){return ;}#safety check 1\n \
   elsif (($dir=~/\\*/)){return ;}#safety check 2\\
n    else\n      {\n	`rm -rf $dir`;\n      }\n    \
return;\n  }\nsub cleanup\n  {\n    #print stderr \
\"\\n----tc: $$ Kills $PIDCHILD\\n\";\n    #kill (\
SIGTERM,$PIDCHILD);\n    my $p=getppid();\n    $CL\
EAN_EXIT_STARTED=1;\n\n\n\n    if (&lock4tc($$,\"L\
ERROR\", \"LCHECK\", \"\"))\n      {\n	my $ppid=ge\
tppid();\n	if (!$ERROR_DONE)\n	  {\n	    &lock4tc(\
$$,\"LERROR\", \"LSET\", \"$$ -- STACK: $p -> $$\\\
n\");\n	    &lock4tc($$,\"LERROR\", \"LSET\", \"$$\
 -- COM: $CL\\n\");\n	  }\n      }\n    my $warnin\
g=&lock4tc($$, \"LWARNING\", \"LREAD\", \"\");\n  \
  my $error=&lock4tc($$,  \"LERROR\", \"LREAD\", \\
"\");\n    #release error and warning lock if root\
\n\n    if (isrootpid() && ($warning || $error) )\\
n      {\n\n	print STDERR \"**************** Summa\
ry *************\\n$error\\n$warning\\n\";\n\n	&lo\
ck4tc($$,\"LERROR\",\"RELEASE\",\"\");\n	&lock4tc(\
$$,\"LWARNING\",\"RELEASE\",\"\");\n      }\n\n\n \
   foreach my $f (@TMPFILE_LIST)\n      {\n	if (-e\
 $f){unlink ($f);}\n      }\n    foreach my $d (@T\
MPDIR_LIST)\n      {\n	clean_dir ($d);\n      }\n \
   #No More Lock Release\n    #&lock4tc($$,\"LLOCK\
\",\"LRELEASE\",\"\"); #release lock\n\n    if ( $\
debug_lock ){print STDERR \"\\n\\n\\n********** EN\
D PG: $PROGRAM ($$) *************\\n\";}\n    if (\
 $debug_lock ){print STDERR \"\\n\\n\\n**********(\
tcg) LOCKDIR: $LOCKDIR $$ *************\\n\";}\n  \
}\nEND\n  {\n\n    &cleanup();\n  }\n\nsub blast_c\
om2new_blast_com\n    {\n      my $com=shift;\n	  \
if ($com=~/formatdb/)\n	    {\n	      $com=~s/form\
atdb/makeblastdb/;\n	      $com=~s/\\-i/\\-in/;\n	\
      if ($com =~/pF/){$com=~s/\\-pF/\\-dbtype nuc\
l/;}\n	      if ($com =~/p F/){$com=~s/\\-p F/\\-d\
btype nucl/;}\n	      $com=\"$com -logfile /dev/nu\
ll\";\n	      return $com;\n	    }\n	  else {retur\
n $com;}\n\n    }\nsub safe_system\n{\n  my $com=s\
hift;\n  my $ntry=shift;\n  my $ctry=shift;\n  my \
$pid;\n  my $status;\n  my $ppid=getppid();\n  if \
($com eq \"\"){return 1;}\n\n  if ( ($com=~/^blast\
/) ||($com=~/^formatdb/)){$com=&blast_com2new_blas\
t_com($com);}\n\n  if (($pid = fork ()) < 0){retur\
n (-1);}\n  if ($pid == 0)\n    {\n      set_lock(\
$$, \" -SHELL- $com (tcg)\");\n      if( $debug_ge\
neric_method ) { printf \"~ exec: %s\\n\", $com; }\
\n      exec ($com);\n      if( $debug_generic_met\
hod ) { printf \"~ exitcode: %s\\n\", $?; }\n    }\
\n  else\n    {\n      lock4tc ($$, \"LLOCK\", \"L\
SET\", \"$pid\\n\");#update parent\n      $PIDCHIL\
D=$pid;\n    }\n  if ($debug_lock){printf STDERR \\
"\\n\\t .... safe_system (fasta_seq2hmm)  p: $$ c:\
 $pid COM: $com\\n\";}\n\n  waitpid ($pid,WTERMSIG\
);\n\n  shift_lock ($pid,$$, \"LWARNING\",\"LWARNI\
NG\", \"LSET\");\n\n  if ($? == $EXIT_FAILURE || l\
ock4tc($pid, \"LERROR\", \"LCHECK\", \"\"))\n    {\
\n      if ($ntry && $ctry <$ntry)\n	{\n\n	  add_w\
arning ($$,$$,\"$com failed [retry: $ctry out of $\
ntry]\");\n	  lock4tc ($pid, \"LRELEASE\", \"LERRO\
R\", \"\");\n	  #if ($com=~/EBI/){$com=~s/EBI/NCBI\
/;}\n	  #elsif ($com=~/NCBI/){$com=~s/NCBI/EBI/;}\\
n\n	  return safe_system ($com, $ntry, ++$ctry);\n\
	}\n      elsif ($ntry == -1)\n	{\n	  if (!shift_l\
ock ($pid, $$, \"LERROR\", \"LWARNING\", \"LSET\")\
)\n	    {\n	      add_warning ($$,$$,\"$com failed\
\");\n	    }\n	  else\n	    {\n	      lock4tc ($pi\
d, \"LRELEASE\", \"LERROR\", \"\");\n	    }\n	  re\
turn $?;}\n      else\n	{\n	  if (!shift_lock ($pi\
d,$$, \"LERROR\",\"LERROR\", \"LSET\"))\n	    {\n	\
      myexit(add_error ($EXIT_FAILURE,$$,$pid,getp\
pid(), \"UNSPECIFIED system\", $com));\n	    }\n	}\
\n    }\n  return $?;\n}\n\nsub check_configuratio\
n\n    {\n      my @l=@_;\n      my $v;\n      for\
each my $p (@l)\n	{\n\n	  if   ( $p eq \"EMAIL\")\\
n	    {\n	      if ( !($EMAIL=~/@/))\n		{\n		add_w\
arning($$,$$,\"Could Not Use EMAIL\");\n		myexit(a\
dd_error ($EXIT_FAILURE,$$,$$,getppid(),\"EMAIL\",\
\"$CL\"));\n	      }\n	    }\n	  elsif( $p eq \"IN\
TERNET\")\n	    {\n	      if ( !&check_internet_co\
nnection())\n		{\n		  myexit(add_error ($EXIT_FAIL\
URE,$$,$$,getppid(),\"INTERNET\",\"$CL\"));\n		}\n\
	    }\n	  elsif( $p eq \"wget\")\n	    {\n	      \
if (!&pg_is_installed (\"wget\") && !&pg_is_instal\
led (\"curl\"))\n		{\n		  myexit(add_error ($EXIT_\
FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:wget\",\
\"$CL\"));\n		}\n	    }\n	  elsif( !(&pg_is_instal\
led ($p)))\n	    {\n	      myexit(add_error ($EXIT\
_FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:$p\",\\
"$CL\"));\n	    }\n	}\n      return 1;\n    }\nsub\
 nw\n      {\n	my($A,$B)=@_;\n	my ($i,$j, $s);\n	m\
y $gep=-2;\n	my $match=+2;\n	my $mismatch=0;\n	my \
($sub, $ins, $del);\n\n\n	if ($A eq $B){return ($A\
,$B);}\n	\n	$A=~s/[\\s\\d]//g;	\n	$B=~s/[\\s\\d]//\
g;	\n\n\n	my @rA=split ('',$A);\n	my @rB=split (''\
,$B);\n	\n	#evaluate substitutions\n	my $lenA=@rA;\
\n	my $lenB=@rB;\n	\n	for ($i=0; $i<=$lenA; $i++){\
$smat[$i][0]=$i*$gep;$tb[$i][0 ]= 1;}\n	for ($j=0;\
 $j<=$lenB; $j++){$smat[0][$j]=$j*$gep;$tb[0 ][$j]\
=-1;}\n	\n	for ($i=1; $i<=$lenA; $i++)\n	  {\n	   \
 for ($j=1; $j<=$lenB; $j++)\n	      {\n		if ($rA[\
$i-1] eq $rB[$j-1]){$s=$match;}\n		else {$s=$misma\
tch;}\n		\n		$sub=$smat[$i-1][$j-1]+$s;\n		$del=$s\
mat[$i  ][$j-1]+$gep;\n		$ins=$smat[$i-1][$j  ]+$g\
ep;\n		\n		if   ($sub>=$del && $sub>=$ins){$smat[$\
i][$j]=$sub;$tb[$i][$j]=0;}\n		elsif($del>$ins){$s\
mat[$i][$j]=$del;$tb[$i][$j]=-1;}\n		else {$smat[$\
i][$j]=$ins;$tb[$i][$j]=1;}\n		}\n	  }\n	#print \"\
\\n---- SCORE=$smat[$lenA][$lenB]\\n\";\n	\n	$i=$l\
enA;\n	$j=$lenB;\n	my $aln_len=0;\n\n	while (!($i=\
=0 && $j==0))\n	  {\n	    if ($tb[$i][$j]==0)\n	  \
  {\n	      $aA[$aln_len]=$rA[--$i];\n	      $aB[$\
aln_len]=$rB[--$j];\n	    }\n	  elsif ($tb[$i][$j]\
==-1)\n	    {\n	      $aA[$aln_len]='-';\n	      $\
aB[$aln_len]=$rB[--$j];\n	    }\n	  elsif ($tb[$i]\
[$j]==1)\n	    {\n	      $aA[$aln_len]=$rA[--$i];\\
n	      $aB[$aln_len]='-';\n	    }\n	  $aln_len++;\
\n	  }\n	\n	\n	@aA=reverse (@aA);\n	@aB=reverse (@\
aB);\n	my $sA=join('',@aA);\n	my $sB=join('',@aB);\
\n	return ($sA,$sB);\n      }\n      \nsub fasta2n\
seq\n	{\n	  \n	  my $f=@_[0];\n	  my $nseq;\n\n	  \
open (F, \"$f\") or return 0;\n	  while (<F>)\n	  \
  {\n	      if ($_=~/\\>/){$nseq++;}\n	    }\n	  c\
lose (F);\n	  return $nseq;\n	}\n	\n$program=\"T-C\
OFFEE (dev_shivashamloo@20190813_14:24)\";\n\n","u\
se Env;\nuse FileHandle;\nuse Cwd;\nuse File::Path\
;\nuse Sys::Hostname;\nmy $f = new FileHandle;\n\n\
open ($f, $ARGV[1]);\n$atom=$ARGV[0];\n\n$atom=~s/\
PRIME/\\'/;\nwhile (<$f>)\n  {\n    my $l=$_;\n\n \
   $l=~s/$atom/CA /;\n    \n    \n    $l=~s/  G /G\
LY /g;\n    $l=~s/  C /CYS /g;\n    $l=~s/  T /THR\
 /g;\n    $l=~s/  A /ALA /g;\n    $l=~s/  U /THR /\
g;\n    \n    $l=~s/ DG /GLY /g;\n    $l=~s/ DC /C\
YS /g;\n    $l=~s/ DT /THR /g;\n    $l=~s/ DA /ALA\
 /g;\n    $l=~s/ DU /THR /g;\n    \n    print $l;\\
n  }\n\n\n\n","*TC_METHOD_FORMAT_01\n*************\
*****generic_method.tc_method*************\n*\n*  \
     Incorporating new methods in T-Coffee\n*     \
  Cedric Notredame 26/08/08\n*\n******************\
*************************************\n*This file \
is a method file\n*Copy it and adapt it to your ne\
ed so that the method \n*you want to use can be in\
corporated within T-Coffee\n**********************\
*********************************\n*              \
    USAGE                              *\n********\
***********************************************\n*\
This file is passed to t_coffee via -in:\n*\n*	t_c\
offee -in Mgeneric_method.method\n*\n*	The method \
is passed to the shell using the following\n*call:\
\n*<EXECUTABLE><PARAM1><IN_FLAG><seq_file><PARAM2>\
<OUT_FLAG><outname><PARAM>\n*\n*Conventions:\n*<FL\
AG_NAME> 	<TYPE>		<VALUE>\n*<VALUE>:	no_name 	<=> \
Replaced with a space\n*<VALUE>:	&nbsp	<=> Replace\
d with a space\n*\n*******************************\
************************\n*                  ALN_M\
ODE                           *\n*****************\
**************************************\n*pairwise \
  ->all Vs all (no self )[(n2-n)/2aln]\n*m_pairwis\
e ->all Vs all (no self)[n^2-n]^2\n*s_pairwise ->a\
ll Vs all (self): [n^2-n]/2 + n\n*multiple   ->All\
 the sequences in one go\n*\nALN_MODE		pairwise\n*\
\n************************************************\
*******\n*                  OUT_MODE              \
             *\n**********************************\
*********************\n* mode for the output:\n*Ex\
ternal methods: \n* aln -> alignmnent File (Fasta \
or ClustalW Format)\n* lib-> Lib file (TC_LIB_FORM\
AT_01)\n*Internal Methods:\n* fL -> Internal Funct\
ion returning a List (Librairie)\n* fA -> Internal\
 Function returning an Alignmnent\n*\nOUT_MODE		al\
n\n***********************************************\
********\n*                  SEQ_TYPE             \
              *\n*********************************\
**********************\n*G: Genomic, S: Sequence, \
P: PDB, R: Profile\n*Examples:\n*SEQTYPE	S	sequenc\
es against sequences (default)\n*SEQTYPE	S_P	seque\
nce against structure\n*SEQTYPE	P_P	structure agai\
nst structure\n*SEQTYPE	PS	mix of sequences and st\
ructure	\n*\nSEQ_TYPE	S\n*\n\n********************\
***********************************\n*            \
    COMMAND LINE                         *\n*EXECU\
TABLE PARAM1 IN_FLAG OUT_FLAG PARAM             *\\
n*************************************************\
******\n******************************************\
*************\n*                  EXECUTABLE      \
                   *\n****************************\
***************************\n*name of the executab\
le\n*passed to the shell: executable\n*	\nEXECUTAB\
LE	tc_generic_method.pl\n*\n**********************\
*********************************\n*              \
    IN_FLAG                             *\n*******\
************************************************\n\
*IN_FLAG\n*flag indicating the name of the in comi\
ng sequences\n*IN_FLAG S no_name ->no flag\n*IN_FL\
AG S &bnsp-in&bnsp -> \" -in \"\n*\nIN_FLAG		-infi\
le=\n*\n******************************************\
*************\n*                  OUT_FLAG        \
                   *\n****************************\
***************************\n*OUT_FLAG\n*flag indi\
cating the name of the out-coming data\n*same conv\
entions as IN_FLAG\n*OUT_FLAG	S no_name ->no flag\\
n*if you want to redirect, pass the parameters via\
 PARAM1\n*set OUT_FLAG to >\n*\nOUT_FLAG		-outfile\
=\n*\n********************************************\
***********\n*                  PARAM_1           \
                   *\n****************************\
***************************\n*<EXECUTABLE><PARAM1>\
<IN_FLAG><seq_file><PARAM2><OUT_FLAG><outname><PAR\
AM>\n*Parameters sent to the EXECUTABLE and specif\
ied *before* IN_FLAG \n*If there is more than 1 PA\
RAM line, the lines are\n*concatenated\n*Command_l\
ine: @EP@PARAM@-gapopen%e10%s-gapext%e20\n*	%s whi\
te space\n*	%e equal sign\n*\n*PARAM1	\n*\n*\n*\n*\
**************************************************\
****\n*                  PARAM_2                  \
            *\n***********************************\
********************\n*<EXECUTABLE><PARAM1><IN_FLA\
G><seq_file><PARAM2><OUT_FLAG><outname><PARAM>\n*P\
arameters sent to the EXECUTABLE and specified \n*\
after* IN_FLAG and *before* OUT_FLAG\n*If there is\
 more than 1 PARAM line, the lines are\n*concatena\
ted\n*\n*PARAM1	\n*\n*\n**************************\
*****************************\n*                  \
PARAM                              *\n************\
*******************************************\n*<EXE\
CUTABLE><PARAM1><IN_FLAG><seq_file><PARAM2><OUT_FL\
AG><outname><PARAM>\n*Parameters sent to the EXECU\
TABLE and specified *after* OUT_FLAG\n*If there is\
 more than 1 PARAM line, the lines are\n*concatena\
ted\n*\nPARAM	-mode=seq_msa -method=clustalw\nPARA\
M   -OUTORDER=INPUT -NEWTREE=core -align -gapopen=\
-15\n*\n******************************************\
*************\n*                  END             \
                   *\n****************************\
***************************\n","*TC_METHOD_FORMAT_\
01\n***************clustalw_method.tc_method******\
***\nEXECUTABLE	clustalw\nALN_MODE		pairwise\nIN_F\
LAG		-INFILE=\nOUT_FLAG		-OUTFILE=\nOUT_MODE		aln\\
nPARAM		-gapopen=-10\nSEQ_TYPE		S\n***************\
**********************************\n","$VersionTag\
 =                                                \
                                                  \
                                 2.43;\nuse Env;\n\
use FileHandle;\nuse Cwd;\nuse File::Path;\nuse Sy\
s::Hostname;\n\nour $PIDCHILD;\nour $ERROR_DONE;\n\
our @TMPFILE_LIST;\nour $EXIT_FAILURE=1;\nour $EXI\
T_SUCCESS=0;\n\nour $REFDIR=getcwd;\nour $EXIT_SUC\
CESS=0;\nour $EXIT_FAILURE=1;\n\nour $PROGRAM=\"ex\
tract_from_pdb\";\nour $CL=$PROGRAM;\n\nour $CLEAN\
_EXIT_STARTED;\nour $debug_lock=$ENV{\"DEBUG_LOCK\\
"};\nour $LOCKDIR=$ENV{\"LOCKDIR_4_TCOFFEE\"};\nif\
 (!$LOCKDIR){$LOCKDIR=getcwd();}\nour $ERRORDIR=$E\
NV{\"ERRORDIR_4_TCOFFEE\"};\nour $ERRORFILE=$ENV{\\
"ERRORFILE_4_TCOFFEE\"};\n&set_lock ($$);\nif (iss\
hellpid(getppid())){lock4tc(getppid(), \"LLOCK\", \
\"LSET\", \"$$\\n\");}\n      \nour $SILENT=\" >/d\
ev/null 2>/dev/null\";\nour $INTERNET=-1;\n\n\n\n\\
n\n\n\nour $BLAST_MAX_NRUNS=2;\nour $EXIT_SUCCESS=\
0;\nour $EXIT_FAILURE=1;\nour $CONFIGURATION=-1;\n\
our $REF_EMAIL=\"\";\nour $PROGRAM=\"extract_from_\
pdb\";\n\n\nmy %onelett_prot=&fill_onelett_prot();\
\nmy %threelett_prot=&fill_threelett_prot();\nmy %\
onelett_RNA=&fill_onelett_RNA();\nmy %threelett_RN\
A=&fill_threelett_RNA();\nmy %onelett_DNA=&fill_on\
elett_DNA();\nmy %threelett_DNA=&fill_threelett_DN\
A();\n\n\n\n\n\nmy %onelett = (\n'P' => \\%onelett\
_prot,\n'D' => \\%onelett_DNA,\n'R' => \\%onelett_\
RNA\n);\n\n\nmy %threelett = (\n'P' => \\%threelet\
t_prot,\n'D' => \\%threelett_DNA,\n'R' => \\%three\
lett_RNA\n);\n\n\n\n\n\n\n\nif($ARGV[0]=~/help/ ||\
$ARGV[0]=~/man/ || $ARGV[0]=~/HELP/ || $ARGV[0]=~/\
Man/ || $ARGV[0] eq \"-h\"  || $ARGV[0] eq \"-H\" \
 )\n{die \"SYNTAX: extract_from_pdb Version $Versi\
onTag	\n	Minimum:             [extract_from_pdb fi\
le] \n			   OR \n			     [... | extract_from_pdb]\\
n 	Flags (Default setting on the first line)\n	   \
-version...................[Returns the Version Nu\
mber]\n           -force.....................[Forc\
es the file to be treated like a PDB file]\n      \
                                [Regenerates the h\
eader and SEQRES fields]\n           -force_name..\
..............[Forces the file to be named after n\
ame]]\n           -infile.....file...........[Flag\
 can be omited]\n			              [File must be pd\
b or fro pgm]\n                                   \
   [File can also be compressed Z or gz]\n        \
                              [In the case of a co\
mpressed file, you can omit the gz|Z extension]\n \
          -netfile...................[File will be\
 fetch from the net using wget]\n                 \
                     [wget or curl must be install\
ed]\n                                      [ftp://\
ftp.gnu.org/pub/gnu/wget/]\n                      \
                [http://curl.haxx.se/]\n          \
                            [Must also be used to \
retrieve the file from a local pdb copy (cf netadd\
ress)]\n           -netaddress................[Add\
ress used for the retrieving the netfile]\n       \
                               [http://www.rcsb.or\
g/pdb/cgi/export.cgi/%%.pdb.gz?format=PDB&pdbId=%%\
&compression=gz]\n                                \
      [http://www.expasy.ch/cgi-bin/get-pdb-entry.\
pl?%%]\n                                      [loc\
al -> will get the file from pdb_dir (see pdb_dir)\
]\n           -netcompression............[Extensio\
n if the netfile comes compressed]\n              \
                        [gz]\n           -pdb_dir.\
..................[address of the repertory where \
the pdb is installed]\n                           \
           [Supports standard ftp style installati\
on OR every stru in DIR]\n                        \
              [Give the ..../pdb/structure/ dir]\n\
                                      [If value om\
itted, the pg gets it from the env variable PDB_DI\
R]\n           -netcompression_pg.........[gunzip]\
\n           -is_pdb_name..........name.[Returns 1\
 if the name is a PDB ID, 0 otherwise]\n          \
 -model_type...........name.[Returns the model typ\
e if valid PDB name]\n           -is_released_pdb_\
name name.[Returns 1 if the name corresponds to a \
released PDB file]\n           -get_pdb_chains....\
.name...[Returns the list of chains corresponding \
to the entry]\n           -get_pdb_id.........name\
...[Returns the PDB id within the provided pdb fil\
e]\n           -get_fugue_name.....name...[Turns a\
 name into a name valid for fugue]\n              \
                        [Uses the netaddress to do\
 so]\n	   -chain......FIRST..........[Extract the \
first chain only]\n		       A B C..........[Extrac\
t Several chains if needed]\n		       ALL.........\
...[Extract all the chains]	\n           -ligand..\
...ALL............[Extract the ligands in the chai\
n (HETATM)]\n                       <name1>,<name2\
>[Extract All the named lignds]\n	   -ligand_only.\
..............[Extract only the ligands]\n        \
   -ligand_list...............[Extract the list of\
 ligands]\n	   -coor.......<start>..<end>.[Coordin\
ates of the fragment to extract]\n			             \
 [Omit end to include the Cter]\n           -num..\
......absolute.......[absolute: relative to the se\
q] \n                       file...........[file: \
relative to file]\n           -num_out....new.....\
.......[new: start 1->L]\n                       o\
ld............[old: keep the file coordinates]\n  \
         -delete.....<start>..<end>.[Delete from r\
esidue start to residue end]\n	   -atom.......CA..\
...........[Atoms to include, ALL for all of them]\
\n		       CA O N.........[Indicate several atoms \
if needed]\n	   -code.......3..............[Use th\
e 1 letter code or the 3 letters code]\n	   -mode.\
......raw............[Output original pdb file]\n \
                      pdb............[Output somet\
hing that looks like pdb]\n		       fasta.........\
.[Output the sequences in fasta format]\n		       \
simple.........[Output a format easy to parse in C\
 ]\n            -seq_field..ATOM...........[Field \
used to extract the sequence]\n		       SEQRES....\
.....[Use the complete sequence]\n	   -seq........\
...............[Equivalent to  -mode fasta]\n	   -\
model......1..............[Chosen Model in an NMR \
file]\n           -nodiagnostic..............[Swit\
ches Error Messages off]\n           -debug.......\
..............[Sets the DEBUG ON]\n           -no_\
remote_pdb_dir.........[Do not look for a remote f\
ile]\n           -cache_pdb.................[Cache\
 Value, default is $HOME/.t_coffee/cache, other va\
lues: NO<=> No cache]\n\n      Environement Variab\
les\n           These variables can be set from th\
e environement\n           Command line values wit\
h the corresponding flag superseed evironement val\
ue\n           NO_REMOTE_PDB_DIR..........[Prevent\
s the program from searching remote file: faster]\\
n           PDB_DIR....................[Indicates \
where PDB file must be fetched (localy)]\n\n	 PROB\
LEMS: please contact cedric.notredame\\@europe.com\
\\n\";\n	 exit ($EXIT_SUCCESS);\n}\n\n$np=0;\n$n_p\
ara=$#ARGV;\n$model=1;\n$pdb_dir=$ENV{'PDB_DIR'};i\
f ($pdb_dir){$pdb_dir.=\"/\";}\n$debug=$ENV{'DEBUG\
_EXTRACT_FROM_PDB'};\n\n$no_remote_pdb_dir=$ENV{NO\
_REMOTE_PDB_DIR};\n$HOME=$ENV{'HOME'};\nif ( $ENV{\
CACHE_4_TCOFFEE})\n{$cache=$ENV{CACHE_4_TCOFFEE};}\
\nelse\n{\n    $cache=\"$HOME/.t_coffee/cache/\";\\
n}\n\n   \n$netaddress=\"http://www.rcsb.org/pdb/f\
iles/%%.pdb.gz\";\n$netcompression_pg=\"gunzip\";\\
n$netcompression=\"gz\";\n\nforeach ($np=0; $np<=$\
n_para; $np++)\n  {        \n    $value=$ARGV[$np]\
;\n   \n    if  ($np==0 && !($value=~/^-.*/))\n   \
   { \n	$pdb_file= $ARGV[$np];\n      }\n    elsif\
 ( !($value=~/^-.*/))\n      {\n	print \"@ARGV\";\\
n	die;\n      } \n    \n    elsif ($value eq \"-no\
diagnostic\"){$nodiagnostic=1;}\n    elsif ($value\
 eq \"-force\")\n      {\n	$force_pdb=1;\n      }\\
n    elsif ($value eq \"-force_name\")\n      {\n	\
$force_name=$ARGV[++$np];\n	$force_pdb=1;\n      }\
\n    \n    elsif ($value eq \"-is_pdb_name\")\n  \
    {\n	$pdb_file= $ARGV[++$np];	\n	$is_pdb_name=1\
;	\n      } \n    elsif ($value eq \"-is_released_\
pdb_name\")\n      {\n	$pdb_file= $ARGV[++$np];	\n\
	$is_released_pdb_name=1;\n      }\n    elsif ($va\
lue eq \"-model_type\")\n      {\n	$pdb_file= $ARG\
V[++$np];	\n	$model_type=1;\n      }\n    elsif ($\
value eq \"-debug\")\n{\n	$debug=1;\n}\n    elsif \
($value eq \"-get_pdb_chains\")\n{\n	$pdb_file= $A\
RGV[++$np];\n	$get_pdb_chains=1;\n}\n    elsif ($v\
alue eq \"-get_pdb_ligands\")\n{\n	$get_pdb_ligand\
s=1;\n}\n    \n    elsif ($value eq \"-get_pdb_id\\
")\n{\n	$pdb_file= $ARGV[++$np];\n	$get_pdb_id=1;\\
n	\n}\n    \n    elsif ( $value eq \"-get_fugue_na\
me\")\n{\n	$pdb_file= $ARGV[++$np];\n	$get_fugue_n\
ame=1;\n}\n    elsif ( $value eq \"-infile\")\n{\n\
       $pdb_file= $ARGV[++$np];\n} \n    elsif ($v\
alue eq \"-netfile\")\n{\n	$netfile=1;\n	if ( !($A\
RGV[$np+1]=~/^-.*/)){$pdb_file= $ARGV[++$np];}\n}\\
n    elsif (  $value eq \"-num\")\n{\n       $numb\
ering= $ARGV[++$np];\n}\n    elsif (  $value eq \"\
-num_out\")\n{\n       $numbering_out= $ARGV[++$np\
];\n}\n    elsif ( $value eq \"-netaddress\")\n{\n\
	$netadress=$ARGV[++$np];\n}\n     \n    elsif ( $\
value eq \"-netcompression\")\n{\n	 $netcompressio\
n=$ARGV[++$np];\n}\n    elsif ( $value eq \"-pdb_d\
ir\")\n{\n	 if ( !($ARGV[$np+1]=~/^-.*/)){$pdb_dir\
= \"$ARGV[++$np]/\";}\n}\n     elsif ( $value eq \\
"-no_remote_pdb_dir\")\n{\n	$no_remote_pdb_dir=1;\\
n	if ( !($ARGV[$np+1]=~/^-.*/)){$pdb_dir= \"$ARGV[\
++$np]/\";}\n}\n    elsif ( $value eq \"-cache\")\\
n{\n	$cache=$ARGV[++$np];\n}\n    \n    elsif ($va\
lue eq \"-netcompression_pg\")\n{\n	  $netcompress\
ion_pg=$ARGV[++$np];\n}\n     elsif ($value eq \"-\
mode\")\n{\n       $MODE=$ARGV[++$np];\n}\n\n    e\
lsif ( $value eq \"-model\")\n{\n       $model= $A\
RGV[++$np];\n}\n    elsif ($value eq \"-seq_field\\
" )\n{\n       $seq_field= $ARGV[++$np];\n}   \n  \
  elsif ($value eq \"-coor\" )\n{\n       $start= \
$ARGV[++$np];\n  \n       if (($ARGV[$np+1] eq \"\\
") ||($ARGV[$np+1]=~/^-.*/)){$end=\"*\";} \n      \
 else {$end=   $ARGV[++$np];}     \n       $coor_s\
et=1;\n}\n    elsif ($value eq \"-delete\" )\n{\n \
      $delete_start= $ARGV[++$np];\n       $delete\
_end= $ARGV[++$np];\n       $delete_set=1;\n}\n   \
 elsif  ($value eq \"-code\")\n{\n       $code= $A\
RGV[++$np];\n}\n    elsif  ($value eq \"-no_hetatm\
\")\n{\n       $no_hetatm=1;\n}\n    elsif ($value\
 eq \"-chain\")\n{\n       while (!($ARGV[$np+1] e\
q \"\") &&!($ARGV[$np+1]=~/^-.*/))\n{\n	      ++$n\
p;\n	      @c_chain=(@chain,  $ARGV[$np]);\n	     \
 $hc_chain{$ARGV[$np]}=$#c_chain+1;\n}           \\
n}\n    elsif ($value eq \"-atom\")\n{\n\n       w\
hile (!($ARGV[$np+1] eq \"\") && !($ARGV[$np+1]=~/\
^-.*/))\n{\n	      ++$np;\n	      $atom[$n_atom++]\
=  $ARGV[$np];\n	      $atom_list{$ARGV[$np]}=1;	 \
     \n} \n       \n}\n    elsif ( $value eq \"-un\
fold\")\n{\n	$unfold=1;\n}\n    elsif ($value eq \\
"-seq\" ||$value eq \"-fasta\" )\n{\n       $MODE=\
\"fasta\";\n}\n    elsif ( $value eq \"-version\")\
\n{\n	print STDERR  \"\\nextract_from_pdb: Version\
 $VersionTag\\n\";\n	&myexit ($EXIT_SUCCESS);\n}\n\
    elsif ( $value eq \"-ligand\")\n{\n	while (!($\
ARGV[$np+1] eq \"\") && !($ARGV[$np+1]=~/^-.*/))\n\
{\n	    ++$np;\n	    $ligand=1;\n	    $ligand_list\
{$ARGV[$np]}=1;	      \n} \n	$hc_chain{'LIGAND'}=1\
;\n}\n    elsif ( $value eq \"-ligand_only\")\n{\n\
	$ligand_only=1;\n}\n}\nif ( $debug)\n{\n    print\
 STDERR \"\\n[DEBUG:extract_from_pdb] NO_REMOTE_PD\
B_DIR: $no_remote_pdb_dir\\n\";\n    print STDERR \
\"\\n[DEBUG:extract_from_pdb] PDB_DIR: $pdb_dir\\n\
\";\n}\n\n\nif ( $is_pdb_name)\n  {\n    if (&remo\
te_is_pdb_name($pdb_file))\n      {\n	print \"1\";\
\n      }\n    else\n      {\n	print \"0\";\n     \
 }\n    exit ($EXIT_SUCCESS);\n  }\n\nif ( $is_rel\
eased_pdb_name)\n  {\n    \n    if (&is_released($\
pdb_file))\n      {\n	print \"1\";\n      }\n    e\
lse\n      {\n	print \"0\";\n      }\n    exit ($E\
XIT_SUCCESS);\n  }\nif ($model_type)\n  {\n   \n  \
  printf \"%s\", &pdb2model_type($pdb_file);\n    \
exit ($EXIT_SUCCESS);\n    \n  }\n    \n\nif (!$fo\
rce_name)\n{\n    $pdb_file=~/([^\\/]*)$/;\n    $f\
orce_name=$1;\n}\n\n$local_pdb_file=$pdb_file;\n\n\
if ( $debug){print STDERR \"\\n[DEBUG: extract_fro\
m_pdb] Scan For $local_pdb_file\\n\";}\n\n$mem=$no\
_remote_pdb_dir;\n$no_remote_pdb_dir=1;\n$tmp_pdb_\
file=get_pdb_file ($local_pdb_file);\n\nif ( !-e $\
tmp_pdb_file || $tmp_pdb_file eq \"\")\n  {\n    $\
local_pdb_file=$pdb_file;\n    ($local_pdb_file, $\
suffix_chain)=&pdb_name2name_and_chain($local_pdb_\
file);\n\n    if ($local_pdb_file)\n      {\n	if (\
 $debug){print STDERR \"\\nSplit $pdb_file into $l\
ocal_pdb_file and $suffix_chain \\n\";}\n	$tmp_pdb\
_file=get_pdb_file ($local_pdb_file);\n	if ( $tmp_\
pdb_file ne \"\")\n	  {\n	    @c_chain=();\n	    @\
c_chain=($suffix_chain);\n	    %hc_chain=();\n	   \
 $hc_chain{$suffix_chain}=1;\n	  }\n      }\n  }\n\
\n$no_remote_pdb_dir=$mem;\nif ($no_remote_pdb_dir\
==0)\n  {\n    \n    if ( !-e $tmp_pdb_file || $tm\
p_pdb_file eq \"\")\n      {\n	\n	$local_pdb_file=\
$pdb_file;\n	($local_pdb_file, $suffix_chain)=&pdb\
_name2name_and_chain($local_pdb_file);\n	if ($loca\
l_pdb_file)\n	  {\n	    \n	    if ( $debug){print \
STDERR \"\\nSplit $pdb_file into $local_pdb_file a\
nd $suffix_chain \\n\";}\n	    \n	    $tmp_pdb_fil\
e=get_pdb_file ($local_pdb_file);    \n	    \n	   \
 if ( $tmp_pdb_file ne \"\")\n	      {\n		@c_chain\
=();\n		@c_chain=($suffix_chain);\n		%hc_chain=();\
\n		$hc_chain{$suffix_chain}=1;\n	      }\n	  }\n \
     }\n  }\n\nif ( $debug){print STDERR \"\\n$pdb\
_file copied into ##$tmp_pdb_file##\\n\";}\n\nif (\
 !-e $tmp_pdb_file || $tmp_pdb_file eq \"\")\n{\n	\
if ($is_pdb_name)\n{\n	    print \"0\\n\"; exit ($\
EXIT_SUCCESS);\n}\n	else\n{\n  \n	    print \"\\nE\
XTRACT_FROM_PDB: NO RESULT for $pdb_file\\n\";\n	 \
   &myexit ($EXIT_SUCCESS);	\n}\n}\n\n\n\n\n%molec\
ule_type=&pdbfile2chaintype($tmp_pdb_file);\nif ( \
$debug){print STDERR \"\\n\\tSequence Type determi\
ned\\n\";}\n\n$pdb_id=&get_pdb_id ($tmp_pdb_file);\
\nif ( $debug){print STDERR \"\\n\\tID: $pdb_id (f\
or $tmp_pdb_file)\\n\";}\n\nif ( $pdb_id eq \"\"){\
$pdb_id=$force_name;}\n\n@f_chain=&get_chain_list \
($tmp_pdb_file);\nif ( $debug){print STDERR \"\\n\\
\tChain_list:@f_chain\\n\";}\n\nif ( $get_pdb_chai\
ns)\n{\n    print \"@f_chain\\n\";\n    &myexit ($\
EXIT_SUCCESS);\n}\nif ( $get_pdb_ligands)\n{\n    \
%complete_ligand_list=&get_ligand_list ($tmp_pdb_f\
ile);\n    print $complete_ligand_list{\"result\"}\
;\n    &myexit ($EXIT_SUCCESS);\n}\n\nelsif ( $get\
_pdb_id ||$get_fugue_name )\n{\n    if    (@c_chai\
n && $c_chain[0] eq \"FIRST\"){$pdb_id=$pdb_id.$f_\
chain[0];}\n    elsif (@c_chain && $c_chain[0] ne \
\" \"){$pdb_id=$pdb_id.$c_chain[0];}\n    \n    pr\
int \"$pdb_id\\n\";\n    &myexit ($EXIT_SUCCESS);\\
n    \n}\nelsif ( $is_pdb_name)\n{\n    printf \"1\
\\n\";\n    &myexit ($EXIT_SUCCESS);\n}\n\n\n\n$st\
ructure_file=vtmpnam();\n\nif ( $debug){print STDE\
RR \"\\n\\tCheck_point #1: $tmp_pdb_file  $structu\
re_file\\n\";}\n\n$INFILE=vfopen (\"$tmp_pdb_file\\
", \"r\"); \nmy $TMP=vfopen (\"$structure_file\", \
\"w\");\n\n$print_model=1;\n$in_model=0;\n\nif ( $\
debug){print STDERR \"\\n\\tCheck_point #2\\n\";}\\
nwhile ( <$INFILE>)\n{\n  my $first_model=0;\n  $l\
ine=$_;\n\n  if ( !$first_model && ($line =~/^MODE\
L\\s*(\\d*)/))\n    {\n      $first_model=$1;\n   \
   if ($model==1){$model=$first_model;}\n    }\n  \
\n  if (($line =~/^MODEL\\s*(\\d*)/))\n    {\n    \
  if ($1==$model)\n	{\n	  $in_model=1;\n	  $print_\
model=1;\n	  $is_nmr=1;\n	}\n      elsif ( $in_mod\
el==0)\n	{\n	  $print_model=0;\n	}\n      elsif ( \
$in_model==1)\n	{\n	  last;\n	}\n    }\n  if ($pri\
nt_model){print $TMP $line;}  \n}\nclose ($TMP);\n\
close ($INFILE);\n\nif ( $debug){print STDERR \"\\\
n\\tCheck_point #3\\n\";}	\n\n  if ($numbering eq \
\"\"){$numbering=\"absolute\";}\n  if ($numbering_\
out eq \"\"){$numbering_out=\"new\";}\n\n  if ( $d\
elete_set && $coor_set) {die \"-delete and -coor a\
re mutually exclusive, sorry\\n\";}\n  if ( $n_ato\
m==0){$atom_list[$n_atom++]=\"ALL\";$atom_list{$at\
om_list[0]}=1;}\n  if ( $seq_field eq \"\"){$seq_f\
ield=\"ATOM\";}\n  \n  if ( $MODE eq \"\"){$MODE=\\
"pdb\";}\n  elsif ( $MODE eq \"simple\" && $code==\
0){$code=1;}\n\n  if ( $code==0){$code=3;}\n\n\nif\
 ($f_chain[0] eq \" \"){$hc_chain{' '}=1;$c_chain[\
0]=\" \";}\nelsif (!@c_chain){$hc_chain{FIRST}=1;$\
c_chain[0]=\"FIRST\";}#make sure the first chain i\
s taken by default\n\nif    ($hc_chain{ALL}) \n{\n\
      @c_chain=@f_chain;\n      foreach $e (@c_cha\
in){$hc_chain{$e}=1;}\n}\nelsif($hc_chain{FIRST})\\
n{\n	@c_chain=($f_chain[0]);\n	$hc_chain{$f_chain[\
0]}=1;\n}\n\n\n$MAIN_HOM_CODE=&get_main_hom_code (\
$structure_file);\n$INFILE=vfopen ($structure_file\
, \"r\");\n\n\nif ( $MODE eq \"raw_pdb\" || $MODE \
eq \"raw\")\n{\n    while (<$INFILE>)\n{	print \"$\
_\";}\n    close ( $INFILE);\n    &myexit($EXIT_SU\
CCESS);\n}    \nif ( $MODE eq \"raw4fugue\" )\n{\n\
    while (<$INFILE>)\n{	\n	$l=$_;\n	if ($l=~/^SEQ\
RES/)\n{\n	    \n	    $c= substr($l,11,1);\n	    i\
f ($hc_chain {$c}){print \"$l\";}\n}\n	elsif ( $l=\
~/^ATOM/)\n{\n	    $c=substr($l,21,1);\n	    if ($\
hc_chain {$c}){print \"$l\";}\n}\n}\n    close ( $\
INFILE);\n    &myexit($EXIT_SUCCESS);\n}    \n\nif\
 ( $MODE eq \"pdb\")\n{\n\n    $read_header=0;\n  \
  while (<$INFILE>) \n{\n	    $line=$_;\n	    if  \
  ($line =~ /^HEADER/){print \"$line\";$read_heade\
r=1;}\n}\n    close ($INFILE);\n\n    if (!$read_h\
eader)\n{\n	print \"HEADER    UNKNOWN             \
                    00-JAN-00   $force_name\\n\";\\
n}\n\n    $INFILE=vfopen ($structure_file, \"r\");\
\n    \n    print \"COMPND   1 CHAIN:\";\n    $las\
t=pop(@c_chain);\n    foreach $c ( @c_chain){ prin\
t \" $c,\";}\n    if ( $last eq \" \"){print \" NU\
LL;\\n\";}\n    else \n{\n      print \" $last;\\n\
\";\n}\n    @c_chain=(@c_chain, $last);\n    \n   \
 print \"REMARK Output of the program extract_from\
_pdb (Version $VersionTag)\\n\";\n    print \"REMA\
RK Legal PDB format not Guaranteed\\n\";\n    prin\
t \"REMARK This format is not meant to be used in \
place of the PDB format\\n\";\n    print \"REMARK \
The header refers to the original entry\\n\";\n   \
 print \"REMARK The sequence from the original fil\
e has been taken in the field: $seq_field\\n\";\n \
   print \"REMARK extract_from_pdb, 2001, 2002, 20\
03, 2004, 2005 2006 (c) CNRS and Cedric Notredame\\
\n\";   \n    if ( $coor_set)\n{\n       print \"R\
EMARK Partial chain: Start $start End $end\\n\";\n\
}\n    if ( $is_nmr)\n{\n       print \"REMARK NMR\
 structure: MODEL $model\\n\";\n}\n   \n    if ( $\
n_atom!=0)\n{\n       print \"REMARK Contains Coor\
dinates of: \";\n       foreach $a (@atom){print \\
"$a \";}\n       print \"\\n\";\n}  \n}\n\n\n\n\nm\
y $residue_index = -999;\nmy $old_c = \"TemporaryC\
hain\";\n\nwhile (<$INFILE>) \n{\n	$line=$_;\n\n\n\
	if ($line =~ /^SEQRES/)\n{\n\n		@field=/(\\S*)\\s\
*/g;\n\n		$c= substr($_,11,1);\n\n		\n		$l=$#field\
;\n		for ($a=4; $a<$#field ;)\n{\n			if (!$onelett\
{$molecule_type{$c}}->{$field[$a]})\n{\n				splice\
 @field, $a, 1;\n}\n			else \n{\n				$a++;\n}\n}\n\
	\n		if ( $c ne $in_chain)\n{\n			$pdb_chain_list[\
$n_pdb_chains]=$c;\n			$pdb_chain_len [$n_pdb_chai\
ns]=$len;\n			$in_chain=$c;\n			$n_pdb_chains++;\n\
}\n	\n		for ( $a=4; $a<$#field;$a++)\n{\n			$compl\
ete_seq{$c}[$complete_seq_len{$c}++]=$field[$a];\n\
}\n}\n    elsif ( $line=~/^ATOM/ || ($line=~/^HETA\
TM/ && &is_aa(substr($line,17,3),substr($line,21,1\
)) && !$no_hetatm))\n{\n\n	 \n    $RAW_AT_ID=$AT_I\
D=substr($line,12,4);\n	$RES_ID=&is_aa(substr($lin\
e,17,3),substr($line,21,1));\n	$CHAIN=substr($line\
,21,1);\n\n    $RES_NO=substr($line,22,4);\n	$HOM_\
CODE=substr ($line, 26, 1);\n	$TEMP=substr($line,6\
0,6);\n	\n	$TEMP=~s/\\s//g;\n        $AT_ID=~s/\\s\
//g;\n	$RES_ID=~s/\\s//g;\n        $RES_NO=~s/\\s/\
/g;\n		\n	if ( $HOM_CODE ne $MAIN_HOM_CODE){next;}\
\n	elsif ( $already_read2{$CHAIN}{$RES_ID}{$AT_ID}\
{$RES_NO}){next;}\n	else{$already_read2{$CHAIN}{$R\
ES_ID}{$AT_ID}{$RES_NO}=1;}\n	\n	\n	if ($coor_set \
&& $numbering eq \"file\" && $residue_index ne $RE\
S_NO)\n{\n	    \n	    if ( $RES_NO<=$start){$real_\
start{$CHAIN}++;}\n	    if ( $RES_NO<=$end){$real_\
end{$CHAIN}++;}\n}\n	elsif ($numbering eq \"absolu\
te\")\n{\n	    $real_start{$CHAIN}=$start;\n	    $\
real_end{$CHAIN}=$end;\n}\n\n        $KEY=\"ALL\";\
\n        if ( $CHAIN ne $in_atom_chain)\n{\n	    \
\n	  $pdb_atom_chain_list[$n_pdb_atom_chains]=$c;\\
n	  $pdb_atom_chain_len [$n_pdb_atom_chains]=$len;\
\n	  $in_atom_chain=$c;\n	  $n_pdb_atom_chains++;\\
n}\n	\n	if ( $residue_index ne $RES_NO)\n{\n	     \
$residue_index = $RES_NO;\n	     $atom_seq{$CHAIN}\
[$atom_seq_len{$CHAIN}++]=$RES_ID;;\n}\n}\n\n}\ncl\
ose ($INFILE);\n\n\n\n\n\n\n$INFILE=vfopen ($struc\
ture_file, \"r\");\nforeach $c (@c_chain)\n{\n\n	i\
f    ( $seq_field eq \"SEQRES\"){@pdb_seq=@{$compl\
ete_seq{$c}};}\n	elsif ( $seq_field eq \"ATOM\")  \
{@pdb_seq=@{$atom_seq{$c}};}\n	\n\n	$full_length=$\
l=$#pdb_seq+1;\n		\n	if ( $real_end{$c}==\"*\"){$r\
eal_end{$c}=$full_length;}\n	if ( $coor_set)\n{	  \
 \n\n	   if ( $real_end{$c} < $l){splice @pdb_seq,\
 $real_end{$c}, $l;}\n	   if ( $real_start{$c} < $\
l){splice @pdb_seq, 0, $real_start{$c}-1;}	  	   \\
n	   $l=$#pdb_seq;\n}\n\n	elsif ( $delete_set)\n{\\
n	   splice @pdb_seq, $delete_start, $delete_end-$\
delete_start+1;\n	   $l=$#pdb_seq;\n}\n	\n	$new_fa\
sta_name=\"$pdb_id$c\";\n	if ( $coor_set)\n{\n	   \
if ( $n_pdb_chains==0){$new_fasta_name=\"$new_fast\
a_name$c\";}\n	   $new_fasta_name= $new_fasta_name\
.\"\\_$start\\_$end\";\n}\n	   \n	if ( $MODE eq \"\
pdb\")\n{\n	   $nl=1;\n	   $n=0;\n	   \n	   foreac\
h $res ( @pdb_seq)\n		{\n		if ( !$n)\n		{\n		\n		 \
printf \"SEQRES %3d %1s %4d  \", $nl,$c, $l;\n		 $\
nl++;\n	}\n	     $res=~s/\\s//g;\n	     \n	     if\
 ($code==1){ printf \"%3s \",$onelett{$molecule_ty\
pe{$c}}->{$res};}\n	     elsif  ($code==3){ printf\
 \"%3s \",$res};\n	     \n	     $n++;		  \n	     i\
f ( $n==13){$n=0;print \"\\n\";}\n}\n	  if ( $n!=0\
){print \"\\n\"; $n=0;}\n}\n	elsif ( $MODE eq \"si\
mple\")\n{\n	  print \"# SIMPLE_PDB_FORMAT\\n\";\n\
	  if ( $new_fasta_name eq \" \"){$new_fasta_name=\
\"dummy_name\";}\n	  print \">$new_fasta_name\\n\"\
;\n\n	  foreach $res ( @pdb_seq)\n{\n	      print \
\"$onelett{$molecule_type{$c}}->{$res}\";\n}\n	  p\
rint \"\\n\";\n}\n	elsif ( $MODE eq \"fasta\")\n{\\
n	  $n=0;\n	  print \">$new_fasta_name\\n\";\n	  \\
n	  foreach $res ( @pdb_seq)\n{\n\n	      print \"\
$onelett{$molecule_type{$c}}->{$res}\";\n         \
     $n++;\n	      if ( $n==60){print \"\\n\"; $n=\
0;}\n}\n	  print \"\\n\"; \n}\n}\n\nif ( $MODE eq \
\"fasta\")\n{\n     &myexit($EXIT_SUCCESS);\n  \n}\
\n\n  \n  $charcount=0;\n  $inchain=\"BEGIN\";\n  \
$n=0;\n  while (<$INFILE>) \n{\n    $line=$_;\n   \
  \n    if ($line =~/^ATOM/  ||  ($line=~/^HETATM/\
))\n{\n	$line_header=\"UNKNWN\";\n	$RES_ID=substr(\
$line,17,3);\n	$chain = substr($line,21,1);\n\n	if\
 ($line =~/^ATOM/)\n{\n	    $line_header=\"ATOM\";\
\n	    $RES_ID=(&is_aa($RES_ID,$chain))?&is_aa($RE\
S_ID,$chain):$RES_ID;\n}\n	elsif ($line=~/^HETATM/\
 && ($ligand_list {$RES_ID} ||$ligand_list {'ALL'}\
 || !&is_aa($RES_ID,$chain)))\n{\n	    $line_heade\
r=\"HETATM\";\n}\n	elsif ($line=~/^HETATM/ && (&is\
_aa($RES_ID,$chain) && !$no_hetatm))\n{\n	    $lin\
e_header=\"ATOM\";\n	    $RES_ID=&is_aa($RES_ID,$c\
hain);\n}\n	else\n{\n	    next;\n}\n\n	\n\n	$X=sub\
str($line,30,8);     \n	$Y=substr($line,38,8);\n	$\
Z=substr($line,46,8);\n	$TEMP=substr($line,60,6);\\
n	\n	$RAW_AT_ID=$AT_ID=substr($line,12,4);\n	$CHAI\
N=substr($line,21,1);\n	$RES_NO=substr($line,22,4)\
;\n	$HOM_CODE=substr ($line, 26, 1);\n	\n	$X=~s/\\\
s//g;\n	$Y=~s/\\s//g;\n	$Z=~s/\\s//g;\n	$TEMP=~s/\\
\s//g;\n	\n	$AT_ID=~s/\\s//g;\n	$RES_ID=~s/\\s//g;\
\n	$RES_NO=~s/\\s//g;\n\n	\n	if ( $HOM_CODE ne $MA\
IN_HOM_CODE){next;}\n	elsif ( $already_read{$CHAIN\
}{$RES_ID}{$AT_ID}{$RES_NO}){next;}\n	else{$alread\
y_read{$CHAIN}{$RES_ID}{$AT_ID}{$RES_NO}=1;}\n	\n	\
$KEY=\"ALL\";\n\n      	if ( $RES_NO ==0){$start_a\
t_zero=1;}\n\n	$RES_NO+=$start_at_zero;    \n	\n	i\
f ( $current_chain ne $CHAIN)\n{\n	    $current_ch\
ain=$CHAIN;\n	    $pos=$current_residue=0;\n	    $\
offset=($coor_set)?($real_start{$CHAIN}-1):0;\n	  \
  if    ( $seq_field eq \"SEQRES\"){@ref_seq=@{$co\
mplete_seq{$CHAIN}};}\n	    elsif ( $seq_field eq \
\"ATOM\")  {@ref_seq=@{$atom_seq{$CHAIN}};}\n}\n	\\
n	if ($current_residue != $RES_NO)\n{\n	    $curre\
nt_residue=$RES_NO;\n	    if    ( $seq_field eq \"\
SEQRES\"){$pos=$current_residue;}\n	    elsif ( $s\
eq_field eq \"ATOM\"){$pos++;}\n}\n	\n	\n	if ($n_a\
tom==0 || $atom_list{$AT_ID}==1 || $atom_list{$KEY\
}==1)\n{ 	\n	    \n	    $do_it=(!@c_chain || $hc_c\
hain{$CHAIN} ||$hc_chain{'LIGAND'} );\n	    \n	   \
 $do_it= ($do_it==1) && ($coor_set==0 ||($pos>=$re\
al_start{$CHAIN} && $pos<=$real_end{$CHAIN}));\n	 \
   $do_it= ($do_it==1) && ($delete_set==0 || $pos<\
$delete_start ||$pos>$delete_end );\n	    if ($lig\
and==0 && $line_header eq \"HETATM\" ){$do_it=0;}\\
n	    if ($ligand_only==1 && $line_header eq \"ATO\
M\" ){$do_it=0;}\n	    if ($ligand==1 && $line_hea\
der eq \"HETATM\" && $ligand_list{$RES_ID}==0 && $\
ligand_list{\"ALL\"}==0){$do_it=0;} \n	    \n	    \
\n	    if ( $do_it)\n{\n		$n++;\n		$out_pos=$pos;\\
n		\n	       if ( $delete_set)\n{\n		  if ( $out_p\
os< $delete_start){;}\n		  else {$offset=$delete_e\
nd-$delete_start;}\n}       \n	       \n	       if\
 ( $numbering_out eq \"new\"){$out_pos-=$offset;}\\
n	       elsif ( $numbering_out eq \"old\"){$out_p\
os=$RES_NO;}\n	       \n       \n	       \n	      \
 if ( $code==1){$RES_ID=$onelett{$molecule_type{$c\
}}->{$RES_ID};}\n	    \n	       if ($unfold)\n{\n	\
	   $unfolded_x+=5;\n		   $X=$unfolded_x;\n		   $Y\
=0;\n		   $Z=0;\n		   $float=1;\n}\n	       else\n\
{\n		   $float=3;\n}\n\n	       if ( $MODE eq \"pd\
b\")\n{\n		   printf \"%-6s%5d %-4s %3s %s%4d    %\
8.3f%8.3f%8.3f  1.00 %5.2f\\n\",$line_header, $n, \
$RAW_AT_ID,$RES_ID,$CHAIN,$out_pos, $X, $Y, $Z,$TE\
MP;		  \n}\n	       elsif ( $MODE eq \"simple\")\n\
{\n		    if ( $RES_ID eq \"\"){$RES_ID=\"X\";}\n		\
  printf \"%-6s %5s %s %2s %4d    %8.3f %8.3f %8.3\
f\\n\",$line_header, $AT_ID, $RES_ID,($CHAIN eq\"\\
" || $CHAIN eq \" \")?\"A\":$CHAIN,$out_pos, $X, $\
Y, $Z,$TEMP;\n}\n\n}\n}\n}\n}\nprint \"\\n\";\nclo\
se($INFILE);\n\n\nif ( $error ne \"\") \n{$error=$\
error.\"\\nDiagnostic:    SEQRES and the residues \
in ATOM are probably Incompatible\\n\";\n    $erro\
r=$error.  \"Recomendation: Rerun with '-fix 1' in\
 order to ignore the SEQRES sequences\\n\";\n}\nif\
 (!$nodiagnostic){print STDERR $error;}\n&myexit (\
 $EXIT_SUCCESS);\n\nsub get_pdb_entry_type_file\n \
 {\n    my $cache_file=\"$cache/pdb_entry_type.txt\
\";\n    my $env_file  = $ENV{\"PDB_ENTRY_TYPE_FIL\
E\"};\n    my $pdb_file  =\"$ENV{'PDB_DIR'}/derive\
d_data/pdb_entry_type.txt\";\n    \n    \n    if (\
-z $cache_file){unlink ($cache_file);}#will get up\
dated\n    if (-z $env_file){$env_file=\"\";}    #\
cannot update\n    if (-z $pdb_file){$pdb_file=\"\\
";}    #cannot update\n    \n    if    (-e $env_fi\
le){return $env_file;} #env wins: user decides\n  \
  elsif (-e $pdb_file){return $pdb_file;} #local d\
atabase wins: network file may be out of sync\n   \
 elsif ($no_remote_pdb_dir==1)\n      {\n	if (-e $\
cache_file){return $cache_file;}\n	else\n	  {add_w\
arning($$,$$,\"PDB_ENTRY_TYPE_FILE must be set to \
the location of <pdb>/derived_data/pdb_entry_type.\
txt when using NO_REMOTE_PDB_DIR=1\");\n	   return\
 \"\";\n	 }\n      }\n    else #update can only ta\
ke place if the file lives in cache\n      {\n	my \
$new_file;\n	if (!-e $cache_file || (-M $cache_fil\
e)>1)\n	  {\n	    $new_file=vtmpnam();\n	    &url2\
file(\"ftp://ftp.wwpdb.org/pub/pdb/derived_data/pd\
b_entry_type.txt\", $new_file);\n	    if ( !-z $ne\
w_file){system (\"mv $new_file $cache_file\"); unl\
ink ($new_file); $new_file=$cache_file;}\n	    els\
e {unlink($new_file);}\n	  }\n	else\n	  {\n	    $n\
ew_file=$cache_file;\n	  }\n	\n	if (!-e $cache_fil\
e && !-e $new_file)\n	  {\n	    add_warning($$,$$,\
\"Could not download ftp://ftp.wwpdb.org/pub/pdb/d\
erived_data/pdb_entry_type.txt\");\n	    return \"\
\";\n	  }\n	elsif (-e $cache_file && !-e $new_file\
)\n	  {\n	    my $m=(-M $cache_file);\n	    add_wa\
rning($$,$$,\"Could not update file ftp://ftp.wwpd\
b.org/pub/pdb/derived_data/pdb_entry_type.txt. Old\
er Version [$cache_file]($m Month(s) old) will be \
used instead\");\n	    return $cache_file;\n	  }\n\
	else\n	  {\n	    return $new_file;\n	  }\n      }\
\n  }\n\n\n\nsub get_unrealeased_file\n  {\n    my\
 $cache_file=\"$cache/unrealeased.xml\";\n    my $\
env_file  = $ENV{\"PDB_UNREALEASED_FILE\"};\n    m\
y $pdb_file  =\"$ENV{'PDB_DIR'}/derived_data/unrea\
leased.xml\";\n    \n    \n    if ($env_file eq \"\
NO\" || $env_file eq \"No\" || $env_file eq \"no\"\
 || $env_file eq \"0\"){return \"NO\";}\n\n    if \
(-z $cache_file){unlink ($cache_file);}#will get u\
pdated\n    if (-z $env_file){unlink($env_file);} \
    #will update\n    if (-z $pdb_file){$pdb_file=\
\"\";}          #cannot update\n    \n    if    (-\
e $env_file){return $env_file;} #env wins: user de\
cides\n    elsif (-e $pdb_file){return $pdb_file;}\
 #local database wins: network file may be out of \
sync\n    elsif ($no_remote_pdb_dir==1)        \n \
     {\n	if (-e $cache_file){return $cache_file;}\\
n	elsif ( $env_file && ! -e $env_file)\n	  {\n	   \
 &url2file(\"http://www.rcsb.org/pdb/rest/getUnrel\
eased\",$env_file);\n	    if ( -e $env_file && !-z\
 $env_file){return $env_file;}\n	  }\n	else\n	  {\\
n	    add_warning($$,$$,\"UNREALEASED_FILE must be\
 set to the location of your unrealeased.xml file \
as downloaded from http://www.rcsb.org/pdb/rest/ge\
tUnreleased when using NO_REMOTE_PDB_DIR=1\");\n	 \
   return \"\";\n	  }\n      }\n    else #update c\
an only take place if the file lives in cache\n   \
   {\n	my $new_file=vtmpnam ();\n	if (!-e $cache_f\
ile || (-M $cache_file)>1)\n	  {\n	    &url2file(\\
"http://www.rcsb.org/pdb/rest/getUnreleased\",$new\
_file);\n	    if ( !-z $new_file){system (\"mv $ne\
w_file $cache_file\"); unlink ($new_file); $new_fi\
le=$cache_file;}\n	    else {unlink($new_file);}\n\
	  }\n	else\n	  {\n	    $new_file=$cache_file;\n	 \
 }\n	\n	if (!-e $cache_file && !-e $new_file)\n	  \
{\n	    add_warning($$,$$,\"Could not download htt\
p://www.rcsb.org/pdb/rest/getUnreleased\");\n	    \
return \"\";\n	  }\n	elsif (-e $cache_file && !-e \
$new_file)\n	  {\n	    my $m=(-M $cache_file);\n	 \
   add_warning($$,$$,\"Could not update file http:\
//www.rcsb.org/pdb/rest/getUnreleased. Older Versi\
on [$cache_file]($m Month(s) ) will be used\");\n	\
    return $cache_file;\n	  }\n	else\n	  {\n	    r\
eturn $new_file;\n	  }\n      }\n  }\n\nsub is_rel\
eased \n  {\n    my ($r);\n    my $in=@_[0];\n    \
my $name=&remote_is_pdb_name ($in);\n    my $hold=\
&remote_is_on_hold($in);\n    \n    $r=($name && !\
$hold)?1:0;\n    return $r;\n  }\n\nsub remote_is_\
pdb_name \n  {\n    my $in=@_[0];\n    my ($pdb);\\
n    my ($value,$value1,$value2);\n    my $max=2;\\
n    \n    \n    \n    my $ref_file=&get_pdb_entry\
_type_file();\n    \n    if ( $in=~/[^\\w\\d\\:\\_\
]/){return 0;}\n    elsif (!-e $ref_file)\n      {\
\n	add_warning ($$,$$,\"Cannot find pdb_entry_type\
.txt;  $in is assumed to be valid; add ftp://ftp.w\
wpdb.org/pub/pdb/derived_data/pdb_entry_type.txt i\
n $cache to automatically check name status\");\n	\
return 1;\n      }\n    else\n      {\n	$pdb=subst\
r ($in,0, 4);\n	chomp(($value1=`grep -c $pdb $ref_\
file`));\n	$pdb=lc($pdb);\n	chomp(($value2=`grep -\
c $pdb $ref_file`));\n	$value=($value1 || $value2)\
?1:0;\n	$value=($value>0)?1:0;\n	\n	return $value;\
\n      }\n  }\n\n\n\nsub pdb2model_type\n{\n    m\
y $in=@_[0];\n    my ($ref_file, $pdb);\n    my ($\
value, $ret);\n\n    if ( $in=~/[^\\w\\d\\:\\_]/){\
return 0;}\n    $ref_file=&get_pdb_entry_type_file\
();\n    if (!-e $ref_file)\n      {\n	add_warning\
 ($$,$$,\"Cannot find pdb_entry_type.txt;  $in is \
assumed to be diffraction; add ftp://ftp.wwpdb.org\
/pub/pdb/derived_data/pdb_entry_type.txt in $cache\
 to check name status\");\n	return \"diffraction\"\
;\n      }\n    else\n      {\n	$pdb=substr ($in,0\
, 4);\n	$pdb=lc($pdb);\n	\n	chomp(($value=`grep $p\
db $ref_file`));\n	\n	$value=~/^\\S+\\s+\\S+\\s+(\\
\S+)/;\n	$ret=$1;\n	if ( $ret eq\"\"){return \"UNK\
NOWN\";}\n	\n	return $ret;\n      }\n  }\nsub remo\
te_is_on_hold\n  {\n    my $in=@_[0];\n    my ($re\
f_file, $pdb);\n    my ($value1, $value2,$value);\\
n    \n\n\n    \n    $ref_file=&get_unrealeased_fi\
le();\n    if ($ref_file eq \"NO\"){return 0;}\n\n\
\n    if ($no_remote_pdb==1){return 0;}\n    if ( \
$in=~/[^\\w\\d\\:\\_]/){return 0;}\n    \n    $ref\
_file=&get_unrealeased_file();\n    if (!-e $ref_f\
ile)\n      {\n	add_warning ($$,$$,\"Cannot find u\
nrealeased.xml;  $in is assumed to be released;\")\
;\n	return 1;\n      }\n    \n    $pdb=substr ($in\
,0, 4);\n    chomp(($value1=`grep -c $pdb $ref_fil\
e`));\n    $pdb=lc($pdb);\n    chomp(($value2=`gre\
p -c $pdb $ref_file`));\n    $value=($value1 || $v\
alue2)?1:0;\n    $value=($value>0)?1:0;\n    retur\
n $value;\n  }\n\nsub is_pdb_file\n  {\n    my @ar\
g=@_;\n    \n    if ( !-e $arg[0]){return 0;}\n   \
 \n    $F=vfopen ($arg[0], \"r\");\n    while ( <$\
F>)\n      {\n	if (/^HEADER/)\n	  {\n	    close $F\
;\n	    return 1;\n	  }\n	elsif ( /^SEQRES/)\n	  {\
\n	    close $F;\n	    return 1;\n	  }\n	elsif ( /\
^ATOM/)\n	  {\n	    close $F;\n	    return 1;\n	  \
}\n      }\n    return 0;\n  }\nsub get_pdb_id\n{\\
n    my $header_file=@_[0];\n    my $id;\n    my $\
F= new FileHandle;\n    \n    \n    $F=vfopen (\"$\
header_file\", \"r\");\n\n    while ( <$F>)\n     \
 {\n	if ( /HEADER/)\n	  {\n	    if ($debug){print \
\"$_\";}\n	    $id=substr($_,62,4 );\n	    return \
$id;\n	  }\n      }\n    close ($F);\n    \n    re\
turn \"\";\n}\n\nsub get_ligand_list\n{\n    my $p\
db_file=@_[0];\n    my $chain;\n    my $ligand;\n \
   my %complete_ligand_list;\n    \n\n    $F=vfope\
n ($pdb_file, \"r\");\n    while ( <$F>)\n{\n	if (\
 /^HETATM/)\n{\n	    $line=$_;\n	    $chain=substr\
($line,21,1);\n	    $ligand=substr($line,17,3);\n	\
    \n	    if (!$complete_ligand_list{$chain}{$lig\
and})\n{\n		\n		$complete_ligand_list{\"result\"}.\
=\"CHAIN $chain LIGAND $ligand\\n\";\n		$complete_\
ligand_list{$chain}{$ligand}=1;\n}\n}\n}\n    clos\
e ($F);\n    return %complete_ligand_list;\n}\n\ns\
ub get_chain_list \n{\n    my $header_file;\n    m\
y @chain_list;\n    my @list;\n    my $n_chains;\n\
    my %chain_hasch;\n    my $pdb_file=@_[0];\n   \
 my $c;\n    my %hasch;\n    my $chain;\n  \n    \\
n    $F=vfopen ($pdb_file, \"r\");\n    while ( <$\
F>)\n{\n\n\n	if (/SEQRES\\s+\\d+\\s+(\\S+)/)\n	  {\
\n	    $chain = substr($_,11,1);$chain=~s/\\s//g;i\
f ( $chain eq \"\"){$chain=\" \";}\n	    if (!$has\
ch{$chain}){$hasch{$chain}=1;push @chain_list, $ch\
ain;}\n	  }\n	if (/^ATOM/ || /^HETATM/)\n	  {\n	  \
  $chain = substr($_,21,1); $chain=~s/\\s//g;if ( \
$chain eq \"\"){$chain=\" \";}\n	    if (!$hasch{$\
chain}){$hasch{$chain}=1;push @chain_list, $chain;\
}\n	  }\n      }\n\n\nclose ($F);\nif (!@chain_lis\
t)\n  {\n    @chain_list=(\"A\");\n  }\n\n\nreturn\
 @chain_list;\n}\n\nsub token_is_in_list\n{\n\n   \
 my @list=@_;\n    my $a;\n    \n    for ($a=1; $a\
<=$#list; $a++)\n{\n	if ( $list[$a] eq $list[0]){r\
eturn $a;}\n}\n}\n\nsub pdb_name2name_and_chain \n\
{\n    my $pdb_file=@_[0];\n    my $pdb_file_in;\n\
    my @array;\n    my $chain;\n    my $c;\n\n    \
$pdb_file_in=$pdb_file;\n\n    $pdb_file=~/^(.{4})\
/;$pdb_id=$1;\n    @array=($pdb_file=~/([\\w])/g);\
\n  \n  \n    $chain=uc ($array[4]);\n    $chain=(\
$chain eq \"\")?\"FIRST\":$chain;\n    \n    retur\
n ( $pdb_id, $chain);\n\n    if ( $#array==3){retu\
rn ($pdb_id, \"FIRST\");}\n    elsif ( $#array<4){\
 return ($pdb_id, \"\");}\n    else {return ( $pdb\
_id, $chain);}\n      \n    \n    \n}\nsub get_mai\
n_hom_code \n{\n    my $pdb_file=@_[0];\n    my %h\
om, $n, $best, $best_h;\n    open (F, $pdb_file);\\
n    while (<F>)\n{\n	if ( $_=~/^ATOM/)\n{\n	    $\
h=substr ($_,26, 1);\n	    $n=++$hom{$h};\n	    if\
 ($n>$best)\n{\n		$best=$n;\n		$best_h=$h;\n}\n}\n\
}\n    close (F);\n    return $best_h;\n}\n\n\nsub\
 get_pdb_file \n{\n    my ($pdb_file_in)=(@_);\n  \
  my $result;\n    my @letter;\n    my @chain;\n  \
  my $v;\n    my $pdb_file=$pdb_file_in;\n\n    $p\
db_file=($pdb_file_in=~/\\S+_S_(\\S+)/)?$1:$pdb_fi\
le_in;\n    \n    if ($no_remote_pdb_dir==0)\n    \
  {\n	$no_remote_pdb_dir=1;\n	$result=get_pdb_file\
3 ($pdb_file);\n	$no_remote_pdb_dir=0;\n	if ( $res\
ult){return $result;}\n	else\n	  {\n	    \n	    lc\
 ($pdb_file);\n	    $result=get_pdb_file3($pdb_fil\
e);\n	    return  $result;\n	  }\n      }\n    els\
e\n      {\n	return get_pdb_file3 ($pdb_file);\n  \
    }\n    \n  }\n\nsub get_pdb_file3 \n{\n    my \
$pdb_file_in=@_[0];\n    my $result;\n    my @lett\
er;\n    my @chain;\n    my $lcfile;\n    my $ucfi\
le;\n    my $pdb_file=$pdb_file_in;\n    \n    $lc\
file=lc $pdb_file;\n    $ucfile=uc $pdb_file;\n\n \
   if ( ($result=get_pdb_file2 ($pdb_file))){retur\
n $result;}\n    \n\n    if ($lcfile ne $pdb_file \
&& ($result=get_pdb_file2 ($lcfile))){return $resu\
lt;}\n    if ($ucfile ne $pdb_file && ($result=get\
_pdb_file2 ($ucfile))){return $result;}\n    \n   \
\n    \n    return \"\";\n}\nsub get_pdb_file2\n{\\
n    my $pdb_file=@_[0];\n    my $return_value;\n \
   \n    $return_value=\"\";\n    \n    if ( ($res\
ult=get_pdb_file1 ($pdb_file))){$return_value=$res\
ult;}\n    elsif ( !($pdb_file=~/\\.pdb/) && !($pd\
b_file=~/\\.PDB/))\n{\n	if ( ($result=get_pdb_file\
1 (\"$pdb_file.pdb\"))){$return_value=$result;}\n	\
elsif ( ($result=get_pdb_file1 (\"$pdb_file.PDB\")\
)){$return_value=$result;}\n\n	elsif ( ($result=ge\
t_pdb_file1 (\"pdb$pdb_file.pdb\"))){$return_value\
=$result;}	\n	elsif ( ($result=get_pdb_file1 (\"pd\
b$pdb_file.PDB\"))){$return_value=$result;}\n	elsi\
f ( ($result=get_pdb_file1 (\"PDB$pdb_file.PDB\"))\
){$return_value=$result;}\n	elsif ( ($result=get_p\
db_file1 (\"PDB$pdb_file.pdb\"))){$return_value=$r\
esult;}\n	\n	\n	elsif ( ($result=get_pdb_file1 (\"\
$pdb_file.ent\"))){$return_value=$result;}\n	elsif\
 ( ($result=get_pdb_file1 (\"pdb$pdb_file.ent\")))\
{$return_value=$result;}\n	elsif ( ($result=get_pd\
b_file1 (\"PDB$pdb_file.ent\"))){$return_value=$re\
sult;}\n\n	elsif ( ($result=get_pdb_file1 (\"$pdb_\
file.ENT\"))){$return_value=$result;}\n	elsif ( ($\
result=get_pdb_file1 (\"pdb$pdb_file.ENT\"))){$ret\
urn_value=$result;}\n	elsif ( ($result=get_pdb_fil\
e1 (\"PDB$pdb_file.ENT\"))){$return_value=$result;\
}\n	\n	\n	\n}\n    return $return_value;\n}\n    \\
nsub get_pdb_file1\n{\n    my ($pdb_file)=(@_);\n \
   my $return_value;\n    \n\n    $return_value=\"\
\";\n    if ( ($result=get_pdb_file0 ($pdb_file)))\
{$return_value=$result;}\n    elsif ( ($result=get\
_pdb_file0 (\"$pdb_file.Z\"))){$return_value=$resu\
lt;}\n    elsif ( ($result=get_pdb_file0 (\"$pdb_f\
ile.gz\"))){$return_value=$result;}\n    elsif ( (\
$result=get_pdb_file0 (\"$pdb_file.GZ\"))){$return\
_value=$result;}\n    return $return_value;\n}\nsu\
b get_pdb_file0 \n{ \n    my ($pdb_file, $attempt)\
=(@_);\n    my $pdb_file=@_[0];\n    my $tmp_pdb_f\
ile;    \n    my $return_value;\n\n    if ( !$atte\
mpt){$attempt=1;}\n    \n    $local_pdb_file=\"$pd\
b_file\";\n    if ( $local_pdb_file eq \"\")\n{\n	\
$tmp_pdb_file=vtmpnam();\n	open F, \">$tmp_pdb_fil\
e\";\n	\n	while (<STDIN>){print F \"$_\";}\n	close\
 (F);\n	\n	if (-e $tmp_pdb_file && &is_pdb_file ( \
$local_pdb_file))\n{return $tmp_pdb_file;}\n}\n\n \
   $local_pdb_file=\"$pdb_file\";\n    &debug_prin\
t (\"\\nTry access local file: $local_pdb_file\");\
\n    \n    $local_pdb_file=&check_pdb_file4compre\
ssion ($local_pdb_file);\n    if ( -e $local_pdb_f\
ile && (&is_pdb_file ($local_pdb_file) || $force_p\
db))\n{\n	&debug_print ( \"\\n\\tIs in Current Dir\
\");\n	$tmp_pdb_file=vtmpnam();\n	`cp $local_pdb_f\
ile $tmp_pdb_file`;\n	return $tmp_pdb_file;\n}\n  \
  else\n{\n	&debug_print (\"\\n\\tFile Not in Curr\
ent Dir\");\n}\n\n    if ($pdb_file=~/^pdb/||$pdb_\
file=~/^PDB/){$pdb_div=substr ($pdb_file, 4, 2);}\\
n    else\n{\n	  $pdb_div=substr ($pdb_file, 1, 2)\
;\n}\n    $local_pdb_file=\"$pdb_dir/$pdb_div/$pdb\
_file\";\n    $local_pdb_file=&check_pdb_file4comp\
ression ( $local_pdb_file);\n    &debug_print (\"\\
\nTry access file From PDB_DIR: $local_pdb_file\")\
;\n    if ($pdb_dir && -e $local_pdb_file && &is_p\
db_file ($local_pdb_file))\n{\n	&debug_print ( \"\\
\n\\tIs in Local PDB DIR\");\n	$tmp_pdb_file=vtmpn\
am();\n	`cp $local_pdb_file $tmp_pdb_file`;\n	retu\
rn $tmp_pdb_file;\n}\n\n    $local_pdb_file=\"$pdb\
_dir/$pdb_file\";\n    $local_pdb_file=&check_pdb_\
file4compression ( $local_pdb_file);\n    &debug_p\
rint (\"\\nTry access file From PDB_DIR: local_pdb\
_file\");\n    if ($pdb_dir && -e $local_pdb_file \
&& &is_pdb_file ($local_pdb_file))\n{\n	&debug_pri\
nt ( \"\\n\\tIs in Local PDB DIR\");\n	$tmp_pdb_fi\
le=vtmpnam();\n	`cp $local_pdb_file $tmp_pdb_file`\
;\n	return $tmp_pdb_file;\n}\n\n    $local_pdb_fil\
e=\"$pdb_dir$pdb_file\";\n    $local_pdb_file=&che\
ck_pdb_file4compression ( $local_pdb_file);\n    &\
debug_print (\"\\nTry access file From PDB_DIR: $l\
ocal_pdb_file\");\n    if ($pdb_dir && -e $local_p\
db_file && &is_pdb_file ($local_pdb_file))\n{\n	&d\
ebug_print ( \"\\n\\tIs in Local PDB DIR\");\n	$tm\
p_pdb_file=vtmpnam();\n	`cp $local_pdb_file $tmp_p\
db_file`;\n	return $tmp_pdb_file;\n}\n    else\n{&\
debug_print ( \"\\n\\tNot In Local Pdb Dir\");}\n\\
n    if ($cache ne \"NO\" && $cache ne \"no\")\n{\\
n\n	$local_pdb_file=\"$cache/$pdb_file\";\n	$local\
_pdb_file=&check_pdb_file4compression ( $local_pdb\
_file);\n	&debug_print(\"\\nTry access file From C\
ache: $local_pdb_file\");\n	if (-e $local_pdb_file\
 && &is_pdb_file ($local_pdb_file))\n{\n	    &debu\
g_print ( \"\\n\\tIs in T-Coffee Cache\");\n	    $\
tmp_pdb_file=vtmpnam();\n	    `cp $local_pdb_file \
$tmp_pdb_file`;\n	    return $tmp_pdb_file;\n}\n	e\
lse{&debug_print ( \"\\n\\tNot in Cache Dir\");}\n\
}\n\nif (!$no_remote_pdb_dir) \n  {\n    my $value\
=&is_released ($pdb_file);\n    my $return_value=\\
"\";\n    if ($value==1)\n      {\n	\n	&debug_prin\
t (\"\\n******************************************\
***********\\nTry Remote Access for $pdb_file\");\\
n	$tmp_pdb_file=vtmpnam();\n	$netcommand=$netaddre\
ss;\n	$netcommand=~s/%%/$pdb_file/g;\n	&url2file(\\
"$netcommand\", \"$tmp_pdb_file.$netcompression\")\
;\n	&debug_print(\"\\nREMOTE: $netcommand\\n\");\n\
	\n	$compressed_tmp_file_name=\"$tmp_pdb_file.$net\
compression\";\n	\n	if ($netcompression && -B $com\
pressed_tmp_file_name && $attempt<5)\n	  {\n	    m\
y $r;\n	    &debug_print (\"\\n\\tFile Found Remot\
ely\");\n	    if (($r=safe_system ( \"$netcompress\
ion_pg $compressed_tmp_file_name\")!=$EXIT_SUCCESS\
) && $attempts<5)\n	      {\n		&debug_print (\"\\n\
\\tProper Download Failed Try again\");\n		unlink \
$compressed_tmp_file_name;\n		print \"\\nFailed to\
 Download $compressed_tmp_file_name. New Attempt $\
attempt/5\\n\";\n		return &get_pdb_file0($pdb_file\
, $attempt+1);\n	      }\n	    elsif ($r== $EXIT_S\
UCCESS)\n	      {\n		&debug_print (\"\\n\\tProper \
Download Succeeded \");\n		$return_value=$tmp_pdb_\
file;\n	      }\n	    else\n	      {\n		&debug_pri\
nt (\"\\n\\tProper Download Failed \");\n		&debug_\
print (\"\\nFile Not Found Remotely\");\n		unlink \
$compressed_tmp_file_name;\n	      }\n	  }\n	else\\
n	  {\n\n	    &debug_print (\"\\nFile Not Found Re\
motely\");\n	    unlink $compressed_tmp_file_name;\
\n	  }\n	#Update cache if required\n	if ($cache ne\
 \"no\" && $cache ne \"update\" && -e $return_valu\
e)\n	  {\n	    `cp $return_value $cache/$pdb_file.\
pdb`;\n	    #`t_coffee -other_pg clean_cache.pl -f\
ile $pdb_file.pdb -dir $cache`;\n	  }\n      }\n  \
  &debug_print (\"\\nRemote Download Finished\");\\
n    return $return_value;\n  }\nreturn \"\";\n}\n\
\nsub check_pdb_file4compression \n{\n    my $file\
=@_[0];\n    my $tmp;\n    my $r;\n    \n    $tmp=\
&vtmpnam();\n    if (-e $tmp){unlink $tmp;}\n    \\
n    $file=~s/\\/\\//\\//g;\n    if    (-B $file &\
& ($file=~/\\.Z/)) {`cp $file $tmp.Z`;`rm $tmp`;`g\
unzip $tmp.Z $SILENT`;$r=$tmp;}\n    elsif (-B $fi\
le && ($file=~/\\.gz/)){`cp $file $tmp.gz`;`gunzip\
 $tmp.gz $SILENT`;return $r=$tmp;}\n    elsif (-B \
$file ){`cp $file $tmp.gz`;`gunzip $tmp.gz $SILENT\
`;$r=$tmp;}\n    elsif ( -e $file ) {$r= $file;}\n\
    elsif ( -e \"$file.gz\" ){ `cp $file.gz $tmp.g\
z`;`gunzip     $tmp.gz $SILENT`;$r=$tmp;}    \n   \
 elsif ( -e \"$file.Z\") {`cp $file.Z  $tmp.Z`; `g\
unzip $tmp.Z $SILENT`;$r=$tmp;}\n    else  {$r= $f\
ile;}\n\n    if ( -e \"$tmp.Z\"){unlink \"$tmp.Z\"\
;}\n    if ( -e \"$tmp.gz\"){unlink \"$tmp.gz\";}\\
n    \n    return $r;\n    \n}\n\n\n\n\n\n    \n\n\
\n\n\n\n\n\nsub vfopen \n{\n    my $file=@_[0];\n \
   my $mode=@_[1];\n    my $tmp;\n    my $F = new \
FileHandle;\n    \n    \n    $tmp=$file;\n	\n    \\
n    if ( $mode eq \"r\" && !-e $file){ myexit(flu\
sh_error (\"Cannot open file $file\"));}\n    elsi\
f ($mode eq \"w\"){$tmp=\">$file\";}\n    elsif ($\
mode eq \"a\"){$tmp=\">>$file\";}\n    \n    \n   \
 open ($F,$tmp);\n    return $F;\n}\nsub debug_pri\
nt\n{\n    my $message =@_[0];\n    if ($debug){pr\
int STDERR \"NO_REMOTE_PDB_DIR: $no_remote_pdb_dir\
 - $message [DEBUG:extract_from_pdb]\";}\n    retu\
rn;\n}\nsub is_aa \n{\n    my ($aa, $chain) =@_;\n\
\n    my $one;\n    my $trhee;\n    \n    if ( $on\
elett{$molecule_type{$chain}}->{$aa} eq 'X' || !$o\
nelett{$molecule_type{$chain}}->{$aa} ){return '';\
}\n    else\n      {\n	$one=$onelett{$molecule_typ\
e{$chain}}->{$aa};\n\n	$three=$threelett{$molecule\
_type{$chain}}->{$one};\n	\n\n	return $three;\n   \
   }\n  }\n\n\n\n\n\nsub url2file\n{\n    my ($add\
ress, $out, $wget_arg, $curl_arg)=(@_);\n    my ($\
pg, $flag, $r, $arg, $count);\n    \n    if (!$CON\
FIGURATION){&check_configuration (\"wget\", \"INTE\
RNET\", \"gzip\");$CONFIGURATION=1;}\n    \n    if\
 (&pg_is_installed (\"wget\"))   {$pg=\"wget\"; $f\
lag=\"-O\";$arg=$wget_arg;}\n    elsif (&pg_is_ins\
talled (\"curl\")){$pg=\"curl\"; $flag=\"-o\";$arg\
=$curl_arg;}\n    return safe_system (\"$pg $flag$\
out $address >/dev/null 2>/dev/null\");\n\n}\n\n\n\
\n\nsub pdbfile2chaintype\n  {\n    my $file=@_[0]\
;\n    my %ct;\n    my $F;\n    \n    $F=vfopen ($\
file, \"r\");\n    while (<$F>)\n      {\n	my $lin\
e=$_;\n	if ($line =~/^ATOM/)\n	  {\n	    my $C=sub\
str($line,21,1);\n	    if (!$ct{$C})\n	      {	\n	\
	my $r=substr($line,17,3);\n		$r=~s/\\s+//;\n		if \
(length ($r)==1){$ct{$C}=\"R\";}\n		elsif (length \
($r)==2){$ct{$C}=\"D\";}\n		elsif (length ($r)==3)\
{$ct{$C}=\"P\";}\n		else \n		  {\n		    myexit(flu\
sh_error(\"ERROR: Could not read RES_ID field in f\
ile $file\"));\n		  }\n	      }\n	  }\n      }\n  \
  close ($F);\n    return %ct;\n  }\n   \n    \n\n\
\n\nsub fill_threelett_RNA\n{\n\n	my %threelett=(\\
n	'A', '  A',\n	'T', '  T',\n	'U', '  U',\n	'C', '\
  C',\n	'G', '  G',\n	'I', '  I', #Inosine\n	);\n	\
\n	return %threelett;\n\n}\n\n\nsub fill_onelett_R\
NA\n{\n	my   %onelett=(\n	'  A' => 'A',\n	'  T' =>\
 'T',\n	'  U' => 'U',\n	'  C' => 'C',\n	'  G' => '\
G',\n	'CSL' => 'X',\n	'UMS' => 'X',\n	'  I' => 'I'\
,\n	'A' => 'A',\n	'T' => 'T',\n	'U' => 'U',\n	'C' \
=> 'C',\n	'G' => 'G',\n	'I' => 'I',\n	);\n\n	retur\
n %onelett;\n\n}\n\n\nsub fill_onelett_DNA\n{\n	my\
   %onelett=(\n	' DA', 'A',\n	' DT', 'T',\n	' DC',\
 'C',\n	' DG', 'G',\n	'DA', 'A',\n	'DT', 'T',\n	'D\
C', 'C',\n	'DG', 'G',\n	);\n\n	return %onelett;\n\\
n}\n\nsub fill_threelett_DNA\n{\n\n	my %threelett=\
(\n	'A', ' DA',\n	'T', ' DT',\n	'C', ' DC',\n	'G',\
 ' DG',\n	);\n\n	return %threelett;\n\n}\n\n\n\n\n\
sub fill_threelett_prot\n{  \n  my %threelett;\n\n\
  %threelett=(\n'A', 'ALA',\n'C', 'CYS',\n'D', 'AS\
P',\n'E', 'GLU',\n'F', 'PHE',\n'G', 'GLY',\n'H', '\
HIS',\n'I', 'ILE',\n'K', 'LYS',\n'L', 'LEU',\n'N',\
 'ASN',\n'M', 'MET',\n'P', 'PRO',\n'Q', 'GLN',\n'R\
', 'ARG',\n'S', 'SER',\n'T', 'THR',\n'V', 'VAL',\n\
'W', 'TRP',\n'Y', 'TYR',\n);\n\nreturn %threelett;\
\n\n\n}\n\nsub fill_onelett_prot\n{\n    my %onele\
tt;\n    \n    %onelett=(\n\n'10A', 'X',\n'11O', '\
X',\n'12A', 'X',\n'13P', 'X',\n'13R', 'X',\n'13S',\
 'X',\n'14W', 'X',\n'15P', 'X',\n'16A', 'X',\n'16G\
', 'X',\n'1AN', 'X',\n'1AP', 'X',\n'1AR', 'X',\n'1\
BH', 'X',\n'1BO', 'X',\n'1C5', 'X',\n'1CU', 'X',\n\
'1DA', 'X',\n'1GL', 'X',\n'1GN', 'X',\n'1IN', 'X',\
\n'1LU', 'L',\n'1MA', 'X',\n'1MC', 'X',\n'1MG', 'X\
',\n'1MZ', 'X',\n'1NA', 'X',\n'1NB', 'X',\n'1NI', \
'X',\n'1PA', 'A',\n'1PC', 'X',\n'1PE', 'X',\n'1PG'\
, 'X',\n'1PI', 'A',\n'1PM', 'X',\n'1PN', 'X',\n'1P\
U', 'X',\n'1PY', 'X',\n'1UN', 'X',\n'24T', 'X',\n'\
25T', 'X',\n'26P', 'X',\n'2AB', 'X',\n'2AM', 'X',\\
n'2AN', 'X',\n'2AP', 'X',\n'2AR', 'X',\n'2AS', 'D'\
,\n'2BL', 'X',\n'2BM', 'X',\n'2CP', 'X',\n'2DA', '\
X',\n'2DG', 'X',\n'2DP', 'X',\n'2DT', 'X',\n'2EP',\
 'X',\n'2EZ', 'X',\n'2FG', 'X',\n'2FL', 'X',\n'2FP\
', 'X',\n'2FU', 'X',\n'2GL', 'X',\n'2GP', 'X',\n'2\
HP', 'X',\n'2IB', 'X',\n'2IP', 'X',\n'2LU', 'L',\n\
'2MA', 'X',\n'2MD', 'X',\n'2ME', 'X',\n'2MG', 'X',\
\n'2ML', 'L',\n'2MO', 'X',\n'2MR', 'R',\n'2MU', 'X\
',\n'2MZ', 'X',\n'2NO', 'X',\n'2NP', 'X',\n'2OG', \
'X',\n'2PA', 'X',\n'2PC', 'X',\n'2PE', 'X',\n'2PG'\
, 'X',\n'2PH', 'X',\n'2PI', 'X',\n'2PL', 'X',\n'2P\
P', 'X',\n'2PU', 'X',\n'2SI', 'X',\n'2TB', 'X',\n'\
34C', 'X',\n'35G', 'X',\n'3AA', 'X',\n'3AD', 'X',\\
n'3AH', 'H',\n'3AN', 'X',\n'3AP', 'X',\n'3AT', 'X'\
,\n'3BT', 'X',\n'3CH', 'X',\n'3CN', 'X',\n'3CO', '\
X',\n'3CP', 'X',\n'3DR', 'X',\n'3EP', 'X',\n'3FM',\
 'X',\n'3GA', 'X',\n'3GP', 'X',\n'3HB', 'X',\n'3HC\
', 'X',\n'3HP', 'X',\n'3IB', 'X',\n'3ID', 'X',\n'3\
IN', 'X',\n'3MA', 'X',\n'3MB', 'X',\n'3MC', 'X',\n\
'3MD', 'D',\n'3MF', 'X',\n'3MP', 'X',\n'3MT', 'X',\
\n'3OL', 'X',\n'3PA', 'X',\n'3PG', 'X',\n'3PO', 'X\
',\n'3PP', 'X',\n'3PY', 'X',\n'49A', 'X',\n'4AB', \
'X',\n'4AM', 'X',\n'4AN', 'X',\n'4AP', 'X',\n'4BA'\
, 'X',\n'4BT', 'X',\n'4CA', 'X',\n'4CO', 'X',\n'4H\
P', 'X',\n'4IP', 'X',\n'4MO', 'X',\n'4MV', 'X',\n'\
4MZ', 'X',\n'4NC', 'X',\n'4NP', 'X',\n'4OX', 'X',\\
n'4PB', 'X',\n'4PN', 'X',\n'4PP', 'X',\n'4SC', 'X'\
,\n'4SU', 'X',\n'4TB', 'X',\n'55C', 'X',\n'5AD', '\
X',\n'5AN', 'X',\n'5AT', 'X',\n'5CM', 'X',\n'5GP',\
 'X',\n'5HP', 'E',\n'5HT', 'X',\n'5IT', 'X',\n'5IU\
', 'X',\n'5MB', 'X',\n'5MC', 'X',\n'5MD', 'X',\n'5\
MP', 'X',\n'5MU', 'X',\n'5NC', 'X',\n'5OB', 'X',\n\
'5PA', 'X',\n'5PV', 'X',\n'6AB', 'X',\n'6CT', 'X',\
\n'6HA', 'X',\n'6HC', 'X',\n'6HG', 'X',\n'6HT', 'X\
',\n'6IN', 'X',\n'6MO', 'X',\n'6MP', 'X',\n'6PG', \
'X',\n'6WO', 'X',\n'70U', 'X',\n'7DG', 'X',\n'7HP'\
, 'X',\n'7I2', 'X',\n'7MG', 'X',\n'7MQ', 'X',\n'7N\
I', 'X',\n'87Y', 'X',\n'8AD', 'X',\n'8BR', 'X',\n'\
8IG', 'X',\n'8IN', 'X',\n'8OG', 'X',\n'95A', 'X',\\
n'9AD', 'X',\n'9AM', 'X',\n'9AP', 'X',\n'9DG', 'X'\
,\n'9DI', 'X',\n'9HX', 'X',\n'9OH', 'X',\n'9TA', '\
X',\n'A12', 'X',\n'A15', 'X',\n'A23', 'X',\n'A24',\
 'X',\n'A26', 'X',\n'A2G', 'X',\n'A2P', 'X',\n'A32\
', 'X',\n'A3P', 'X',\n'A4P', 'X',\n'A5P', 'X',\n'A\
70', 'X',\n'A76', 'X',\n'A77', 'X',\n'A78', 'X',\n\
'A79', 'X',\n'A80', 'X',\n'A85', 'X',\n'A88', 'X',\
\n'A9A', 'X',\n'AA3', 'X',\n'AA4', 'X',\n'AA6', 'X\
',\n'AAA', 'X',\n'AAB', 'X',\n'AAC', 'X',\n'AAE', \
'X',\n'AAG', 'R',\n'AAH', 'X',\n'AAM', 'X',\n'AAN'\
, 'X',\n'AAP', 'X',\n'AAR', 'R',\n'AAS', 'X',\n'AA\
T', 'X',\n'ABA', 'X',\n'ABC', 'X',\n'ABD', 'X',\n'\
ABE', 'X',\n'ABH', 'X',\n'ABI', 'X',\n'ABK', 'X',\\
n'ABM', 'X',\n'ABN', 'X',\n'ABP', 'X',\n'ABR', 'X'\
,\n'ABS', 'X',\n'ABU', 'X',\n'AC1', 'X',\n'AC2', '\
X',\n'ACA', 'X',\n'ACB', 'D',\n'ACC', 'C',\n'ACD',\
 'X',\n'ACE', 'X',\n'ACH', 'X',\n'ACI', 'X',\n'ACL\
', 'R',\n'ACM', 'X',\n'ACN', 'X',\n'ACO', 'X',\n'A\
CP', 'X',\n'ACQ', 'X',\n'ACR', 'X',\n'ACS', 'X',\n\
'ACT', 'X',\n'ACV', 'V',\n'ACX', 'X',\n'ACY', 'X',\
\n'AD2', 'X',\n'AD3', 'X',\n'ADC', 'X',\n'ADD', 'X\
',\n'ADE', 'X',\n'ADH', 'X',\n'ADI', 'X',\n'ADM', \
'X',\n'ADN', 'X',\n'ADP', 'X',\n'ADQ', 'X',\n'ADR'\
, 'X',\n'ADS', 'X',\n'ADT', 'X',\n'ADU', 'X',\n'AD\
W', 'X',\n'ADX', 'X',\n'AE2', 'X',\n'AEA', 'X',\n'\
AEB', 'X',\n'AEI', 'D',\n'AEN', 'X',\n'AET', 'T',\\
n'AF1', 'X',\n'AF3', 'X',\n'AFA', 'D',\n'AFP', 'X'\
,\n'AG7', 'X',\n'AGB', 'X',\n'AGF', 'X',\n'AGL', '\
X',\n'AGM', 'R',\n'AGN', 'X',\n'AGP', 'X',\n'AGS',\
 'X',\n'AGU', 'X',\n'AH0', 'X',\n'AH1', 'X',\n'AHA\
', 'X',\n'AHB', 'D',\n'AHC', 'X',\n'AHF', 'X',\n'A\
HG', 'X',\n'AHH', 'X',\n'AHM', 'X',\n'AHO', 'X',\n\
'AHP', 'X',\n'AHS', 'X',\n'AHT', 'Y',\n'AHU', 'X',\
\n'AHX', 'X',\n'AI1', 'X',\n'AI2', 'X',\n'AIB', 'X\
',\n'AIC', 'X',\n'AIM', 'X',\n'AIP', 'X',\n'AIQ', \
'X',\n'AIR', 'X',\n'AJ3', 'X',\n'AKB', 'X',\n'AKG'\
, 'X',\n'AKR', 'X',\n'AL1', 'X',\n'AL2', 'X',\n'AL\
3', 'X',\n'AL4', 'X',\n'AL5', 'X',\n'AL6', 'X',\n'\
AL7', 'X',\n'AL8', 'X',\n'AL9', 'X',\n'ALA', 'A',\\
n'ALB', 'X',\n'ALC', 'X',\n'ALD', 'L',\n'ALE', 'X'\
,\n'ALF', 'X',\n'ALG', 'X',\n'ALL', 'X',\n'ALM', '\
A',\n'ALN', 'A',\n'ALO', 'T',\n'ALP', 'X',\n'ALQ',\
 'X',\n'ALR', 'X',\n'ALS', 'X',\n'ALT', 'A',\n'ALY\
', 'K',\n'ALZ', 'X',\n'AMA', 'X',\n'AMB', 'X',\n'A\
MC', 'X',\n'AMD', 'X',\n'AMG', 'X',\n'AMH', 'X',\n\
'AMI', 'X',\n'AML', 'X',\n'AMN', 'X',\n'AMO', 'X',\
\n'AMP', 'X',\n'AMQ', 'X',\n'AMR', 'X',\n'AMS', 'X\
',\n'AMT', 'X',\n'AMU', 'X',\n'AMW', 'X',\n'AMX', \
'X',\n'AMY', 'X',\n'ANA', 'X',\n'ANB', 'X',\n'ANC'\
, 'X',\n'AND', 'X',\n'ANE', 'X',\n'ANI', 'X',\n'AN\
L', 'X',\n'ANO', 'X',\n'ANP', 'X',\n'ANS', 'X',\n'\
ANT', 'X',\n'AOE', 'X',\n'AOP', 'X',\n'AP1', 'X',\\
n'AP2', 'X',\n'AP3', 'X',\n'AP4', 'X',\n'AP5', 'X'\
,\n'AP6', 'X',\n'APA', 'X',\n'APB', 'X',\n'APC', '\
X',\n'APE', 'F',\n'APF', 'X',\n'APG', 'X',\n'APH',\
 'A',\n'API', 'X',\n'APL', 'X',\n'APM', 'X',\n'APN\
', 'G',\n'APP', 'X',\n'APQ', 'X',\n'APR', 'X',\n'A\
PS', 'X',\n'APT', 'X',\n'APU', 'X',\n'APX', 'X',\n\
'APY', 'X',\n'APZ', 'X',\n'AQS', 'X',\n'AR1', 'X',\
\n'AR2', 'X',\n'ARA', 'X',\n'ARB', 'X',\n'ARC', 'X\
',\n'ARD', 'X',\n'ARG', 'R',\n'ARH', 'X',\n'ARI', \
'X',\n'ARM', 'R',\n'ARN', 'X',\n'ARO', 'R',\n'ARP'\
, 'X',\n'ARQ', 'X',\n'ARS', 'X',\n'AS1', 'R',\n'AS\
2', 'X',\n'ASA', 'D',\n'ASB', 'D',\n'ASC', 'X',\n'\
ASD', 'X',\n'ASE', 'X',\n'ASF', 'X',\n'ASI', 'X',\\
n'ASK', 'D',\n'ASL', 'X',\n'ASM', 'N',\n'ASO', 'X'\
,\n'ASP', 'D',\n'ASQ', 'X',\n'ASU', 'X',\n'ATA', '\
X',\n'ATC', 'X',\n'ATD', 'X',\n'ATF', 'X',\n'ATG',\
 'X',\n'ATH', 'X',\n'ATM', 'X',\n'ATO', 'X',\n'ATP\
', 'X',\n'ATQ', 'X',\n'ATR', 'X',\n'ATT', 'X',\n'A\
TY', 'X',\n'ATZ', 'X',\n'AUC', 'X',\n'AUR', 'X',\n\
'AVG', 'X',\n'AXP', 'X',\n'AYA', 'A',\n'AZ2', 'X',\
\n'AZA', 'X',\n'AZC', 'X',\n'AZD', 'X',\n'AZE', 'X\
',\n'AZI', 'X',\n'AZL', 'X',\n'AZM', 'X',\n'AZR', \
'X',\n'AZT', 'X',\n'B12', 'X',\n'B1F', 'F',\n'B2A'\
, 'A',\n'B2F', 'F',\n'B2I', 'I',\n'B2V', 'V',\n'B3\
I', 'X',\n'B3P', 'X',\n'B7G', 'X',\n'B96', 'X',\n'\
B9A', 'X',\n'BA1', 'X',\n'BAA', 'X',\n'BAB', 'X',\\
n'BAC', 'X',\n'BAF', 'X',\n'BAH', 'X',\n'BAI', 'X'\
,\n'BAK', 'X',\n'BAL', 'A',\n'BAM', 'X',\n'BAO', '\
X',\n'BAP', 'X',\n'BAR', 'X',\n'BAS', 'X',\n'BAT',\
 'F',\n'BAY', 'X',\n'BAZ', 'X',\n'BB1', 'X',\n'BB2\
', 'X',\n'BBA', 'X',\n'BBH', 'X',\n'BBS', 'X',\n'B\
BT', 'X',\n'BBZ', 'X',\n'BCA', 'X',\n'BCB', 'X',\n\
'BCC', 'X',\n'BCD', 'X',\n'BCL', 'X',\n'BCN', 'X',\
\n'BCR', 'X',\n'BCS', 'C',\n'BCT', 'X',\n'BCY', 'X\
',\n'BCZ', 'X',\n'BDA', 'X',\n'BDG', 'X',\n'BDK', \
'X',\n'BDM', 'X',\n'BDN', 'X',\n'BDS', 'X',\n'BE1'\
, 'X',\n'BE2', 'X',\n'BEA', 'X',\n'BEF', 'X',\n'BE\
N', 'X',\n'BEO', 'X',\n'BEP', 'X',\n'BER', 'X',\n'\
BES', 'X',\n'BET', 'X',\n'BEZ', 'X',\n'BF2', 'X',\\
n'BFA', 'X',\n'BFD', 'X',\n'BFP', 'X',\n'BFS', 'X'\
,\n'BFU', 'X',\n'BG6', 'X',\n'BGF', 'X',\n'BGG', '\
X',\n'BGL', 'X',\n'BGN', 'X',\n'BGP', 'X',\n'BGX',\
 'X',\n'BH4', 'X',\n'BHA', 'X',\n'BHC', 'X',\n'BHD\
', 'D',\n'BHO', 'X',\n'BHS', 'X',\n'BIC', 'X',\n'B\
IN', 'X',\n'BIO', 'X',\n'BIP', 'X',\n'BIS', 'X',\n\
'BIZ', 'X',\n'BJH', 'X',\n'BJI', 'X',\n'BJP', 'X',\
\n'BLA', 'X',\n'BLB', 'X',\n'BLE', 'L',\n'BLG', 'P\
',\n'BLI', 'X',\n'BLM', 'X',\n'BLV', 'X',\n'BLY', \
'K',\n'BM1', 'X',\n'BM2', 'X',\n'BM5', 'X',\n'BM9'\
, 'X',\n'BMA', 'X',\n'BMD', 'X',\n'BME', 'X',\n'BM\
P', 'X',\n'BMQ', 'X',\n'BMS', 'X',\n'BMT', 'T',\n'\
BMU', 'X',\n'BMY', 'X',\n'BMZ', 'X',\n'BNA', 'X',\\
n'BNG', 'X',\n'BNI', 'X',\n'BNN', 'F',\n'BNO', 'L'\
,\n'BNS', 'X',\n'BNZ', 'X',\n'BO3', 'X',\n'BO4', '\
X',\n'BOC', 'X',\n'BOG', 'X',\n'BOM', 'X',\n'BOT',\
 'X',\n'BOX', 'X',\n'BOZ', 'X',\n'BPA', 'X',\n'BPB\
', 'X',\n'BPD', 'X',\n'BPG', 'X',\n'BPH', 'X',\n'B\
PI', 'X',\n'BPJ', 'X',\n'BPM', 'X',\n'BPN', 'X',\n\
'BPO', 'X',\n'BPP', 'X',\n'BPT', 'X',\n'BPY', 'X',\
\n'BRB', 'X',\n'BRC', 'X',\n'BRE', 'X',\n'BRI', 'X\
',\n'BRL', 'X',\n'BRM', 'X',\n'BRN', 'X',\n'BRO', \
'X',\n'BRS', 'X',\n'BRU', 'X',\n'BRZ', 'X',\n'BSB'\
, 'X',\n'BSI', 'X',\n'BSP', 'X',\n'BT1', 'X',\n'BT\
2', 'X',\n'BT3', 'X',\n'BTA', 'L',\n'BTB', 'X',\n'\
BTC', 'C',\n'BTD', 'X',\n'BTN', 'X',\n'BTP', 'X',\\
n'BTR', 'W',\n'BU1', 'X',\n'BUA', 'X',\n'BUB', 'X'\
,\n'BUC', 'X',\n'BUG', 'X',\n'BUL', 'X',\n'BUM', '\
X',\n'BUQ', 'X',\n'BUT', 'X',\n'BVD', 'X',\n'BX3',\
 'X',\n'BYS', 'X',\n'BZ1', 'X',\n'BZA', 'X',\n'BZB\
', 'X',\n'BZC', 'X',\n'BZD', 'X',\n'BZF', 'X',\n'B\
ZI', 'X',\n'BZM', 'X',\n'BZO', 'X',\n'BZP', 'X',\n\
'BZQ', 'X',\n'BZS', 'X',\n'BZT', 'X',\n'C02', 'X',\
\n'C11', 'X',\n'C1O', 'X',\n'C20', 'X',\n'C24', 'X\
',\n'C2F', 'X',\n'C2O', 'X',\n'C2P', 'X',\n'C3M', \
'X',\n'C3P', 'X',\n'C3X', 'X',\n'C48', 'X',\n'C4M'\
, 'X',\n'C4X', 'X',\n'C5C', 'X',\n'C5M', 'X',\n'C5\
P', 'X',\n'C5X', 'X',\n'C60', 'X',\n'C6C', 'X',\n'\
C6M', 'X',\n'C78', 'X',\n'C8E', 'X',\n'CA3', 'X',\\
n'CA5', 'X',\n'CAA', 'X',\n'CAB', 'X',\n'CAC', 'X'\
,\n'CAD', 'X',\n'CAF', 'C',\n'CAG', 'X',\n'CAH', '\
X',\n'CAL', 'X',\n'CAM', 'X',\n'CAN', 'X',\n'CAO',\
 'X',\n'CAP', 'X',\n'CAQ', 'X',\n'CAR', 'X',\n'CAS\
', 'C',\n'CAT', 'X',\n'CAV', 'X',\n'CAY', 'C',\n'C\
AZ', 'X',\n'CB3', 'X',\n'CB4', 'X',\n'CBA', 'X',\n\
'CBD', 'X',\n'CBG', 'X',\n'CBI', 'X',\n'CBL', 'X',\
\n'CBM', 'X',\n'CBN', 'X',\n'CBO', 'X',\n'CBP', 'X\
',\n'CBS', 'X',\n'CBX', 'X',\n'CBZ', 'X',\n'CC0', \
'X',\n'CC1', 'X',\n'CCC', 'X',\n'CCH', 'X',\n'CCI'\
, 'X',\n'CCM', 'X',\n'CCN', 'X',\n'CCO', 'X',\n'CC\
P', 'X',\n'CCR', 'X',\n'CCS', 'C',\n'CCV', 'X',\n'\
CCY', 'X',\n'CD1', 'X',\n'CDC', 'X',\n'CDE', 'X',\\
n'CDF', 'X',\n'CDI', 'X',\n'CDL', 'X',\n'CDM', 'X'\
,\n'CDP', 'X',\n'CDR', 'X',\n'CDU', 'X',\n'CE1', '\
X',\n'CEA', 'C',\n'CEB', 'X',\n'CEC', 'X',\n'CED',\
 'X',\n'CEF', 'X',\n'CEH', 'X',\n'CEM', 'X',\n'CEO\
', 'X',\n'CEP', 'X',\n'CEQ', 'X',\n'CER', 'X',\n'C\
ES', 'G',\n'CET', 'X',\n'CFC', 'X',\n'CFF', 'X',\n\
'CFM', 'X',\n'CFO', 'X',\n'CFP', 'X',\n'CFS', 'X',\
\n'CFX', 'X',\n'CGN', 'X',\n'CGP', 'X',\n'CGS', 'X\
',\n'CGU', 'E',\n'CH2', 'X',\n'CH3', 'X',\n'CHA', \
'X',\n'CHB', 'X',\n'CHD', 'X',\n'CHF', 'X',\n'CHG'\
, 'G',\n'CHI', 'X',\n'CHN', 'X',\n'CHO', 'X',\n'CH\
P', 'G',\n'CHR', 'X',\n'CHS', 'F',\n'CHT', 'X',\n'\
CHX', 'X',\n'CIC', 'X',\n'CIN', 'X',\n'CIP', 'X',\\
n'CIR', 'X',\n'CIT', 'X',\n'CIU', 'X',\n'CKI', 'X'\
,\n'CL1', 'X',\n'CL2', 'X',\n'CLA', 'X',\n'CLB', '\
A',\n'CLC', 'S',\n'CLD', 'A',\n'CLE', 'L',\n'CLF',\
 'X',\n'CLK', 'S',\n'CLL', 'X',\n'CLM', 'X',\n'CLN\
', 'X',\n'CLO', 'X',\n'CLP', 'X',\n'CLQ', 'X',\n'C\
LR', 'X',\n'CLS', 'X',\n'CLT', 'X',\n'CLX', 'X',\n\
'CLY', 'X',\n'CMA', 'R',\n'CMC', 'X',\n'CMD', 'X',\
\n'CME', 'C',\n'CMG', 'X',\n'CMK', 'X',\n'CMN', 'X\
',\n'CMO', 'X',\n'CMP', 'X',\n'CMR', 'X',\n'CMS', \
'X',\n'CMT', 'C',\n'CMX', 'X',\n'CNA', 'X',\n'CNC'\
, 'X',\n'CND', 'X',\n'CNH', 'X',\n'CNM', 'X',\n'CN\
N', 'X',\n'CNO', 'X',\n'CNP', 'X',\n'CO2', 'X',\n'\
CO3', 'X',\n'CO5', 'X',\n'CO8', 'X',\n'COA', 'X',\\
n'COB', 'X',\n'COC', 'X',\n'COD', 'X',\n'COE', 'X'\
,\n'COF', 'X',\n'COH', 'X',\n'COI', 'X',\n'COJ', '\
X',\n'COL', 'X',\n'COM', 'X',\n'CON', 'X',\n'COP',\
 'X',\n'COR', 'X',\n'COS', 'X',\n'COT', 'X',\n'COY\
', 'X',\n'CP1', 'G',\n'CP2', 'X',\n'CP4', 'X',\n'C\
PA', 'X',\n'CPB', 'X',\n'CPC', 'X',\n'CPD', 'X',\n\
'CPG', 'X',\n'CPH', 'X',\n'CPI', 'X',\n'CPM', 'X',\
\n'CPN', 'G',\n'CPO', 'X',\n'CPP', 'X',\n'CPQ', 'X\
',\n'CPR', 'X',\n'CPS', 'X',\n'CPT', 'X',\n'CPU', \
'X',\n'CPV', 'X',\n'CPY', 'X',\n'CR1', 'X',\n'CR6'\
, 'X',\n'CRA', 'X',\n'CRB', 'X',\n'CRC', 'X',\n'CR\
G', 'X',\n'CRH', 'X',\n'CRO', 'T',\n'CRP', 'X',\n'\
CRQ', 'X',\n'CRS', 'X',\n'CRT', 'X',\n'CRY', 'X',\\
n'CSA', 'C',\n'CSB', 'X',\n'CSD', 'C',\n'CSE', 'C'\
,\n'CSH', 'X',\n'CSI', 'X',\n'CSN', 'X',\n'CSO', '\
C',\n'CSP', 'C',\n'CSR', 'C',\n'CSS', 'C',\n'CST',\
 'X',\n'CSW', 'C',\n'CSX', 'C',\n'CSY', 'X',\n'CSZ\
', 'C',\n'CT3', 'X',\n'CTA', 'X',\n'CTB', 'X',\n'C\
TC', 'X',\n'CTD', 'X',\n'CTH', 'T',\n'CTO', 'X',\n\
'CTP', 'X',\n'CTR', 'X',\n'CTS', 'X',\n'CTT', 'X',\
\n'CTY', 'X',\n'CTZ', 'X',\n'CU1', 'X',\n'CUA', 'X\
',\n'CUC', 'X',\n'CUL', 'X',\n'CUO', 'X',\n'CUZ', \
'X',\n'CVI', 'X',\n'CXF', 'X',\n'CXL', 'X',\n'CXM'\
, 'M',\n'CXN', 'X',\n'CXP', 'X',\n'CXS', 'X',\n'CY\
1', 'C',\n'CY3', 'X',\n'CYB', 'X',\n'CYC', 'X',\n'\
CYF', 'C',\n'CYG', 'C',\n'CYH', 'X',\n'CYL', 'X',\\
n'CYM', 'C',\n'CYN', 'X',\n'CYO', 'X',\n'CYP', 'X'\
,\n'CYQ', 'C',\n'CYS', 'C',\n'CYU', 'X',\n'CYY', '\
X',\n'CYZ', 'X',\n'CZH', 'X',\n'CZZ', 'C',\n'D12',\
 'X',\n'D13', 'X',\n'D16', 'X',\n'D18', 'X',\n'D19\
', 'X',\n'D1P', 'X',\n'D24', 'X',\n'D34', 'X',\n'D\
35', 'X',\n'D4D', 'X',\n'D4T', 'X',\n'D6G', 'X',\n\
'DA2', 'R',\n'DA3', 'X',\n'DA6', 'X',\n'DA7', 'X',\
\n'DAA', 'X',\n'DAB', 'X',\n'DAC', 'X',\n'DAD', 'X\
',\n'DAE', 'X',\n'DAF', 'X',\n'DAG', 'X',\n'DAH', \
'A',\n'DAJ', 'X',\n'DAK', 'X',\n'DAL', 'A',\n'DAM'\
, 'A',\n'DAN', 'X',\n'DAO', 'X',\n'DAP', 'X',\n'DA\
Q', 'X',\n'DAR', 'R',\n'DAS', 'D',\n'DAT', 'X',\n'\
DAU', 'X',\n'DAV', 'X',\n'DBA', 'X',\n'DBD', 'X',\\
n'DBF', 'X',\n'DBG', 'X',\n'DBI', 'X',\n'DBV', 'X'\
,\n'DBY', 'Y',\n'DCA', 'X',\n'DCB', 'X',\n'DCE', '\
X',\n'DCF', 'X',\n'DCG', 'X',\n'DCH', 'X',\n'DCI',\
 'I',\n'DCL', 'X',\n'DCM', 'X',\n'DCP', 'X',\n'DCS\
', 'X',\n'DCT', 'X',\n'DCY', 'C',\n'DCZ', 'X',\n'D\
DA', 'X',\n'DDB', 'X',\n'DDC', 'X',\n'DDF', 'X',\n\
'DDG', 'X',\n'DDH', 'X',\n'DDL', 'X',\n'DDM', 'X',\
\n'DDO', 'L',\n'DDP', 'X',\n'DDQ', 'X',\n'DDT', 'Y\
',\n'DDU', 'X',\n'DEA', 'X',\n'DEB', 'X',\n'DEC', \
'X',\n'DEF', 'X',\n'DEL', 'X',\n'DEM', 'X',\n'DEN'\
, 'X',\n'DEP', 'X',\n'DEQ', 'X',\n'DES', 'X',\n'DE\
T', 'X',\n'DFC', 'X',\n'DFG', 'X',\n'DFI', 'X',\n'\
DFL', 'X',\n'DFO', 'X',\n'DFP', 'X',\n'DFR', 'X',\\
n'DFT', 'X',\n'DFV', 'X',\n'DFX', 'X',\n'DG2', 'X'\
,\n'DG3', 'X',\n'DG6', 'X',\n'DGA', 'X',\n'DGD', '\
X',\n'DGG', 'X',\n'DGL', 'E',\n'DGN', 'Q',\n'DGP',\
 'X',\n'DGT', 'X',\n'DGX', 'X',\n'DH2', 'X',\n'DHA\
', 'A',\n'DHB', 'X',\n'DHC', 'X',\n'DHD', 'X',\n'D\
HE', 'X',\n'DHF', 'X',\n'DHG', 'X',\n'DHI', 'H',\n\
'DHL', 'X',\n'DHM', 'X',\n'DHN', 'V',\n'DHP', 'X',\
\n'DHQ', 'X',\n'DHR', 'X',\n'DHS', 'X',\n'DHT', 'X\
',\n'DHU', 'X',\n'DHY', 'X',\n'DHZ', 'X',\n'DI2', \
'X',\n'DI3', 'G',\n'DI4', 'X',\n'DI5', 'X',\n'DIA'\
, 'X',\n'DIC', 'X',\n'DIF', 'X',\n'DIG', 'X',\n'DI\
I', 'X',\n'DIL', 'I',\n'DIM', 'X',\n'DIO', 'X',\n'\
DIP', 'X',\n'DIQ', 'X',\n'DIS', 'X',\n'DIT', 'X',\\
n'DIV', 'V',\n'DIX', 'X',\n'DIY', 'X',\n'DKA', 'X'\
,\n'DLA', 'X',\n'DLE', 'L',\n'DLF', 'X',\n'DLS', '\
K',\n'DLY', 'K',\n'DM1', 'X',\n'DM2', 'X',\n'DM3',\
 'X',\n'DM4', 'X',\n'DM5', 'X',\n'DM6', 'X',\n'DM7\
', 'X',\n'DM8', 'X',\n'DM9', 'X',\n'DMA', 'X',\n'D\
MB', 'X',\n'DMC', 'X',\n'DMD', 'X',\n'DME', 'X',\n\
'DMF', 'X',\n'DMG', 'G',\n'DMH', 'N',\n'DMI', 'X',\
\n'DMJ', 'X',\n'DML', 'X',\n'DMM', 'X',\n'DMN', 'X\
',\n'DMO', 'X',\n'DMP', 'X',\n'DMQ', 'X',\n'DMR', \
'X',\n'DMS', 'X',\n'DMT', 'X',\n'DMV', 'X',\n'DMY'\
, 'X',\n'DNC', 'X',\n'DND', 'X',\n'DNH', 'X',\n'DN\
J', 'X',\n'DNN', 'X',\n'DNP', 'X',\n'DNQ', 'X',\n'\
DNR', 'X',\n'DO2', 'X',\n'DO3', 'X',\n'DOA', 'X',\\
n'DOB', 'X',\n'DOC', 'X',\n'DOH', 'D',\n'DOM', 'X'\
,\n'DOS', 'X',\n'DOX', 'X',\n'DP5', 'X',\n'DP7', '\
X',\n'DPA', 'X',\n'DPC', 'X',\n'DPD', 'X',\n'DPE',\
 'X',\n'DPG', 'X',\n'DPH', 'F',\n'DPM', 'X',\n'DPN\
', 'F',\n'DPO', 'X',\n'DPP', 'X',\n'DPR', 'P',\n'D\
PS', 'X',\n'DPT', 'X',\n'DPX', 'X',\n'DPY', 'X',\n\
'DPZ', 'X',\n'DQH', 'X',\n'DQN', 'X',\n'DR1', 'X',\
\n'DRB', 'X',\n'DRC', 'X',\n'DRI', 'X',\n'DRP', 'X\
',\n'DRT', 'X',\n'DRU', 'X',\n'DSA', 'X',\n'DSB', \
'X',\n'DSC', 'X',\n'DSD', 'X',\n'DSE', 'S',\n'DSI'\
, 'X',\n'DSN', 'S',\n'DSP', 'D',\n'DSR', 'X',\n'DS\
S', 'X',\n'DSX', 'X',\n'DSY', 'X',\n'DTB', 'X',\n'\
DTD', 'X',\n'DTH', 'T',\n'DTN', 'X',\n'DTO', 'X',\\
n'DTP', 'X',\n'DTQ', 'X',\n'DTR', 'W',\n'DTT', 'X'\
,\n'DTY', 'Y',\n'DUD', 'X',\n'DUO', 'X',\n'DUR', '\
X',\n'DUT', 'X',\n'DVA', 'V',\n'DVR', 'X',\n'DX9',\
 'X',\n'DXA', 'X',\n'DXB', 'X',\n'DXC', 'X',\n'DXG\
', 'X',\n'DXX', 'X',\n'DZF', 'X',\n'E09', 'X',\n'E\
20', 'X',\n'E2P', 'X',\n'E3G', 'X',\n'E4N', 'X',\n\
'E4P', 'X',\n'E64', 'X',\n'E6C', 'X',\n'E96', 'X',\
\n'E97', 'X',\n'EA2', 'X',\n'EAA', 'X',\n'EAP', 'X\
',\n'EBP', 'X',\n'EBW', 'X',\n'ECO', 'X',\n'EDA', \
'X',\n'EDC', 'X',\n'EDE', 'X',\n'EDO', 'X',\n'EDR'\
, 'X',\n'EEB', 'X',\n'EEE', 'X',\n'EFC', 'X',\n'EF\
Z', 'X',\n'EG1', 'X',\n'EG2', 'X',\n'EG3', 'X',\n'\
EGC', 'X',\n'EGL', 'X',\n'EHP', 'A',\n'EIC', 'X',\\
n'EJT', 'X',\n'ELA', 'X',\n'EMB', 'X',\n'EMC', 'X'\
,\n'EMD', 'X',\n'EMM', 'X',\n'EMO', 'X',\n'EMP', '\
X',\n'EMR', 'X',\n'ENA', 'X',\n'ENC', 'X',\n'ENH',\
 'X',\n'ENO', 'X',\n'ENP', 'X',\n'EOA', 'X',\n'EOH\
', 'X',\n'EOT', 'X',\n'EOX', 'X',\n'EPA', 'X',\n'E\
PE', 'X',\n'EPH', 'X',\n'EPI', 'X',\n'EPN', 'X',\n\
'EPO', 'X',\n'EPT', 'X',\n'EPU', 'X',\n'EPX', 'X',\
\n'EPY', 'X',\n'EQI', 'X',\n'EQP', 'X',\n'EQU', 'X\
',\n'ERG', 'X',\n'ERI', 'X',\n'ERY', 'X',\n'ESC', \
'X',\n'ESD', 'X',\n'ESI', 'X',\n'ESO', 'X',\n'ESP'\
, 'X',\n'EST', 'X',\n'ESX', 'X',\n'ETA', 'X',\n'ET\
C', 'X',\n'ETD', 'X',\n'ETF', 'X',\n'ETH', 'X',\n'\
ETI', 'X',\n'ETN', 'X',\n'ETO', 'X',\n'ETP', 'X',\\
n'ETR', 'X',\n'ETS', 'X',\n'ETY', 'X',\n'EU3', 'X'\
,\n'EUG', 'X',\n'EYS', 'C',\n'F09', 'X',\n'F2B', '\
X',\n'F3S', 'X',\n'F42', 'X',\n'F43', 'X',\n'F4S',\
 'X',\n'F6B', 'X',\n'F6P', 'X',\n'F89', 'X',\n'FA1\
', 'X',\n'FA5', 'F',\n'FAA', 'X',\n'FAB', 'X',\n'F\
AC', 'X',\n'FAD', 'X',\n'FAF', 'X',\n'FAG', 'X',\n\
'FAM', 'X',\n'FAR', 'X',\n'FAS', 'X',\n'FAT', 'X',\
\n'FBA', 'X',\n'FBE', 'X',\n'FBI', 'X',\n'FBP', 'X\
',\n'FBQ', 'X',\n'FBS', 'X',\n'FBT', 'X',\n'FBU', \
'X',\n'FCA', 'X',\n'FCB', 'X',\n'FCI', 'X',\n'FCN'\
, 'X',\n'FCO', 'X',\n'FCR', 'X',\n'FCT', 'X',\n'FC\
X', 'X',\n'FCY', 'C',\n'FD1', 'F',\n'FD2', 'F',\n'\
FD3', 'F',\n'FD4', 'F',\n'FDA', 'X',\n'FDC', 'X',\\
n'FDI', 'X',\n'FDP', 'X',\n'FDS', 'X',\n'FE2', 'X'\
,\n'FEA', 'X',\n'FEL', 'X',\n'FEM', 'X',\n'FEN', '\
X',\n'FEO', 'X',\n'FEP', 'X',\n'FER', 'X',\n'FES',\
 'X',\n'FFB', 'X',\n'FFC', 'X',\n'FFF', 'X',\n'FFO\
', 'X',\n'FGL', 'G',\n'FHB', 'X',\n'FHC', 'X',\n'F\
HP', 'X',\n'FHU', 'X',\n'FID', 'X',\n'FII', 'X',\n\
'FIP', 'X',\n'FK5', 'X',\n'FKA', 'X',\n'FKI', 'X',\
\n'FKP', 'X',\n'FL2', 'X',\n'FL9', 'X',\n'FLA', 'A\
',\n'FLC', 'X',\n'FLD', 'X',\n'FLE', 'L',\n'FLF', \
'X',\n'FLO', 'X',\n'FLP', 'X',\n'FLT', 'Y',\n'FLU'\
, 'X',\n'FLX', 'X',\n'FM1', 'X',\n'FM2', 'X',\n'FM\
A', 'X',\n'FMB', 'X',\n'FMC', 'X',\n'FME', 'M',\n'\
FMN', 'X',\n'FMP', 'X',\n'FMR', 'X',\n'FMS', 'X',\\
n'FMT', 'X',\n'FNE', 'X',\n'FNP', 'X',\n'FNS', 'X'\
,\n'FOC', 'X',\n'FOE', 'X',\n'FOG', 'F',\n'FOH', '\
X',\n'FOK', 'X',\n'FOL', 'X',\n'FON', 'X',\n'FOP',\
 'X',\n'FOR', 'X',\n'FOS', 'X',\n'FPA', 'X',\n'FPC\
', 'X',\n'FPI', 'X',\n'FPO', 'X',\n'FPP', 'X',\n'F\
PT', 'X',\n'FQP', 'X',\n'FRA', 'X',\n'FRD', 'F',\n\
'FRU', 'X',\n'FS3', 'X',\n'FS4', 'X',\n'FSB', 'X',\
\n'FSO', 'X',\n'FSX', 'X',\n'FTC', 'X',\n'FTP', 'X\
',\n'FTR', 'W',\n'FTT', 'X',\n'FTY', 'Y',\n'FUA', \
'X',\n'FUC', 'X',\n'FUM', 'X',\n'FUP', 'X',\n'FVF'\
, 'X',\n'FXP', 'X',\n'FXV', 'X',\n'FYA', 'F',\n'G1\
6', 'X',\n'G1P', 'X',\n'G20', 'X',\n'G21', 'X',\n'\
G23', 'X',\n'G26', 'X',\n'G28', 'X',\n'G2F', 'X',\\
n'G37', 'X',\n'G39', 'X',\n'G3H', 'X',\n'G3P', 'X'\
,\n'G4D', 'X',\n'G6D', 'X',\n'G6P', 'X',\n'G6Q', '\
X',\n'G7M', 'X',\n'GA2', 'X',\n'GAA', 'X',\n'GAB',\
 'X',\n'GAC', 'X',\n'GAI', 'X',\n'GAL', 'X',\n'GAM\
', 'X',\n'GAN', 'X',\n'GAO', 'X',\n'GAP', 'X',\n'G\
AR', 'G',\n'GAS', 'X',\n'GAT', 'X',\n'GBC', 'X',\n\
'GBI', 'X',\n'GBP', 'X',\n'GBS', 'X',\n'GBX', 'X',\
\n'GC4', 'X',\n'GCA', 'X',\n'GCD', 'X',\n'GCG', 'G\
',\n'GCH', 'G',\n'GCK', 'X',\n'GCL', 'X',\n'GCM', \
'X',\n'GCN', 'X',\n'GCO', 'X',\n'GCP', 'X',\n'GCR'\
, 'X',\n'GCS', 'X',\n'GCU', 'X',\n'GD3', 'X',\n'GD\
B', 'X',\n'GDM', 'X',\n'GDN', 'X',\n'GDP', 'X',\n'\
GDS', 'X',\n'GDU', 'X',\n'GE1', 'X',\n'GE2', 'X',\\
n'GE3', 'X',\n'GEA', 'X',\n'GEL', 'X',\n'GEM', 'X'\
,\n'GEN', 'X',\n'GEP', 'X',\n'GER', 'X',\n'GFP', '\
X',\n'GGB', 'X',\n'GGL', 'E',\n'GGP', 'X',\n'GHP',\
 'G',\n'GIP', 'X',\n'GIS', 'X',\n'GKR', 'X',\n'GL2\
', 'X',\n'GL3', 'G',\n'GL4', 'X',\n'GL5', 'X',\n'G\
L7', 'X',\n'GL9', 'X',\n'GLA', 'X',\n'GLB', 'X',\n\
'GLC', 'X',\n'GLD', 'X',\n'GLE', 'X',\n'GLF', 'X',\
\n'GLG', 'X',\n'GLH', 'Q',\n'GLI', 'X',\n'GLL', 'X\
',\n'GLM', 'G',\n'GLN', 'Q',\n'GLO', 'X',\n'GLP', \
'X',\n'GLR', 'X',\n'GLS', 'X',\n'GLT', 'X',\n'GLU'\
, 'E',\n'GLV', 'X',\n'GLW', 'X',\n'GLY', 'G',\n'GL\
Z', 'X',\n'GM1', 'X',\n'GMA', 'X',\n'GMC', 'X',\n'\
GMH', 'X',\n'GMP', 'X',\n'GMY', 'X',\n'GN7', 'X',\\
n'GNA', 'X',\n'GNB', 'X',\n'GNH', 'X',\n'GNP', 'X'\
,\n'GNT', 'X',\n'GOA', 'X',\n'GOL', 'X',\n'GOX', '\
X',\n'GP1', 'X',\n'GP3', 'X',\n'GP4', 'X',\n'GP6',\
 'X',\n'GP8', 'X',\n'GPB', 'E',\n'GPC', 'X',\n'GPE\
', 'X',\n'GPG', 'X',\n'GPI', 'X',\n'GPJ', 'X',\n'G\
PL', 'K',\n'GPM', 'X',\n'GPN', 'G',\n'GPP', 'X',\n\
'GPR', 'X',\n'GPS', 'X',\n'GPX', 'X',\n'GR1', 'X',\
\n'GR3', 'X',\n'GR4', 'X',\n'GSA', 'X',\n'GSB', 'X\
',\n'GSC', 'G',\n'GSE', 'S',\n'GSH', 'X',\n'GSP', \
'X',\n'GSR', 'X',\n'GSS', 'X',\n'GT9', 'C',\n'GTA'\
, 'X',\n'GTB', 'X',\n'GTD', 'X',\n'GTE', 'X',\n'GT\
H', 'T',\n'GTN', 'X',\n'GTO', 'X',\n'GTP', 'X',\n'\
GTR', 'X',\n'GTS', 'X',\n'GTT', 'X',\n'GTX', 'X',\\
n'GTZ', 'X',\n'GU7', 'X',\n'GUA', 'X',\n'GUD', 'X'\
,\n'GUM', 'X',\n'GUN', 'X',\n'GUP', 'X',\n'GUR', '\
X',\n'GW3', 'X',\n'GZZ', 'X',\n'H2B', 'X',\n'H2P',\
 'H',\n'H2S', 'X',\n'H2U', 'X',\n'H4B', 'X',\n'H5M\
', 'P',\n'H5P', 'X',\n'HAA', 'X',\n'HAB', 'X',\n'H\
AC', 'A',\n'HAD', 'X',\n'HAE', 'X',\n'HAG', 'X',\n\
'HAI', 'X',\n'HAM', 'X',\n'HAP', 'X',\n'HAQ', 'X',\
\n'HAR', 'R',\n'HAS', 'X',\n'HAV', 'V',\n'HAX', 'X\
',\n'HAZ', 'X',\n'HBA', 'X',\n'HBC', 'X',\n'HBD', \
'X',\n'HBI', 'X',\n'HBO', 'X',\n'HBU', 'X',\n'HBY'\
, 'X',\n'HC0', 'X',\n'HC1', 'X',\n'HC4', 'X',\n'HC\
A', 'X',\n'HCC', 'X',\n'HCI', 'X',\n'HCS', 'X',\n'\
HDA', 'X',\n'HDD', 'X',\n'HDF', 'X',\n'HDN', 'X',\\
n'HDS', 'X',\n'HDZ', 'X',\n'HE1', 'X',\n'HE6', 'X'\
,\n'HEA', 'X',\n'HEB', 'X',\n'HEC', 'X',\n'HED', '\
X',\n'HEE', 'X',\n'HEF', 'X',\n'HEG', 'X',\n'HEM',\
 'X',\n'HEN', 'X',\n'HEO', 'X',\n'HEP', 'X',\n'HEU\
', 'X',\n'HEV', 'X',\n'HEX', 'X',\n'HEZ', 'X',\n'H\
F1', 'X',\n'HFA', 'X',\n'HFP', 'X',\n'HGA', 'Q',\n\
'HGB', 'X',\n'HGC', 'X',\n'HGI', 'X',\n'HGU', 'X',\
\n'HHO', 'X',\n'HHP', 'X',\n'HIB', 'X',\n'HIC', 'H\
',\n'HII', 'X',\n'HIN', 'X',\n'HIO', 'X',\n'HIP', \
'H',\n'HIS', 'H',\n'HLE', 'X',\n'HLT', 'X',\n'HMA'\
, 'A',\n'HMB', 'X',\n'HMC', 'X',\n'HMD', 'X',\n'HM\
F', 'A',\n'HMG', 'X',\n'HMH', 'X',\n'HMI', 'L',\n'\
HMM', 'X',\n'HMN', 'X',\n'HMO', 'X',\n'HMP', 'X',\\
n'HMR', 'R',\n'HNI', 'X',\n'HNP', 'X',\n'HOA', 'X'\
,\n'HOE', 'X',\n'HOH', 'X',\n'HOM', 'X',\n'HOP', '\
X',\n'HOQ', 'X',\n'HP1', 'A',\n'HP2', 'A',\n'HP3',\
 'X',\n'HPA', 'X',\n'HPB', 'X',\n'HPC', 'X',\n'HPD\
', 'X',\n'HPE', 'A',\n'HPG', 'X',\n'HPH', 'F',\n'H\
PP', 'X',\n'HPQ', 'F',\n'HPR', 'X',\n'HPT', 'X',\n\
'HPY', 'X',\n'HQO', 'X',\n'HQQ', 'X',\n'HQU', 'X',\
\n'HRG', 'R',\n'HRI', 'X',\n'HSA', 'X',\n'HSE', 'S\
',\n'HSF', 'X',\n'HSM', 'X',\n'HSO', 'H',\n'HSP', \
'X',\n'HT1', 'X',\n'HT2', 'X',\n'HTA', 'X',\n'HTL'\
, 'X',\n'HTO', 'X',\n'HTP', 'X',\n'HTR', 'W',\n'HU\
P', 'X',\n'HUX', 'X',\n'HV5', 'A',\n'HV7', 'X',\n'\
HV8', 'X',\n'HXA', 'X',\n'HXC', 'X',\n'HXP', 'X',\\
n'HY1', 'X',\n'HYA', 'X',\n'HYB', 'X',\n'HYD', 'X'\
,\n'HYG', 'X',\n'HYP', 'P',\n'I06', 'X',\n'I10', '\
X',\n'I11', 'X',\n'I17', 'X',\n'I2P', 'X',\n'I3N',\
 'X',\n'I3P', 'X',\n'I40', 'X',\n'I48', 'X',\n'I4B\
', 'X',\n'I52', 'X',\n'I5P', 'X',\n'I84', 'G',\n'I\
AG', 'G',\n'IAS', 'X',\n'IB2', 'X',\n'IBB', 'X',\n\
'IBP', 'X',\n'IBR', 'X',\n'IBS', 'X',\n'IBZ', 'X',\
\n'IC1', 'X',\n'ICA', 'X',\n'ICI', 'X',\n'ICL', 'X\
',\n'ICP', 'X',\n'ICT', 'X',\n'ICU', 'X',\n'ID2', \
'X',\n'IDC', 'X',\n'IDG', 'X',\n'IDH', 'X',\n'IDM'\
, 'X',\n'IDO', 'X',\n'IDP', 'X',\n'IDR', 'X',\n'ID\
S', 'X',\n'IDT', 'X',\n'IDU', 'X',\n'IFG', 'X',\n'\
IFP', 'X',\n'IGL', 'X',\n'IGN', 'X',\n'IGP', 'X',\\
n'IGU', 'X',\n'IH1', 'X',\n'IH2', 'X',\n'IH3', 'X'\
,\n'IHB', 'X',\n'IHN', 'X',\n'IHP', 'X',\n'IIC', '\
X',\n'IIL', 'I',\n'IIP', 'X',\n'IK2', 'X',\n'IKT',\
 'X',\n'ILA', 'I',\n'ILE', 'I',\n'ILG', 'X',\n'ILO\
', 'X',\n'ILX', 'I',\n'IM1', 'X',\n'IM2', 'X',\n'I\
MC', 'X',\n'IMD', 'X',\n'IME', 'X',\n'IMF', 'X',\n\
'IMG', 'X',\n'IMH', 'X',\n'IMI', 'X',\n'IML', 'I',\
\n'IMM', 'X',\n'IMN', 'X',\n'IMO', 'X',\n'IMP', 'X\
',\n'IMR', 'X',\n'IMU', 'X',\n'IN0', 'D',\n'IN1', \
'R',\n'IN2', 'K',\n'IN3', 'L',\n'IN4', 'X',\n'IN5'\
, 'A',\n'IN6', 'L',\n'IN7', 'X',\n'IN8', 'X',\n'IN\
9', 'X',\n'INA', 'L',\n'INB', 'X',\n'INC', 'X',\n'\
IND', 'X',\n'INE', 'X',\n'INF', 'F',\n'ING', 'F',\\
n'INH', 'R',\n'INI', 'X',\n'INJ', 'X',\n'INK', 'X'\
,\n'INL', 'X',\n'INM', 'X',\n'INN', 'A',\n'INO', '\
X',\n'INP', 'X',\n'INQ', 'X',\n'INR', 'X',\n'INS',\
 'X',\n'INT', 'V',\n'INU', 'X',\n'INV', 'X',\n'INW\
', 'X',\n'INX', 'X',\n'INY', 'X',\n'INZ', 'X',\n'I\
OA', 'X',\n'IOB', 'X',\n'IOC', 'X',\n'IOD', 'X',\n\
'IOE', 'X',\n'IOF', 'X',\n'IOH', 'X',\n'IOL', 'X',\
\n'IOP', 'X',\n'IP1', 'X',\n'IP2', 'X',\n'IP3', 'X\
',\n'IP4', 'X',\n'IPA', 'X',\n'IPB', 'X',\n'IPD', \
'X',\n'IPG', 'G',\n'IPH', 'X',\n'IPL', 'X',\n'IPM'\
, 'X',\n'IPN', 'X',\n'IPO', 'F',\n'IPP', 'X',\n'IP\
S', 'X',\n'IPT', 'X',\n'IPU', 'X',\n'IPY', 'A',\n'\
IQB', 'X',\n'IQP', 'X',\n'IQS', 'X',\n'IR3', 'X',\\
n'IRI', 'X',\n'IRP', 'X',\n'ISA', 'X',\n'ISF', 'X'\
,\n'ISO', 'X',\n'ISP', 'X',\n'ISQ', 'X',\n'ISU', '\
X',\n'ITM', 'X',\n'ITP', 'X',\n'ITR', 'W',\n'ITS',\
 'X',\n'ITU', 'X',\n'IU5', 'X',\n'IUM', 'X',\n'IUR\
', 'X',\n'IVA', 'X',\n'IYG', 'G',\n'IYR', 'Y',\n'J\
77', 'X',\n'J78', 'X',\n'J80', 'X',\n'JE2', 'X',\n\
'JEN', 'X',\n'JST', 'X',\n'K21', 'X',\n'KAH', 'X',\
\n'KAI', 'X',\n'KAM', 'X',\n'KAN', 'X',\n'KAP', 'X\
',\n'KCP', 'X',\n'KCX', 'K',\n'KDO', 'X',\n'KEF', \
'X',\n'KET', 'X',\n'KGR', 'X',\n'KH1', 'X',\n'KIF'\
, 'X',\n'KIV', 'V',\n'KNI', 'X',\n'KPH', 'K',\n'KT\
H', 'X',\n'KTN', 'X',\n'KTP', 'X',\n'KWT', 'X',\n'\
L04', 'X',\n'L1P', 'X',\n'L24', 'E',\n'L2P', 'X',\\
n'L34', 'E',\n'L37', 'E',\n'L3P', 'X',\n'L4P', 'X'\
,\n'L75', 'X',\n'LAC', 'X',\n'LAD', 'X',\n'LAK', '\
X',\n'LAM', 'X',\n'LAR', 'X',\n'LAT', 'X',\n'LAX',\
 'X',\n'LCO', 'X',\n'LCP', 'X',\n'LCS', 'X',\n'LDA\
', 'X',\n'LDO', 'L',\n'LDP', 'X',\n'LEA', 'X',\n'L\
EO', 'X',\n'LEU', 'L',\n'LG2', 'X',\n'LG6', 'X',\n\
'LGC', 'X',\n'LGP', 'X',\n'LHG', 'X',\n'LHY', 'F',\
\n'LI1', 'X',\n'LIG', 'X',\n'LIL', 'X',\n'LIM', 'X\
',\n'LIN', 'X',\n'LIO', 'X',\n'LIP', 'X',\n'LLA', \
'X',\n'LLP', 'K',\n'LLY', 'K',\n'LMG', 'X',\n'LML'\
, 'X',\n'LMT', 'X',\n'LMU', 'X',\n'LMZ', 'X',\n'LN\
K', 'X',\n'LNL', 'X',\n'LNO', 'X',\n'LOF', 'X',\n'\
LOL', 'L',\n'LOM', 'X',\n'LOR', 'X',\n'LOS', 'X',\\
n'LOV', 'L',\n'LOX', 'X',\n'LP1', 'X',\n'LP2', 'R'\
,\n'LPA', 'X',\n'LPC', 'X',\n'LPF', 'X',\n'LPL', '\
X',\n'LPM', 'X',\n'LPP', 'X',\n'LRB', 'X',\n'LRU',\
 'X',\n'LS1', 'X',\n'LS2', 'X',\n'LS3', 'X',\n'LS4\
', 'X',\n'LS5', 'X',\n'LTA', 'X',\n'LTL', 'X',\n'L\
TR', 'W',\n'LUM', 'X',\n'LVS', 'L',\n'LXC', 'X',\n\
'LY2', 'X',\n'LY3', 'X',\n'LYA', 'X',\n'LYB', 'X',\
\n'LYC', 'X',\n'LYD', 'X',\n'LYM', 'K',\n'LYN', 'X\
',\n'LYS', 'K',\n'LYT', 'X',\n'LYW', 'X',\n'LYZ', \
'K',\n'M1A', 'X',\n'M1G', 'X',\n'M2G', 'X',\n'M3L'\
, 'K',\n'M6P', 'X',\n'M6T', 'X',\n'M7G', 'X',\n'MA\
1', 'X',\n'MA2', 'X',\n'MA3', 'X',\n'MA4', 'X',\n'\
MA6', 'X',\n'MAA', 'A',\n'MAB', 'X',\n'MAC', 'X',\\
n'MAE', 'X',\n'MAG', 'X',\n'MAH', 'X',\n'MAI', 'R'\
,\n'MAK', 'X',\n'MAL', 'X',\n'MAM', 'X',\n'MAN', '\
X',\n'MAO', 'X',\n'MAP', 'X',\n'MAR', 'X',\n'MAS',\
 'X',\n'MAT', 'X',\n'MAU', 'X',\n'MAZ', 'X',\n'MBA\
', 'X',\n'MBD', 'X',\n'MBG', 'X',\n'MBH', 'X',\n'M\
BN', 'X',\n'MBO', 'X',\n'MBR', 'X',\n'MBS', 'X',\n\
'MBV', 'X',\n'MBZ', 'X',\n'MCA', 'X',\n'MCD', 'X',\
\n'MCE', 'X',\n'MCG', 'G',\n'MCI', 'X',\n'MCN', 'X\
',\n'MCP', 'X',\n'MCT', 'X',\n'MCY', 'X',\n'MD2', \
'X',\n'MDA', 'X',\n'MDC', 'X',\n'MDG', 'X',\n'MDH'\
, 'X',\n'MDL', 'X',\n'MDM', 'X',\n'MDN', 'X',\n'MD\
P', 'X',\n'ME6', 'X',\n'MEB', 'X',\n'MEC', 'X',\n'\
MEL', 'X',\n'MEN', 'N',\n'MEP', 'X',\n'MER', 'X',\\
n'MES', 'X',\n'MET', 'M',\n'MEV', 'X',\n'MF2', 'X'\
,\n'MF3', 'M',\n'MFB', 'X',\n'MFD', 'X',\n'MFU', '\
X',\n'MG7', 'X',\n'MGA', 'X',\n'MGB', 'X',\n'MGD',\
 'X',\n'MGG', 'R',\n'MGL', 'X',\n'MGN', 'Q',\n'MGO\
', 'X',\n'MGP', 'X',\n'MGR', 'X',\n'MGS', 'X',\n'M\
GT', 'X',\n'MGU', 'X',\n'MGY', 'G',\n'MHB', 'X',\n\
'MHF', 'X',\n'MHL', 'L',\n'MHM', 'X',\n'MHO', 'M',\
\n'MHS', 'H',\n'MHZ', 'X',\n'MIA', 'X',\n'MIC', 'X\
',\n'MID', 'X',\n'MIL', 'X',\n'MIM', 'X',\n'MIN', \
'G',\n'MIP', 'X',\n'MIS', 'S',\n'MIT', 'X',\n'MJI'\
, 'X',\n'MK1', 'X',\n'MKC', 'X',\n'MLA', 'X',\n'ML\
C', 'X',\n'MLE', 'L',\n'MLN', 'X',\n'MLT', 'X',\n'\
MLY', 'K',\n'MLZ', 'K',\n'MM3', 'X',\n'MM4', 'X',\\
n'MMA', 'X',\n'MMC', 'X',\n'MME', 'M',\n'MMO', 'R'\
,\n'MMP', 'X',\n'MMQ', 'X',\n'MMT', 'X',\n'MN1', '\
X',\n'MN2', 'X',\n'MN3', 'X',\n'MN5', 'X',\n'MN7',\
 'X',\n'MN8', 'X',\n'MNA', 'X',\n'MNB', 'X',\n'MNC\
', 'X',\n'MNG', 'X',\n'MNL', 'L',\n'MNO', 'X',\n'M\
NP', 'X',\n'MNQ', 'X',\n'MNS', 'X',\n'MNT', 'X',\n\
'MNV', 'V',\n'MO1', 'X',\n'MO2', 'X',\n'MO3', 'X',\
\n'MO4', 'X',\n'MO5', 'X',\n'MO6', 'X',\n'MOA', 'X\
',\n'MOB', 'X',\n'MOC', 'X',\n'MOE', 'X',\n'MOG', \
'X',\n'MOH', 'X',\n'MOL', 'X',\n'MOO', 'X',\n'MOP'\
, 'X',\n'MOR', 'X',\n'MOS', 'X',\n'MOT', 'X',\n'MO\
X', 'X',\n'MP1', 'X',\n'MP3', 'X',\n'MPA', 'X',\n'\
MPB', 'X',\n'MPC', 'X',\n'MPD', 'X',\n'MPG', 'X',\\
n'MPH', 'M',\n'MPI', 'X',\n'MPJ', 'M',\n'MPL', 'X'\
,\n'MPN', 'X',\n'MPO', 'X',\n'MPP', 'X',\n'MPQ', '\
G',\n'MPR', 'X',\n'MPS', 'X',\n'MQ0', 'X',\n'MQ7',\
 'X',\n'MQ8', 'X',\n'MQ9', 'X',\n'MQI', 'X',\n'MR2\
', 'X',\n'MRC', 'X',\n'MRM', 'X',\n'MRP', 'X',\n'M\
S2', 'X',\n'MSA', 'X',\n'MSB', 'X',\n'MSD', 'X',\n\
'MSE', 'M',\n'MSF', 'X',\n'MSI', 'X',\n'MSO', 'M',\
\n'MSQ', 'X',\n'MST', 'X',\n'MSU', 'X',\n'MTA', 'X\
',\n'MTB', 'X',\n'MTC', 'X',\n'MTD', 'X',\n'MTE', \
'X',\n'MTF', 'X',\n'MTG', 'X',\n'MTO', 'X',\n'MTS'\
, 'X',\n'MTT', 'X',\n'MTX', 'X',\n'MTY', 'Y',\n'MU\
G', 'X',\n'MUP', 'X',\n'MUR', 'X',\n'MVA', 'V',\n'\
MW1', 'X',\n'MW2', 'X',\n'MXA', 'X',\n'MXY', 'X',\\
n'MYA', 'X',\n'MYC', 'X',\n'MYG', 'X',\n'MYR', 'X'\
,\n'MYS', 'X',\n'MYT', 'X',\n'MZM', 'X',\n'N1T', '\
X',\n'N25', 'X',\n'N2B', 'X',\n'N3T', 'X',\n'N4B',\
 'X',\n'NA2', 'X',\n'NA5', 'X',\n'NA6', 'X',\n'NAA\
', 'X',\n'NAB', 'X',\n'NAC', 'X',\n'NAD', 'X',\n'N\
AE', 'X',\n'NAF', 'X',\n'NAG', 'X',\n'NAH', 'X',\n\
'NAI', 'X',\n'NAL', 'A',\n'NAM', 'A',\n'NAN', 'X',\
\n'NAO', 'X',\n'NAP', 'X',\n'NAQ', 'X',\n'NAR', 'X\
',\n'NAS', 'X',\n'NAU', 'X',\n'NAV', 'X',\n'NAW', \
'X',\n'NAX', 'X',\n'NAY', 'X',\n'NBA', 'X',\n'NBD'\
, 'X',\n'NBE', 'X',\n'NBG', 'X',\n'NBN', 'X',\n'NB\
P', 'X',\n'NBS', 'X',\n'NBU', 'X',\n'NCA', 'X',\n'\
NCB', 'A',\n'NCD', 'X',\n'NCH', 'X',\n'NCM', 'X',\\
n'NCN', 'X',\n'NCO', 'X',\n'NCR', 'X',\n'NCS', 'X'\
,\n'ND4', 'X',\n'NDA', 'X',\n'NDC', 'X',\n'NDD', '\
X',\n'NDO', 'X',\n'NDP', 'X',\n'NDT', 'X',\n'NEA',\
 'X',\n'NEB', 'X',\n'NED', 'X',\n'NEM', 'H',\n'NEN\
', 'X',\n'NEO', 'X',\n'NEP', 'H',\n'NEQ', 'X',\n'N\
ES', 'X',\n'NET', 'X',\n'NEV', 'X',\n'NFA', 'F',\n\
'NFE', 'X',\n'NFG', 'X',\n'NFP', 'X',\n'NFS', 'X',\
\n'NG6', 'X',\n'NGA', 'X',\n'NGL', 'X',\n'NGM', 'X\
',\n'NGO', 'X',\n'NGP', 'X',\n'NGT', 'X',\n'NGU', \
'X',\n'NH2', 'X',\n'NH3', 'X',\n'NH4', 'X',\n'NHD'\
, 'X',\n'NHE', 'X',\n'NHM', 'X',\n'NHP', 'X',\n'NH\
R', 'X',\n'NHS', 'X',\n'NI1', 'X',\n'NI2', 'X',\n'\
NIC', 'X',\n'NID', 'X',\n'NIK', 'X',\n'NIO', 'X',\\
n'NIP', 'X',\n'NIT', 'X',\n'NIU', 'X',\n'NIY', 'Y'\
,\n'NLA', 'X',\n'NLE', 'L',\n'NLG', 'X',\n'NLN', '\
L',\n'NLP', 'L',\n'NM1', 'X',\n'NMA', 'A',\n'NMB',\
 'X',\n'NMC', 'G',\n'NMD', 'X',\n'NME', 'X',\n'NMN\
', 'X',\n'NMO', 'X',\n'NMQ', 'X',\n'NMX', 'X',\n'N\
MY', 'X',\n'NNH', 'R',\n'NNO', 'X',\n'NO2', 'X',\n\
'NO3', 'X',\n'NOA', 'X',\n'NOD', 'X',\n'NOJ', 'X',\
\n'NON', 'X',\n'NOP', 'X',\n'NOR', 'X',\n'NOS', 'X\
',\n'NOV', 'X',\n'NOX', 'X',\n'NP3', 'X',\n'NPA', \
'X',\n'NPC', 'X',\n'NPD', 'X',\n'NPE', 'X',\n'NPF'\
, 'X',\n'NPH', 'C',\n'NPI', 'X',\n'NPL', 'X',\n'NP\
N', 'X',\n'NPO', 'X',\n'NPP', 'X',\n'NPT', 'X',\n'\
NPY', 'X',\n'NRG', 'R',\n'NRI', 'X',\n'NS1', 'X',\\
n'NS5', 'X',\n'NSP', 'X',\n'NTA', 'X',\n'NTB', 'X'\
,\n'NTC', 'X',\n'NTH', 'X',\n'NTM', 'X',\n'NTP', '\
X',\n'NTS', 'X',\n'NTU', 'X',\n'NTZ', 'X',\n'NU1',\
 'X',\n'NVA', 'V',\n'NVI', 'X',\n'NVP', 'X',\n'NW1\
', 'X',\n'NYP', 'X',\n'O4M', 'X',\n'OAA', 'X',\n'O\
AI', 'X',\n'OAP', 'X',\n'OAR', 'X',\n'OAS', 'S',\n\
'OBA', 'X',\n'OBN', 'X',\n'OC1', 'X',\n'OC2', 'X',\
\n'OC3', 'X',\n'OC4', 'X',\n'OC5', 'X',\n'OC6', 'X\
',\n'OC7', 'X',\n'OCL', 'X',\n'OCM', 'X',\n'OCN', \
'X',\n'OCO', 'X',\n'OCP', 'X',\n'OCS', 'C',\n'OCT'\
, 'X',\n'OCV', 'K',\n'OCY', 'C',\n'ODA', 'X',\n'OD\
S', 'X',\n'OES', 'X',\n'OET', 'X',\n'OF1', 'X',\n'\
OF2', 'X',\n'OF3', 'X',\n'OFL', 'X',\n'OFO', 'X',\\
n'OHE', 'X',\n'OHO', 'X',\n'OHT', 'X',\n'OIC', 'X'\
,\n'OIP', 'X',\n'OKA', 'X',\n'OLA', 'X',\n'OLE', '\
X',\n'OLI', 'X',\n'OLO', 'X',\n'OMB', 'X',\n'OMC',\
 'X',\n'OMD', 'X',\n'OME', 'X',\n'OMG', 'X',\n'OMP\
', 'X',\n'OMT', 'M',\n'OMU', 'X',\n'ONE', 'X',\n'O\
NL', 'L',\n'ONP', 'X',\n'OPA', 'X',\n'OPD', 'X',\n\
'OPE', 'X',\n'OPG', 'X',\n'OPH', 'X',\n'OPN', 'X',\
\n'OPP', 'X',\n'OPR', 'R',\n'ORN', 'X',\n'ORO', 'X\
',\n'ORP', 'X',\n'OSB', 'X',\n'OSS', 'X',\n'OTA', \
'X',\n'OTB', 'X',\n'OTE', 'X',\n'OTG', 'X',\n'OUT'\
, 'X',\n'OVA', 'X',\n'OWQ', 'X',\n'OXA', 'X',\n'OX\
E', 'X',\n'OXI', 'X',\n'OXL', 'X',\n'OXM', 'X',\n'\
OXN', 'X',\n'OXO', 'X',\n'OXP', 'X',\n'OXS', 'X',\\
n'OXY', 'X',\n'P11', 'A',\n'P24', 'X',\n'P28', 'X'\
,\n'P2P', 'X',\n'P2U', 'X',\n'P3M', 'X',\n'P4C', '\
X',\n'P4P', 'X',\n'P5P', 'X',\n'P6G', 'X',\n'PA1',\
 'X',\n'PA2', 'X',\n'PA3', 'X',\n'PA4', 'X',\n'PA5\
', 'X',\n'PAA', 'X',\n'PAB', 'X',\n'PAC', 'X',\n'P\
AD', 'X',\n'PAE', 'X',\n'PAG', 'X',\n'PAH', 'X',\n\
'PAI', 'X',\n'PAL', 'D',\n'PAM', 'X',\n'PAN', 'X',\
\n'PAO', 'X',\n'PAP', 'A',\n'PAQ', 'F',\n'PAR', 'X\
',\n'PAS', 'X',\n'PAT', 'W',\n'PBA', 'X',\n'PBB', \
'X',\n'PBC', 'X',\n'PBF', 'F',\n'PBG', 'X',\n'PBI'\
, 'X',\n'PBM', 'X',\n'PBN', 'X',\n'PBP', 'X',\n'PB\
R', 'X',\n'PBZ', 'X',\n'PC2', 'X',\n'PCA', 'E',\n'\
PCB', 'X',\n'PCD', 'X',\n'PCE', 'X',\n'PCG', 'X',\\
n'PCH', 'X',\n'PCL', 'X',\n'PCM', 'X',\n'PCP', 'X'\
,\n'PCR', 'X',\n'PCS', 'X',\n'PCU', 'X',\n'PCV', '\
X',\n'PCY', 'X',\n'PD1', 'X',\n'PDA', 'X',\n'PDC',\
 'X',\n'PDD', 'A',\n'PDE', 'A',\n'PDI', 'X',\n'PDL\
', 'A',\n'PDN', 'X',\n'PDO', 'X',\n'PDP', 'X',\n'P\
DT', 'X',\n'PDU', 'X',\n'PE2', 'X',\n'PE6', 'X',\n\
'PEA', 'X',\n'PEB', 'X',\n'PEC', 'X',\n'PED', 'X',\
\n'PEE', 'X',\n'PEF', 'X',\n'PEG', 'X',\n'PEL', 'X\
',\n'PEO', 'X',\n'PEP', 'X',\n'PEQ', 'X',\n'PER', \
'X',\n'PET', 'X',\n'PFB', 'X',\n'PFC', 'X',\n'PFG'\
, 'X',\n'PFL', 'X',\n'PFM', 'X',\n'PFZ', 'X',\n'PG\
4', 'X',\n'PG5', 'X',\n'PG6', 'X',\n'PGA', 'X',\n'\
PGC', 'X',\n'PGD', 'X',\n'PGE', 'X',\n'PGG', 'G',\\
n'PGH', 'X',\n'PGL', 'X',\n'PGO', 'X',\n'PGP', 'X'\
,\n'PGQ', 'X',\n'PGR', 'X',\n'PGS', 'X',\n'PGU', '\
X',\n'PGX', 'X',\n'PGY', 'G',\n'PH1', 'X',\n'PH2',\
 'X',\n'PH3', 'X',\n'PHA', 'F',\n'PHB', 'X',\n'PHC\
', 'X',\n'PHD', 'X',\n'PHE', 'F',\n'PHG', 'X',\n'P\
HH', 'X',\n'PHI', 'F',\n'PHL', 'F',\n'PHM', 'X',\n\
'PHN', 'X',\n'PHO', 'X',\n'PHP', 'X',\n'PHQ', 'X',\
\n'PHS', 'H',\n'PHT', 'X',\n'PHW', 'P',\n'PHY', 'X\
',\n'PI1', 'X',\n'PI2', 'X',\n'PI3', 'X',\n'PI4', \
'X',\n'PI5', 'X',\n'PI6', 'X',\n'PI7', 'X',\n'PI8'\
, 'X',\n'PI9', 'X',\n'PIA', 'X',\n'PIB', 'X',\n'PI\
C', 'X',\n'PID', 'X',\n'PIG', 'X',\n'PIH', 'X',\n'\
PIM', 'X',\n'PIN', 'X',\n'PIO', 'X',\n'PIP', 'X',\\
n'PIQ', 'X',\n'PIR', 'X',\n'PIV', 'X',\n'PKF', 'X'\
,\n'PL1', 'X',\n'PL9', 'X',\n'PLA', 'D',\n'PLC', '\
X',\n'PLE', 'L',\n'PLG', 'G',\n'PLH', 'X',\n'PLM',\
 'X',\n'PLP', 'X',\n'PLS', 'S',\n'PLT', 'W',\n'PLU\
', 'L',\n'PLY', 'X',\n'PMA', 'X',\n'PMB', 'X',\n'P\
MC', 'X',\n'PME', 'F',\n'PML', 'X',\n'PMM', 'X',\n\
'PMO', 'X',\n'PMP', 'X',\n'PMS', 'X',\n'PMY', 'X',\
\n'PN2', 'X',\n'PNA', 'X',\n'PNB', 'X',\n'PNC', 'G\
',\n'PND', 'X',\n'PNE', 'A',\n'PNF', 'X',\n'PNG', \
'X',\n'PNI', 'X',\n'PNL', 'X',\n'PNM', 'X',\n'PNN'\
, 'X',\n'PNO', 'X',\n'PNP', 'X',\n'PNQ', 'X',\n'PN\
S', 'X',\n'PNT', 'X',\n'PNU', 'X',\n'PO2', 'X',\n'\
PO4', 'X',\n'POB', 'X',\n'POC', 'X',\n'POL', 'X',\\
n'POM', 'P',\n'PON', 'X',\n'POP', 'X',\n'POR', 'X'\
,\n'POS', 'X',\n'PP1', 'X',\n'PP2', 'X',\n'PP3', '\
A',\n'PP4', 'X',\n'PP5', 'X',\n'PP6', 'X',\n'PP7',\
 'X',\n'PP8', 'N',\n'PP9', 'X',\n'PPB', 'X',\n'PPC\
', 'X',\n'PPD', 'X',\n'PPE', 'E',\n'PPG', 'X',\n'P\
PH', 'F',\n'PPI', 'X',\n'PPJ', 'V',\n'PPL', 'X',\n\
'PPM', 'X',\n'PPN', 'A',\n'PPO', 'X',\n'PPP', 'X',\
\n'PPQ', 'X',\n'PPR', 'X',\n'PPS', 'X',\n'PPT', 'X\
',\n'PPU', 'X',\n'PPX', 'F',\n'PPY', 'X',\n'PPZ', \
'X',\n'PQ0', 'X',\n'PQN', 'X',\n'PQQ', 'X',\n'PR1'\
, 'X',\n'PR2', 'X',\n'PR3', 'X',\n'PRA', 'X',\n'PR\
B', 'X',\n'PRC', 'X',\n'PRD', 'X',\n'PRE', 'X',\n'\
PRF', 'X',\n'PRH', 'X',\n'PRI', 'P',\n'PRL', 'X',\\
n'PRN', 'X',\n'PRO', 'P',\n'PRP', 'X',\n'PRR', 'A'\
,\n'PRS', 'P',\n'PRZ', 'X',\n'PS0', 'X',\n'PSA', '\
X',\n'PSD', 'X',\n'PSE', 'X',\n'PSF', 'S',\n'PSG',\
 'X',\n'PSI', 'X',\n'PSO', 'X',\n'PSQ', 'X',\n'PSS\
', 'X',\n'PST', 'X',\n'PSU', 'X',\n'PT1', 'X',\n'P\
T3', 'X',\n'PTA', 'X',\n'PTC', 'X',\n'PTD', 'X',\n\
'PTE', 'X',\n'PTH', 'Y',\n'PTL', 'X',\n'PTM', 'Y',\
\n'PTN', 'X',\n'PTO', 'X',\n'PTP', 'X',\n'PTR', 'Y\
',\n'PTS', 'X',\n'PTT', 'X',\n'PTU', 'X',\n'PTY', \
'X',\n'PUA', 'X',\n'PUB', 'X',\n'PUR', 'X',\n'PUT'\
, 'X',\n'PVA', 'X',\n'PVB', 'X',\n'PVH', 'H',\n'PV\
L', 'X',\n'PXA', 'X',\n'PXF', 'X',\n'PXG', 'X',\n'\
PXP', 'X',\n'PXY', 'X',\n'PXZ', 'X',\n'PY2', 'X',\\
n'PY4', 'X',\n'PY5', 'X',\n'PY6', 'X',\n'PYA', 'A'\
,\n'PYC', 'X',\n'PYD', 'X',\n'PYE', 'X',\n'PYL', '\
X',\n'PYM', 'X',\n'PYO', 'X',\n'PYP', 'X',\n'PYQ',\
 'X',\n'PYR', 'X',\n'PYS', 'X',\n'PYT', 'X',\n'PYX\
', 'X',\n'PYY', 'X',\n'PYZ', 'X',\n'PZQ', 'X',\n'Q\
82', 'X',\n'QNC', 'X',\n'QND', 'X',\n'QSI', 'Q',\n\
'QTR', 'X',\n'QUA', 'X',\n'QUE', 'X',\n'QUI', 'X',\
\n'QUO', 'X',\n'R11', 'X',\n'R12', 'X',\n'R13', 'X\
',\n'R18', 'X',\n'R1P', 'X',\n'R56', 'X',\n'R5P', \
'X',\n'RA2', 'X',\n'RAD', 'X',\n'RAI', 'X',\n'RAL'\
, 'X',\n'RAM', 'X',\n'RAN', 'X',\n'RAP', 'X',\n'RB\
F', 'X',\n'RBU', 'X',\n'RCA', 'X',\n'RCL', 'X',\n'\
RCO', 'X',\n'RDC', 'X',\n'RDF', 'W',\n'RE9', 'X',\\
n'REA', 'X',\n'RED', 'K',\n'REO', 'X',\n'REP', 'X'\
,\n'RET', 'X',\n'RFA', 'X',\n'RFB', 'X',\n'RFL', '\
X',\n'RFP', 'X',\n'RG1', 'X',\n'RGS', 'X',\n'RH1',\
 'X',\n'RHA', 'X',\n'RHC', 'X',\n'RHD', 'X',\n'RHM\
', 'X',\n'RHO', 'X',\n'RHQ', 'X',\n'RHS', 'X',\n'R\
IA', 'X',\n'RIB', 'X',\n'RIC', 'X',\n'RIF', 'X',\n\
'RIN', 'X',\n'RIP', 'X',\n'RIT', 'X',\n'RMB', 'X',\
\n'RMN', 'X',\n'RMP', 'X',\n'RNG', 'X',\n'RNS', 'X\
',\n'RNT', 'X',\n'RO2', 'X',\n'RO4', 'X',\n'ROC', \
'N',\n'ROI', 'X',\n'ROM', 'X',\n'RON', 'V',\n'ROP'\
, 'X',\n'ROS', 'X',\n'ROX', 'X',\n'RPA', 'X',\n'RP\
D', 'X',\n'RPH', 'X',\n'RPL', 'X',\n'RPP', 'X',\n'\
RPR', 'X',\n'RPX', 'X',\n'RQ3', 'X',\n'RR1', 'X',\\
n'RR6', 'X',\n'RRS', 'X',\n'RS1', 'X',\n'RS2', 'X'\
,\n'RS7', 'X',\n'RSS', 'X',\n'RTA', 'X',\n'RTB', '\
X',\n'RTC', 'X',\n'RTL', 'X',\n'RUB', 'X',\n'RUN',\
 'X',\n'RWJ', 'X',\n'RXP', 'X',\n'S02', 'X',\n'S11\
', 'X',\n'S1H', 'S',\n'S27', 'X',\n'S2C', 'C',\n'S\
3P', 'X',\n'S4U', 'X',\n'S57', 'X',\n'S58', 'X',\n\
'S5H', 'X',\n'S6G', 'X',\n'S80', 'X',\n'SAA', 'X',\
\n'SAB', 'X',\n'SAC', 'S',\n'SAD', 'X',\n'SAE', 'X\
',\n'SAF', 'X',\n'SAH', 'C',\n'SAI', 'C',\n'SAL', \
'X',\n'SAM', 'M',\n'SAN', 'X',\n'SAP', 'X',\n'SAR'\
, 'X',\n'SAS', 'X',\n'SB1', 'X',\n'SB2', 'X',\n'SB\
3', 'X',\n'SB4', 'X',\n'SB5', 'X',\n'SB6', 'X',\n'\
SBA', 'L',\n'SBB', 'X',\n'SBD', 'A',\n'SBI', 'X',\\
n'SBL', 'A',\n'SBN', 'X',\n'SBO', 'X',\n'SBR', 'X'\
,\n'SBS', 'X',\n'SBT', 'X',\n'SBU', 'X',\n'SBX', '\
X',\n'SC4', 'X',\n'SCA', 'X',\n'SCC', 'X',\n'SCD',\
 'X',\n'SCH', 'C',\n'SCI', 'X',\n'SCL', 'X',\n'SCM\
', 'X',\n'SCN', 'X',\n'SCO', 'X',\n'SCP', 'S',\n'S\
CR', 'X',\n'SCS', 'X',\n'SCV', 'C',\n'SCY', 'C',\n\
'SD8', 'X',\n'SDK', 'X',\n'SDZ', 'X',\n'SE4', 'X',\
\n'SEA', 'X',\n'SEB', 'S',\n'SEC', 'X',\n'SEG', 'A\
',\n'SEI', 'X',\n'SEL', 'S',\n'SEM', 'X',\n'SEO', \
'X',\n'SEP', 'S',\n'SER', 'S',\n'SES', 'X',\n'SET'\
, 'S',\n'SEU', 'X',\n'SF4', 'X',\n'SFG', 'X',\n'SF\
N', 'X',\n'SFO', 'X',\n'SGA', 'X',\n'SGC', 'X',\n'\
SGL', 'X',\n'SGM', 'X',\n'SGN', 'X',\n'SGP', 'X',\\
n'SHA', 'X',\n'SHC', 'X',\n'SHF', 'X',\n'SHH', 'X'\
,\n'SHP', 'G',\n'SHR', 'E',\n'SHT', 'T',\n'SHU', '\
X',\n'SI2', 'X',\n'SIA', 'X',\n'SIF', 'X',\n'SIG',\
 'X',\n'SIH', 'X',\n'SIM', 'X',\n'SIN', 'X',\n'SKD\
', 'X',\n'SKF', 'X',\n'SLB', 'X',\n'SLE', 'X',\n'S\
LZ', 'K',\n'SMA', 'X',\n'SMC', 'C',\n'SME', 'M',\n\
'SML', 'X',\n'SMM', 'M',\n'SMN', 'X',\n'SMP', 'X',\
\n'SMS', 'X',\n'SN1', 'X',\n'SN6', 'X',\n'SN7', 'X\
',\n'SNC', 'C',\n'SNN', 'X',\n'SNP', 'X',\n'SO1', \
'X',\n'SO2', 'X',\n'SO3', 'X',\n'SO4', 'X',\n'SOA'\
, 'X',\n'SOC', 'C',\n'SOM', 'X',\n'SOR', 'X',\n'SO\
T', 'X',\n'SOX', 'X',\n'SPA', 'X',\n'SPB', 'X',\n'\
SPC', 'X',\n'SPD', 'X',\n'SPE', 'X',\n'SPG', 'X',\\
n'SPH', 'X',\n'SPI', 'X',\n'SPK', 'X',\n'SPM', 'X'\
,\n'SPN', 'X',\n'SPO', 'X',\n'SPP', 'X',\n'SPS', '\
X',\n'SPY', 'X',\n'SQU', 'X',\n'SRA', 'X',\n'SRB',\
 'X',\n'SRD', 'X',\n'SRL', 'X',\n'SRM', 'X',\n'SRS\
', 'X',\n'SRY', 'X',\n'SSA', 'X',\n'SSB', 'X',\n'S\
SG', 'X',\n'SSP', 'X',\n'ST1', 'X',\n'ST2', 'X',\n\
'ST3', 'X',\n'ST4', 'X',\n'ST5', 'X',\n'ST6', 'X',\
\n'STA', 'X',\n'STB', 'X',\n'STE', 'X',\n'STG', 'X\
',\n'STI', 'X',\n'STL', 'X',\n'STN', 'X',\n'STO', \
'X',\n'STP', 'X',\n'STR', 'X',\n'STU', 'X',\n'STY'\
, 'Y',\n'SU1', 'X',\n'SU2', 'X',\n'SUC', 'X',\n'SU\
I', 'X',\n'SUL', 'X',\n'SUR', 'X',\n'SVA', 'S',\n'\
SWA', 'X',\n'T16', 'X',\n'T19', 'X',\n'T23', 'X',\\
n'T29', 'X',\n'T33', 'X',\n'T3P', 'X',\n'T42', 'A'\
,\n'T44', 'X',\n'T5A', 'X',\n'T6A', 'T',\n'T6P', '\
X',\n'T80', 'X',\n'T87', 'X',\n'TA1', 'X',\n'TAA',\
 'X',\n'TAB', 'X',\n'TAC', 'X',\n'TAD', 'X',\n'TAF\
', 'X',\n'TAM', 'X',\n'TAP', 'X',\n'TAR', 'X',\n'T\
AS', 'X',\n'TAU', 'X',\n'TAX', 'X',\n'TAZ', 'X',\n\
'TB9', 'X',\n'TBA', 'X',\n'TBD', 'X',\n'TBG', 'G',\
\n'TBH', 'X',\n'TBM', 'T',\n'TBO', 'X',\n'TBP', 'X\
',\n'TBR', 'X',\n'TBS', 'X',\n'TBT', 'X',\n'TBU', \
'X',\n'TBZ', 'X',\n'TC4', 'X',\n'TCA', 'X',\n'TCB'\
, 'X',\n'TCH', 'X',\n'TCK', 'X',\n'TCL', 'X',\n'TC\
M', 'X',\n'TCN', 'X',\n'TCP', 'X',\n'TCR', 'W',\n'\
TCS', 'X',\n'TCZ', 'X',\n'TDA', 'X',\n'TDB', 'X',\\
n'TDG', 'X',\n'TDP', 'X',\n'TDR', 'X',\n'TDX', 'X'\
,\n'TEA', 'X',\n'TEM', 'X',\n'TEN', 'X',\n'TEO', '\
X',\n'TEP', 'X',\n'TER', 'X',\n'TES', 'X',\n'TET',\
 'X',\n'TFA', 'X',\n'TFB', 'X',\n'TFH', 'X',\n'TFI\
', 'X',\n'TFK', 'X',\n'TFP', 'X',\n'THA', 'X',\n'T\
HB', 'X',\n'THC', 'T',\n'THD', 'X',\n'THE', 'X',\n\
'THF', 'X',\n'THJ', 'X',\n'THK', 'X',\n'THM', 'X',\
\n'THN', 'X',\n'THO', 'T',\n'THP', 'X',\n'THQ', 'X\
',\n'THR', 'T',\n'THS', 'X',\n'THT', 'X',\n'THU', \
'X',\n'THX', 'X',\n'THZ', 'X',\n'TI1', 'X',\n'TI2'\
, 'X',\n'TI3', 'P',\n'TIA', 'X',\n'TIH', 'A',\n'TK\
4', 'X',\n'TLA', 'X',\n'TLC', 'X',\n'TLM', 'X',\n'\
TLN', 'X',\n'TLX', 'X',\n'TM5', 'X',\n'TM6', 'X',\\
n'TMA', 'X',\n'TMB', 'T',\n'TMC', 'X',\n'TMD', 'T'\
,\n'TME', 'X',\n'TMF', 'X',\n'TML', 'K',\n'TMM', '\
X',\n'TMN', 'X',\n'TMP', 'X',\n'TMQ', 'X',\n'TMR',\
 'X',\n'TMT', 'X',\n'TMZ', 'X',\n'TNB', 'C',\n'TND\
', 'X',\n'TNK', 'X',\n'TNP', 'X',\n'TNT', 'X',\n'T\
OA', 'X',\n'TOB', 'X',\n'TOC', 'X',\n'TOL', 'X',\n\
'TOP', 'X',\n'TOS', 'X',\n'TOT', 'X',\n'TP1', 'G',\
\n'TP2', 'P',\n'TP3', 'E',\n'TP4', 'E',\n'TP7', 'T\
',\n'TPA', 'X',\n'TPE', 'X',\n'TPF', 'X',\n'TPI', \
'X',\n'TPL', 'W',\n'TPM', 'X',\n'TPN', 'G',\n'TPO'\
, 'T',\n'TPP', 'X',\n'TPQ', 'A',\n'TPR', 'P',\n'TP\
S', 'X',\n'TPT', 'X',\n'TPV', 'X',\n'TPX', 'X',\n'\
TPY', 'X',\n'TQ3', 'X',\n'TQ4', 'X',\n'TQ5', 'X',\\
n'TQ6', 'X',\n'TR1', 'X',\n'TRA', 'X',\n'TRB', 'X'\
,\n'TRC', 'X',\n'TRD', 'X',\n'TRE', 'X',\n'TRF', '\
W',\n'TRG', 'K',\n'TRH', 'X',\n'TRI', 'X',\n'TRJ',\
 'X',\n'TRM', 'X',\n'TRN', 'W',\n'TRO', 'W',\n'TRP\
', 'W',\n'TRQ', 'X',\n'TRS', 'X',\n'TRX', 'W',\n'T\
RZ', 'X',\n'TS2', 'X',\n'TS3', 'X',\n'TS4', 'X',\n\
'TS5', 'X',\n'TSA', 'X',\n'TSB', 'X',\n'TSI', 'X',\
\n'TSM', 'X',\n'TSN', 'X',\n'TSP', 'X',\n'TSU', 'X\
',\n'TTA', 'X',\n'TTE', 'X',\n'TTN', 'X',\n'TTO', \
'X',\n'TTP', 'X',\n'TTX', 'X',\n'TXL', 'X',\n'TYA'\
, 'Y',\n'TYB', 'Y',\n'TYD', 'X',\n'TYI', 'Y',\n'TY\
L', 'X',\n'TYM', 'W',\n'TYN', 'Y',\n'TYQ', 'Y',\n'\
TYR', 'Y',\n'TYS', 'Y',\n'TYV', 'X',\n'TYY', 'A',\\
n'TZB', 'X',\n'TZC', 'X',\n'TZE', 'X',\n'TZL', 'X'\
,\n'TZO', 'X',\n'TZP', 'X',\n'U01', 'X',\n'U02', '\
X',\n'U03', 'X',\n'U04', 'X',\n'U05', 'X',\n'U0E',\
 'X',\n'U10', 'X',\n'U18', 'X',\n'U2G', 'X',\n'U3P\
', 'X',\n'U49', 'X',\n'U55', 'X',\n'U5P', 'X',\n'U\
66', 'X',\n'U89', 'X',\n'U8U', 'X',\n'UAA', 'X',\n\
'UAG', 'A',\n'UAP', 'X',\n'UAR', 'X',\n'UC1', 'X',\
\n'UC2', 'X',\n'UC3', 'X',\n'UC4', 'X',\n'UD1', 'X\
',\n'UD2', 'X',\n'UDP', 'X',\n'UDX', 'X',\n'UFG', \
'X',\n'UFM', 'X',\n'UFP', 'X',\n'UGA', 'X',\n'UIN'\
, 'X',\n'UKP', 'A',\n'UM3', 'X',\n'UMA', 'A',\n'UM\
G', 'X',\n'UMP', 'X',\n'UNA', 'X',\n'UND', 'X',\n'\
UNI', 'X',\n'UNK', 'X',\n'UNN', 'X',\n'UNX', 'X',\\
n'UP5', 'X',\n'UP6', 'X',\n'UPA', 'X',\n'UPF', 'X'\
,\n'UPG', 'X',\n'UPP', 'X',\n'UQ1', 'X',\n'UQ2', '\
X',\n'UQ6', 'X',\n'UR2', 'X',\n'URA', 'X',\n'URE',\
 'X',\n'URF', 'X',\n'URI', 'X',\n'URS', 'X',\n'UTP\
', 'X',\n'UVC', 'X',\n'UVW', 'X',\n'V35', 'X',\n'V\
36', 'X',\n'V4O', 'X',\n'V7O', 'X',\n'VAA', 'V',\n\
'VAC', 'X',\n'VAD', 'V',\n'VAF', 'V',\n'VAG', 'X',\
\n'VAL', 'V',\n'VAN', 'X',\n'VAS', 'X',\n'VAX', 'X\
',\n'VDX', 'X',\n'VDY', 'X',\n'VG1', 'X',\n'VIB', \
'X',\n'VIR', 'X',\n'VIT', 'X',\n'VK3', 'X',\n'VO3'\
, 'X',\n'VO4', 'X',\n'VS1', 'F',\n'VS2', 'F',\n'VS\
3', 'F',\n'VS4', 'F',\n'VXA', 'X',\n'W01', 'X',\n'\
W02', 'X',\n'W03', 'X',\n'W11', 'X',\n'W33', 'X',\\
n'W35', 'X',\n'W42', 'X',\n'W43', 'X',\n'W54', 'X'\
,\n'W56', 'X',\n'W59', 'X',\n'W71', 'X',\n'W84', '\
X',\n'W8R', 'X',\n'W91', 'X',\n'WAY', 'X',\n'WCC',\
 'X',\n'WO2', 'X',\n'WO4', 'X',\n'WRB', 'X',\n'WRR\
', 'X',\n'WRS', 'X',\n'WW7', 'X',\n'X2F', 'X',\n'X\
7O', 'X',\n'XAA', 'X',\n'XAN', 'X',\n'XAO', 'X',\n\
'XBB', 'X',\n'XBP', 'X',\n'XDN', 'X',\n'XDP', 'X',\
\n'XIF', 'X',\n'XIM', 'X',\n'XK2', 'X',\n'XL1', 'X\
',\n'XLS', 'X',\n'XMP', 'X',\n'XN1', 'X',\n'XN2', \
'X',\n'XN3', 'X',\n'XUL', 'X',\n'XV6', 'X',\n'XYD'\
, 'X',\n'XYH', 'X',\n'XYL', 'X',\n'XYP', 'X',\n'XY\
S', 'X',\n'YOF', 'Y',\n'YRR', 'X',\n'YT3', 'X',\n'\
YZ9', 'X',\n'Z34', 'G',\n'Z5A', 'X',\n'ZAF', 'X',\\
n'ZAP', 'X',\n'ZEB', 'X',\n'ZEN', 'X',\n'ZES', 'X'\
,\n'ZID', 'X',\n'ZMR', 'X',\n'ZN3', 'X',\n'ZNH', '\
X',\n'ZNO', 'X',\n'ZO3', 'X',\n'ZPR', 'P',\n'ZRA',\
 'A',\n'ZST', 'X',\n'ZYA', 'A',\n\n\n'ASN','N');\n\
} \n\n\nsub file2head\n      {\n	my $file = shift;\
\n	my $size = shift;\n	my $f= new FileHandle;\n	my\
 $line;\n	open ($f,$file);\n	read ($f,$line, $size\
);\n	close ($f);\n	return $line;\n      }\nsub fil\
e2tail\n      {\n	my $file = shift;\n	my $size = s\
hift;\n	my $f= new FileHandle;\n	my $line;\n	\n	op\
en ($f,$file);\n	seek ($f,$size*-1, 2);\n	read ($f\
,$line, $size);\n	close ($f);\n	return $line;\n   \
   }\n\n\nsub vtmpnam\n      {\n	my $r=rand(100000\
);\n	my $f=\"file.$r.$$\";\n	while (-e $f)\n	  {\n\
	    $f=vtmpnam();\n	  }\n	push (@TMPFILE_LIST, $f\
);\n	return $f;\n      }\n\nsub myexit\n  {\n    m\
y $code=@_[0];\n    if ($CLEAN_EXIT_STARTED==1){re\
turn;}\n    else {$CLEAN_EXIT_STARTED=1;}\n    ###\
 ONLY BARE EXIT\n    exit ($code);\n  }\nsub set_e\
rror_lock\n    {\n      my $name = shift;\n      m\
y $pid=$$;\n\n      \n      &lock4tc ($$,\"LERROR\\
", \"LSET\", \"$$ -- ERROR: $name $PROGRAM\\n\");\\
n      return;\n    }\nsub set_lock\n  {\n    my $\
pid=shift;\n    my $msg= shift;\n    my $p=getppid\
();\n    &lock4tc ($pid,\"LLOCK\",\"LRESET\",\"$p$\
msg\\n\");\n  }\nsub unset_lock\n   {\n     \n    \
my $pid=shift;\n    &lock4tc ($pid,\"LLOCK\",\"LRE\
LEASE\",\"\");\n  }\nsub shift_lock\n  {\n    my $\
from=shift;\n    my $to=shift;\n    my $from_type=\
shift;\n    my $to_type=shift;\n    my $action=shi\
ft;\n    my $msg;\n    \n    if (!&lock4tc($from, \
$from_type, \"LCHECK\", \"\")){return 0;}\n    $ms\
g=&lock4tc ($from, $from_type, \"LREAD\", \"\");\n\
    &lock4tc ($from, $from_type,\"LRELEASE\", $msg\
);\n    &lock4tc ($to, $to_type, $action, $msg);\n\
    return;\n  }\nsub isshellpid\n  {\n    my $p=s\
hift;\n    if (!lock4tc ($p, \"LLOCK\", \"LCHECK\"\
)){return 0;}\n    else\n      {\n	my $c=lock4tc($\
p, \"LLOCK\", \"LREAD\");\n	if ( $c=~/-SHELL-/){re\
turn 1;}\n      }\n    return 0;\n  }\nsub isrootp\
id\n  {\n    if(lock4tc (getppid(), \"LLOCK\", \"L\
CHECK\")){return 0;}\n    else {return 1;}\n  }\ns\
ub lock4tc\n	{\n	  my ($pid,$type,$action,$value)=\
@_;\n	  my $fname;\n	  my $host=hostname;\n	  \n	 \
 if ($type eq \"LLOCK\"){$fname=\"$LOCKDIR/.$pid.$\
host.lock4tcoffee\";}\n	  elsif ( $type eq \"LERRO\
R\"){ $fname=\"$LOCKDIR/.$pid.$host.error4tcoffee\\
";}\n	  elsif ( $type eq \"LWARNING\"){ $fname=\"$\
LOCKDIR/.$pid.$host.warning4tcoffee\";}\n	  \n	  i\
f ($debug_lock)\n	    {\n	      print STDERR \"\\n\
\\t---lock4tc(tcg): $action => $fname =>$value (RD\
: $LOCKDIR)\\n\";\n	    }\n\n	  if    ($action eq \
\"LCHECK\") {return -e $fname;}\n	  elsif ($action\
 eq \"LREAD\"){return file2string($fname);}\n	  el\
sif ($action eq \"LSET\") {return string2file ($va\
lue, $fname, \">>\");}\n	  elsif ($action eq \"LRE\
SET\") {return string2file ($value, $fname, \">\")\
;}\n	  elsif ($action eq \"LRELEASE\") \n	    {\n	\
      if ( $debug_lock)\n		{\n		  my $g=new FileHa\
ndle;\n		  open ($g, \">>$fname\");\n		  print $g \
\"\\nDestroyed by $$\\n\";\n		  close ($g);\n		  s\
afe_system (\"mv $fname $fname.old\");\n		}\n	    \
  else\n		{\n		  unlink ($fname);\n		}\n	    }\n	 \
 return \"\";\n	}\n	\nsub file2string\n	{\n	  my $\
file=@_[0];\n	  my $f=new FileHandle;\n	  my $r;\n\
	  open ($f, \"$file\");\n	  while (<$f>){$r.=$_;}\
\n	  close ($f);\n	  return $r;\n	}\nsub string2fi\
le \n    {\n    my ($s,$file,$mode)=@_;\n    my $f\
=new FileHandle;\n    \n    open ($f, \"$mode$file\
\");\n    print $f  \"$s\";\n    close ($f);\n  }\\
n\nBEGIN\n    {\n      srand;\n    \n      $SIG{'S\
IGUP'}='signal_cleanup';\n      $SIG{'SIGINT'}='si\
gnal_cleanup';\n      $SIG{'SIGQUIT'}='signal_clea\
nup';\n      $SIG{'SIGILL'}='signal_cleanup';\n   \
   $SIG{'SIGTRAP'}='signal_cleanup';\n      $SIG{'\
SIGABRT'}='signal_cleanup';\n      $SIG{'SIGEMT'}=\
'signal_cleanup';\n      $SIG{'SIGFPE'}='signal_cl\
eanup';\n      \n      $SIG{'SIGKILL'}='signal_cle\
anup';\n      $SIG{'SIGPIPE'}='signal_cleanup';\n \
     $SIG{'SIGSTOP'}='signal_cleanup';\n      $SIG\
{'SIGTTIN'}='signal_cleanup';\n      $SIG{'SIGXFSZ\
'}='signal_cleanup';\n      $SIG{'SIGINFO'}='signa\
l_cleanup';\n      \n      $SIG{'SIGBUS'}='signal_\
cleanup';\n      $SIG{'SIGALRM'}='signal_cleanup';\
\n      $SIG{'SIGTSTP'}='signal_cleanup';\n      $\
SIG{'SIGTTOU'}='signal_cleanup';\n      $SIG{'SIGV\
TALRM'}='signal_cleanup';\n      $SIG{'SIGUSR1'}='\
signal_cleanup';\n\n\n      $SIG{'SIGSEGV'}='signa\
l_cleanup';\n      $SIG{'SIGTERM'}='signal_cleanup\
';\n      $SIG{'SIGCONT'}='signal_cleanup';\n     \
 $SIG{'SIGIO'}='signal_cleanup';\n      $SIG{'SIGP\
ROF'}='signal_cleanup';\n      $SIG{'SIGUSR2'}='si\
gnal_cleanup';\n\n      $SIG{'SIGSYS'}='signal_cle\
anup';\n      $SIG{'SIGURG'}='signal_cleanup';\n  \
    $SIG{'SIGCHLD'}='signal_cleanup';\n      $SIG{\
'SIGXCPU'}='signal_cleanup';\n      $SIG{'SIGWINCH\
'}='signal_cleanup';\n      \n      $SIG{'INT'}='s\
ignal_cleanup';\n      $SIG{'TERM'}='signal_cleanu\
p';\n      $SIG{'KILL'}='signal_cleanup';\n      $\
SIG{'QUIT'}='signal_cleanup';\n      \n      our $\
debug_lock=$ENV{\"DEBUG_LOCK\"};\n      \n      \n\
      \n      \n      foreach my $a (@ARGV){$CL.=\\
" $a\";}\n      if ( $debug_lock ){print STDERR \"\
\\n\\n\\n********** START PG: $PROGRAM ***********\
**\\n\";}\n      if ( $debug_lock ){print STDERR \\
"\\n\\n\\n**********(tcg) LOCKDIR: $LOCKDIR $$ ***\
**********\\n\";}\n      if ( $debug_lock ){print \
STDERR \"\\n --- $$ -- $CL\\n\";}\n      \n	     \\
n      \n      \n    }\nsub flush_error\n  {\n    \
my $msg=shift;\n    return add_error ($EXIT_FAILUR\
E,$$, $$,getppid(), $msg, $CL);\n  }\nsub add_erro\
r \n  {\n    my $code=shift;\n    my $rpid=shift;\\
n    my $pid=shift;\n    my $ppid=shift;\n    my $\
type=shift;\n    my $com=shift;\n    \n    $ERROR_\
DONE=1;\n    lock4tc ($rpid, \"LERROR\",\"LSET\",\\
"$pid -- ERROR: $type\\n\");\n    lock4tc ($$, \"L\
ERROR\",\"LSET\", \"$pid -- COM: $com\\n\");\n    \
lock4tc ($$, \"LERROR\",\"LSET\", \"$pid -- STACK:\
 $ppid -> $pid\\n\");\n   \n    return $code;\n  }\
\nsub add_warning \n  {\n    my $rpid=shift;\n    \
my $pid =shift;\n    my $command=shift;\n    my $m\
sg=\"$$ -- WARNING: $command\\n\";\n    print STDE\
RR \"$msg\";\n    lock4tc ($$, \"LWARNING\", \"LSE\
T\", $msg);\n  }\n\nsub signal_cleanup\n  {\n    p\
rint dtderr \"\\n**** $$ (tcg) was killed\\n\";\n \
   &cleanup;\n    exit ($EXIT_FAILURE);\n  }\nsub \
clean_dir\n  {\n    my $dir=@_[0];\n    if ( !-d $\
dir){return ;}\n    elsif (!($dir=~/tmp/)){return \
;}#safety check 1\n    elsif (($dir=~/\\*/)){retur\
n ;}#safety check 2\n    else\n      {\n	`rm -rf $\
dir`;\n      }\n    return;\n  }\nsub cleanup\n  {\
\n    #print stderr \"\\n----tc: $$ Kills $PIDCHIL\
D\\n\";\n    #kill (SIGTERM,$PIDCHILD);\n    my $p\
=getppid();\n    $CLEAN_EXIT_STARTED=1;\n    \n   \
 \n    \n    if (&lock4tc($$,\"LERROR\", \"LCHECK\\
", \"\"))\n      {\n	my $ppid=getppid();\n	if (!$E\
RROR_DONE) \n	  {\n	    &lock4tc($$,\"LERROR\", \"\
LSET\", \"$$ -- STACK: $p -> $$\\n\");\n	    &lock\
4tc($$,\"LERROR\", \"LSET\", \"$$ -- COM: $CL\\n\"\
);\n	  }\n      }\n    my $warning=&lock4tc($$, \"\
LWARNING\", \"LREAD\", \"\");\n    my $error=&lock\
4tc($$,  \"LERROR\", \"LREAD\", \"\");\n    #relea\
se error and warning lock if root\n    \n    if (i\
srootpid() && ($warning || $error) )\n      {\n	\n\
	print STDERR \"**************** Summary *********\
****\\n$error\\n$warning\\n\";\n\n	&lock4tc($$,\"L\
ERROR\",\"RELEASE\",\"\");\n	&lock4tc($$,\"LWARNIN\
G\",\"RELEASE\",\"\");\n      } \n    \n    \n    \
foreach my $f (@TMPFILE_LIST)\n      {\n	if (-e $f\
){unlink ($f);} \n      }\n    foreach my $d (@TMP\
DIR_LIST)\n      {\n	clean_dir ($d);\n      }\n   \
 #No More Lock Release\n    #&lock4tc($$,\"LLOCK\"\
,\"LRELEASE\",\"\"); #release lock \n\n    if ( $d\
ebug_lock ){print STDERR \"\\n\\n\\n********** END\
 PG: $PROGRAM ($$) *************\\n\";}\n    if ( \
$debug_lock ){print STDERR \"\\n\\n\\n**********(t\
cg) LOCKDIR: $LOCKDIR $$ *************\\n\";}\n  }\
\nEND \n  {\n    \n    &cleanup();\n  }\n   \n\nsu\
b safe_system \n{\n  my $com=shift;\n  my $ntry=sh\
ift;\n  my $ctry=shift;\n  my $pid;\n  my $status;\
\n  my $ppid=getppid();\n  if ($com eq \"\"){retur\
n 1;}\n  \n  \n\n  if (($pid = fork ()) < 0){retur\
n (-1);}\n  if ($pid == 0)\n    {\n      set_lock(\
$$, \" -SHELL- $com (tcg)\");\n      exec ($com);\\
n    }\n  else\n    {\n      lock4tc ($$, \"LLOCK\\
", \"LSET\", \"$pid\\n\");#update parent\n      $P\
IDCHILD=$pid;\n    }\n  if ($debug_lock){printf ST\
DERR \"\\n\\t .... safe_system (fasta_seq2hmm)  p:\
 $$ c: $pid COM: $com\\n\";}\n\n  waitpid ($pid,WT\
ERMSIG);\n\n  shift_lock ($pid,$$, \"LWARNING\",\"\
LWARNING\", \"LSET\");\n\n  if ($? == $EXIT_FAILUR\
E || lock4tc($pid, \"LERROR\", \"LCHECK\", \"\"))\\
n    {\n      if ($ntry && $ctry <$ntry)\n	{\n	  a\
dd_warning ($$,$$,\"$com failed [retry: $ctry]\");\
\n	  lock4tc ($pid, \"LRELEASE\", \"LERROR\", \"\"\
);\n	  return safe_system ($com, $ntry, ++$ctry);\\
n	}\n      elsif ($ntry == -1)\n	{\n	  if (!shift_\
lock ($pid, $$, \"LERROR\", \"LWARNING\", \"LSET\"\
))\n	    {\n	      add_warning ($$,$$,\"$com faile\
d\");\n	    }\n	  else\n	    {\n	      lock4tc ($p\
id, \"LRELEASE\", \"LERROR\", \"\");\n	    }\n	  r\
eturn $?;}\n      else\n	{\n	  if (!shift_lock ($p\
id,$$, \"LERROR\",\"LERROR\", \"LSET\"))\n	    {\n\
	      myexit(add_error ($EXIT_FAILURE,$$,$pid,get\
ppid(), \"UNSPECIFIED system\", $com));\n	    }\n	\
}\n    }\n  return $?;\n}\n\nsub check_configurati\
on \n    {\n      my @l=@_;\n      my $v;\n      f\
oreach my $p (@l)\n	{\n	  \n	  if   ( $p eq \"EMAI\
L\")\n	    { \n	      if ( !($EMAIL=~/@/))\n		{\n	\
	add_warning($$,$$,\"Could Not Use EMAIL\");\n		my\
exit(add_error ($EXIT_FAILURE,$$,$$,getppid(),\"EM\
AIL\",\"$CL\"));\n	      }\n	    }\n	  elsif( $p e\
q \"INTERNET\")\n	    {\n	      if ( !&check_inter\
net_connection())\n		{\n		  myexit(add_error ($EXI\
T_FAILURE,$$,$$,getppid(),\"INTERNET\",\"$CL\"));\\
n		}\n	    }\n	  elsif( $p eq \"wget\")\n	    {\n	\
      if (!&pg_is_installed (\"wget\") && !&pg_is_\
installed (\"curl\"))\n		{\n		  myexit(add_error (\
$EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:w\
get\",\"$CL\"));\n		}\n	    }\n	  elsif( !(&pg_is_\
installed ($p)))\n	    {\n	      myexit(add_error \
($EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:\
$p\",\"$CL\"));\n	    }\n	}\n      return 1;\n    \
}\nsub pg_is_installed\n  {\n    my @ml=@_;\n    m\
y $r, $p, $m;\n    my $supported=0;\n    \n    my \
$p=shift (@ml);\n    if ($p=~/::/)\n      {\n	if (\
safe_system (\"perl -M$p -e 1\")==$EXIT_SUCCESS){r\
eturn 1;}\n	else {return 0;}\n      }\n    else\n \
     {\n	$r=`which $p 2>/dev/null`;\n	if ($r eq \"\
\"){return 0;}\n	else {return 1;}\n      }\n  }\n\\
n\n\nsub check_internet_connection\n  {\n    my $i\
nternet;\n    my $tmp;\n    &check_configuration (\
 \"wget\"); \n    \n    $tmp=&vtmpnam ();\n    \n \
   if     (&pg_is_installed    (\"wget\")){`wget w\
ww.google.com -O$tmp >/dev/null 2>/dev/null`;}\n  \
  elsif  (&pg_is_installed    (\"curl\")){`curl ww\
w.google.com -o$tmp >/dev/null 2>/dev/null`;}\n   \
 \n    if ( !-e $tmp || -s $tmp < 10){$internet=0;\
}\n    else {$internet=1;}\n    if (-e $tmp){unlin\
k $tmp;}\n\n    return $internet;\n  }\nsub check_\
pg_is_installed\n  {\n    my @ml=@_;\n    my $r=&p\
g_is_installed (@ml);\n    if (!$r && $p=~/::/)\n \
     {\n	print STDERR \"\\nYou Must Install the pe\
rl package $p on your system.\\nRUN:\\n\\tsudo per\
l -MCPAN -e 'install $pg'\\n\";\n      }\n    elsi\
f (!$r)\n      {\n	myexit(flush_error(\"\\nProgram\
 $p Supported but Not Installed on your system\"))\
;\n      }\n    else\n      {\n	return 1;\n      }\
\n  }\n\n\nsub remote_is_pdb_name_deprecated\n{\n \
   my $in=@_[0];\n    my ($ref_file, $pdb);\n    m\
y ($value,$value1,$value2);\n    my $max=2;\n    \\
n    \n    \n    $ref_file=\"$cache/pdb_entry_type\
.txt\";\n    \n    if ( $in=~/[^\\w\\d\\:\\_]/){re\
turn 0;}\n    elsif ($no_remote_pdb_dir==1)\n     \
 {\n	my $pdbdir=$ENV{'PDB_DIR'};\n	\n	my $r1=\"$pd\
bdir/derived_data/pdb_entry_type.txt\";\n	my $r2=$\
ref_file;\n	if    (-e $r1){$ref_file=$r1;}\n	elsif\
 (-e $r2){$ref_file=$r2;}\n	else\n	  {\n	    my $p\
=substr ($in,0, 4);\n	    add_warning ($$, $$, \"C\
annot find pdb_entry_type.txt;  $p is assumed to b\
e valid; add ftp://ftp.wwpdb.org/pub/pdb/derived_d\
ata/pdb_entry_type.txt in $cache to check name sta\
tus\");\n	  }\n      }\n    elsif ( !-e $ref_file \
|| (-M $ref_file)>$max || -z $ref_file)\n      {\n\
	&url2file(\"ftp://ftp.wwpdb.org/pub/pdb/derived_d\
ata/pdb_entry_type.txt\", $ref_file);\n      }\n  \
  $pdb=substr ($in,0, 4);\n    chomp(($value1=`gre\
p -c $pdb $ref_file`));\n    $pdb=lc($pdb);\n    c\
homp(($value2=`grep -c $pdb $ref_file`));\n    $va\
lue=($value1 || $value2)?1:0;\n    $value=($value>\
0)?1:0;\n    \n    return $value;\n  }\n","use Cwd\
;\nuse Env;\nuse File::Path;\nuse FileHandle;\nuse\
 strict;\n\n\nour (%MODE, %PG, %ENV_SET, %SUPPORTE\
D_OS);\n\n\nour $EXIT_SUCCESS=0;\nour $EXIT_FAILUR\
E=1;\nour $INTERNET=0;\n\nour $CP=\"cp \"; #was ca\
using a crash on MacOSX\nour $SILENT=\">/dev/null \
2>/dev/null\";\nour $WEB_BASE=\"http://www.tcoffee\
.org\";\nour $TCLINKDB_ADDRESS=\"$WEB_BASE/Resourc\
es/tclinkdb.txt\";\nour $OS=get_os();\nour $ROOT=&\
get_root();\nour $CD=cwd();\nour $CDIR=$CD;\nour $\
HOME=$ENV{'HOME'};\n\nour $OSNAME=$ENV{'OSNAME'};\\
nour $OSARCH=$ENV{'OSARCH'};\nour $REPO_ROOT=\"\";\
\n\nour $TCDIR;\nour $TCCACHE;\nour $TCTMP;\nour $\
TCM;\nour $TCMETHODS;\nour $TCPLUGINS;\nour $PLUGI\
NS_DIR=\"\";\nour $INSTALL_DIR=\"\";\n\nour $CXX=\\
"g++\";\nour $CXXFLAGS=\"\";\n\nour $CPP=\"g++\";\\
nour $CPPFLAGS=\"\";\n\nour $CC=\"gcc\";\nour $CFL\
AGS=$ENV{'CFLAGS'};\n\nour $FC=\"f77\";\nour $FFLA\
GS=\"\";\n\nmy $install=\"all\";\nmy $default_upda\
te_action=\"no_update\";\nmy @required_application\
s=(\"wget_OR_curl\");\nmy @smode=(\"all\", \"clean\
\", \"install\");\n\n&initialize_PG();\n\nmy $cl=j\
oin( \" \", @ARGV);\nif ($#ARGV==-1 || ($cl=~/-h/)\
 ||($cl=~/-H/) )\n  {\n     print \"\\n!!!!!!! ./i\
nstall  t_coffee             --> installs t_coffee\
 only\";\n     print \"\\n!!!!!!! ./install  all  \
                --> installs all the modes [mcoffe\
e, expresso, psicoffee,rcoffee..]\";\n     print \\
"\\n!!!!!!! ./install  [mcoffee|rcoffee|..] --> in\
stalls the specified mode\";\n     print \"\\n!!!!\
!!! ./install  -h                   --> print usag\
e\\n\\n\";\n     if ( $#ARGV==-1){exit ($EXIT_FAIL\
URE);}\n   }\n     \nif (($cl=~/-h/) ||($cl=~/-H/)\
 )\n  {\n    my $m;\n    print \"\\n\\n!!!!!!! adv\
anced mode\\n\";\n    foreach $m ((keys (%MODE)),@\
smode)\n      {\n	print \"!!!!!!!       ./install \
$m\\n\";\n      }\n    \n    print \"!!!!!!! ./ins\
tall [target:package|mode|] [-update|-force|-exec=\
dir|-dis=dir|-root|-tclinkdb=file|-] [CC=|FCC=|CXX\
=|CFLAGS=|CXXFLAGS=]\\n\";\n    print \"!!!!!!! ./\
install clean    [removes all executables]\\n\";\n\
    print \"!!!!!!! ./install [optional:target] -u\
pdate               [updates package already insta\
lled]\\n\";\n    print \"!!!!!!! ./install [option\
al:target] -force                [Forces recompila\
tion over everything]\\n\";\n    \n    print \"!!!\
!!!! ./install [optional:target] -root            \
     [You are running as root]\\n\";\n    print \"\
!!!!!!! ./install [optional:target] -exec=/foo/bar\
/       [address for the T-Coffee executable]\\n\"\
;\n    print \"!!!!!!! ./install [optional:target]\
 -dis=/foo/bar/        [Address where distribution\
s should be stored]\\n\";\n    print \"!!!!!!! ./i\
nstall [optional:target] -tclinkdb=foo|update  [fi\
le containing all the packages to be installed]\\n\
\";\n    print \"!!!!!!! ./install [optional:targe\
t] -clean                [clean everything]\\n\";\\
n    print \"!!!!!!! ./install [optional:target] -\
plugins              [plugins directory]\\n\";\n  \
  print \"!!!!!!! ./install [optional:target] -tcd\
ir=/foor/bar      [base path where T-Coffee will b\
e installed]\\n\";\n    print \"!!!!!!! ./install \
[optional:target] -repo=/path/to/repo   [binaries \
repository root directory]\\n\";\n    print \"!!!!\
!!! mode:\";\n    foreach $m (keys(%MODE)){print \\
"$m \";}\n    print \"\\n\";\n    print \"!!!!!!! \
Packages:\";\n    foreach $m (keys (%PG)){print \"\
$m \";}\n    print \"\\n\";\n    \n    print \"\\n\
\\n\";\n    exit ($EXIT_FAILURE);\n  }\n\n\n\nmy (\
@argl)=($cl=~/(\\S+=[^=]+)\\s\\w+=/g);\npush (@arg\
l, ($cl=~/(\\S+=[^=]+\\S)\\s*$/g));\n\nforeach $a \
(@argl)\n  {\n    if ( ($cl=~/CXX=(.*)/)){$CXX=$1;\
}\n    if ( ($cl=~/-CC=(.*)/    )){$CC=$1;}\n    i\
f ( ($cl=~/-FC=(.*)/    )){$FC=$1;}\n    if ( ($cl\
=~/-CFLAGS=(.*)/)){$CFLAGS=$1;}\n    if ( ($cl=~/-\
CXXFLAGS=(.*)/)){$CXXFLAGS=$1;}\n  }\nour ($ROOT_I\
NSTALL, $NO_QUESTION, $default_update_action,$BINA\
RIES_ONLY,$force, $default_update_action, $INSTALL\
_DIR, $PLUGINS_DIR, $DISTRIBUTIONS,$tclinkdb, $pro\
xy, $clean);\nif ( ($cl=~/-root/)){$ROOT_INSTALL=1\
;}\nif ( ($cl=~/-no_question/)){$NO_QUESTION=1;}\n\
if ( ($cl=~/-update/)){$default_update_action=\"up\
date\";}\n$BINARIES_ONLY=1;\nif ( ($cl=~/-nobinari\
es/)){$BINARIES_ONLY=0;}\nif ( ($cl=~/-force/)){$f\
orce=1;$default_update_action=\"update\"}\nif ( ($\
cl=~/-exec=\\s*(\\S+)/)){$INSTALL_DIR=$1;}\nif ( (\
$cl=~/-plugins=\\s*(\\S+)/)){$PLUGINS_DIR=$1;}\nif\
 ( ($cl=~/-dis=\\s*(\\S+)/)){$DISTRIBUTIONS=$1;}\n\
\nif ( ($cl=~/-tclinkdb=\\s*(\\S+)/)){$tclinkdb=$1\
;}\nif ( ($cl=~/-proxy=\\s*(\\S+)/)){$proxy=$1;}\n\
if ( ($cl=~/-clean/)){$clean=1;}\nif ( ($cl=~/-rep\
o=\\s*(\\S+)/)){ $REPO_ROOT=$1; }\nif ( ($cl=~/-tc\
dir=\\s*(\\S+)/)){ $TCDIR=$1; }\nif ($tclinkdb){&u\
pdate_tclinkdb ($tclinkdb);}\n\n\nif( $REPO_ROOT n\
e \"\" ) {\n	if( $OSNAME eq \"\" ) { print \"You h\
ave specified the repository folder but the requir\
ed \\\"OSNAME\\\" enviroment variable is missing. \
\\n\"; exit 1; } \n	if( $OSARCH eq \"\" ) { print \
\"You have specified the repository folder but the\
 required \\\"OSARCH\\\" enviroment variable is mi\
ssing. \\n\"; exit 1; } \n}\n\n\nif(!$TCDIR) { $TC\
DIR=\"$HOME/.t_coffee\"; }\n&add_dir ($TCDIR);\n&a\
dd_dir ($TCCACHE=\"$TCDIR/cache\");\n&add_dir ($TC\
TMP=\"$CDIR/tmp\");\n&add_dir ($TCM=\"$TCDIR/mcoff\
ee\");\n&add_dir ($TCMETHODS=\"$TCDIR/methods\");\\
n&add_dir ($TCPLUGINS=\"$TCDIR/plugins/$OS\");\n\n\
\nour $BASE=\"$CD/bin\";\nour $BIN=\"$BASE/binarie\
s/$OS\";\nour $DOWNLOAD_DIR=\"$BASE/download\";\no\
ur $DOWNLOAD_FILE=\"$DOWNLOAD_DIR/files\";\nour $T\
MP=\"$BASE/tmp\";\n\n&add_dir($BASE);\n&add_dir($B\
IN);\n&add_dir($DOWNLOAD_DIR);\n&add_dir($DOWNLOAD\
_FILE);\nif (!$DISTRIBUTIONS){$DISTRIBUTIONS=\"$DO\
WNLOAD_DIR/distributions\";}\n&add_dir ($DISTRIBUT\
IONS);\n&add_dir ($TMP);\n\n\nif    (!$PLUGINS_DIR\
 && !$ROOT_INSTALL){$PLUGINS_DIR=$TCPLUGINS;}\nels\
if (!$PLUGINS_DIR &&  $ROOT_INSTALL){$PLUGINS_DIR=\
\"/usr/local/bin/\";}\n\nif    (!$INSTALL_DIR && !\
$ROOT_INSTALL){$INSTALL_DIR=\"$HOME/bin/\";mkpath \
($INSTALL_DIR);}\nelsif (!$INSTALL_DIR &&  $ROOT_I\
NSTALL){$INSTALL_DIR=\"/usr/local/bin/\";}\n\nif (\
-d \"mcoffee\"){`cp mcoffee/* $TCM`;}\n\n\nour $EN\
V_FILE=\"$TCDIR/t_coffee_env\";\n&env_file2putenv \
($ENV_FILE);\n&set_proxy($proxy);\nmy ($target, $p\
, $r);\n$target=$p;\n\nforeach $p (  ((keys (%PG))\
,(keys(%MODE)),(@smode)) )\n  {\n    if ($ARGV[0] \
eq $p && $target eq \"\"){$target=$p;}\n  }\nif ($\
target eq \"\"){exit ($EXIT_FAILURE);}\n\n\nforeac\
h $r (@required_applications)\n  {\n    my @app_li\
st;\n    my $i;\n    $i=0;\n    \n    @app_list=sp\
lit (/_OR_/, $r);\n    foreach my $pg (@app_list)\\
n      {\n	$i+=&pg_is_installed ($pg);\n      }\n \
   if ($i==0)\n      {\n      print \"One of the f\
ollowing packages must be installed to proceed: \"\
;\n      foreach my $pg (@app_list)\n	{\n	  print \
(\"$pg \");\n	}\n      die;\n    }\n  }\n\n\n\n\n\\
n\n&sign_license_ni();\n\n\n$PG{C}{compiler}=get_C\
_compiler($CC);\n$PG{Fortran}{compiler}=get_F_comp\
iler($FC);\n$PG{CXX}{compiler}=$PG{CPP}{compiler}=\
$PG{GPP}{compiler}=get_CXX_compiler($CXX);\nif ($C\
XXFLAGS){$PG{CPP}{options}=$PG{GPP}{options}=$PG{C\
XX}{options}=$CXXFLAGS;}\nif ($CFLAGS ne \"\" ){$P\
G{C}{options}=$CFLAGS;}\nforeach my $c (keys(%PG))\
\n  {\n    my $arguments;\n    if ($PG{$c}{compile\
r})\n      {\n	$arguments=\"$PG{$c}{compiler_flag}\
=$PG{$c}{compiler} \";\n	if ($PG{$c}{options})\n	 \
 {\n	    $arguments.=\"$PG{$c}{options_flag}='\" .\
 $PG{$c}{options} . \"' \";\n	  }\n	$PG{$c}{argume\
nts}=$arguments;\n      }\n  }\n\nif ($PG{$target}\
){$PG{$target}{install}=1;}\nelse\n  {\n    foreac\
h my $pg (keys(%PG))\n      {\n	if ( $target eq \"\
all\" || ($PG{$pg}{mode}=~/$target/))\n	  {\n	    \
$PG{$pg} {install}=1;\n	  }\n      }\n  }\n\nforea\
ch my $pg (keys(%PG))\n  {\n    if (!$PG{$pg}{upda\
te_action}){$PG{$pg}{update_action}=$default_updat\
e_action;}\n    elsif ($PG{$pg}{update_action} eq \
\"never\"){$PG{$pg}{install}=0;}\n    if ( $force \
&& $PG{$pg}{install})\n      {\n	`rm $BIN/$pg $BIN\
/$pg.exe $SILENT`;\n      }\n    if ($PG{$pg}{upda\
te_action} eq \"update\" && $PG{$pg}{install}){$PG\
{$pg}{update}=1;}\n  }\n\nif (($target=~/clean/))\\
n  {\n    print \"------- cleaning executables ---\
--\\n\";\n    `rm bin/* $SILENT`;\n    exit ($EXIT\
_SUCCESS);\n  }\n\nif ( !$PG{$target}){print \"---\
---- Installing T-Coffee Modes\\n\";}\n\nforeach m\
y $m (keys(%MODE))\n  {\n    if ( $target eq \"all\
\" || $target eq $m)\n      {\n	print \"\\n-------\
 The installer will now install the $m components \
$MODE{$m}{description}\\n\";\n	foreach my $pg (key\
s(%PG))\n	  {\n	    if ( $PG{$pg}{mode} =~/$m/ && \
$PG{$pg}{install})\n	      {\n		if ($PG{$pg}{touch\
ed}){print \"------- $PG{$pg}{dname}: already proc\
essed\\n\";}\n		else {$PG{$pg}{success}=&install_p\
g($pg);$PG{$pg}{touched}=1;}\n	      }\n	  }\n    \
  }\n  }\n\nif ( $PG{$target}){print \"------- Ins\
talling Individual Package\\n\";}\nforeach my $pg \
(keys (%PG))\n  {\n    \n    if ( $PG{$pg}{install\
} && !$PG{$pg}{touched})\n      {\n	print \"\\n---\
---- Install $pg\\n\";\n	$PG{$pg}{success}=&instal\
l_pg($pg);$PG{$pg}{touched}=1;\n      }\n  }\nprin\
t \"------- Finishing The installation\\n\";\nmy $\
final_report=&install ($INSTALL_DIR);\n\nprint \"\\
\n\";\nprint \"***********************************\
**********************************\\n\";\nprint \"\
********              INSTALLATION SUMMARY        \
  *****************\\n\";\nprint \"***************\
**************************************************\
****\\n\";\nprint \"------- SUMMARY package Instal\
lation:\\n\";\nprint \"-------   Executable Instal\
led in: $PLUGINS_DIR\\n\";\n\nforeach my $pg (keys\
(%PG))\n  {\n    if ( $PG{$pg}{install})\n      {\\
n	my $bin_status=($PG{$pg}{from_binary} && $PG{$pg\
}{success})?\"[from binary]\":\"\";\n	if     ( $PG\
{$pg}{new} && !$PG{$pg}{old})                     \
{print \"*------        $PG{$pg}{dname}: installed\
 $bin_status\\n\"; $PG{$pg}{status}=1;}\n	elsif  (\
 $PG{$pg}{new} &&  $PG{$pg}{old})                 \
    {print \"*------        $PG{$pg}{dname}: updat\
ed $bin_status\\n\"  ; $PG{$pg}{status}=1;} \n	els\
if  (!$PG{$pg}{new} &&  $PG{$pg}{old} && !$PG{$pg}\
{update}){print \"*------        $PG{$pg}{dname}: \
previous\\n\" ; $PG{$pg}{status}=1;}\n	elsif  (!$P\
G{$pg}{new} &&  $PG{$pg}{old} &&  $PG{$pg}{update}\
){print \"*------        $PG{$pg}{dname}: failed u\
pdate (previous installation available)\\n\";$PG{$\
pg}{status}=0;}\n	else                            \
                              {print \"*------    \
    $PG{$pg}{dname}: failed installation\\n\";$PG{\
$pg}{status}=0;}\n      }\n  }\nmy $failure;\n\nif\
 ( !$PG{$target}){print \"*------ SUMMARY mode Ins\
tallation:\\n\";}\nforeach my $m (keys(%MODE))\n  \
{\n  \n    if ( $target eq \"all\" || $target eq $\
m)\n      {\n	my $succesful=1;\n	foreach my $pg (k\
eys(%PG))\n	  {\n	    if (($PG{$pg}{mode}=~/$m/) &\
& $PG{$pg}{install} && $PG{$pg}{status}==0)\n	    \
  {\n		$succesful=0;\n		print \"*!!!!!!       $PG{\
$pg}{dname}: Missing\\n\";\n	      }\n	  }\n	if ( \
$succesful)\n	  {\n	    $MODE{$m}{status}=1;\n	   \
 print \"*------       MODE $MODE{$m}{dname} SUCCE\
SSFULLY installed\\n\";\n	  }\n	else\n	  {\n	    $\
failure++;\n	    $MODE{$m}{status}=0;\n	    print \
\"*!!!!!!       MODE $MODE{$m}{dname} UNSUCCESSFUL\
LY installed\\n\";\n	  }\n      }\n  }\n\n    \n  \
    \nif ($clean==1 && ($BASE=~/install4tcoffee/) \
){print \"*------ Clean Installation Directory: $B\
ASE\\n\";`rm -rf $BASE`;}\nforeach my $pg (keys(%P\
G)){if ($PG{$pg}{install} && $PG{$pg}{status}==0){\
exit ($EXIT_FAILURE);}}\n\nif ($failure)\n  {\n   \
 print \"*****************************************\
****************************\\n\";\n    print \"**\
******     SOME PACKAGES FAILED TO INSTALL        \
*****************\\n\";\n    print \"*************\
**************************************************\
******\\n\";\n    print \"\\nSome of the reported \
failures may be due to connectivity problems\";\n \
   print \"\\nRerun the installation and the insta\
ller will specifically try to install the missing \
packages\";\n    print \"\\nIf this Fails, go to t\
he original website and install the package manual\
ly\";\n  }\n\nprint \"****************************\
*****************************************\\n\";\np\
rint \"********              FINALIZE YOUR INSTALL\
ATION    *****************\\n\";\nprint \"********\
**************************************************\
***********\\n\";\nprint \"------- Your executable\
s are in:\\n\"; \nprint \"-------       $PLUGINS_D\
IR:\\n\";\nprint \"------- Add this directory to y\
our path with the following command:\\n\";\nprint \
\"-------       export PATH=$PLUGINS_DIR:\\$PATH\\\
n\";\nprint \"------- Make this permanent by addin\
g this line to the file:\\n\";\nprint \"-------   \
    $HOME/.bashrc\\n\";\nexit ($EXIT_SUCCESS);  \n\
  \nsub get_CXX_compiler\n  {\n    my $c=@_[0];\n \
   my (@clist)=(\"g++\");\n    \n    return get_co\
mpil ($c, @clist);\n }\nsub get_C_compiler\n  {\n \
   my $c=@_[0];\n    my (@clist)=(\"gcc\", \"cc\",\
 \"icc\");\n    \n    return get_compil ($c, @clis\
t);\n }\n\nsub get_F_compiler\n  {\n    my ($c)=@_\
[0];\n    my @clist=(\"f77\", \"g77\",\"g95\", \"g\
fortran\", \"ifort\");\n    return get_compil ($c,\
 @clist);\n  } \n       \nsub get_compil\n  {\n   \
 my ($fav,@clist)=(@_);\n    \n    #return the fir\
st compiler found installed in the system. Check f\
irst the favorite\n    foreach my $c ($fav,@clist)\
\n      {\n	if  (&pg_is_installed ($c)){return $c;\
}\n      }\n    return \"\";\n  }\nsub exit_if_pg_\
not_installed\n  {\n    my (@arg)=(@_);\n    \n   \
 foreach my $p (@arg)\n      {\n	if ( !&pg_is_inst\
alled ($p))\n	  {\n	    print \"!!!!!!!! The $p ut\
ility must be installed for this installation to p\
roceed [FATAL]\\n\";\n	    die;\n	  }\n      }\n  \
  return 1;\n  }\nsub set_proxy\n  {\n    my ($pro\
xy)=(@_);\n    my (@list,$p);\n    \n    @list= (\\
"HTTP_proxy\", \"http_proxy\", \"HTTP_PROXY\", \"A\
LL_proxy\", \"all_proxy\",\"HTTP_proxy_4_TCOFFEE\"\
,\"http_proxy_4_TCOFFEE\");\n    \n    if (!$proxy\
)\n      {\n	foreach my $p (@list)\n	  {\n	    if \
( ($ENV_SET{$p}) || $ENV{$p}){$proxy=$ENV{$p};}\n	\
  }\n      }\n    foreach my $p(@list){$ENV{$p}=$p\
roxy;}\n  }\n	\nsub check_internet_connection\n  {\
\n    my $internet;\n    \n    if ( -e \"x\"){unli\
nk (\"x\");}\n    if     (&pg_is_installed    (\"w\
get\")){`wget www.google.com -Ox >/dev/null 2>/dev\
/null`;}\n    elsif  (&pg_is_installed    (\"curl\\
")){`curl www.google.com -ox >/dev/null 2>/dev/nul\
l`;}\n    else\n      {\n	printf stderr \"\\nERROR\
: No pg for remote file fetching [wget or curl][FA\
TAL]\\n\";\n	exit ($EXIT_FAILURE);\n      }\n    \\
n    if ( !-e \"x\" || -s \"x\" < 10){$internet=0;\
}\n    else {$internet=1;}\n    if (-e \"x\"){unli\
nk \"x\";}\n    return $internet;\n  }\nsub url2fi\
le\n  {\n    my ($cmd, $file,$wget_arg, $curl_arg)\
=(@_);\n    my ($exit,$flag, $pg, $arg);\n    \n  \
  if ($INTERNET || check_internet_connection ()){$\
INTERNET=1;}\n    else\n      {\n	print STDERR \"E\
RROR: No Internet Connection [FATAL:install.pl]\\n\
\";\n	exit ($EXIT_FAILURE);\n      }\n    \n    if\
     (&pg_is_installed    (\"wget\")){$pg=\"wget\"\
; $flag=\"-O\";$arg=\"--tries=2 --connect-timeout=\
10 --no-check-certificate $wget_arg\";}\n    elsif\
  (&pg_is_installed    (\"curl\")){$pg=\"curl\"; $\
flag=\"-o\";$arg=$curl_arg;}\n    else\n      {\n	\
printf stderr \"\\nERROR: No pg for remote file fe\
tching [wget or curl][FATAL]\\n\";\n	exit ($EXIT_F\
AILURE);\n      }\n    \n    \n    if (-e $file){u\
nlink($file);}\n    $exit=system \"$pg $cmd $flag$\
file $arg\";\n    return $exit;\n  }\n\nsub pg_is_\
installed\n  {\n    my ($p, $dir)=(@_);\n    my ($\
r,$m, $ret);\n    my ($supported, $language, $comp\
il);\n    \n  \n    if ( $PG{$p})\n      {\n	$lang\
uage=$PG{$p}{language2};\n	$compil=$PG{$language}{\
compiler};\n      }\n    \n    if ( $compil eq \"C\
PAN\")\n      {\n	if ( system (\"perl -M$p -e 1\")\
==$EXIT_SUCCESS){$ret=1;}\n	else {$ret=0;}\n      \
}\n    elsif ($dir)\n      {\n	if (-e \"$dir/$p\" \
|| -e \"$dir/$p\\.exe\"){$ret=1;}\n	else {$ret=0;}\
\n      }\n    elsif (-e \"$PLUGINS_DIR/$p\" || -e\
 \"$PLUGINS_DIR/$p.exe\"){$ret=1;}\n    else\n    \
  {\n	$r=`which $p 2>/dev/null`;\n	if ($r eq \"\")\
{$ret=0;}\n	else {$ret=1;}\n      }\n   \n    retu\
rn $ret;\n  }\nsub install\n  {\n    my ($new_bin)\
=(@_);\n    my ($copied, $report);\n\n    \n    if\
 (!$ROOT_INSTALL)\n      {\n	\n	if (-e \"$BIN/t_co\
ffee\"){`$CP $BIN/t_coffee $INSTALL_DIR`};\n	`cp $\
BIN/* $PLUGINS_DIR`;\n	$copied=1;\n      }\n    el\
se\n      {\n	$copied=&root_run (\"You must be roo\
t to finalize the installation\", \"$CP $BIN/* $IN\
STALL_DIR $SILENT\");\n      }\n    \n     \n  if \
( !$copied)\n    {\n      $report=\"*!!!!!! Instal\
lation unsuccesful. The executables have been left\
 in $BASE/bin\\n\";\n    }\n  elsif ( $copied && $\
ROOT)\n    {\n      $report=\"*------ Installation\
 succesful. Your executables have been copied in $\
new_bin and are on your PATH\\n\";\n    }\n  elsif\
 ( $copied && !$ROOT)\n    {\n      $report= \"*!!\
!!!! T-Coffee and associated packages have been co\
pied in: $new_bin\\n\";\n      $report.=\"*!!!!!! \
This address is NOT in your PATH sytem variable\\n\
\";\n      $report.=\"*!!!!!! You can do so by add\
ing the following line in your ~/.bashrc file:\\n\\
";\n      $report.=\"*!!!!!! export PATH=$new_bin:\
\\$PATH\\n\";\n    }\n  return $report;\n}\n\nsub \
sign_license_ni\n  {\n    my $F=new FileHandle;\n \
   open ($F, \"license.txt\");\n    while (<$F>)\n\
      {\n	print \"$_\";\n      }\n    close ($F);\\
n    \n    return;\n  }\n\nsub install_pg\n  {\n  \
  my ($pg)=(@_);\n    my ($report, $previous, $lan\
guage, $compiler, $return);\n    \n    if (!$PG{$p\
g}{install}){return 1;}\n    \n    $previous=&pg_i\
s_installed ($pg);\n    \n    if ($PG{$pg}{update_\
action} eq \"no_update\" && $previous)\n      {\n	\
$PG{$pg}{old}=1;\n	$PG{$pg}{new}=0;\n	$return=1;\n\
      }\n    else\n      {\n	$PG{$pg}{old}=$previo\
us;\n	\n	if ($PG{$pg} {language2} eq \"Perl\"){&in\
stall_perl_package ($pg);}\n	elsif ($pg ne \"t_cof\
fee\" && $BINARIES_ONLY && &install_binary_package\
 ($pg)){$PG{$pg}{from_binary}=1;}\n	elsif (&instal\
l_source_package ($pg)){;}\n	else \n	  {\n	    \n	\
    if (!&supported_os($OS))\n	      {\n		print \"\
!!!!!!!! $pg compilation failed, binary unsupporte\
d for $OS\\n\"; \n	      }\n	    elsif (!($PG{$pg}\
{from_binary}=&install_binary_package ($pg)))\n	  \
    {\n		print \"!!!!!!!! $pg compilation and  bin\
ary installation failed\\n\";\n	      }\n	  }\n	$P\
G{$pg}{new}=$return=&pg_is_installed ($pg,$BIN);\n\
      }\n\n    \n    return $return;\n  }\nsub ins\
tall_perl_package\n  {\n    my ($pg)=(@_);\n    my\
 ($report, $language, $compiler);\n    \n    $lang\
uage=$PG{$pg} {language2};\n    $compiler=$PG{$lan\
guage}{compiler};\n    \n    if (!&pg_is_installed\
 ($pg))\n      {\n	if ( $OS eq \"windows\"){`perl \
-M$compiler -e 'install $pg'`;}\n	elsif ( $ROOT eq\
 \"sudo\"){system (\"sudo perl -M$compiler -e 'ins\
tall $pg'\");}\n	else {system (\"su root -c perl -\
M$compiler -e 'install $pg'\");}\n      }\n    ret\
urn &pg_is_installed ($pg);\n  }\n\n\n\nsub instal\
l_source_package\n  {\n    my ($pg)=(@_);\n    my \
($report, $download, $arguments, $language, $addre\
ss, $name, $ext, $main_dir, $distrib);\n    my $wg\
et_tmp=\"$TMP/wget.tmp\";\n    my (@fl);\n    if (\
 -e \"$BIN/$pg\" || -e \"$BIN/$pg.exe\"){return 1;\
}\n    \n    #\n    # check if the module exists i\
n the repository cache \n    #\n	if( repo_load($pg\
) ) {\n		return 1;\n	}\n    \n    if ($pg eq \"t_c\
offee\")  {return   &install_t_coffee ($pg);}\n   \
 elsif ($pg eq \"TMalign\"){return   &install_TMal\
ign ($pg);}\n    \n    chdir $DISTRIBUTIONS;\n    \
\n    $download=$PG{$pg}{source};\n    \n    if ((\
$download =~/tgz/))\n      {\n	($address,$name,$ex\
t)=($download=~/(.+\\/)([^\\/]+)(\\.tgz).*/);\n   \
   }\n    elsif (($download=~/tar\\.gz/))\n      {\
\n	($address,$name,$ext)=($download=~/(.+\\/)([^\\\
/]+)(\\.tar\\.gz).*/);\n      }\n    elsif (($down\
load=~/tar/))\n      {\n	($address,$name,$ext)=($d\
ownload=~/(.+\\/)([^\\/]+)(\\.tar).*/);\n      }\n\
    else\n      {\n	($address,$name)=($download=~/\
(.+\\/)([^\\/]+)/);\n	$ext=\"\";\n      }\n    $di\
strib=\"$name$ext\";\n    \n    if ( !-d $pg){mkdi\
r $pg;}\n    chdir $pg;\n   \n    #get the distrib\
ution if available\n    if ( -e \"$DOWNLOAD_DIR/$d\
istrib\")\n      {\n	`$CP $DOWNLOAD_DIR/$distrib .\
`;\n      }\n    #UNTAR and Prepare everything\n  \
  if (!-e \"$name.tar\" && !-e \"$name\")\n      {\
\n	&check_rm ($wget_tmp);\n	print \"\\n------- Dow\
nloading/Installing $pg\\n\";\n	\n	if (!-e $distri\
b && &url2file (\"$download\", \"$wget_tmp\")==$EX\
IT_SUCCESS)\n	  {\n	    \n	    `mv $wget_tmp $dist\
rib`;\n	    `$CP $distrib $DOWNLOAD_DIR/`;\n	  }\n\
\n	if (!-e $distrib)\n	  {\n	    print \"!!!!!!! D\
ownload of $pg distribution failed\\n\";\n	    pri\
nt \"!!!!!!! Check Address: $PG{$pg}{source}\\n\";\
\n	    return 0;\n	  }\n	print \"\\n------- unzipp\
ing/untaring $name\\n\";\n	if (($ext =~/z/))\n	  {\
 \n	    &flush_command (\"gunzip $name$ext\");\n	 \
   \n	  }\n	if (($ext =~/tar/) || ($ext =~/tgz/))\\
n	  {\n	    &flush_command(\"tar -xvf $name.tar\")\
;\n	  }\n      }\n    #Guess and enter the distrib\
ution directory\n    @fl=ls($p);\n    foreach my $\
f (@fl)\n      {\n	if (-d $f)\n	  {\n	    $main_di\
r=$f;\n	  }\n      }\n    if (-d $main_dir)\n	  \n\
      {\n	chdir $main_dir;}\n    else\n      {\n	p\
rint \"Error: $main_dir does not exist\";\n      }\
\n    print \"\\n------- Compiling/Installing $pg\\
\n\";\n    `make clean $SILENT`;\n    \n    \n    \
#\n    # SAP module\n    #\n    if ($pg eq \"sap\"\
)\n      {\n	if (-e \"./configure\")\n	  {\n	    #\
new sap distribution\n	    \n	    &flush_command (\
\"./configure\");\n	    &flush_command (\"make cle\
an\");\n	    &flush_command (\"make\");\n	    &che\
ck_cp (\"./src/$pg\", \"$BIN\");\n	    repo_store(\
\"./src/$pg\");\n	  }\n	else\n	  {\n	    #old styl\
e distribution\n	    `rm *.o sap  sap.exe ./util/a\
a/*.o  ./util/wt/.o $SILENT`;\n	    &flush_command\
 (\"make $arguments sap\");\n	    &check_cp ($pg, \
\"$BIN\");\n	    repo_store($pg);\n	  }\n      }\n\
    \n    #\n    # CLUSTALW2 module\n    #\n    el\
sif ($pg eq \"clustalw2\")\n      {\n	&flush_comma\
nd(\"./configure\");\n	&flush_command(\"make $argu\
ments\");\n	&check_cp (\"./src/$pg\", \"$BIN\");\n\
	repo_store(\"./src/$pg\");\n      }\n\n    #\n   \
 # CLUSTAL-OMEGA module\n    #\n    elsif ($pg eq \
\"clustalo\")\n      {\n	&flush_command(\"./config\
ure\");\n	&flush_command(\"make $arguments\");\n	&\
check_cp (\"./src/$pg\", \"$BIN\");\n	repo_store(\\
"./src/$pg\");\n      }\n\n    #\n    # STRIKE mod\
ule\n    #\n    elsif ($pg eq \"strike\")\n      {\
\n	&flush_command(\"make $arguments\");\n	&check_c\
p (\"./bin/$pg\", \"$BIN\");\n	repo_store(\"./bin/\
$pg\");\n      }\n    \n    #\n    # FSA module\n \
   # \n    elsif ($pg eq \"fsa\")\n      {\n	&flus\
h_command(\"./configure --prefix=$BIN\");\n	&flush\
_command(\"make $arguments\");\n	&flush_command (\\
"make install\");\n\n	repo_store(\"fsa\", \"$BIN/b\
in\");\n	`mv $BIN/bin/* $BIN`;\n	`rmdir $BIN/bin`;\
\n      }\n    \n    #\n    # CLUSTALW module\n   \
 #\n    elsif ($pg eq \"clustalw\")\n      {\n	&fl\
ush_command(\"make $arguments clustalw\");\n	`$CP \
$pg $BIN $SILENT`;\n	repo_store($pg);\n      }\n  \
  \n    #\n    # MAFFT module\n    #\n    elsif ($\
pg eq \"mafft\")\n      {\n	my $base=cwd();\n	my $\
c;\n	\n	#compile core\n	mkpath (\"./mafft/bin\");\\
n	mkpath (\"./mafft/lib\");\n	chdir \"$base/core\"\
;\n	`make clean $SILENT`;\n	&flush_command (\"make\
 $arguments\");\n	&flush_command (\"make install L\
IBDIR=../mafft/lib BINDIR=../mafft/bin\");\n	\n	#c\
ompile extension\n	chdir \"$base/extensions\";\n	`\
make clean $SILENT`;\n	&flush_command (\"make $arg\
uments\");\n	&flush_command (\"make install LIBDIR\
=../mafft/lib BINDIR=../mafft/bin\");\n	\n	#put ev\
erything in mafft and copy the compiled stuff in b\
in\n	chdir \"$base\";\n	if ($ROOT_INSTALL)\n	  {\n\
	    &root_run (\"You Must be Root to Install MAFF\
T\\n\", \"mkdir /usr/local/mafft/;$CP mafft/lib/* \
/usr/local/mafft;$CP mafft/lib/mafft* /usr/local/b\
in ;$CP mafft/bin/mafft /usr/local/bin/; \");\n	  \
}\n	else\n	  {\n	    `$CP mafft/lib/*  $BIN`;\n	  \
  `$CP mafft/bin/mafft  $BIN`;\n	  }\n	`tar -cvf m\
afft.tar mafft`;\n	`gzip mafft.tar`;\n	`mv mafft.t\
ar.gz $BIN`;\n	\n	repo_store(\"mafft/bin/mafft\", \
\"mafft/lib/\", \"$BIN/mafft.tar.gz\");\n      }\n\
      \n    #\n    # DIALIGN-TX module\n    #\n   \
 elsif ( $pg eq \"dialign-tx\" )\n      {\n	my $f;\
\n	my $base=cwd();\n\n	chdir \"./source\";\n	if ($\
OS eq \"macosx\"){&flush_command (\"cp makefile.MA\
C_OS makefile\");}\n\n	&flush_command (\" make CPP\
FLAGS='-O3 -funroll-loops' all\");\n	\n	chdir \"..\
\";\n	&check_cp (\"./source/$pg\", \"$BIN\");\n	re\
po_store(\"./source/$pg\");\n      }\n      \n    \
#\n    # DIALIGN-T module \n    # (is the same as \
dialign-tx, but it is mantained for backward name \
compatibility with tcoffee)\n    #\n    elsif ( $p\
g eq \"dialign-t\" )\n      {\n	my $f;\n	my $base=\
cwd();\n\n	chdir \"./source\";\n	if ($OS eq \"maco\
sx\"){&flush_command (\"cp makefile.MAC_OS makefil\
e\");}\n\n	&flush_command (\" make CPPFLAGS='-O3 -\
funroll-loops' all\");\n	\n	chdir \"..\";\n	&check\
_cp (\"./source/dialign-tx\", \"$BIN/dialign-t\");\
\n	repo_store(\"$BIN/dialign-t\");	\n      }      \
\n      \n    #\n    # POA module\n    #\n    elsi\
f ($pg eq \"poa\")\n      {\n	&flush_command (\"ma\
ke $arguments poa\");\n	&check_cp (\"$pg\", \"$BIN\
\");\n	repo_store(\"$pg\");\n      }\n     \n     \
\n    #\n    # PROBCONS module\n    #\n    elsif (\
 $pg eq \"probcons\")\n      {\n	&add_C_libraries(\
\"./ProbabilisticModel.h\", \"list\", \"cstring\")\
;\n	\n	`rm *.exe $SILENT`;\n	&flush_command (\"mak\
e $arguments probcons\");\n	&check_cp(\"$pg\", \"$\
BIN/$pg\");\n	repo_store(\"$pg\");\n      }\n     \
 \n    #\n    # PROBCONS RNA module\n    #\n    el\
sif ( $pg eq \"probconsRNA\")\n      {\n	&add_C_li\
braries(\"./ProbabilisticModel.h\", \"list\", \"cs\
tring\");\n	&add_C_libraries(\"./Main.cc\", \"ioma\
nip\", \"cstring\",\"climits\");\n	`rm *.exe $SILE\
NT`;\n	&flush_command (\"make $arguments probcons\\
");\n	&check_cp(\"probcons\", \"$BIN/$pg\");\n	rep\
o_store(\"$BIN/$pg\");\n      }\n\n	#\n	# MUSCLE m\
odule\n	#\n    elsif (  $pg eq \"muscle\")\n      \
{	\n	`rm *.o muscle muscle.exe $SILENT`;\n	if ($OS\
 eq \"macosx\" || $OS eq \"linux\")\n	  {\n	    &r\
eplace_line_in_file (\"./Makefile\", \"LDLIBS = -l\
m -static\",  \"LDLIBS = -lm\");\n	  }\n	elsif ($O\
S eq \"windows\")\n	  {\n	    &replace_line_in_fil\
e (\"./intmath.cpp\",  \"double log2e\",      \"do\
uble cedric_log\");\n	    &replace_line_in_file (\\
"./intmath.cpp\",  \"double log2\",       \"double\
 log_notuse\");\n	    &replace_line_in_file (\"./i\
ntmath.cpp\",  \"double cedric_log\", \"double log\
2e\");\n	  }\n	&flush_command (\"make $arguments a\
ll\");\n	&check_cp(\"$pg\", \"$BIN\");\n	repo_stor\
e(\"$pg\");	\n      }\n      \n     #\n     # MUS4\
 module\n     #\n     elsif (  $pg eq \"mus4\")\n \
     {\n	`rm *.o muscle muscle.exe $SILENT`;\n	&fl\
ush_command (\"./mk\");\n	&check_cp(\"$pg\", \"$BI\
N\");\n	repo_store(\"$pg\");	\n      }\n      \n  \
  #\n    # PCMA module\n    #\n    elsif ( $pg eq \
\"pcma\")\n      {\n	if ($OS eq \"macosx\")\n	  {\\
n	    &replace_line_in_file (\"./alcomp2.c\", \"ma\
lloc.h\",  \"\");\n	  }\n	&flush_command (\"make $\
arguments pcma\");\n	&check_cp(\"$pg\", \"$BIN\");\
\n	repo_store(\"$pg\");	\n      }\n      \n    #\n\
    # KALIGN module\n    #\n    elsif ($pg eq \"ka\
lign\")\n      {\n	&flush_command (\"./configure\"\
);\n	&flush_command(\"make $arguments\");\n	&check\
_cp (\"$pg\",$BIN);\n	repo_store(\"$pg\");	\n     \
 }\n      \n    #\n    # AMAP module\n    #\n    e\
lsif ( $pg eq \"amap\")\n      {\n	&add_C_librarie\
s(\"./Amap.cc\", \"iomanip\", \"cstring\",\"climit\
s\");	\n	`make clean $SILENT`;\n	&flush_command (\\
"make $arguments all\");\n	&check_cp (\"$pg\", $BI\
N);\n	repo_store(\"$pg\");	\n      }\n      \n    \
#\n    # PRODA module\n    #\n    elsif ( $pg eq \\
"proda\")\n      {\n	`sed -i '' 's/int errno = 0;/\
int errno; errno = 0;/' Main.cc`;\n	&add_C_librari\
es(\"AlignedFragment.h\", \"vector\", \"iostream\"\
, \"cstring\",\"cstdlib\");\n	&add_C_libraries(\"M\
ain.cc\", \"vector\", \"climits\");	\n	&add_C_libr\
aries(\"Sequence.cc\", \"stdlib.h\", \"cstdio\");	\
\n	&flush_command (\"make $arguments all\");\n	&ch\
eck_cp (\"$pg\", $BIN);\n	repo_store(\"$pg\");	\n \
     }\n      \n    #\n    # PRANK module\n    #\n\
    elsif ( $pg eq \"prank\")\n      {\n	&flush_co\
mmand (\"make $arguments all\");\n	&check_cp (\"$p\
g\", $BIN);\n	repo_store(\"$pg\");	\n      }\n    \
  \n    #\n    # !!!! MUSTANG module\n    #\n     \
elsif ( $pg eq \"mustang\")\n      {\n	&flush_comm\
and (\"rm ./bin/*\");\n	&flush_command (\"make $ar\
guments all\");\n\n	if ( $OS=~/windows/){&flush_co\
mmand(\"cp ./bin/* $BIN/mustang.exe\");}\n	else {&\
flush_command(\"cp ./bin/* $BIN/mustang\");}\n	\n	\
repo_store(\"$BIN/mustang\");\n      }\n\n	#\n	# R\
NAplfold module\n	#\n    elsif ( $pg eq \"RNAplfol\
d\")\n      {\n	&flush_command(\"./configure\");\n\
	&flush_command (\"make $arguments all\");\n	&chec\
k_cp(\"./Progs/RNAplfold\", \"$BIN\");\n	&check_cp\
(\"./Progs/RNAalifold\", \"$BIN\");\n	&check_cp(\"\
./Progs/RNAfold\", \"$BIN\");\n	\n	repo_store(\"./\
Progs/RNAplfold\", \"./Progs/RNAalifold\", \"./Pro\
gs/RNAfold\");\n      }\n      \n    #\n    # !!! \
RETREE module\n    #\n    elsif ( $pg eq \"retree\\
")\n      {\n	chdir \"src\";\n	&flush_command (\"c\
p Makefile.unx Makefile\");\n	&flush_command (\"ma\
ke $arguments all\");\n	&flush_command (\"make put\
\");\n	system \"cp ../exe/* $BIN\";\n	\n	repo_stor\
e(\"retree\", \"../exe\");\n      }\n	\n    chdir \
$CDIR;\n    return &pg_is_installed ($pg, $BIN);\n\
  }\n\nsub install_t_coffee\n  {\n    my ($pg)=(@_\
);\n    my ($report,$cflags, $arguments, $language\
, $compiler) ;\n    #1-Install T-Coffee\n    chdir\
 \"t_coffee_source\";\n    &flush_command (\"make \
clean\");\n    print \"\\n------- Compiling T-Coff\
ee\\n\";\n    $language=$PG{$pg} {language2};\n   \
 $arguments=$PG{$language}{arguments};\n\n    if (\
 $CC ne \"\"){\n      print \"make -i $arguments t\
_coffee \\n\";\n      &flush_command (\"make -i $a\
rguments t_coffee\");\n    }\n    &check_cp ($pg, \
$BIN);\n    \n    chdir $CDIR;\n    return &pg_is_\
installed ($pg, $BIN);\n  }\nsub install_TMalign\n\
  {\n    my ($pg)=(@_);\n    my $report;\n    chdi\
r \"t_coffee_source\";\n    print \"\\n------- Com\
piling TMalign\\n\";\n    `rm TMalign TMalign.exe \
$SILENT`;\n    if ( $FC ne \"\"){&flush_command (\\
"make -i $PG{Fortran}{arguments} TMalign\");}\n   \
 &check_cp ($pg, $BIN);\n    repo_store($pg);\n\n \
   if ( !-e \"$BIN/$pg\" && pg_has_binary_distrib \
($pg))\n      {\n	print \"!!!!!!! Compilation of $\
pg impossible. Will try to install from binary\\n\\
";\n	return &install_binary_package ($pg);\n      \
}\n    chdir $CDIR;\n    return &pg_is_installed (\
$pg, $BIN);\n  }\n\nsub pg_has_binary_distrib\n  {\
\n    my ($pg)=(@_);\n    if ($PG{$pg}{windows}){r\
eturn 1;}\n    elsif ($PG{$pg}{osx}){return 1;}\n \
   elsif ($PG{$pg}{linux}){return 1;}\n    return \
0;\n  }\nsub install_binary_package\n  {\n    my (\
$pg)=(@_);\n    my ($base,$report,$name, $download\
, $arguments, $language, $dir);\n    my $isdir;\n \
   &input_os();\n    \n    #\n    # - paolodt - Ch\
eck if the module exists in the repository cache \\
n    #\n	if( repo_load($pg) ) {\n	    $PG{$pg}{fro\
m_binary}=1;\n		return 1;\n	}\n    # - paolodt - e\
nd \n    \n    if (!&supported_os($OS)){return 0;}\
\n    if ( $PG{$pg}{binary}){$name=$PG{$pg}{binary\
};}\n    else {$name=$pg;}\n    \n    $download=\"\
$WEB_BASE/Packages/Binaries/$OS/$name\";\n    \n  \
  $base=cwd();\n    chdir $TMP;\n    \n    if (!-e\
 $name)\n      {\n	`rm x $SILENT`;\n	if ( url2file\
(\"$download\",\"x\")==$EXIT_SUCCESS)\n	  {\n	    \
`mv x $name`;\n	  }\n      }\n    \n    if (!-e $n\
ame)\n      {\n	print \"!!!!!!! $PG{$pg}{dname}: D\
ownload of $pg binary failed\\n\";\n	print \"!!!!!\
!! $PG{$pg}{dname}: Check Address: $download\\n\";\
\n	return 0;\n      }\n    print \"\\n------- Inst\
alling $pg\\n\";\n    \n    if ($name =~/tar\\.gz/\
)\n      {\n	`gunzip  $name`;\n	`tar -xvf $pg.tar`\
;\n	chdir $pg;\n	`chmod u+x *`;\n 	`mv * $BIN`;\n	\
#if (!($pg=~/\\*/)){`rm -rf $pg`;}\n      }\n    e\
lse\n      {\n	&check_cp (\"$pg\", \"$BIN\");\n	`c\
hmod u+x $BIN/$pg`; \n	unlink ($pg);\n      }\n   \
 chdir $base;\n    $PG{$pg}{from_binary}=1;\n    r\
eturn &pg_is_installed ($pg, $BIN);\n  }\n\nsub ad\
d_dir \n  {\n    my $dir=@_[0];\n    \n    if (!-e\
 $dir && !-d $dir)\n      {\n	my @l;\n	umask (0000\
);\n	@l=mkpath ($dir,{mode => 0777});\n	\n      }\\
n    else\n      {\n	return 0;\n      }\n  }\nsub \
check_rm \n  {\n    my ($file)=(@_);\n    \n    if\
 ( -e $file)\n      {\n	return unlink($file);\n   \
   }\n    return 0;\n  }\nsub check_cp\n  {\n    m\
y ($from, $to)=(@_);\n    if ( !-e $from && -e \"$\
from\\.exe\"){$from=\"$from\\.exe\";}\n    if ( !-\
e $from){return 0;}\n        \n    `$CP $from $to`\
;\n    return 1;\n  }\n\nsub repo_store \n{\n   # \
check that all required data are available\n   if(\
 $REPO_ROOT eq \"\" ) { return; }\n\n\n    # extra\
ct the package name from the specified path\n    m\
y $pg =`basename $_[0]`;\n    chomp($pg);\n	\n    \
my $VER = $PG{$pg}{version};\n    my $CACHE = \"$R\
EPO_ROOT/$pg/$VER/$OSNAME-$OSARCH\"; \n    \n    p\
rint \"-------- Storing package: \\\"$pg\\\" to pa\
th: $CACHE\\n\";\n    \n    # clean the cache path\
 if exists and create it again\n    `rm -rf $CACHE\
`;\n    `mkdir -p $CACHE`;\n    \n 	for my $path (\
@_) {\n\n	    # check if it is a single file \n	 	\
if( -f $path ) {\n	    	`cp $path $CACHE`;\n		}\n	\
	# .. or a directory, in this case copy all the co\
ntent \n		elsif( -d $path ) {\n			opendir(IMD, $pa\
th);\n			my @thefiles= readdir(IMD);\n			closedir(\
IMD);\n			\n			for my $_file (@thefiles) {\n				if\
( $_file ne \".\" && $_file ne \"..\") {\n	    			\
`cp $path/$_file $CACHE`;\n				}\n			}\n		} \n	}	 \
  \n    \n	\n}   \n\nsub repo_load \n{\n    my ($p\
g)=(@_);\n\n    #Bypass the Repository Cache\n    \
return 0;\n    # check that all required data are \
available\n    if( $REPO_ROOT eq \"\" ) { return 0\
; }\n\n    my $VER = $PG{$pg}{version};\n    my $C\
ACHE = \"$REPO_ROOT/$pg/$VER/$OSNAME-$OSARCH\"; \n\
    if( !-e \"$CACHE/$pg\" ) {\n   	 	print \"----\
---- Module \\\"$pg\\\" NOT found on repository ca\
che.\\n\";\n    	return 0;\n    }\n    \n    print\
 \"-------- Module \\\"$pg\\\" found on repository\
 cache. Using copy on path: $CACHE\\n\";\n    `cp \
$CACHE/* $BIN`;\n    return 1;\n}\n\nsub check_fil\
e_list_exists \n  {\n    my ($base, @flist)=(@_);\\
n    my $f;\n\n    foreach $f (@flist)\n      {\n	\
if ( !-e \"$base/$f\"){return 0;}\n      }\n    re\
turn 1;\n  }\nsub ls\n  {\n    my $f=@_[0];\n    m\
y @fl;\n    chomp(@fl=`ls -1 $f`);\n    return @fl\
;\n  }\nsub flush_command\n  {\n    my $command=@_\
[0];\n    my $F=new FileHandle;\n    open ($F, \"$\
command|\");\n    while (<$F>){print \"    --- $_\\
";}\n    close ($F);\n  }    \n\nsub input_install\
ation_directory\n  {\n    my $dir=@_[0];\n    my $\
new;\n    \n    print \"------- The current instal\
lation directory is: [$dir]\\n\";\n    print \"???\
???? Return to keep the default or new value:\";\n\
   \n    if ($NO_QUESTION==0)\n      {\n	chomp ($n\
ew=<stdin>);\n	while ( $new ne \"\" && !input_yes \
(\"You have entered $new. Is this correct? ([y]/n)\
:\"))\n	  {\n	    print \"???????New installation \
directory:\";\n	    chomp ($new=<stdin>);\n	  }\n	\
$dir=($new eq \"\")?$dir:$new;\n	$dir=~s/\\/$//;\n\
      }\n    \n    if ( -d $dir){return $dir;}\n  \
  elsif (&root_run (\"You must be root to create $\
dir\",\"mkdir $dir\")==$EXIT_SUCCESS){return $dir;\
}\n    else\n      {\n	print \"!!!!!!! $dir could \
not be created\\n\";\n	if ( $NO_QUESTION)\n	  {\n	\
    return \"\";\n	  }\n	elsif ( &input_yes (\"???\
???? Do you want to provide a new directory([y]/n)\
?:\"))\n	  {\n	    return input_installation_direc\
tory ($dir);\n	  }\n	else\n	  {\n	    return \"\";\
\n	  }\n      }\n    \n  }\nsub input_yes\n  {\n  \
  my $question =@_[0];\n    my $answer;\n\n    if \
($NO_QUESTION==1){return 1;}\n    \n    if ($quest\
ion eq \"\"){$question=\"??????? Do you wish to pr\
oceed ([y]/n)?:\";}\n    print $question;\n    cho\
mp($answer=lc(<STDIN>));\n    if (($answer=~/^y/) \
|| $answer eq \"\"){return 1;}\n    elsif ( ($answ\
er=~/^n/)){return 0;}\n    else\n      {\n	return \
input_yes($question);\n      }\n  }\nsub root_run\\
n  {\n    my ($txt, $cmd)=(@_);\n    \n    if ( sy\
stem ($cmd)==$EXIT_SUCCESS){return $EXIT_SUCCESS;}\
\n    else \n      {\n	print \"------- $txt\\n\";\\
n	if ( $ROOT eq \"sudo\"){return system (\"sudo $c\
md\");}\n	else {return system (\"su root -c \\\"$c\
md\\\"\");}\n      }\n  }\nsub get_root\n  {\n    \
if (&pg_is_installed (\"sudo\")){return \"sudo\";}\
\n    else {return \"su\";}\n  }\n\nsub get_os\n  \
{\n    my $raw_os=`uname`;\n    my $os;\n\n    $ra\
w_os=lc ($raw_os);\n    \n    if ($raw_os =~/cygwi\
n/){$os=\"windows\";}\n    elsif ($raw_os =~/linux\
/){$os=\"linux\";}\n    elsif ($raw_os =~/osx/){$o\
s=\"macosx\";}\n    elsif ($raw_os =~/darwin/){$os\
=\"macosx\";}\n    else\n      {\n	$os=$raw_os;\n \
     }\n    return $os;\n  }\nsub input_os\n  {\n \
   my $answer;\n    if ($OS) {return $OS;}\n    \n\
    print \"??????? which os do you use: [w]indows\
, [l]inux, [m]acosx:?\";\n    $answer=lc(<STDIN>);\
\n\n    if (($answer=~/^m/)){$OS=\"macosx\";}\n   \
 elsif ( ($answer=~/^w/)){$OS=\"windows\";}\n    e\
lsif ( ($answer=~/^linux/)){$OS=\"linux\";}\n    \\
n    else\n      {\n	return &input_os();\n      }\\
n    return $OS;\n  }\n\nsub supported_os\n  {\n  \
  my ($os)=(@_[0]);\n    return $SUPPORTED_OS{$os}\
;\n  }\n    \n    \n\n\nsub update_tclinkdb \n  {\\
n    my $file =@_[0];\n    my $name;\n    my $F=ne\
w FileHandle;\n    my ($download, $address, $name,\
 $l, $db);\n    \n    if ( $file eq \"update\"){$f\
ile=$TCLINKDB_ADDRESS;}\n    \n    if ( $file =~/h\
ttp:\\/\\// || $file =~/ftp:\\/\\//)\n      {\n	($\
address, $name)=($download=~/(.*)\\/([^\\/]+)$/);\\
n	`rm x $SILENT`;\n	if (&url2file ($file,\"x\")==$\
EXIT_SUCCESS)\n	  {\n	    print \"------- Susscess\
ful upload of $name\";\n	    `mv x $name`;\n	    $\
file=$name;\n	  }\n      }\n    open ($F, \"$file\\
");\n    while (<$F>)\n      {\n	my $l=$_;\n	if ((\
$l =~/^\\/\\//) || ($db=~/^#/)){;}\n	elsif ( !($l \
=~/\\w/)){;}\n	else\n	  {\n	    my @v=split (/\\s+\
/, $l);\n	    if ( $l=~/^MODE/)\n	      {\n		$MODE\
{$v[1]}{$v[2]}=$v[3];\n	      }\n	    elsif ($l=~/\
^PG/)\n	      {\n		$PG{$v[1]}{$v[2]}=$v[3];\n	    \
  }\n	  }\n      }\n    close ($F);\n    &post_pro\
cess_PG();\n    return;\n  }\n\n\n\nsub initialize\
_PG\n  {\n\n$PG{\"t_coffee\"}{\"4_TCOFFEE\"}=\"TCO\
FFEE\";\n$PG{\"t_coffee\"}{\"type\"}=\"sequence_mu\
ltiple_aligner\";\n$PG{\"t_coffee\"}{\"ADDRESS\"}=\
\"http://www.tcoffee.org\";\n$PG{\"t_coffee\"}{\"l\
anguage\"}=\"C++\";\n$PG{\"t_coffee\"}{\"language2\
\"}=\"CXX\";\n$PG{\"t_coffee\"}{\"source\"}=\"http\
://www.tcoffee.org/Packages/Stable/Latest/T-COFFEE\
_distribution.tar.gz\";\n$PG{\"t_coffee\"}{\"updat\
e_action\"}=\"always\";\n$PG{\"t_coffee\"}{\"mode\\
"}=\"tcoffee,mcoffee,rcoffee,expresso,3dcoffee\";\\
n$PG{\"clustalo\"}{\"4_TCOFFEE\"}=\"CLUSTALO\";\n$\
PG{\"clustalo\"}{\"type\"}=\"sequence_multiple_ali\
gner\";\n$PG{\"clustalo\"}{\"ADDRESS\"}=\"http://w\
ww.clustal.org/omega/\";\n$PG{\"clustalo\"}{\"lang\
uage\"}=\"C++\";\n$PG{\"clustalo\"}{\"language2\"}\
=\"C++\";\n$PG{\"clustalo\"}{\"source\"}=\"http://\
www.clustal.org/omega/clustal-omega-1.2.4.tar.gz\"\
;\n$PG{\"clustalo\"}{\"mode\"}=\"mcoffee\";\n$PG{\\
"clustalo\"}{\"binary\"}=\"clustalo\";\n$PG{\"clus\
talo\"}{\"version\"}=\"1.2.4\";\n$PG{\"strike\"}{\\
"4_TCOFFEE\"}=\"STRIKE\";\n$PG{\"strike\"}{\"type\\
"}=\"sequence_alignment_scoring\";\n$PG{\"strike\"\
}{\"ADDRESS\"}=\"http://www.tcoffee.org/Projects/s\
trike/index.html\";\n$PG{\"strike\"}{\"language\"}\
=\"C++\";\n$PG{\"strike\"}{\"language2\"}=\"CXX\";\
\n$PG{\"strike\"}{\"source\"}=\"http://www.tcoffee\
.org/Projects/strike/strike_v1.2.tar.bz2\";\n$PG{\\
"strike\"}{\"mode\"}=\"tcoffee,expresso\";\n$PG{\"\
strike\"}{\"version\"}=\"1.2\";\n$PG{\"strike\"}{\\
"binary\"}=\"strike\";\n$PG{\"clustalw2\"}{\"4_TCO\
FFEE\"}=\"CLUSTALW2\";\n$PG{\"clustalw2\"}{\"type\\
"}=\"sequence_multiple_aligner\";\n$PG{\"clustalw2\
\"}{\"ADDRESS\"}=\"http://www.clustal.org\";\n$PG{\
\"clustalw2\"}{\"language\"}=\"C++\";\n$PG{\"clust\
alw2\"}{\"language2\"}=\"CXX\";\n$PG{\"clustalw2\"\
}{\"source\"}=\"http://www.clustal.org/download/2.\
0.10/clustalw-2.0.10-src.tar.gz\";\n$PG{\"clustalw\
2\"}{\"mode\"}=\"mcoffee,rcoffee\";\n$PG{\"clustal\
w2\"}{\"binary\"}=\"clustalw2\";\n$PG{\"clustalw2\\
"}{\"version\"}=\"2.0.10\";\n$PG{\"clustalw\"}{\"4\
_TCOFFEE\"}=\"CLUSTALW\";\n$PG{\"clustalw\"}{\"typ\
e\"}=\"sequence_multiple_aligner\";\n$PG{\"clustal\
w\"}{\"ADDRESS\"}=\"http://www.clustal.org\";\n$PG\
{\"clustalw\"}{\"language\"}=\"C\";\n$PG{\"clustal\
w\"}{\"language2\"}=\"C\";\n$PG{\"clustalw\"}{\"so\
urce\"}=\"http://www.clustal.org/download/1.X/ftp-\
igbmc.u-strasbg.fr/pub/ClustalW/clustalw1.82.UNIX.\
tar.gz\";\n$PG{\"clustalw\"}{\"mode\"}=\"mcoffee,r\
coffee\";\n$PG{\"clustalw\"}{\"version\"}=\"1.82\"\
;\n$PG{\"clustalw\"}{\"binary\"}=\"clustalw\";\n$P\
G{\"dialign-t\"}{\"4_TCOFFEE\"}=\"DIALIGNT\";\n$PG\
{\"dialign-t\"}{\"type\"}=\"sequence_multiple_alig\
ner\";\n$PG{\"dialign-t\"}{\"ADDRESS\"}=\"http://d\
ialign-tx.gobics.de/\";\n$PG{\"dialign-t\"}{\"DIR\\
"}=\"/usr/share/dialign-tx/\";\n$PG{\"dialign-t\"}\
{\"language\"}=\"C\";\n$PG{\"dialign-t\"}{\"langua\
ge2\"}=\"C\";\n$PG{\"dialign-t\"}{\"source\"}=\"ht\
tp://dialign-tx.gobics.de/DIALIGN-TX_1.0.2.tar.gz\\
";\n$PG{\"dialign-t\"}{\"mode\"}=\"mcoffee\";\n$PG\
{\"dialign-t\"}{\"binary\"}=\"dialign-t\";\n$PG{\"\
dialign-t\"}{\"version\"}=\"1.0.2\";\n$PG{\"dialig\
n-tx\"}{\"4_TCOFFEE\"}=\"DIALIGNTX\";\n$PG{\"diali\
gn-tx\"}{\"type\"}=\"sequence_multiple_aligner\";\\
n$PG{\"dialign-tx\"}{\"ADDRESS\"}=\"http://dialign\
-tx.gobics.de/\";\n$PG{\"dialign-tx\"}{\"DIR\"}=\"\
/usr/share/dialign-tx/\";\n$PG{\"dialign-tx\"}{\"l\
anguage\"}=\"C\";\n$PG{\"dialign-tx\"}{\"language2\
\"}=\"C\";\n$PG{\"dialign-tx\"}{\"source\"}=\"http\
://dialign-tx.gobics.de/DIALIGN-TX_1.0.2.tar.gz\";\
\n$PG{\"dialign-tx\"}{\"mode\"}=\"mcoffee\";\n$PG{\
\"dialign-tx\"}{\"binary\"}=\"dialign-tx\";\n$PG{\\
"dialign-tx\"}{\"version\"}=\"1.0.2\";\n$PG{\"poa\\
"}{\"4_TCOFFEE\"}=\"POA\";\n$PG{\"poa\"}{\"type\"}\
=\"sequence_multiple_aligner\";\n$PG{\"poa\"}{\"AD\
DRESS\"}=\"http://www.bioinformatics.ucla.edu/poa/\
\";\n$PG{\"poa\"}{\"language\"}=\"C\";\n$PG{\"poa\\
"}{\"language2\"}=\"C\";\n$PG{\"poa\"}{\"source\"}\
=\"http://downloads.sourceforge.net/poamsa/poaV2.t\
ar.gz\";\n$PG{\"poa\"}{\"DIR\"}=\"/usr/share/\";\n\
$PG{\"poa\"}{\"FILE1\"}=\"blosum80.mat\";\n$PG{\"p\
oa\"}{\"mode\"}=\"mcoffee\";\n$PG{\"poa\"}{\"binar\
y\"}=\"poa\";\n$PG{\"poa\"}{\"version\"}=\"2.0\";\\
n$PG{\"probcons\"}{\"4_TCOFFEE\"}=\"PROBCONS\";\n$\
PG{\"probcons\"}{\"type\"}=\"sequence_multiple_ali\
gner\";\n$PG{\"probcons\"}{\"ADDRESS\"}=\"http://p\
robcons.stanford.edu/\";\n$PG{\"probcons\"}{\"lang\
uage2\"}=\"CXX\";\n$PG{\"probcons\"}{\"language\"}\
=\"C++\";\n$PG{\"probcons\"}{\"source\"}=\"http://\
probcons.stanford.edu/probcons_v1_12.tar.gz\";\n$P\
G{\"probcons\"}{\"mode\"}=\"mcoffee\";\n$PG{\"prob\
cons\"}{\"binary\"}=\"probcons\";\n$PG{\"probcons\\
"}{\"version\"}=\"1.12\";\n$PG{\"msaprobs\"}{\"4_T\
COFFEE\"}=\"MSAPROBS\";\n$PG{\"msaprobs\"}{\"type\\
"}=\"sequence_multiple_aligner\";\n$PG{\"msaprobs\\
"}{\"ADDRESS\"}=\"http://msaprobs.sourceforge.net/\
homepage.htm#latest\";\n$PG{\"msaprobs\"}{\"langua\
ge2\"}=\"CXX\";\n$PG{\"msaprobs\"}{\"language\"}=\\
"C++\";\n$PG{\"msaprobs\"}{\"source\"}=\"https://s\
ourceforge.net/projects/msaprobs/files/MSAProbs-MP\
I/MSAProbs-MPI_rel1.0.5.tar.gz\";\n$PG{\"msaprobs\\
"}{\"mode\"}=\"mcoffee\";\n$PG{\"msaprobs\"}{\"bin\
ary\"}=\"msaprobs\";\n$PG{\"msaprobs\"}{\"version\\
"}=\"1.05\";\n$PG{\"msaprobs\"}{\"update_action\"}\
=\"never\";\n$PG{\"upp\"}{\"4_TCOFFEE\"}=\"UPP\";\\
n$PG{\"upp\"}{\"type\"}=\"sequence_multiple_aligne\
r\";\n$PG{\"upp\"}{\"ADDRESS\"}=\"http://www.cs.ut\
exas.edu/users/phylo/software/upp/\";\n$PG{\"upp\"\
}{\"language2\"}=\"CXX\";\n$PG{\"upp\"}{\"language\
\"}=\"C++\";\n$PG{\"upp\"}{\"source\"}=\"https://g\
ithub.com/smirarab/pasta/archive/upp.zip\";\n$PG{\\
"upp\"}{\"mode\"}=\"mcoffee\";\n$PG{\"upp\"}{\"bin\
ary\"}=\"upp\";\n$PG{\"upp\"}{\"version\"}=\"1\";\\
n$PG{\"upp\"}{\"update_action\"}=\"never\";\n$PG{\\
"mafft\"}{\"4_TCOFFEE\"}=\"MAFFT\";\n$PG{\"mafft\"\
}{\"type\"}=\"sequence_multiple_aligner\";\n$PG{\"\
mafft\"}{\"ADDRESS\"}=\"http://align.bmr.kyushu-u.\
ac.jp/mafft/online/server/\";\n$PG{\"mafft\"}{\"la\
nguage\"}=\"C\";\n$PG{\"mafft\"}{\"language\"}=\"C\
\";\n$PG{\"mafft\"}{\"source\"}=\"http://mafft.cbr\
c.jp/alignment/software/mafft-7.310-with-extension\
s-src.tgz\";\n$PG{\"mafft\"}{\"mode\"}=\"mcoffee,r\
coffee\";\n$PG{\"mafft\"}{\"binary\"}=\"mafft.tar.\
gz\";\n$PG{\"mafft\"}{\"version\"}=\"7.310\";\n$PG\
{\"msa\"}{\"4_TCOFFEE\"}=\"MSA\";\n$PG{\"msa\"}{\"\
type\"}=\"sequence_multiple_aligner\";\n$PG{\"msa\\
"}{\"ADDRESS\"}=\"https://www.ncbi.nlm.nih.gov/CBB\
research/Schaffer/msa.html\";\n$PG{\"msa\"}{\"lang\
uage\"}=\"C\";\n$PG{\"msa\"}{\"language\"}=\"C\";\\
n$PG{\"msa\"}{\"source\"}=\"ftp://ftp.ncbi.nih.gov\
/pub/msa/msa.tar.Z\";\n$PG{\"msa\"}{\"mode\"}=\"mc\
offee\";\n$PG{\"msa\"}{\"binary\"}=\"msa.pl\";\n$P\
G{\"msa\"}{\"version\"}=\"1.0\";\n$PG{\"msa\"}{\"u\
pdate_action\"}=\"never\";\n$PG{\"dca\"}{\"4_TCOFF\
EE\"}=\"DCA\";\n$PG{\"dca\"}{\"type\"}=\"sequence_\
multiple_aligner\";\n$PG{\"dca\"}{\"ADDRESS\"}=\"h\
ttps://bibiserv2.cebitec.uni-bielefeld.de/dca\";\n\
$PG{\"dca\"}{\"language\"}=\"C\";\n$PG{\"dca\"}{\"\
language\"}=\"C\";\n$PG{\"dca\"}{\"source\"}=\"htt\
ps://bibiserv2.cebitec.uni-bielefeld.de/applicatio\
ns/dca/resources/downloads/dca-1.1-src.tar.gz\";\n\
$PG{\"dca\"}{\"mode\"}=\"mcoffee\";\n$PG{\"dca\"}{\
\"binary\"}=\"dca.pl\";\n$PG{\"dca\"}{\"version\"}\
=\"1.1\";\n$PG{\"dca\"}{\"update_action\"}=\"never\
\";\n$PG{\"muscle\"}{\"4_TCOFFEE\"}=\"MUSCLE\";\n$\
PG{\"muscle\"}{\"type\"}=\"sequence_multiple_align\
er\";\n$PG{\"muscle\"}{\"ADDRESS\"}=\"http://www.d\
rive5.com/muscle/\";\n$PG{\"muscle\"}{\"language\"\
}=\"C++\";\n$PG{\"muscle\"}{\"language2\"}=\"GPP\"\
;\n$PG{\"muscle\"}{\"source\"}=\"http://www.drive5\
.com/muscle/downloads3.7/muscle3.7_src.tar.gz\";\n\
$PG{\"muscle\"}{\"windows\"}=\"http://www.drive5.c\
om/muscle/downloads3.7/muscle3.7_win32.zip\";\n$PG\
{\"muscle\"}{\"linux\"}=\"http://www.drive5.com/mu\
scle/downloads3.7/muscle3.7_linux_ia32.tar.gz\";\n\
$PG{\"muscle\"}{\"mode\"}=\"mcoffee,rcoffee\";\n$P\
G{\"muscle\"}{\"version\"}=\"3.7\";\n$PG{\"pcma\"}\
{\"4_TCOFFEE\"}=\"PCMA\";\n$PG{\"pcma\"}{\"type\"}\
=\"sequence_multiple_aligner\";\n$PG{\"pcma\"}{\"A\
DDRESS\"}=\"http://prodata.swmed.edu/pcma/pcma.php\
\";\n$PG{\"pcma\"}{\"language\"}=\"C\";\n$PG{\"pcm\
a\"}{\"language2\"}=\"C\";\n$PG{\"pcma\"}{\"source\
\"}=\"http://prodata.swmed.edu/download/pub/PCMA/p\
cma.tar.gz\";\n$PG{\"pcma\"}{\"mode\"}=\"mcoffee\"\
;\n$PG{\"pcma\"}{\"version\"}=\"1.0\";\n$PG{\"kali\
gn\"}{\"4_TCOFFEE\"}=\"KALIGN\";\n$PG{\"kalign\"}{\
\"type\"}=\"sequence_multiple_aligner\";\n$PG{\"ka\
lign\"}{\"ADDRESS\"}=\"http://msa.cgb.ki.se\";\n$P\
G{\"kalign\"}{\"language\"}=\"C\";\n$PG{\"kalign\"\
}{\"language2\"}=\"C\";\n$PG{\"kalign\"}{\"source\\
"}=\"http://msa.cgb.ki.se/downloads/kalign/current\
.tar.gz\";\n$PG{\"kalign\"}{\"mode\"}=\"mcoffee\";\
\n$PG{\"kalign\"}{\"version\"}=\"1.0\";\n$PG{\"ama\
p\"}{\"4_TCOFFEE\"}=\"AMAP\";\n$PG{\"amap\"}{\"typ\
e\"}=\"sequence_multiple_aligner\";\n$PG{\"amap\"}\
{\"ADDRESS\"}=\"http://bio.math.berkeley.edu/amap/\
\";\n$PG{\"amap\"}{\"language\"}=\"C++\";\n$PG{\"a\
map\"}{\"language2\"}=\"CXX\";\n$PG{\"amap\"}{\"so\
urce\"}=\"https://github.com/mes5k/amap-align/arch\
ive/amap.zip\";\n$PG{\"amap\"}{\"mode\"}=\"mcoffee\
\";\n$PG{\"amap\"}{\"version\"}=\"2.0\";\n$PG{\"am\
ap\"}{\"update_action\"}=\"never\";\n$PG{\"proda\"\
}{\"4_TCOFFEE\"}=\"PRODA\";\n$PG{\"proda\"}{\"type\
\"}=\"sequence_multiple_aligner\";\n$PG{\"proda\"}\
{\"ADDRESS\"}=\"http://proda.stanford.edu\";\n$PG{\
\"proda\"}{\"language\"}=\"C++\";\n$PG{\"proda\"}{\
\"language2\"}=\"CXX\";\n$PG{\"proda\"}{\"source\"\
}=\"http://proda.stanford.edu/proda_1_0.tar.gz\";\\
n$PG{\"proda\"}{\"mode\"}=\"mcoffee\";\n$PG{\"prod\
a\"}{\"version\"}=\"1.0\";\n$PG{\"prank\"}{\"4_TCO\
FFEE\"}=\"PRANK\";\n$PG{\"prank\"}{\"type\"}=\"seq\
uence_multiple_aligner\";\n$PG{\"prank\"}{\"ADDRES\
S\"}=\"http://www.ebi.ac.uk/goldman-srv/prank/\";\\
n$PG{\"prank\"}{\"language\"}=\"C++\";\n$PG{\"pran\
k\"}{\"language2\"}=\"CXX\";\n$PG{\"prank\"}{\"sou\
rce\"}=\"http://www.ebi.ac.uk/goldman-srv/prank/sr\
c/prank/prank.src.100802.tgz\";\n$PG{\"prank\"}{\"\
mode\"}=\"mcoffee\";\n$PG{\"prank\"}{\"version\"}=\
\"100303\";\n$PG{\"sap\"}{\"4_TCOFFEE\"}=\"SAP\";\\
n$PG{\"sap\"}{\"type\"}=\"structure_pairwise_align\
er\";\n$PG{\"sap\"}{\"ADDRESS\"}=\"https://mathbio\
.crick.ac.uk/wiki/Software#SAP\";\n$PG{\"sap\"}{\"\
language\"}=\"C\";\n$PG{\"sap\"}{\"language2\"}=\"\
C\";\n$PG{\"sap\"}{\"source\"}=\"https://github.co\
m/jkleinj/SAP/archive/v.1.1.3.tar.gz\";\n$PG{\"sap\
\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\"sap\"}\
{\"version\"}=\"1.1.3\";\n$PG{\"sap\"}{\"binary\"}\
=\"sap\";\n$PG{\"TMalign\"}{\"4_TCOFFEE\"}=\"TMALI\
GN\";\n$PG{\"TMalign\"}{\"type\"}=\"structure_pair\
wise_aligner\";\n$PG{\"TMalign\"}{\"ADDRESS\"}=\"h\
ttp://zhanglab.ccmb.med.umich.edu/TM-align/TMalign\
.f\";\n$PG{\"TMalign\"}{\"language\"}=\"Fortran\";\
\n$PG{\"TMalign\"}{\"language2\"}=\"Fortran\";\n$P\
G{\"TMalign\"}{\"source\"}=\"http://zhanglab.ccmb.\
med.umich.edu/TM-align/TMalign.f\";\n$PG{\"TMalign\
\"}{\"linux\"}=\"http://zhanglab.ccmb.med.umich.ed\
u/TM-align/TMalign_32.gz\";\n$PG{\"TMalign\"}{\"mo\
de\"}=\"expresso,3dcoffee\";\n$PG{\"TMalign\"}{\"v\
ersion\"}=\"2013.05.11\";\n$PG{\"mustang\"}{\"4_TC\
OFFEE\"}=\"MUSTANG\";\n$PG{\"mustang\"}{\"type\"}=\
\"structure_pairwise_aligner\";\n$PG{\"mustang\"}{\
\"ADDRESS\"}=\"http://lcb.infotech.monash.edu.au/m\
ustang/\";\n$PG{\"mustang\"}{\"language\"}=\"C++\"\
;\n$PG{\"mustang\"}{\"language2\"}=\"CXX\";\n$PG{\\
"mustang\"}{\"source\"}=\"http://lcb.infotech.mona\
sh.edu.au/mustang/mustang_v3.2.3.tgz\";\n$PG{\"mus\
tang\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\"mu\
stang\"}{\"version\"}=\"3.2.3\";\n$PG{\"lsqman\"}{\
\"4_TCOFFEE\"}=\"LSQMAN\";\n$PG{\"lsqman\"}{\"type\
\"}=\"structure_pairwise_aligner\";\n$PG{\"lsqman\\
"}{\"ADDRESS\"}=\"empty\";\n$PG{\"lsqman\"}{\"lang\
uage\"}=\"empty\";\n$PG{\"lsqman\"}{\"language2\"}\
=\"empty\";\n$PG{\"lsqman\"}{\"source\"}=\"empty\"\
;\n$PG{\"lsqman\"}{\"update_action\"}=\"never\";\n\
$PG{\"lsqman\"}{\"mode\"}=\"expresso,3dcoffee\";\n\
$PG{\"align_pdb\"}{\"4_TCOFFEE\"}=\"ALIGN_PDB\";\n\
$PG{\"align_pdb\"}{\"type\"}=\"structure_pairwise_\
aligner\";\n$PG{\"align_pdb\"}{\"ADDRESS\"}=\"empt\
y\";\n$PG{\"align_pdb\"}{\"language\"}=\"empty\";\\
n$PG{\"align_pdb\"}{\"language2\"}=\"empty\";\n$PG\
{\"align_pdb\"}{\"source\"}=\"empty\";\n$PG{\"alig\
n_pdb\"}{\"update_action\"}=\"never\";\n$PG{\"alig\
n_pdb\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\"f\
ugueali\"}{\"4_TCOFFEE\"}=\"FUGUE\";\n$PG{\"fuguea\
li\"}{\"type\"}=\"structure_pairwise_aligner\";\n$\
PG{\"fugueali\"}{\"ADDRESS\"}=\"http://mizuguchila\
b.org/fugue/\";\n$PG{\"fugueali\"}{\"language\"}=\\
"empty\";\n$PG{\"fugueali\"}{\"language2\"}=\"empt\
y\";\n$PG{\"fugueali\"}{\"source\"}=\"empty\";\n$P\
G{\"fugueali\"}{\"update_action\"}=\"never\";\n$PG\
{\"fugueali\"}{\"mode\"}=\"expresso,3dcoffee\";\n$\
PG{\"dalilite.pl\"}{\"4_TCOFFEE\"}=\"DALILITEc\";\\
n$PG{\"dalilite.pl\"}{\"type\"}=\"structure_pairwi\
se_aligner\";\n$PG{\"dalilite.pl\"}{\"ADDRESS\"}=\\
"built_in\";\n$PG{\"dalilite.pl\"}{\"ADDRESS2\"}=\\
"http://www.ebi.ac.uk/Tools/webservices/services/d\
alilite\";\n$PG{\"dalilite.pl\"}{\"language\"}=\"P\
erl\";\n$PG{\"dalilite.pl\"}{\"language2\"}=\"Perl\
\";\n$PG{\"dalilite.pl\"}{\"source\"}=\"empty\";\n\
$PG{\"dalilite.pl\"}{\"update_action\"}=\"never\";\
\n$PG{\"dalilite.pl\"}{\"mode\"}=\"expresso,3dcoff\
ee\";\n$PG{\"probconsRNA\"}{\"4_TCOFFEE\"}=\"PROBC\
ONSRNA\";\n$PG{\"probconsRNA\"}{\"type\"}=\"RNA_mu\
ltiple_aligner\";\n$PG{\"probconsRNA\"}{\"ADDRESS\\
"}=\"http://probcons.stanford.edu/\";\n$PG{\"probc\
onsRNA\"}{\"language\"}=\"C++\";\n$PG{\"probconsRN\
A\"}{\"language2\"}=\"CXX\";\n$PG{\"probconsRNA\"}\
{\"source\"}=\"http://probcons.stanford.edu/probco\
nsRNA.tar.gz\";\n$PG{\"probconsRNA\"}{\"mode\"}=\"\
mcoffee,rcoffee\";\n$PG{\"probconsRNA\"}{\"version\
\"}=\"1.0\";\n$PG{\"sfold\"}{\"4_TCOFFEE\"}=\"CONS\
AN\";\n$PG{\"sfold\"}{\"type\"}=\"RNA_pairwise_ali\
gner\";\n$PG{\"sfold\"}{\"ADDRESS\"}=\"http://sela\
b.janelia.org/software/consan/\";\n$PG{\"sfold\"}{\
\"language\"}=\"empty\";\n$PG{\"sfold\"}{\"languag\
e2\"}=\"empty\";\n$PG{\"sfold\"}{\"source\"}=\"emp\
ty\";\n$PG{\"sfold\"}{\"update_action\"}=\"never\"\
;\n$PG{\"sfold\"}{\"mode\"}=\"rcoffee\";\n$PG{\"RN\
Aplfold\"}{\"4_TCOFFEE\"}=\"RNAPLFOLD\";\n$PG{\"RN\
Aplfold\"}{\"type\"}=\"RNA_secondarystructure_pred\
ictor\";\n$PG{\"RNAplfold\"}{\"ADDRESS\"}=\"http:/\
/www.tbi.univie.ac.at/RNA/\";\n$PG{\"RNAplfold\"}{\
\"language\"}=\"C\";\n$PG{\"RNAplfold\"}{\"languag\
e2\"}=\"C\";\n$PG{\"RNAplfold\"}{\"source\"}=\"htt\
p://www.tbi.univie.ac.at/RNA/packages/source/Vienn\
aRNA-2.1.9.tar.gz\";\n$PG{\"RNAplfold\"}{\"mode\"}\
=\"rcoffee,\";\n$PG{\"RNAplfold\"}{\"binary\"}=\"R\
NAplfold.tar.gz\";\n$PG{\"RNAplfold\"}{\"version\"\
}=\"2.1.9\";\n$PG{\"retree\"}{\"4_TCOFFEE\"}=\"PHY\
LIP\";\n$PG{\"retree\"}{\"type\"}=\"Phylogeny\";\n\
$PG{\"retree\"}{\"ADDRESS\"}=\"http://evolution.gs\
.washington.edu/phylip/\";\n$PG{\"retree\"}{\"lang\
uage\"}=\"C\";\n$PG{\"retree\"}{\"language2\"}=\"C\
\";\n$PG{\"retree\"}{\"source\"}=\"http://www.tcof\
fee.org/Packages/mirrors/source/phylip-3.66.tar.gz\
\";\n$PG{\"retree\"}{\"mode\"}=\"trmsd,\";\n$PG{\"\
retree\"}{\"binary\"}=\"retree.tar.gz\";\n$PG{\"re\
tree\"}{\"version\"}=\"3.66\";\n$PG{\"hmmtop\"}{\"\
4_TCOFFEE\"}=\"HMMTOP\";\n$PG{\"hmmtop\"}{\"type\"\
}=\"protein_secondarystructure_predictor\";\n$PG{\\
"hmmtop\"}{\"ADDRESS\"}=\"www.enzim.hu/hmmtop/\";\\
n$PG{\"hmmtop\"}{\"language\"}=\"C\";\n$PG{\"hmmto\
p\"}{\"language2\"}=\"C\";\n$PG{\"hmmtop\"}{\"sour\
ce\"}=\"http://www.tcoffee.org/Packages/mirrors/hm\
mtop2.1.tgz\";\n$PG{\"hmmtop\"}{\"binary\"}=\"hmmt\
op\";\n$PG{\"hmmtop\"}{\"update_action\"}=\"never\\
";\n$PG{\"hmmtop\"}{\"mode\"}=\"tcoffee\";\n$PG{\"\
hmmtop\"}{\"version\"}=\"2.1\";\n$PG{\"gorIV\"}{\"\
4_TCOFFEE\"}=\"GOR4\";\n$PG{\"gorIV\"}{\"type\"}=\\
"protein_secondarystructure_predictor\";\n$PG{\"go\
rIV\"}{\"ADDRESS\"}=\"http://mig.jouy.inra.fr/logi\
ciels/gorIV/\";\n$PG{\"gorIV\"}{\"language\"}=\"C\\
";\n$PG{\"gorIV\"}{\"language2\"}=\"C\";\n$PG{\"go\
rIV\"}{\"source\"}=\"http://www.tcoffee.org/Packag\
es/mirrors/GOR_IV.tar.gz\";\n$PG{\"gorIV\"}{\"upda\
te_action\"}=\"never\";\n$PG{\"gorIV\"}{\"mode\"}=\
\"tcoffee\";\n$PG{\"wublast.pl\"}{\"4_TCOFFEE\"}=\\
"EBIWUBLASTc\";\n$PG{\"wublast.pl\"}{\"type\"}=\"p\
rotein_homology_predictor\";\n$PG{\"wublast.pl\"}{\
\"ADDRESS\"}=\"built_in\";\n$PG{\"wublast.pl\"}{\"\
ADDRESS2\"}=\"http://www.ebi.ac.uk/Tools/webservic\
es/services/wublast\";\n$PG{\"wublast.pl\"}{\"lang\
uage\"}=\"Perl\";\n$PG{\"wublast.pl\"}{\"language2\
\"}=\"Perl\";\n$PG{\"wublast.pl\"}{\"source\"}=\"e\
mpty\";\n$PG{\"wublast.pl\"}{\"update_action\"}=\"\
never\";\n$PG{\"wublast.pl\"}{\"mode\"}=\"psicoffe\
e,expresso,accurate\";\n$PG{\"blastpgp.pl\"}{\"4_T\
COFFEE\"}=\"EBIBLASTPGPc\";\n$PG{\"blastpgp.pl\"}{\
\"type\"}=\"protein_homology_predictor\";\n$PG{\"b\
lastpgp.pl\"}{\"ADDRESS\"}=\"built_in\";\n$PG{\"bl\
astpgp.pl\"}{\"ADDRESS2\"}=\"http://www.ebi.ac.uk/\
Tools/webservices/services/blastpgp\";\n$PG{\"blas\
tpgp.pl\"}{\"language\"}=\"Perl\";\n$PG{\"blastpgp\
.pl\"}{\"language2\"}=\"Perl\";\n$PG{\"blastpgp.pl\
\"}{\"source\"}=\"empty\";\n$PG{\"blastpgp.pl\"}{\\
"update_action\"}=\"never\";\n$PG{\"blastpgp.pl\"}\
{\"mode\"}=\"psicoffee,expresso,accurate\";\n$PG{\\
"blastall\"}{\"4_TCOFFEE\"}=\"blastall\";\n$PG{\"b\
lastall\"}{\"type\"}=\"protein_homology_predictor\\
";\n$PG{\"blastall\"}{\"ADDRESS\"}=\"ftp://ftp.ncb\
i.nih.gov/blast/executables/LATEST\";\n$PG{\"blast\
all\"}{\"language\"}=\"C\";\n$PG{\"blastall\"}{\"l\
anguage2\"}=\"C\";\n$PG{\"blastall\"}{\"source\"}=\
\"ftp://ftp.ncbi.nlm.nih.gov/blast/executables/bla\
st+/2.6.0/ncbi-blast-2.6.0+-src.tar.gz\";\n$PG{\"b\
lastall\"}{\"update_action\"}=\"never\";\n$PG{\"bl\
astall\"}{\"mode\"}=\"psicoffee,expresso,3dcoffee\\
";\n$PG{\"legacy_blast.pl\"}{\"4_TCOFFEE\"}=\"NCBI\
BLAST\";\n$PG{\"legacy_blast.pl\"}{\"type\"}=\"pro\
tein_homology_predictor\";\n$PG{\"legacy_blast.pl\\
"}{\"ADDRESS\"}=\"ftp://ftp.ncbi.nih.gov/blast/exe\
cutables/LATEST\";\n$PG{\"legacy_blast.pl\"}{\"lan\
guage\"}=\"C\";\n$PG{\"legacy_blast.pl\"}{\"langua\
ge2\"}=\"C\";\n$PG{\"legacy_blast.pl\"}{\"source\"\
}=\"ftp://ftp.ncbi.nlm.nih.gov/blast/executables/b\
last+/2.6.0/ncbi-blast-2.6.0+-src.tar.gz\";\n$PG{\\
"legacy_blast.pl\"}{\"update_action\"}=\"never\";\\
n$PG{\"legacy_blast.pl\"}{\"mode\"}=\"psicoffee,ex\
presso,3dcoffee\";\n$PG{\"SOAP::Lite\"}{\"4_TCOFFE\
E\"}=\"SOAPLITE\";\n$PG{\"SOAP::Lite\"}{\"type\"}=\
\"library\";\n$PG{\"SOAP::Lite\"}{\"ADDRESS\"}=\"h\
ttp://cpansearch.perl.org/src/MKUTTER/SOAP-Lite-0.\
710.08/Makefile.PL\";\n$PG{\"SOAP::Lite\"}{\"langu\
age\"}=\"Perl\";\n$PG{\"SOAP::Lite\"}{\"language2\\
"}=\"Perl\";\n$PG{\"SOAP::Lite\"}{\"source\"}=\"em\
pty\";\n$PG{\"SOAP::Lite\"}{\"update_action\"}=\"n\
ever\";\n$PG{\"SOAP::Lite\"}{\"mode\"}=\"none\";\n\
$PG{\"XML::Simple\"}{\"4_TCOFFEE\"}=\"XMLSIMPLE\";\
\n$PG{\"XML::Simple\"}{\"type\"}=\"library\";\n$PG\
{\"XML::Simple\"}{\"ADDRESS\"}=\"http://search.cpa\
n.org/~grantm/XML-Simple-2.18/lib/XML/Simple.pm\";\
\n$PG{\"XML::Simple\"}{\"language\"}=\"Perl\";\n$P\
G{\"XML::Simple\"}{\"language2\"}=\"Perl\";\n$PG{\\
"XML::Simple\"}{\"source\"}=\"empty\";\n$PG{\"XML:\
:Simple\"}{\"mode\"}=\"psicoffee,expresso,accurate\
\";\n$PG{\"x3dna\"}{\"4_TCOFFEE\"}=\"x3dna\";\n$PG\
{\"x3dna\"}{\"type\"}=\"RNA_secondarystructure_pre\
dictor\";\n$PG{\"x3dna\"}{\"ADDRESS\"}=\"http://x3\
dna.bio.columbia.edu/\";\n$PG{\"x3dna\"}{\"source\\
"}=\"http://www.tcoffee.org/Packages/mirrors/sourc\
e/x3dna-v2.3-linux-64bit.tar.gz\";\n$PG{\"x3dna\"}\
{\"mode\"}=\"saracoffee\";\n$PG{\"x3dna\"}{\"updat\
e_action\"}=\"never\";\n$PG{\"fsa\"}{\"4_TCOFFEE\"\
}=\"FSA\";\n$PG{\"fsa\"}{\"type\"}=\"sequence_mult\
iple_aligner\";\n$PG{\"fsa\"}{\"ADDRESS\"}=\"http:\
//fsa.sourceforge.net/\";\n$PG{\"fsa\"}{\"language\
\"}=\"C++\";\n$PG{\"fsa\"}{\"language2\"}=\"CXX\";\
\n$PG{\"fsa\"}{\"source\"}=\"http://sourceforge.ne\
t/projects/fsa/files/fsa-1.15.3.tar.gz/download/\"\
;\n$PG{\"fsa\"}{\"mode\"}=\"mcoffee\";\n$PG{\"fsa\\
"}{\"version\"}=\"1.15.3\";\n$PG{\"fsa\"}{\"update\
_action\"}=\"never\";\n$PG{\"mus4\"}{\"4_TCOFFEE\"\
}=\"MUS4\";\n$PG{\"mus4\"}{\"type\"}=\"sequence_mu\
ltiple_aligner\";\n$PG{\"mus4\"}{\"ADDRESS\"}=\"ht\
tp://www.drive5.com/muscle/\";\n$PG{\"mus4\"}{\"la\
nguage\"}=\"C++\";\n$PG{\"mus4\"}{\"language2\"}=\\
"GPP\";\n$PG{\"mus4\"}{\"source\"}=\"http://www.dr\
ive5.com/muscle/muscle4.0_src.tar.gz\";\n$PG{\"mus\
4\"}{\"mode\"}=\"mcoffee,rcoffee\";\n$PG{\"mus4\"}\
{\"version\"}=\"4.0\";\n$PG{\"mus4\"}{\"update_act\
ion\"}=\"never\";\n$MODE{\"tcoffee\"}{\"name\"}=\"\
tcoffee\";\n$MODE{\"rcoffee\"}{\"name\"}=\"rcoffee\
\";\n$MODE{\"3dcoffee\"}{\"name\"}=\"3dcoffee\";\n\
$MODE{\"mcoffee\"}{\"name\"}=\"mcoffee\";\n$MODE{\\
"expresso\"}{\"name\"}=\"expresso\";\n$MODE{\"trms\
d\"}{\"name\"}=\"trmsd\";\n$MODE{\"accurate\"}{\"n\
ame\"}=\"accurate\";\n$MODE{\"seq_reformat\"}{\"na\
me\"}=\"seq_reformat\";\n\n\n$PG{C}{compiler}=\"gc\
c\";\n$PG{C}{compiler_flag}=\"CC\";\n$PG{C}{option\
s}=\"\";\n$PG{C}{options_flag}=\"CFLAGS\";\n$PG{C}\
{type}=\"compiler\";\n\n$PG{\"CXX\"}{compiler}=\"g\
++\";\n$PG{\"CXX\"}{compiler_flag}=\"CXX\";\n$PG{\\
"CXX\"}{options}=\"\";\n$PG{\"CXX\"}{options_flag}\
=\"CXXFLAGS\";\n$PG{CXX}{type}=\"compiler\";\n\n$P\
G{\"CPP\"}{compiler}=\"g++\";\n$PG{\"CPP\"}{compil\
er_flag}=\"CPP\";\n$PG{\"CPP\"}{options}=\"\";\n$P\
G{\"CPP\"}{options_flag}=\"CPPFLAGS\";\n$PG{CPP}{t\
ype}=\"compiler\";\n\n$PG{\"GPP\"}{compiler}=\"g++\
\";\n$PG{\"GPP\"}{compiler_flag}=\"GPP\";\n$PG{\"G\
PP\"}{options}=\"\";\n$PG{\"GPP\"}{options_flag}=\\
"CFLAGS\";\n$PG{GPP}{type}=\"compiler\";\n\n$PG{Fo\
rtran}{compiler}=\"g77\";\n$PG{Fortran}{compiler_f\
lag}=\"FCC\";\n$PG{Fortran}{type}=\"compiler\";\n\\
n$PG{Perl}{compiler}=\"CPAN\";\n$PG{Perl}{type}=\"\
compiler\";\n\n$SUPPORTED_OS{macosx}=\"Macintosh\"\
;\n$SUPPORTED_OS{linux}=\"Linux\";\n$SUPPORTED_OS{\
windows}=\"Cygwin\";\n\n\n\n$MODE{t_coffee}{descri\
ption}=\" for regular multiple sequence alignments\
\";\n$MODE{rcoffee} {description}=\" for RNA multi\
ple sequence alignments\";\n\n$MODE{psicoffee} {de\
scription}=\" for Homology Extended multiple seque\
nce alignments\";\n$MODE{expresso}{description}=\"\
 for very accurate structure based multiple sequen\
ce alignments\";\n$MODE{\"3dcoffee\"}{description}\
=\" for multiple structure alignments\";\n$MODE{mc\
offee} {description}=\" for combining alternative \
multiple sequence alignment packages\\n------- int\
o a unique meta-package. The installer will upload\
 several MSA packages and compile them\\n\n\";\n\n\
\n&post_process_PG();\nreturn;\n}\n\nsub post_proc\
ess_PG\n  {\n    my $p;\n    \n    %PG=&name2dname\
 (%PG);\n    %MODE=&name2dname(%MODE);\n    foreac\
h $p (keys(%PG)){if ( $PG{$p}{type} eq \"compiler\\
"){$PG{$p}{update_action}=\"never\";}}\n    \n  }\\
n\nsub name2dname\n  {\n    my (%L)=(@_);\n    my \
($l, $ml);\n    \n    foreach my $pg (keys(%L))\n \
     {\n	$l=length ($pg);\n	if ( $l>$ml){$ml=$l;}\\
n      }\n    $ml+=1;\n    foreach my $pg (keys(%L\
))\n      {\n	my $name;\n	$l=$ml-length ($pg);\n	$\
name=$pg;\n	for ( $b=0; $b<$l; $b++)\n	  {\n	    $\
name .=\" \";\n	  }\n	$L{$pg}{dname}=$name;\n     \
 }\n    return %L;\n  }\n\nsub env_file2putenv\n  \
{\n    my $f=@_[0];\n    my $F=new FileHandle;\n  \
  my $n;\n    \n    open ($F, \"$f\");\n    while \
(<$F>)\n      {\n	my $line=$_;\n	my($var, $value)=\
($_=~/(\\S+)\\=(\\S*)/);\n	$ENV{$var}=$value;\n	$E\
NV_SET{$var}=1;\n	$n++;\n      }\n    close ($F);\\
n    return $n;\n  }\n\nsub replace_line_in_file\n\
  {\n    my ($file, $wordin, $wordout)=@_;\n    my\
 $O=new FileHandle;\n    my $I=new FileHandle;\n  \
  my $l;\n    if (!-e $file){return;}\n    \n    s\
ystem (\"mv $file $file.old\");\n    open ($O, \">\
$file\");\n    open ($I, \"$file.old\");\n    whil\
e (<$I>)\n      {\n	$l=$_;\n	if (!($l=~/$wordin/))\
{print $O \"$l\";}\n	elsif ( $wordout ne \"\"){$l=\
~s/$wordin/$wordout/g;print $O \"$l\";}\n      }\n\
    close ($O);\n    close ($I);\n    return;\n  }\
\n\nsub add_C_libraries\n  {\n   my ($file,$first,\
@list)=@_;\n   \n    my $O=new FileHandle;\n    my\
 $I=new FileHandle;\n    my ($l,$anchor);\n    if \
(!-e $file){return;}\n   \n    $anchor=\"#include \
<$first>\";\n	 \n    system (\"mv $file $file.old\\
");\n    open ($O, \">$file\");\n    open ($I, \"$\
file.old\");\n    while (<$I>)\n      {\n	$l=$_;\n\
	print $O \"$l\";\n	if (!($l=~/$anchor/))\n	   {\n\
	    \n	    foreach my $lib (@list)\n	       {\n  \
                print $O \"#include <$lib>\\n\";\n\
	       }\n           }\n      }\n    close ($O);\\
n    close ($I);\n    return;\n    }\n","use Env;\\
nuse Cwd;\n@suffix=(\"tmp\", \"temp\", \"cache\", \
\"t_coffee\", \"core\", \"tcoffee\");\n\nif ($#ARG\
V==-1)\n  {\n    print \"clean_cache.pl -file <fil\
e to add in -dir> -dir=<dir> -size=<value in Mb>\\\
n0: unlimited -1 always.\\nWill only clean directo\
ries matching:[\";\n    foreach $k(@suffix){print \
\"*$k* \";}\n    print \"]\\n\";\n    exit (EXIT_F\
AILURE);\n  }\n\n$cl=join (\" \",@ARGV);\nif (($cl\
=~/\\-no_action/))\n  {\n    exit (EXIT_SUCCESS);\\
n  }\n\nif (($cl=~/\\-debug/))\n  {\n    $DEBUG=1;\
\n  }\nelse\n  {\n    $DEBUG=0;\n  }\n\nif (($cl=~\
/\\-dir=(\\S+)/))\n  {\n    $dir=$1;\n  }\nelse\n \
 {\n    $dir=\"./\";\n  }\n\nif ($cl=~/\\-file=(\\\
S+)/)\n  {\n    $file=$1;\n  }\nelse\n  {\n    $fi\
le=0;\n  }\n\nif ($cl=~/\\-size=(\\S+)/)\n  {\n   \
 $max_size=$1;\n  }\nelse\n  {\n    $max_size=0;#u\
nlimited\n  }\nif ($cl=~/\\-force/)\n  {\n    $for\
ce=1;\n  }\nelse\n  {\n    $force=0;\n  }\n\nif ($\
cl=~/\\-age=(\\S+)/)\n  {\n    $max_age=$1;\n  }\n\
else\n  {\n    $max_age=0;#unlimited\n  }\n\n$max_\
size*=1000000;\nif ( ! -d $dir)\n  {\n    print ST\
DERR \"\\nCannot process $dir: does not exist \\n\\
";\n    exit (EXIT_FAILURE);\n  }\n\nif ( !($dir=~\
/^\\//))\n  {\n    $base=cwd();\n    $dir=\"$base/\
$dir\";\n  }\n\n$proceed=0;\nforeach $s (@suffix)\\
n  {\n    \n    if (($dir=~/$s/)){$proceed=1;}\n  \
  $s=uc ($s);\n    if (($dir=~/$s/)){$proceed=1;}\\
n  }\nif ( $proceed==0)\n  {\n    print STDERR \"C\
lean_cache.pl can only clean directories whose abs\
olute path name contains the following strings:\";\
\n    foreach $w (@suffix) {print STDERR \"$w \";$\
w=lc($w); print STDERR \"$w \";}\n    print STDERR\
 \"\\nCannot process $dir\\n\";\n    exit (EXIT_FA\
ILURE);\n  }\n\n$name_file=\"$dir/name_file.txt\";\
\n$size_file=\"$dir/size_file.txt\";\nif ( $force)\
{&create_ref_file ($dir,$name_file,$size_file);}\n\
if ($file){&add_file ($dir, $name_file, $size_file\
, $file);}\n&clean_dir ($dir, $name_file, $size_fi\
le, $max_size,$max_age);\nexit (EXIT_SUCCESS);\n\n\
sub clean_dir \n  {\n    my ($dir, $name_file, $si\
ze_file, $max_size, $max_age)=@_;\n    my ($tot_si\
ze, $size, $f, $s);\n\n  \n    $tot_size=&get_tot_\
size ($dir, $name_file, $size_file);\n\n    if ( $\
tot_size<=$max_size){return ;}\n    else {$max_siz\
e/=2;}\n    \n    #recreate the name file in case \
some temprary files have not been properly registe\
red\n    &create_ref_file ($dir, $name_file, $size\
_file, $max_age);\n  \n    $new_name_file=&vtmpnam\
();\n    open (R, \"$name_file\");\n    open (W, \\
">$new_name_file\");\n    while (<R>)\n      {\n	m\
y $line=$_;\n	\n	($f, $s)=($line=~/(\\S+) (\\S+)/)\
;\n	if ( !($f=~/\\S/)){next;}\n	\n	elsif ($max_siz\
e && $tot_size>=$max_size && !($f=~/name_file/))\n\
	  {\n	    remove ( \"$dir/$f\");\n	    $tot_size-\
=$s;\n	  }\n	elsif ( $max_age && -M(\"$dir/$f\")>=\
$max_age)\n	  {\n	    remove ( \"$dir/$f\");\n	   \
 $tot_size-=$s;\n	  }\n	else\n	  {\n	    print W \\
"$f $s\\n\";\n	  }\n      }\n    close (R);\n    c\
lose (W);\n    open (F, \">$size_file\");\n    pri\
nt F \"$tot_size\";\n    if ( -e $new_name_file){`\
mv $new_name_file $name_file`;}\n    close (F);\n \
 }\nsub get_tot_size\n  {\n    my ($dir, $name_fil\
e, $size_file)=@_;\n    my $size;\n    \n    if ( \
!-d $dir){return 0;}\n    if ( !-e $name_file)\n  \
    {\n	\n	&create_ref_file ($dir, $name_file, $si\
ze_file);\n      }\n    open (F, \"$size_file\");\\
n    $size=<F>;\n    close (F);\n    chomp ($size)\
;\n    return $size;\n  }\nsub size \n  {\n    my \
$f=@_[0];\n\n    if ( !-d $f){return -s($f);}\n   \
 else {return &dir2size($f);}\n  }\nsub dir2size\n\
  {\n    my $d=@_[0];\n    my ($s, $f);\n    \n   \
 if ( !-d $d) {return 0;}\n    \n    foreach $f (&\
dir2list ($d))\n      {\n	if ( -d $f){$s+=&dir2siz\
e (\"$d/$f\");}\n	else {$s+= -s \"$dir/$f\";}\n   \
   }\n    return $s;\n  }\n\nsub remove \n  {\n   \
 my $file=@_[0];\n    my ($f);\n    \n    debug_pr\
int( \"--- $file ---\\n\");\n    if (($file eq \".\
\") || ($file eq \"..\") || ($file=~/\\*/)){return\
 EXIT_FAILURE;}\n    elsif ( !-d $file)\n      {\n\
	debug_print (\"unlink $file\\n\");\n	if (-e $file\
){unlink ($file);}\n      }\n    elsif ( -d $file)\
\n      {\n	debug_print (\"++++++++ $file +++++++\\
\n\");\n	foreach $f (&dir2list($file))\n	  {\n	   \
 &remove (\"$file/$f\");\n	  }\n	debug_print (\"rm\
dir $file\\n\");\n	rmdir $file;\n      }\n    else\
\n      {\n	debug_print (\"????????? $file ???????\
?\\n\");\n      }\n    return EXIT_SUCCESS;\n  }\n\
\nsub dir2list\n  {\n    my $dir=@_[0];\n    my (@\
list1, @list2,@list3, $l);\n\n    opendir (DIR,$di\
r);\n    @list1=readdir (DIR);\n    closedir (DIR)\
;\n    \n    foreach $l (@list1)\n      {\n	if ( $\
l ne \".\" && $l ne \"..\"){@list2=(@list2, $l);}\\
n      }\n    @list3 = sort { (-M \"$dir/$list2[$b\
]\") <=> (-M \"$dir/$list2[$a]\")} @list2;\n    re\
turn @list3;\n    \n  }\n\nsub debug_print\n  {\n \
   \n    if ($DEBUG==1){print @_;}\n    \n  }\nsub\
 create_ref_file\n  {\n    my ($dir,$name_file,$si\
ze_file)=@_;\n    my ($f, $s, $tot_size, @l);\n   \
 \n    if ( !-d $dir){return;}\n    \n    @l=&dir2\
list ($dir);\n    open (F, \">$name_file\");\n    \
foreach $f (@l)\n      {\n	$s=&size(\"$dir/$f\");\\
n	$tot_size+=$s;\n	print F \"$f $s\\n\";\n      }\\
n    &myecho ($tot_size, \">$size_file\");\n    cl\
ose (F);\n  }\nsub add_file \n  {\n    my ($dir,$n\
ame_file,$size_file,$file)=@_;\n    my ($s, $tot_s\
ize);\n    \n    if ( !-d $dir)   {return;}\n    i\
f ( !-e \"$dir/$file\" ) {return;}\n    if ( !-e $\
name_file){&create_ref_file ($dir,$name_file,$size\
_file);}\n					    \n    $s=&size(\"$dir/$file\");\
\n    open (F, \">>$name_file\");\n    print F \"$\
file\\n\";\n    close (F);\n\n    $tot_size=&get_t\
ot_size ($dir,$name_file,$size_file);\n    $tot_si\
ze+=$s;\n    &myecho ($tot_size, \">$size_file\");\
\n    \n  }\n	\nsub myecho\n  {\n    my ($string, \
$file)=@_;\n    open (ECHO, $file) || die;\n    pr\
int ECHO \"$string\";\n    close (ECHO);\n  }\n   \
 \n		\n	\nsub vtmpnam\n  {\n    my $tmp_file_name;\
\n    $tmp_name_counter++;\n    $tmp_file_name=\"t\
mp_file_for_clean_cache_pdb$$.$tmp_name_counter\";\
\n    $tmp_file_list[$ntmp_file++]=$tmp_file_name;\
\n    if ( -e $tmp_file_name) {return &vtmpnam ();\
}\n    else {return $tmp_file_name;}\n  }\n","\nmy\
 $address=\"http://www.tcoffee.org/Data/Datasets/N\
atureProtocolsDataset.tar.gz\";\nmy $out=\"NatureP\
rotocolsDataset.tar.gz\";\n&url2file ($address,$ou\
t);\n\nif ( -e $out)\n  {\n    \n    system (\"gun\
zip NatureProtocolsDataset.tar.gz\");\n    system \
(\"tar -xvf NatureProtocolsDataset.tar\");\n  	sys\
tem (\"rm -rf NatureProtocolsDataset.tar\");  \n  \
  print \"Your Data Set is in the Folder 'NaturePr\
otocolsDataset'\\n\";\n  }\nelse \n  {\n    print \
\"Could not Download Dataset --- Web site may be d\
own -- Try again later\\n\";\n  }\n\n\n\n\nsub url\
2file\n{\n    my ($address, $out, $wget_arg, $curl\
_arg)=(@_);\n    my ($pg, $flag, $r, $arg, $count)\
;\n    \n    if (!$CONFIGURATION){&check_configura\
tion (\"wget\", \"INTERNET\", \"gzip\");$CONFIGURA\
TION=1;}\n    \n    if (&pg_is_installed (\"wget\"\
))   {$pg=\"wget\"; $flag=\"-O\";$arg=$wget_arg;}\\
n    elsif (&pg_is_installed (\"curl\")){$pg=\"cur\
l\"; $flag=\"-o\";$arg=$curl_arg;}\n    return sys\
tem (\"$pg $address $flag $out>/dev/null 2>/dev/nu\
ll\");\n\n}\n\nsub pg_is_installed\n  {\n    my @m\
l=@_;\n    my $r, $p, $m;\n    my $supported=0;\n \
   \n    my $p=shift (@ml);\n    if ($p=~/::/)\n  \
    {\n	if (system (\"perl -M$p -e 1\")==$EXIT_SUC\
CESS){return 1;}\n	else {return 0;}\n      }\n    \
else\n      {\n	$r=`which $p 2>/dev/null`;\n	if ($\
r eq \"\"){return 0;}\n	else {return 1;}\n      }\\
n  }\nsub check_configuration \n    {\n      my @l\
=@_;\n      my $v;\n      foreach my $p (@l)\n	{\n\
	  \n	  if   ( $p eq \"EMAIL\")\n	    { \n	      i\
f ( !($EMAIL=~/@/))\n		{\n		  exit (EXIT_FAILURE);\
\n		}\n	    }\n	  elsif( $p eq \"INTERNET\")\n	   \
 {\n	      if ( !&check_internet_connection())\n		\
{\n		  exit (EXIT_FAILURE);\n		}\n	    }\n	  elsif\
( $p eq \"wget\")\n	    {\n	      if (!&pg_is_inst\
alled (\"wget\") && !&pg_is_installed (\"curl\"))\\
n		{\n		  exit (EXIT_FAILURE);\n		}\n	    }\n	  el\
sif( !(&pg_is_installed ($p)))\n	    {\n	      exi\
t (EXIT_FAILURE);\n	    }\n	}\n      return 1;\n  \
  }\nsub check_internet_connection\n  {\n    my $i\
nternet;\n    my $tmp;\n    &check_configuration (\
 \"wget\"); \n    \n    $tmp=&vtmpnam ();\n    \n \
   if     (&pg_is_installed    (\"wget\")){`wget w\
ww.google.com -O$tmp >/dev/null 2>/dev/null`;}\n  \
  elsif  (&pg_is_installed    (\"curl\")){`curl ww\
w.google.com -o$tmp >/dev/null 2>/dev/null`;}\n   \
 \n    if ( !-e $tmp || -s $tmp < 10){$internet=0;\
}\n    else {$internet=1;}\n    if (-e $tmp){unlin\
k $tmp;}\n\n    return $internet;\n  }\n\nsub vtmp\
nam\n      {\n	my $r=rand(100000);\n	my $f=\"file.\
$r.$$\";\n	while (-e $f)\n	  {\n	    $f=vtmpnam();\
\n	  }\n	push (@TMPFILE_LIST, $f);\n	return $f;\n \
     }\n\n","\n$t_coffee=\"t_coffee\";\n\nforeach \
$value ( @ARGV)\n  {\n    $seq_file=$seq_file.\" \\
".$value;\n  }\n\n$name=$ARGV[0];\n$name=~s/\\.[^\\
\.]*$//;\n$lib_name=\"$name.mocca_lib\";\n$type=`t\
_coffee $seq_file -get_type -quiet`;\nchop ($type)\
;\n\nif ( $type eq \"PROTEIN\"){$lib_mode=\"lalign\
_rs_s_pair -lalign_n_top 20\";}\nelsif ( $type eq\\
"DNA\"){$lib_mode=\"lalign_rs_s_dna_pair -lalign_n\
_top 40\";}\n\nif ( !(-e $lib_name))\n  {\n	  \n  \
$command=\"$t_coffee -mocca -seq_weight=no -cosmet\
ic_penalty=0 -mocca_interactive -in $lib_mode -out\
_lib $lib_name -infile $seq_file\";\n  \n  }\nelsi\
f ( (-e $lib_name))\n  {\n  $command=\"$t_coffee -\
mocca -seq_weight=no -cosmetic_penalty=0 -mocca_in\
teractive -in $lib_name -infile $seq_file\";\n  \n\
  }\n\nsystem ($command);\n\nexit;\n\n","my $WSDL \
= 'http://www.ebi.ac.uk/Tools/webservices/wsdl/WSD\
aliLite.wsdl';\n\nuse SOAP::Lite;\nuse Data::Dumpe\
r;\nuse Getopt::Long qw(:config no_ignore_case bun\
dling);\nuse File::Basename;\n\nmy $checkInterval \
= 5;\n\nmy %params=(\n	    'async' => '1', # Use a\
sync mode and simulate sync mode in client\n	    )\
;\nGetOptions(\n    'pdb1=s'     => \\$params{'seq\
uence1'},\n    'chainid1=s' => \\$params{'chainid1\
'},\n    'pdb2=s'     => \\$params{'sequence2'},\n\
    'chainid2=s' => \\$params{'chainid2'},\n    \"\
help|h\"	 => \\$help, # Usage info\n    \"async|a\\
"	 => \\$async, # Asynchronous submission\n    \"p\
olljob\"	 => \\$polljob, # Get results\n    \"stat\
us\"	 => \\$status, # Get status\n    \"jobid|j=s\\
"  => \\$jobid, # JobId\n    \"email|S=s\"  => \\$\
params{email}, # E-mail address\n    \"trace\"    \
  => \\$trace, # SOAP messages\n    \"sequence=s\"\
 => \\$sequence, # Input PDB\n    );\n\nmy $script\
Name = basename($0, ());\nif($help) {\n    &usage(\
);\n    exit(0);\n}\n\nif($trace) {\n    print \"T\
racing active\\n\";\n    SOAP::Lite->import(+trace\
 => 'debug');\n}\n\nmy $soap = SOAP::Lite\n    ->s\
ervice($WSDL)\n    ->on_fault(sub {\n        my $s\
oap = shift;\n        my $res = shift;\n        # \
Throw an exception for all faults\n        if(ref(\
$res) eq '') {\n            die($res);\n        } \
else {\n            die($res->faultstring);\n     \
   }\n        return new SOAP::SOM;\n    }\n      \
         );\n\nif( !($polljob || $status) &&\n    \
!( defined($params{'sequence1'}) && defined($param\
s{'sequence2'}) )\n    ) {\n    print STDERR 'Erro\
r: bad option combination', \"\\n\";\n    &usage()\
;\n    exit(1);\n}\nelsif($polljob && defined($job\
id)) {\n    print \"Getting results for job $jobid\
\\n\";\n    getResults($jobid);\n}\nelsif($status \
&& defined($jobid)) {\n    print STDERR \"Getting \
status for job $jobid\\n\";\n    my $result = $soa\
p->checkStatus($jobid);\n    print STDOUT \"$resul\
t\", \"\\n\";\n    if($result eq 'DONE') {\n	print\
 STDERR \"To get results: $scriptName --polljob --\
jobid $jobid\\n\";\n    }\n}\nelse {\n    if(-f $p\
arams{'sequence1'}) {\n	$params{'sequence1'} = rea\
d_file($params{'sequence1'});\n    }\n    if(-f $p\
arams{'sequence2'}) {\n	$params{'sequence2'} = rea\
d_file($params{'sequence2'});\n    }\n\n    my $jo\
bid;\n    my $paramsData = SOAP::Data->name('param\
s')->type(map=>\\%params);\n    # For SOAP::Lite 0\
.60 and earlier parameters are passed directly\n  \
  if($SOAP::Lite::VERSION eq '0.60' || $SOAP::Lite\
::VERSION =~ /0\\.[1-5]/) {\n        $jobid = $soa\
p->runDaliLite($paramsData);\n    }\n    # For SOA\
P::Lite 0.69 and later parameter handling is diffe\
rent, so pass\n    # undef's for templated params,\
 and then pass the formatted args.\n    else {\n  \
      $jobid = $soap->runDaliLite(undef,\n				    \
 $paramsData);\n    }\n\n    if (defined($async)) \
{\n	print STDOUT $jobid, \"\\n\";\n        print S\
TDERR \"To check status: $scriptName --status --jo\
bid $jobid\\n\";\n    } else { # Synchronous mode\\
n        print STDERR \"JobId: $jobid\\n\";\n     \
   sleep 1;\n        getResults($jobid);\n    }\n}\
\n\nsub clientPoll($) {\n    my $jobid = shift;\n \
   my $result = 'PENDING';\n    # Check status and\
 wait if not finished\n    #print STDERR \"Checkin\
g status: $jobid\\n\";\n    while($result eq 'RUNN\
ING' || $result eq 'PENDING') {\n        $result =\
 $soap->checkStatus($jobid);\n        print STDERR\
 \"$result\\n\";\n        if($result eq 'RUNNING' \
|| $result eq 'PENDING') {\n            # Wait bef\
ore polling again.\n            sleep $checkInterv\
al;\n        }\n    }\n}\n\nsub getResults($) {\n \
   $jobid = shift;\n    # Check status, and wait i\
f not finished\n    clientPoll($jobid);\n    # Use\
 JobId if output file name is not defined\n    unl\
ess(defined($outfile)) {\n        $outfile=$jobid;\
\n    }\n    # Get list of data types\n    my $res\
ultTypes = $soap->getResults($jobid);\n    # Get t\
he data and write it to a file\n    if(defined($ou\
tformat)) { # Specified data type\n        my $sel\
ResultType;\n        foreach my $resultType (@$res\
ultTypes) {\n            if($resultType->{type} eq\
 $outformat) {\n                $selResultType = $\
resultType;\n            }\n        }\n        $re\
s=$soap->poll($jobid, $selResultType->{type});\n  \
      write_file($outfile.'.'.$selResultType->{ext\
}, $res);\n    } else { # Data types available\n  \
      # Write a file for each output type\n       \
 for my $resultType (@$resultTypes){\n            \
#print \"Getting $resultType->{type}\\n\";\n      \
      $res=$soap->poll($jobid, $resultType->{type}\
);\n            write_file($outfile.'.'.$resultTyp\
e->{ext}, $res);\n        }\n    }\n}\n\nsub read_\
file($) {\n    my $filename = shift;\n    open(FIL\
E, $filename);\n    my $content;\n    my $buffer;\\
n    while(sysread(FILE, $buffer, 1024)) {\n	$cont\
ent.= $buffer;\n    }\n    close(FILE);\n    retur\
n $content;\n}\n\nsub write_file($$) {\n    my ($t\
mp,$entity) = @_;\n    print STDERR \"Creating res\
ult file: \".$tmp.\"\\n\";\n    unless(open (FILE,\
 \">$tmp\")) {\n	return 0;\n    }\n    syswrite(FI\
LE, $entity);\n    close (FILE);\n    return 1;\n}\
\n\nsub usage {\n    print STDERR <<EOF\nDaliLite\\
n========\n\nPairwise comparison of protein struct\
ures\n\n[Required]\n\n  --pdb1                : st\
r  : PDB ID for structure 1\n  --pdb2             \
   : str  : PDB ID for structure 2\n\n[Optional]\n\
\n  --chain1              : str  : Chain identifer\
 in structure 1\n  --chain2              : str  : \
Chain identifer in structure 2\n\n[General]\n\n  -\
h, --help            :      : prints this help tex\
t\n  -S, --email           : str  : user email add\
ress\n  -a, --async           :      : asynchronou\
s submission\n      --status          :      : pol\
l for the status of a job\n      --polljob        \
 :      : poll for the results of a job\n  -j, --j\
obid           : str  : jobid for an asynchronous \
job\n  -O, --outfile         : str  : file name fo\
r results (default is jobid)\n      --trace	      \
  :      : show SOAP messages being interchanged \\
n\nSynchronous job:\n\n  The results/errors are re\
turned as soon as the job is finished.\n  Usage: $\
scriptName --email <your\\@email> [options] pdbFil\
e [--outfile string]\n  Returns: saves the results\
 to disk\n\nAsynchronous job:\n\n  Use this if you\
 want to retrieve the results at a later time. The\
 results \n  are stored for up to 24 hours. \n  Th\
e asynchronous submission mode is recommended when\
 users are submitting \n  batch jobs or large data\
base searches	\n  Usage: $scriptName --email <your\
\\@email> --async [options] pdbFile\n  Returns: jo\
bid\n\n  Use the jobid to query for the status of \
the job. \n  Usage: $scriptName --status --jobid <\
jobId>\n  Returns: string indicating the status of\
 the job:\n    DONE - job has finished\n    RUNNIN\
G - job is running\n    NOT_FOUND - job cannot be \
found\n    ERROR - the jobs has encountered an err\
or\n\n  When done, use the jobid to retrieve the s\
tatus of the job. \n  Usage: $scriptName --polljob\
 --jobid <jobId> [--outfile string]\n\n[Help]\n\n \
 For more detailed help information refer to\n  ht\
tp://www.ebi.ac.uk/DaliLite/\nEOF\n;\n}\n","my $WS\
DL = 'http://www.ebi.ac.uk/Tools/webservices/wsdl/\
WSWUBlast.wsdl';\n\nuse strict;\nuse SOAP::Lite;\n\
use Getopt::Long qw(:config no_ignore_case bundlin\
g);\nuse File::Basename;\n\nmy $checkInterval = 15\
;\n\nmy $numOpts = scalar(@ARGV);\nmy ($outfile, $\
outformat, $help, $async, $polljob, $status, $ids,\
 $jobid, $trace, $sequence);\nmy %params= ( # Defa\
ults\n	      'async' => 1, # Force into async mode\
\n	      'exp' => 10.0, # E-value threshold\n	    \
  'numal' => 50, # Maximum number of alignments\n	\
      'scores' => 100, # Maximum number of scores\\
n            );\nGetOptions( # Map the options int\
o variables\n    \"program|p=s\"     => \\$params{\
program}, # BLAST program\n    \"database|D=s\"   \
 => \\$params{database}, # Search database\n    \"\
matrix|m=s\"      => \\$params{matrix}, # Scoring \
matrix\n    \"exp|E=f\"         => \\$params{exp},\
 # E-value threshold\n    \"echofilter|e\"    => \\
\$params{echofilter}, # Display filtered sequence\\
n    \"filter|f=s\"      => \\$params{filter}, # L\
ow complexity filter name\n    \"alignments|b=i\" \
 => \\$params{numal}, # Number of alignments\n    \
\"scores|s=i\"      => \\$params{scores}, # Number\
 of scores\n    \"sensitivity|S=s\" => \\$params{s\
ensitivity}, # Search sensitivity\n    \"sort|t=s\\
"	      => \\$params{sort}, # Sort hits by...\n   \
 \"stats|T=s\"       => \\$params{stats}, # Scorin\
g statistic to use\n    \"strand|d=s\"      => \\$\
params{strand}, # Strand to use in DNA vs. DNA sea\
rch\n    \"topcombon|c=i\"   => \\$params{topcombo\
n}, # Consistent sets of HSPs\n    \"outfile=s\"  \
     => \\$outfile, # Output file\n    \"outformat\
|o=s\"   => \\$outformat, # Output format\n    \"h\
elp|h\"	      => \\$help, # Usage info\n    \"asyn\
c|a\"	      => \\$async, # Asynchronous mode\n    \
\"polljob\"	      => \\$polljob, # Get results\n  \
  \"status\"	      => \\$status, # Get job status\\
n    \"ids\"             => \\$ids, # Get ids from\
 result\n    \"jobid|j=s\"       => \\$jobid, # Jo\
bId\n    \"email=s\"         => \\$params{email}, \
# E-mail address\n    \"trace\"           => \\$tr\
ace, # SOAP trace\n    \"sequence=s\"      => \\$s\
equence, # Query sequence\n    );\n\nmy $scriptNam\
e = basename($0, ());\nif($help || $numOpts == 0) \
{\n    &usage();\n    exit(0);\n}\n\nif($trace){\n\
    print STDERR \"Tracing active\\n\";\n    SOAP:\
:Lite->import(+trace => 'debug');\n}\n\nmy $soap =\
 SOAP::Lite\n    ->service($WSDL)\n    ->proxy('ht\
tp://localhost/',\n    #proxy => ['http' => 'http:\
//your.proxy.server/'], # HTTP proxy\n    timeout \
=> 600, # HTTP connection timeout\n    )\n    ->on\
_fault(sub { # SOAP fault handler\n        my $soa\
p = shift;\n        my $res = shift;\n        # Th\
row an exception for all faults\n        if(ref($r\
es) eq '') {\n            die($res);\n        } el\
se {\n            die($res->faultstring);\n       \
 }\n        return new SOAP::SOM;\n    }\n        \
       );\n\nif( !($polljob || $status || $ids) &&\
\n    !( defined($ARGV[0]) || defined($sequence) )\
\n    ) {\n    print STDERR 'Error: bad option com\
bination', \"\\n\";\n    &usage();\n    exit(1);\n\
}\nelsif($polljob && defined($jobid)) {\n    print\
 \"Getting results for job $jobid\\n\";\n    getRe\
sults($jobid);\n}\nelsif($status && defined($jobid\
)) {\n    print STDERR \"Getting status for job $j\
obid\\n\";\n    my $result = $soap->checkStatus($j\
obid);\n    print STDOUT \"$result\\n\";\n    if($\
result eq 'DONE') {\n	print STDERR \"To get result\
s: $scriptName --polljob --jobid $jobid\\n\";\n   \
 }\n}  \nelsif($ids && defined($jobid)) {\n    pri\
nt STDERR \"Getting ids from job $jobid\\n\";\n   \
 getIds($jobid);\n}\nelse {\n    # Prepare input d\
ata\n    my $content;\n    my (@contents) = ();\n \
   if(-f $ARGV[0] || $ARGV[0] eq '-') {	\n	$conten\
t={type=>'sequence',content=>read_file($ARGV[0])};\
	\n    }\n    if($sequence) {	\n	if(-f $sequence |\
| $sequence eq '-') {	\n	    $content={type=>'sequ\
ence',content=>read_file($ARGV[0])};	\n	} else {\n\
	    $content={type=>'sequence',content=>$sequence\
};\n	}\n    }\n    push @contents, $content;\n\n  \
  # Submit the job\n    my $paramsData = SOAP::Dat\
a->name('params')->type(map=>\\%params);\n    my $\
contentData = SOAP::Data->name('content')->value(\\
\@contents);\n    # For SOAP::Lite 0.60 and earlie\
r parameters are passed directly\n    if($SOAP::Li\
te::VERSION eq '0.60' || $SOAP::Lite::VERSION =~ /\
0\\.[1-5]/) {\n        $jobid = $soap->runWUBlast(\
$paramsData, $contentData);\n    }\n    # For SOAP\
::Lite 0.69 and later parameter handling is differ\
ent, so pass\n    # undef's for templated params, \
and then pass the formatted args.\n    else {\n   \
     $jobid = $soap->runWUBlast(undef, undef,\n			\
	   $paramsData, $contentData);\n    }\n\n    # As\
ynchronous mode: output jobid and exit.\n    if (d\
efined($async)) {\n	print STDOUT $jobid, \"\\n\";\\
n        print STDERR \"To check status: $scriptNa\
me --status --jobid $jobid\\n\";\n    }\n    # Syn\
chronous mode: try to get results\n    else {\n   \
     print STDERR \"JobId: $jobid\\n\";\n        s\
leep 1;\n        getResults($jobid);\n    }\n}\n\n\
sub getIds($) {\n    my $jobid = shift;\n    my $r\
esults = $soap->getIds($jobid);\n    for my $resul\
t (@$results){\n	print \"$result\\n\";\n    }\n}\n\
\nsub clientPoll($) {\n    my $jobid = shift;\n   \
 my $result = 'PENDING';\n    # Check status and w\
ait if not finished\n    while($result eq 'RUNNING\
' || $result eq 'PENDING') {\n        $result = $s\
oap->checkStatus($jobid);\n        print STDERR \"\
$result\\n\";\n        if($result eq 'RUNNING' || \
$result eq 'PENDING') {\n            # Wait before\
 polling again.\n            sleep $checkInterval;\
\n        }\n    }\n}\n\nsub getResults($) {\n    \
my $jobid = shift;\n    my $res;\n    # Check stat\
us, and wait if not finished\n    clientPoll($jobi\
d);\n    # Use JobId if output file name is not de\
fined\n    unless(defined($outfile)) {\n        $o\
utfile=$jobid;\n    }\n    # Get list of data type\
s\n    my $resultTypes = $soap->getResults($jobid)\
;\n    # Get the data and write it to a file\n    \
if(defined($outformat)) { # Specified data type\n	\
if($outformat eq 'xml') {$outformat = 'toolxml';}\\
n	if($outformat eq 'txt') {$outformat = 'tooloutpu\
t';}\n        my $selResultType;\n        foreach \
my $resultType (@$resultTypes) {\n            if($\
resultType->{type} eq $outformat) {\n             \
   $selResultType = $resultType;\n            }\n \
       }\n        $res=$soap->poll($jobid, $selRes\
ultType->{type});\n	if($outfile eq '-') {\n	     w\
rite_file($outfile, $res);\n	} else {\n	    write_\
file($outfile.'.'.$selResultType->{ext}, $res);\n	\
}\n    } else { # Data types available\n        # \
Write a file for each output type\n        for my \
$resultType (@$resultTypes){\n            #print S\
TDERR \"Getting $resultType->{type}\\n\";\n       \
     $res=$soap->poll($jobid, $resultType->{type})\
;\n	    if($outfile eq '-') {\n		write_file($outfi\
le, $res);\n	    } else {\n		write_file($outfile.'\
.'.$resultType->{ext}, $res);\n	    }\n        }\n\
    }\n}\n\nsub read_file($) {\n    my $filename =\
 shift;\n    my ($content, $buffer);\n    if($file\
name eq '-') {\n	while(sysread(STDIN, $buffer, 102\
4)) {\n	    $content .= $buffer;\n	}\n    }\n    e\
lse { # File\n	open(FILE, $filename) or die \"Erro\
r: unable to open input file\";\n	while(sysread(FI\
LE, $buffer, 1024)) {\n	    $content .= $buffer;\n\
	}\n	close(FILE);\n    }\n    return $content;\n}\\
n\nsub write_file($$) {\n    my ($filename, $data)\
 = @_;\n    print STDERR 'Creating result file: ' \
. $filename . \"\\n\";\n    if($filename eq '-') {\
\n	print STDOUT $data;\n    }\n    else {\n	open(F\
ILE, \">$filename\") or die \"Error: unable to ope\
n output file\";\n	syswrite(FILE, $data);\n	close(\
FILE);\n    }\n}\n\nsub usage {\n    print STDERR \
<<EOF\nWU-BLAST\n========\n\nRapid sequence databa\
se search programs utilizing the BLAST algorithm.\\
n   \n[Required]\n\n      --email       : str  : u\
ser email address \n  -p, --program	    : str  : B\
LAST program to use: blastn, blastp, blastx, \n   \
                          tblastn or tblastx\n  -D\
, --database    : str  : database to search\n  seq\
File           : file : query sequence data file (\
\"-\" for STDIN)\n\n[Optional]\n\n  -m, --matrix	 \
   : str  : scoring matrix\n  -E, --exp	    : real\
 : 0<E<= 1000. Statistical significance threshold\\
n                             for reporting databa\
se sequence matches.\n  -e, --echofilter  :      :\
 display the filtered query sequence in the output\
\n  -f, --filter	    : str  : activates filtering \
of the query sequence\n  -b, --alignments  : int  \
: number of alignments to be reported\n  -s, --sco\
res	    : int  : number of scores to be reported\n\
  -S, --sensitivity : str  :\n  -t, --sort	    : s\
tr  :\n  -T, --stats       : str  :\n  -d, --stran\
d      : str  : DNA strand to search with in DNA v\
s. DNA searches \n  -c, --topcombon   :      :\n\n\
[General]	\n\n  -h, --help       :      : prints t\
his help text\n  -a, --async      :      : forces \
to make an asynchronous query\n      --status     \
:      : poll for the status of a job\n      --pol\
ljob    :      : poll for the results of a job\n  \
-j, --jobid      : str  : jobid that was returned \
when an asynchronous job \n                       \
     was submitted.\n  -O, --outfile    : str  : n\
ame of the file results should be written to \n   \
                         (default is based on the \
jobid; \"-\" for STDOUT)\n  -o, --outformat  : str\
  : txt or xml output (no file is written)\n      \
--trace	   :      : show SOAP messages being inter\
changed \n\nSynchronous job:\n\n  The results/erro\
rs are returned as soon as the job is finished.\n \
 Usage: $scriptName --email <your\\@email> [option\
s...] seqFile\n  Returns: saves the results to dis\
k\n\nAsynchronous job:\n\n  Use this if you want t\
o retrieve the results at a later time. The result\
s \n  are stored for up to 24 hours. \n  The async\
hronous submission mode is recommended when users \
are submitting \n  batch jobs or large database se\
arches	\n  Usage: $scriptName --async --email <you\
r\\@email> [options...] seqFile\n  Returns : jobid\
\n\n  Use the jobid to query for the status of the\
 job. \n  Usage: $scriptName --status --jobid <job\
Id>\n  Returns : string indicating the status of t\
he job:\n    DONE - job has finished\n    RUNNING \
- job is running\n    NOT_FOUND - job cannot be fo\
und\n    ERROR - the jobs has encountered an error\
\n\n  When done, use the jobid to retrieve the sta\
tus of the job. \n  Usage: $scriptName --polljob -\
-jobid <jobId> [--outfile string]\n  Returns: save\
s the results to disk\n\n[Help]\n\nFor more detail\
ed help information refer to \nhttp://www.ebi.ac.u\
k/blast2/WU-Blast2_Help_frame.html\n \nEOF\n;\n}\n\
","\nmy $WSDL = 'http://www.ebi.ac.uk/Tools/webser\
vices/wsdl/WSBlastpgp.wsdl';\n\nuse SOAP::Lite;\nu\
se Getopt::Long qw(:config no_ignore_case bundling\
);\nuse File::Basename;\n\nmy $checkInterval = 15;\
\n\nmy %params=(\n	    'async' => '1', # Use async\
 mode and simulate sync mode in client\n	    );\nG\
etOptions(\n    \"mode=s\"           => \\$params{\
mode}, # Search mode: PSI-Blast or PHI-Blast\n    \
\"database|d=s\"     => \\$params{database}, # Dat\
abase to search\n    \"matrix|M=s\"       => \\$pa\
rams{matrix},# Scoring maxtrix\n    \"exp|e=f\"   \
       => \\$params{exp}, # E-value\n    \"expmult\
i|h=f\"     => \\$params{expmulti}, # E-value\n   \
 \"filter|F=s\"       => \\$params{filter}, # Low \
complexity filter\n    \"dropoff|X=i\"      => \\$\
params{dropoff}, # Dropoff score\n    \"finaldropo\
ff|Z=i\" => \\$params{finaldropoff}, # Final dropo\
ff score\n    \"scores|v=i\"       => \\$params{sc\
ores}, # Max number of scores\n    \"align=i\"    \
      => \\$params{align}, # Alignment view\n    \\
"startregion|S=i\"  => \\$params{startregion}, # S\
tart of region in query\n    \"endregion|H=i\"    \
=> \\$params{endregion}, # End of region in query\\
n    \"maxpasses|j=i\"    => \\$params{maxpasses},\
 # Number of PSI iterations\n    \"opengap|G=i\"  \
    => \\$params{opengap}, # Gap open penalty\n   \
 \"extendgap|E=i\"    => \\$params{extendgap}, # G\
ap extension penalty\n    \"pattern=s\"        => \
\\$params{pattern}, # PHI-BLAST pattern\n    \"usa\
gemode|p=s\"    => \\$params{usagemode}, # PHI-BLA\
ST program\n    \"appxml=s\"         => \\$params{\
appxml}, # Application XML\n    \"sequence=s\"    \
   => \\$sequence, # Query sequence\n    \"help\"	\
       => \\$help, # Usage info\n    \"polljob\"	 \
      => \\$polljob, # Get results\n    \"status\"\
	       => \\$status, # Get status\n    \"ids\"   \
   	       => \\$ids, # Get ids from result\n    \\
"jobid=s\"          => \\$jobid, # JobId\n    \"ou\
tfile=s\"        => \\$outfile, # Output filename\\
n    \"outformat|o=s\"    => \\$outformat, # Outpu\
t file format\n    \"async|a\"	       => \\$async,\
 # Async submission\n    \"email=s\"          => \\
\$params{email}, # User e-mail address\n    \"trac\
e\"            => \\$trace, # Show SOAP messages\n\
    );\n\nmy $scriptName = basename($0, ());\nif($\
help) {\n    &usage();\n    exit(0);\n}\n\nif ($tr\
ace){\n    print \"Tracing active\\n\";\n    SOAP:\
:Lite->import(+trace => 'debug');\n}\n\nmy $soap =\
 SOAP::Lite\n    ->service($WSDL)\n    ->on_fault(\
sub {\n        my $soap = shift;\n        my $res \
= shift;\n        # Throw an exception for all fau\
lts\n        if(ref($res) eq '') {\n            di\
e($res);\n        } else {\n            die($res->\
faultstring);\n        }\n        return new SOAP:\
:SOM;\n    }\n               );\n\nif( !($polljob \
|| $status || $ids) &&\n    !( (defined($ARGV[0]) \
&& -f $ARGV[0]) || defined($sequence) )\n    ) {\n\
    print STDERR 'Error: bad option combination', \
\"\\n\";\n    &usage();\n    exit(1);\n}\nelsif($p\
olljob && defined($jobid)) {\n    print \"Getting \
results for job $jobid\\n\";\n    getResults($jobi\
d);\n}\nelsif($status && defined($jobid)) {\n    p\
rint STDERR \"Getting status for job $jobid\\n\";\\
n    my $result = $soap->checkStatus($jobid);\n   \
 print STDOUT $result, \"\\n\";\n    if($result eq\
 'DONE') {\n	print STDERR \"To get results: $scrip\
tName --polljob --jobid $jobid\\n\";\n    }\n}  \n\
elsif($ids && defined($jobid)) {\n    print STDERR\
 \"Getting ids from job $jobid\\n\";\n    getIds($\
jobid);\n}\nelse {\n    if(-f $ARGV[0]) {	\n	$cont\
ent={type=>'sequence', content=>read_file($ARGV[0]\
)};	\n    }\n    if($sequence) {	\n	if(-f $sequenc\
e) {\n	    $content={type=>'sequence', content=>re\
ad_file($sequence)};	\n	} else {\n	    $content={t\
ype=>'sequence', content=>$sequence};\n	}\n    }\n\
    push @content, $content;\n\n    my $jobid;\n  \
  my $paramsData = SOAP::Data->name('params')->typ\
e(map=>\\%params);\n    my $contentData = SOAP::Da\
ta->name('content')->value(\\@content);\n    # For\
 SOAP::Lite 0.60 and earlier parameters are passed\
 directly\n    if($SOAP::Lite::VERSION eq '0.60' |\
| $SOAP::Lite::VERSION =~ /0\\.[1-5]/) {\n        \
$jobid = $soap->runBlastpgp($paramsData, $contentD\
ata);\n    }\n    # For SOAP::Lite 0.69 and later \
parameter handling is different, so pass\n    # un\
def's for templated params, and then pass the form\
atted args.\n    else {\n        $jobid = $soap->r\
unBlastpgp(undef, undef,\n				    $paramsData, $co\
ntentData);\n    }\n\n    if (defined($async)) {\n\
	print STDOUT $jobid, \"\\n\";\n        print STDE\
RR \"To check status: $scriptName --status --jobid\
 $jobid\\n\";\n    } else { # Synchronous mode\n  \
      print STDERR \"JobId: $jobid\\n\";\n        \
sleep 1;\n        getResults($jobid);\n    }\n}\n\\
nsub getIds($) {\n    $jobid = shift;\n    my $res\
ults = $soap->getIds($jobid);\n    for $result (@$\
results){\n	print \"$result\\n\";\n    }\n}\n\nsub\
 clientPoll($) {\n    my $jobid = shift;\n    my $\
result = 'PENDING';\n    # Check status and wait i\
f not finished\n    #print STDERR \"Checking statu\
s: $jobid\\n\";\n    while($result eq 'RUNNING' ||\
 $result eq 'PENDING') {\n        $result = $soap-\
>checkStatus($jobid);\n        print STDERR \"$res\
ult\\n\";\n        if($result eq 'RUNNING' || $res\
ult eq 'PENDING') {\n            # Wait before pol\
ling again.\n            sleep $checkInterval;\n  \
      }\n    }\n}\n\nsub getResults($) {\n    $job\
id = shift;\n    # Check status, and wait if not f\
inished\n    clientPoll($jobid);\n    # Use JobId \
if output file name is not defined\n    unless(def\
ined($outfile)) {\n        $outfile=$jobid;\n    }\
\n    # Get list of data types\n    my $resultType\
s = $soap->getResults($jobid);\n    # Get the data\
 and write it to a file\n    if(defined($outformat\
)) { # Specified data type\n        my $selResultT\
ype;\n        foreach my $resultType (@$resultType\
s) {\n            if($resultType->{type} eq $outfo\
rmat) {\n                $selResultType = $resultT\
ype;\n            }\n        }\n        $res=$soap\
->poll($jobid, $selResultType->{type});\n        w\
rite_file($outfile.'.'.$selResultType->{ext}, $res\
);\n    } else { # Data types available\n        #\
 Write a file for each output type\n        for my\
 $resultType (@$resultTypes){\n            #print \
\"Getting $resultType->{type}\\n\";\n            $\
res=$soap->poll($jobid, $resultType->{type});\n   \
         write_file($outfile.'.'.$resultType->{ext\
}, $res);\n        }\n    }\n}\n\nsub read_file($)\
 {\n    my $filename = shift;\n    open(FILE, $fil\
ename);\n    my $content;\n    my $buffer;\n    wh\
ile(sysread(FILE, $buffer, 1024)) {\n	$content.= $\
buffer;\n    }\n    close(FILE);  \n    return $co\
ntent;\n}\n\nsub write_file($$) {\n    my ($tmp,$e\
ntity) = @_;\n    print STDERR \"Creating result f\
ile: \".$tmp.\"\\n\";\n    unless(open (FILE, \">$\
tmp\")) {\n	return 0;\n    }\n    syswrite(FILE, $\
entity);\n    close (FILE);\n    return 1;\n}\n\ns\
ub usage {\n    print STDERR <<EOF\nBlastpgp\n====\
====\n   \nThe blastpgp program implements the PSI\
-BLAST and PHI-BLAST variations\nof NCBI BLAST.\n\\
nFor more detailed help information refer to\nhttp\
://www.ebi.ac.uk/blastpgp/blastpsi_help_frame.html\
\n \nBlastpgp specific options:\n\n[Required]\n\n \
     --mode            : str  : search mode to use\
: PSI-Blast or PHI-Blast\n  -d, --database        \
: str  : protein database to search\n  seqFile    \
           : file : query sequence\n\n[Optional]\n\
\n  -M, --matrix          : str  : scoring matrix\\
n  -e, --exp             : real : Expectation valu\
e\n  -h, --expmulti        : real : threshold (mul\
tipass model)\n  -F, --filter          : str  : fi\
lter query sequence with SEG [T,F]\n  -m, --align \
          : int  : alignment view option:\n       \
                          0 - pairwise, 1 - M/S id\
entities,\n                                 2 - M/\
S non-identities, 3 - Flat identities,\n          \
                       4 - Flat non-identities\n  \
-G, --opengap         : int  : cost to open a gap\\
n  -E, --extendgap       : int  : cost to extend a\
 gap\n  -g, --gapalign        : str  : Gapped [T,F\
]\n  -v, --scores          : int  : number of scor\
es to be reported\n  -j, --maxpasses       : int  \
: number of iterations\n  -X, --dropoff         : \
int  : Dropoff score\n  -Z, --finaldropoff    : in\
t  : Dropoff for final alignment\n  -S, --startreg\
ion     : int  : Start of required region in query\
\n  -H, --endregion       : int  : End of required\
 region in query\n  -k, --pattern         : str  :\
 Hit File (PHI-BLAST only)\n  -p, --usagemode     \
  : str  : Program option (PHI-BLAST only):\n     \
                            blastpgp, patseedp, se\
edp\n\n[General]\n\n      --help            :     \
 : prints this help text\n  -a, --async           \
:      : forces to make an asynchronous query\n   \
   --status          :      : poll for the status \
of a job\n      --polljob         :      : poll fo\
r the results of a job\n      --jobid           : \
str  : jobid of an asynchronous job\n      --ids  \
           :      : get hit identifiers for result\
 \n  -O, --outfile         : str  : name of the fi\
le results should be written to\n                 \
                (default is based on the jobid)\n \
 -o, --outformat       : str  : txt or xml output \
(no file is written)\n      --trace           :   \
   : show SOAP messages being interchanged\n\nSync\
hronous job:\n\n  The results/errors are returned \
as soon as the job is finished.\n  Usage: blastpgp\
.pl --email <your@email> [options...] seqfile\n  R\
eturns: saves the results to disk\n\nAsynchronous \
job:\n\n  Use this if you want to retrieve the res\
ults at a later time. The results\n  are stored fo\
r up to 24 hours.\n  The asynchronous submission m\
ode is recommended when users are submitting\n  ba\
tch jobs or large database searches\n  Usage: blas\
tpgp.pl --email <your@email> --async [options...] \
seqFile\n  Returns: jobid\n\n  Use the jobid to qu\
ery for the status of the job.\n  Usage: blastpgp.\
pl --status --jobid <jobId>\n  Returns: string ind\
icating the status of the job\n    DONE - job has \
finished\n    RUNNING - job is running\n    NOT_FO\
UND - job cannot be found\n    ERROR - the jobs ha\
s encountered an error\n\n  When done, use the job\
id to retrieve the results of the job.\n  Usage: b\
lastpgp.pl --polljob --jobid <jobId> [--outfile <f\
ileName>]\n  Returns: saves the results to disk\nE\
OF\n;\n}\n","\n=head1 NAME\n\nncbiblast_lwp.pl\n\n\
=head1 DESCRIPTION\n\nNCBI BLAST (REST) web servic\
e Perl client using L<LWP>.\n\nTested with:\n\n=ov\
er\n\n=item *\nL<LWP> 5.79, L<XML::Simple> 2.12 an\
d Perl 5.8.3\n\n=item *\nL<LWP> 5.808, L<XML::Simp\
le> 2.18 and Perl 5.8.8 (Ubuntu 8.04 LTS)\n\n=item\
 *\nL<LWP> 5.834, L<XML::Simple> 2.18 and Perl 5.1\
0.1 (Ubuntu 10.04 LTS)\n\n=item *\nL<LWP> 6.03, L<\
XML::Simple> 2.18 and Perl 5.14.2 (Ubuntu 12.04 LT\
S)\n\n=back\n\nFor further information see:\n\n=ov\
er\n\n=item *\nL<http://www.ebi.ac.uk/Tools/webser\
vices/services/sss/ncbi_blast_rest>\n\n=item *\nL<\
http://www.ebi.ac.uk/Tools/webservices/tutorials/p\
erl>\n\n=back\n\n=head1 LICENSE\n\nCopyright 2012-\
2013 EMBL - European Bioinformatics Institute\n\nL\
icensed under the Apache License, Version 2.0 (the\
 \"License\");\nyou may not use this file except i\
n compliance with the License.\nYou may obtain a c\
opy of the License at\n\n    http://www.apache.org\
/licenses/LICENSE-2.0\n\nUnless required by applic\
able law or agreed to in writing, software\ndistri\
buted under the License is distributed on an \"AS \
IS\" BASIS,\nWITHOUT WARRANTIES OR CONDITIONS OF A\
NY KIND, either express or implied.\nSee the Licen\
se for the specific language governing permissions\
 and\nlimitations under the License.\n\n=head1 VER\
SION\n\n$Id: ncbiblast_lwp.pl 2560 2013-03-20 12:5\
6:31Z hpm $\n\n=cut\n\nuse strict;\nuse warnings;\\
n\nuse English;\nuse LWP;\nuse XML::Simple;\nuse G\
etopt::Long qw(:config no_ignore_case bundling);\n\
use File::Basename;\nuse Data::Dumper;\n\nmy $base\
Url = 'http://www.ebi.ac.uk/Tools/services/rest/nc\
biblast';\n\nmy $checkInterval = 3;\n\nmy $outputL\
evel = 1;\n\nmy $numOpts = scalar(@ARGV);\nmy %par\
ams = ( 'debugLevel' => 0 );\n\nmy %tool_params = \
();\nGetOptions(\n\n	# Tool specific options\n	'pr\
ogram|p=s'  => \\$tool_params{'program'},   # blas\
tp, blastn, blastx, etc.\n	'database|D=s' => \\$pa\
rams{'database'},       # Database(s) to search\n	\
'matrix|m=s'   => \\$tool_params{'matrix'},    # S\
coring martix to use\n	'exp|E=f'      => \\$tool_p\
arams{'exp'},       # E-value threshold\n	'filter|\
f=s'   => \\$tool_params{'filter'},    # Low compl\
exity filter\n	'align|A=i'    => \\$tool_params{'a\
lign'},     # Pairwise alignment format\n	'scores|\
s=i'   => \\$tool_params{'scores'},    # Number of\
 scores\n	'alignments|n=i' => \\$tool_params{'alig\
nments'},   # Number of alignments\n	'dropoff|d=i'\
    => \\$tool_params{'dropoff'},      # Dropoff s\
core\n	'match_scores=s' => \\$tool_params{'match_s\
cores'}, # Match/missmatch scores\n	'match|u=i'   \
   => \\$params{'match'},             # Match scor\
e\n	'mismatch|v=i'   => \\$params{'mismatch'},    \
      # Mismatch score\n	'gapopen|o=i'    => \\$to\
ol_params{'gapopen'},      # Open gap penalty\n	'g\
apext|x=i'     => \\$tool_params{'gapext'},       \
# Gap extension penality\n	'gapalign|g'     => \\$\
tool_params{'gapalign'},     # Optimise gap alignm\
ents\n	'stype=s' => \\$tool_params{'stype'},    # \
Sequence type\n	'seqrange=s' => \\$tool_params{'se\
qrange'},    # Query subsequence\n	'sequence=s' =>\
 \\$params{'sequence'},         # Query sequence\n\
	'multifasta' => \\$params{'multifasta'},       # \
Multiple fasta input\n\n	# Compatability options, \
old command-line\n	'numal|n=i'     => \\$params{'n\
umal'},        # Number of alignments\n	'opengap|o\
=i'   => \\$params{'opengap'},      # Open gap pen\
alty\n	'extendgap|x=i' => \\$params{'extendgap'}, \
   # Gap extension penality\n	\n	# Generic options\
\n	'email=s'       => \\$params{'email'},         \
 # User e-mail address\n	'title=s'       => \\$par\
ams{'title'},          # Job title\n	'outfile=s'  \
   => \\$params{'outfile'},        # Output file n\
ame\n	'outformat=s'   => \\$params{'outformat'},  \
    # Output file type\n	'jobid=s'       => \\$par\
ams{'jobid'},          # JobId\n	'help|h'        =\
> \\$params{'help'},           # Usage help\n	'asy\
nc'         => \\$params{'async'},          # Asyn\
chronous submission\n	'polljob'       => \\$params\
{'polljob'},        # Get results\n	'resultTypes' \
  => \\$params{'resultTypes'},    # Get result typ\
es\n	'status'        => \\$params{'status'},      \
   # Get status\n	'params'        => \\$params{'pa\
rams'},         # List input parameters\n	'paramDe\
tail=s' => \\$params{'paramDetail'},    # Get deta\
ils for parameter\n	'quiet'         => \\$params{'\
quiet'},          # Decrease output level\n	'verbo\
se'       => \\$params{'verbose'},        # Increa\
se output level\n	'debugLevel=i'  => \\$params{'de\
bugLevel'},     # Debug output level\n	'baseUrl=s'\
     => \\$baseUrl,                  # Base URL fo\
r service.\n);\nif ( $params{'verbose'} ) { $outpu\
tLevel++ }\nif ( $params{'quiet'} )  { $outputLeve\
l-- }\n\n&print_debug_message( 'MAIN', 'LWP::VERSI\
ON: ' . $LWP::VERSION,\n	1 );\n\n&print_debug_mess\
age( 'MAIN', \"params:\\n\" . Dumper( \\%params ),\
           11 );\n&print_debug_message( 'MAIN', \"\
tool_params:\\n\" . Dumper( \\%tool_params ), 11 )\
;\n\nmy $ua;\n\nmy $scriptName = basename( $0, () \
);\n\nif ( $params{'help'} || $numOpts == 0 ) {\n	\
&usage();\n	exit(0);\n}\n\n&print_debug_message( '\
MAIN', 'baseUrl: ' . $baseUrl, 1 );\n\nif (\n	!(\n\
		   $params{'polljob'}\n		|| $params{'resultTypes\
'}\n		|| $params{'status'}\n		|| $params{'params'}\
\n		|| $params{'paramDetail'}\n	)\n	&& !( defined(\
 $ARGV[0] ) || defined( $params{'sequence'} ) )\n \
 )\n{\n\n	# Bad argument combination, so print err\
or message and usage\n	print STDERR 'Error: bad op\
tion combination', \"\\n\";\n	&usage();\n	exit(1);\
\n}\n\nelsif ( $params{'params'} ) {\n	&print_tool\
_params();\n}\n\nelsif ( $params{'paramDetail'} ) \
{\n	&print_param_details( $params{'paramDetail'} )\
;\n}\n\nelsif ( $params{'status'} && defined( $par\
ams{'jobid'} ) ) {\n	&print_job_status( $params{'j\
obid'} );\n}\n\nelsif ( $params{'resultTypes'} && \
defined( $params{'jobid'} ) ) {\n	&print_result_ty\
pes( $params{'jobid'} );\n}\n\nelsif ( $params{'po\
lljob'} && defined( $params{'jobid'} ) ) {\n	&get_\
results( $params{'jobid'} );\n}\n\nelse {\n\n	# Mu\
ltiple input sequence mode, assume fasta format.\n\
	if ( $params{'multifasta'} ) {\n		&multi_submit_j\
ob();\n	}\n\n	# Entry identifier list file.\n	elsi\
f (( defined( $params{'sequence'} ) && $params{'se\
quence'} =~ m/^\\@/ )\n		|| ( defined( $ARGV[0] ) \
&& $ARGV[0] =~ m/^\\@/ ) )\n	{\n		my $list_filenam\
e = $params{'sequence'} || $ARGV[0];\n		$list_file\
name =~ s/^\\@//;\n		&list_file_submit_job($list_f\
ilename);\n	}\n\n	# Default: single sequence/ident\
ifier.\n	else {\n\n		# Load the sequence data and \
submit.\n		&submit_job( &load_data() );\n	}\n}\n\n\
=head1 FUNCTIONS\n\n=cut\n\n\n=head2 rest_user_age\
nt()\n\nGet a LWP UserAgent to use to perform REST\
 requests.\n\n  my $ua = &rest_user_agent();\n\n=c\
ut\n\nsub rest_user_agent() {\n	print_debug_messag\
e( 'rest_user_agent', 'Begin', 21 );\n	# Create an\
 LWP UserAgent for making HTTP calls.\n	my $ua = L\
WP::UserAgent->new();\n	# Set 'User-Agent' HTTP he\
ader to identifiy the client.\n	'$Revision: 2560 $\
' =~ m/(\\d+)/;\n	$ua->agent(\"EBI-Sample-Client/$\
1 ($scriptName; $OSNAME) \" . $ua->agent());\n	# C\
onfigure HTTP proxy support from environment.\n	$u\
a->env_proxy;\n	print_debug_message( 'rest_user_ag\
ent', 'End', 21 );\n	return $ua;\n}\n\n=head2 rest\
_error()\n\nCheck a REST response for an error con\
dition. An error is mapped to a die.\n\n  &rest_er\
ror($response, $content_data);\n\n=cut\n\nsub rest\
_error() {\n	print_debug_message( 'rest_error', 'B\
egin', 21 );\n	my $response = shift;\n	my $content\
data;\n	if(scalar(@_) > 0) {\n		$contentdata = shi\
ft;\n	}\n	if(!defined($contentdata) || $contentdat\
a eq '') {\n		$contentdata = $response->content();\
\n	}\n	# Check for HTTP error codes\n	if ( $respon\
se->is_error ) {\n		my $error_message = '';\n		# H\
TML response.\n		if(	$contentdata =~ m/<h1>([^<]+)\
<\\/h1>/ ) {\n			$error_message = $1;\n		}\n		#  X\
ML response.\n		elsif($contentdata =~ m/<descripti\
on>([^<]+)<\\/description>/) {\n			$error_message \
= $1;\n		}\n		die 'http status: ' . $response->cod\
e . ' ' . $response->message . '  ' . $error_messa\
ge;\n	}\n	print_debug_message( 'rest_error', 'End'\
, 21 );\n}\n\n=head2 rest_request()\n\nPerform a R\
EST request (HTTP GET).\n\n  my $response_str = &r\
est_request($url);\n\n=cut\n\nsub rest_request {\n\
	print_debug_message( 'rest_request', 'Begin', 11 \
);\n	my $requestUrl = shift;\n	print_debug_message\
( 'rest_request', 'URL: ' . $requestUrl, 11 );\n\n\
	# Get an LWP UserAgent.\n	$ua = &rest_user_agent(\
) unless defined($ua);\n	# Available HTTP compress\
ion methods.\n	my $can_accept;\n	eval {\n	    $can\
_accept = HTTP::Message::decodable();\n	};\n	$can_\
accept = '' unless defined($can_accept);\n	# Perfo\
rm the request\n	my $response = $ua->get($requestU\
rl,\n		'Accept-Encoding' => $can_accept, # HTTP co\
mpression.\n	);\n	print_debug_message( 'rest_reque\
st', 'HTTP status: ' . $response->code,\n		11 );\n\
	print_debug_message( 'rest_request',\n		'response\
 length: ' . length($response->content()), 11 );\n\
	print_debug_message( 'rest_request',\n		'request:\
' .\"\\n\" . $response->request()->as_string(), 32\
 );\n	print_debug_message( 'rest_request',\n		'res\
ponse: ' . \"\\n\" . $response->as_string(), 32 );\
\n	# Unpack possibly compressed response.\n	my $re\
tVal;\n	if ( defined($can_accept) && $can_accept n\
e '') {\n	    $retVal = $response->decoded_content\
();\n	}\n	# If unable to decode use orginal conten\
t.\n	$retVal = $response->content() unless defined\
($retVal);\n	# Check for an error.\n	&rest_error($\
response, $retVal);\n	print_debug_message( 'rest_r\
equest', 'retVal: ' . $retVal, 12 );\n	print_debug\
_message( 'rest_request', 'End', 11 );\n\n	# Retur\
n the response data\n	return $retVal;\n}\n\n=head2\
 rest_get_parameters()\n\nGet list of tool paramet\
er names.\n\n  my (@param_list) = &rest_get_parame\
ters();\n\n=cut\n\nsub rest_get_parameters {\n	pri\
nt_debug_message( 'rest_get_parameters', 'Begin', \
1 );\n	my $url                = $baseUrl . '/param\
eters/';\n	my $param_list_xml_str = rest_request($\
url);\n	my $param_list_xml     = XMLin($param_list\
_xml_str);\n	my (@param_list)       = @{ $param_li\
st_xml->{'id'} };\n	print_debug_message( 'rest_get\
_parameters', 'End', 1 );\n	return (@param_list);\\
n}\n\n=head2 rest_get_parameter_details()\n\nGet d\
etails of a tool parameter.\n\n  my $paramDetail =\
 &rest_get_parameter_details($param_name);\n\n=cut\
\n\nsub rest_get_parameter_details {\n	print_debug\
_message( 'rest_get_parameter_details', 'Begin', 1\
 );\n	my $parameterId = shift;\n	print_debug_messa\
ge( 'rest_get_parameter_details',\n		'parameterId:\
 ' . $parameterId, 1 );\n	my $url                 \
 = $baseUrl . '/parameterdetails/' . $parameterId;\
\n	my $param_detail_xml_str = rest_request($url);\\
n	my $param_detail_xml     = XMLin($param_detail_x\
ml_str);\n	print_debug_message( 'rest_get_paramete\
r_details', 'End', 1 );\n	return ($param_detail_xm\
l);\n}\n\n=head2 rest_run()\n\nSubmit a job.\n\n  \
my $job_id = &rest_run($email, $title, \\%params )\
;\n\n=cut\n\nsub rest_run {\n	print_debug_message(\
 'rest_run', 'Begin', 1 );\n	my $email  = shift;\n\
	my $title  = shift;\n	my $params = shift;\n	print\
_debug_message( 'rest_run', 'email: ' . $email, 1 \
);\n	if ( defined($title) ) {\n		print_debug_messa\
ge( 'rest_run', 'title: ' . $title, 1 );\n	}\n	pri\
nt_debug_message( 'rest_run', 'params: ' . Dumper(\
$params), 1 );\n\n	# Get an LWP UserAgent.\n	$ua =\
 &rest_user_agent() unless defined($ua);\n\n	# Cle\
an up parameters\n	my (%tmp_params) = %{$params};\\
n	$tmp_params{'email'} = $email;\n	$tmp_params{'ti\
tle'} = $title;\n	foreach my $param_name ( keys(%t\
mp_params) ) {\n		if ( !defined( $tmp_params{$para\
m_name} ) ) {\n			delete $tmp_params{$param_name};\
\n		}\n	}\n\n	# Submit the job as a POST\n	my $url\
 = $baseUrl . '/run';\n	my $response = $ua->post( \
$url, \\%tmp_params );\n	print_debug_message( 'res\
t_run', 'HTTP status: ' . $response->code, 11 );\n\
	print_debug_message( 'rest_run',\n		'request:' .\\
"\\n\" . $response->request()->as_string(), 11 );\\
n	print_debug_message( 'rest_run',\n		'response: '\
 . length($response->as_string()) . \"\\n\" . $res\
ponse->as_string(), 11 );\n\n	# Check for an error\
.\n	&rest_error($response);\n\n	# The job id is re\
turned\n	my $job_id = $response->content();\n	prin\
t_debug_message( 'rest_run', 'End', 1 );\n	return \
$job_id;\n}\n\n=head2 rest_get_status()\n\nCheck t\
he status of a job.\n\n  my $status = &rest_get_st\
atus($job_id);\n\n=cut\n\nsub rest_get_status {\n	\
print_debug_message( 'rest_get_status', 'Begin', 1\
 );\n	my $job_id = shift;\n	print_debug_message( '\
rest_get_status', 'jobid: ' . $job_id, 2 );\n	my $\
status_str = 'UNKNOWN';\n	my $url        = $baseUr\
l . '/status/' . $job_id;\n	$status_str = &rest_re\
quest($url);\n	print_debug_message( 'rest_get_stat\
us', 'status_str: ' . $status_str, 2 );\n	print_de\
bug_message( 'rest_get_status', 'End', 1 );\n	retu\
rn $status_str;\n}\n\n=head2 rest_get_result_types\
()\n\nGet list of result types for finished job.\n\
\n  my (@result_types) = &rest_get_result_types($j\
ob_id);\n\n=cut\n\nsub rest_get_result_types {\n	p\
rint_debug_message( 'rest_get_result_types', 'Begi\
n', 1 );\n	my $job_id = shift;\n	print_debug_messa\
ge( 'rest_get_result_types', 'jobid: ' . $job_id, \
2 );\n	my (@resultTypes);\n	my $url               \
       = $baseUrl . '/resulttypes/' . $job_id;\n	m\
y $result_type_list_xml_str = &rest_request($url);\
\n	my $result_type_list_xml     = XMLin($result_ty\
pe_list_xml_str);\n	(@resultTypes) = @{ $result_ty\
pe_list_xml->{'type'} };\n	print_debug_message( 'r\
est_get_result_types',\n		scalar(@resultTypes) . '\
 result types', 2 );\n	print_debug_message( 'rest_\
get_result_types', 'End', 1 );\n	return (@resultTy\
pes);\n}\n\n=head2 rest_get_result()\n\nGet result\
 data of a specified type for a finished job.\n\n \
 my $result = rest_get_result($job_id, $result_typ\
e);\n\n=cut\n\nsub rest_get_result {\n	print_debug\
_message( 'rest_get_result', 'Begin', 1 );\n	my $j\
ob_id = shift;\n	my $type   = shift;\n	print_debug\
_message( 'rest_get_result', 'jobid: ' . $job_id, \
1 );\n	print_debug_message( 'rest_get_result', 'ty\
pe: ' . $type,    1 );\n	my $url    = $baseUrl . '\
/result/' . $job_id . '/' . $type;\n	my $result = \
&rest_request($url);\n	print_debug_message( 'rest_\
get_result', length($result) . ' characters',\n		1\
 );\n	print_debug_message( 'rest_get_result', 'End\
', 1 );\n	return $result;\n}\n\n\n=head2 print_deb\
ug_message()\n\nPrint debug message at specified d\
ebug level.\n\n  &print_debug_message($method_name\
, $message, $level);\n\n=cut\n\nsub print_debug_me\
ssage {\n	my $function_name = shift;\n	my $message\
       = shift;\n	my $level         = shift;\n	if \
( $level <= $params{'debugLevel'} ) {\n		print STD\
ERR '[', $function_name, '()] ', $message, \"\\n\"\
;\n	}\n}\n\n=head2 print_tool_params()\n\nPrint li\
st of tool parameters.\n\n  &print_tool_params();\\
n\n=cut\n\nsub print_tool_params {\n	print_debug_m\
essage( 'print_tool_params', 'Begin', 1 );\n	my (@\
param_list) = &rest_get_parameters();\n	foreach my\
 $param ( sort(@param_list) ) {\n		print $param, \\
"\\n\";\n	}\n	print_debug_message( 'print_tool_par\
ams', 'End', 1 );\n}\n\n=head2 print_param_details\
()\n\nPrint details of a tool parameter.\n\n  &pri\
nt_param_details($param_name);\n\n=cut\n\nsub prin\
t_param_details {\n	print_debug_message( 'print_pa\
ram_details', 'Begin', 1 );\n	my $paramName = shif\
t;\n	print_debug_message( 'print_param_details', '\
paramName: ' . $paramName, 2 );\n	my $paramDetail \
= &rest_get_parameter_details($paramName);\n	print\
 $paramDetail->{'name'}, \"\\t\", $paramDetail->{'\
type'}, \"\\n\";\n	print $paramDetail->{'descripti\
on'}, \"\\n\";\n	if(defined($paramDetail->{'values\
'}->{'value'})) {\n		if(ref($paramDetail->{'values\
'}->{'value'}) eq 'ARRAY') {\n			foreach my $value\
 ( @{ $paramDetail->{'values'}->{'value'} } ) {\n	\
			&print_param_value($value);\n			}\n		}\n		else \
{\n				&print_param_value($paramDetail->{'values'}\
->{'value'});\n		}\n	}\n	print_debug_message( 'pri\
nt_param_details', 'End', 1 );\n}\n\n=head2 print_\
param_value()\n\nPrint details of a tool parameter\
 value.\n\n  &print_param_details($param_value);\n\
\nUsed by print_param_details() to handle both sin\
gluar and array values.\n\n=cut\n\nsub print_param\
_value {\n	my $value = shift;\n	print $value->{'va\
lue'};\n	if ( $value->{'defaultValue'} eq 'true' )\
 {\n		print \"\\t\", 'default';\n	}\n	print \"\\n\\
";\n	print \"\\t\", $value->{'label'}, \"\\n\";\n	\
if ( defined( $value->{'properties'} ) ) {\n		fore\
ach\n		  my $key ( sort( keys( %{ $value->{'proper\
ties'}{'property'} } ) ) )\n		{\n			if ( ref( $val\
ue->{'properties'}{'property'}{$key} ) eq 'HASH'\n\
				&& defined( $value->{'properties'}{'property'}\
{$key}{'value'} )\n			  )\n			{\n				print \"\\t\"\
, $key, \"\\t\",\n				  $value->{'properties'}{'pr\
operty'}{$key}{'value'}, \"\\n\";\n			}\n			else {\
\n				print \"\\t\", $value->{'properties'}{'prope\
rty'}{'key'},\n				  \"\\t\", $value->{'properties\
'}{'property'}{'value'}, \"\\n\";\n				last;\n			}\
\n		}\n	}\n}\n\n=head2 print_job_status()\n\nPrint\
 status of a job.\n\n  &print_job_status($job_id);\
\n\n=cut\n\nsub print_job_status {\n	print_debug_m\
essage( 'print_job_status', 'Begin', 1 );\n	my $jo\
bid = shift;\n	print_debug_message( 'print_job_sta\
tus', 'jobid: ' . $jobid, 1 );\n	if ( $outputLevel\
 > 0 ) {\n		print STDERR 'Getting status for job '\
, $jobid, \"\\n\";\n	}\n	my $result = &rest_get_st\
atus($jobid);\n	print \"$result\\n\";\n	if ( $resu\
lt eq 'FINISHED' && $outputLevel > 0 ) {\n		print \
STDERR \"To get results: $scriptName --polljob --j\
obid \" . $jobid\n		  . \"\\n\";\n	}\n	print_debug\
_message( 'print_job_status', 'End', 1 );\n}\n\n=h\
ead2 print_result_types()\n\nPrint available resul\
t types for a job.\n\n  &print_result_types($job_i\
d);\n\n=cut\n\nsub print_result_types {\n	print_de\
bug_message( 'result_types', 'Begin', 1 );\n	my $j\
obid = shift;\n	print_debug_message( 'result_types\
', 'jobid: ' . $jobid, 1 );\n	if ( $outputLevel > \
0 ) {\n		print STDERR 'Getting result types for jo\
b ', $jobid, \"\\n\";\n	}\n	my $status = &rest_get\
_status($jobid);\n	if ( $status eq 'PENDING' || $s\
tatus eq 'RUNNING' ) {\n		print STDERR 'Error: Job\
 status is ', $status,\n		  '. To get result types\
 the job must be finished.', \"\\n\";\n	}\n	else {\
\n		my (@resultTypes) = &rest_get_result_types($jo\
bid);\n		if ( $outputLevel > 0 ) {\n			print STDOU\
T 'Available result types:', \"\\n\";\n		}\n		fore\
ach my $resultType (@resultTypes) {\n			print STDO\
UT $resultType->{'identifier'}, \"\\n\";\n			if ( \
defined( $resultType->{'label'} ) ) {\n				print S\
TDOUT \"\\t\", $resultType->{'label'}, \"\\n\";\n	\
		}\n			if ( defined( $resultType->{'description'}\
 ) ) {\n				print STDOUT \"\\t\", $resultType->{'d\
escription'}, \"\\n\";\n			}\n			if ( defined( $re\
sultType->{'mediaType'} ) ) {\n				print STDOUT \"\
\\t\", $resultType->{'mediaType'}, \"\\n\";\n			}\\
n			if ( defined( $resultType->{'fileSuffix'} ) ) \
{\n				print STDOUT \"\\t\", $resultType->{'fileSu\
ffix'}, \"\\n\";\n			}\n		}\n		if ( $status eq 'FI\
NISHED' && $outputLevel > 0 ) {\n			print STDERR \\
"\\n\", 'To get results:', \"\\n\",\n			  \"  $scr\
iptName --polljob --jobid \" . $params{'jobid'} . \
\"\\n\",\n			  \"  $scriptName --polljob --outform\
at <type> --jobid \"\n			  . $params{'jobid'} . \"\
\\n\";\n		}\n	}\n	print_debug_message( 'result_typ\
es', 'End', 1 );\n}\n\n=head2 submit_job()\n\nSubm\
it a job to the service.\n\n  &submit_job($seq);\n\
\n=cut\n\nsub submit_job {\n	print_debug_message( \
'submit_job', 'Begin', 1 );\n\n	# Set input sequen\
ce\n	$tool_params{'sequence'} = shift;\n\n	# Load \
parameters\n	&load_params();\n\n	# Submit the job\\
n	my $jobid = &rest_run( $params{'email'}, $params\
{'title'}, \\%tool_params );\n\n	# Simulate sync/a\
sync mode\n	if ( defined( $params{'async'} ) ) {\n\
		print STDOUT $jobid, \"\\n\";\n		if ( $outputLev\
el > 0 ) {\n			print STDERR\n			  \"To check statu\
s: $scriptName --status --jobid $jobid\\n\";\n		}\\
n	}\n	else {\n		if ( $outputLevel > 0 ) {\n			prin\
t STDERR \"JobId: $jobid\\n\";\n		}\n		sleep 1;\n	\
	&get_results($jobid);\n	}\n	print_debug_message( \
'submit_job', 'End', 1 );\n}\n\n=head2 multi_submi\
t_job()\n\nSubmit multiple jobs assuming input is \
a collection of fasta formatted sequences.\n\n  &m\
ulti_submit_job();\n\n=cut\n\nsub multi_submit_job\
 {\n	print_debug_message( 'multi_submit_job', 'Beg\
in', 1 );\n	my $jobIdForFilename = 1;\n	$jobIdForF\
ilename = 0 if ( defined( $params{'outfile'} ) );\\
n	my (@filename_list) = ();\n\n	# Query sequence\n\
	if ( defined( $ARGV[0] ) ) {    # Bare option\n		\
if ( -f $ARGV[0] || $ARGV[0] eq '-' ) {    # File\\
n			push( @filename_list, $ARGV[0] );\n		}\n		else\
 {\n			warn 'Warning: Input file \"' . $ARGV[0] . \
'\" does not exist'\n		}\n	}\n	if ( $params{'seque\
nce'} ) {                   # Via --sequence\n		if\
 ( -f $params{'sequence'} || $params{'sequence'} e\
q '-' ) {    # File\n			push( @filename_list, $par\
ams{'sequence'} );\n		}\n		else {\n			warn 'Warnin\
g: Input file \"' . $params{'sequence'} . '\" does\
 not exist'\n		}\n	}\n\n	$/ = '>';\n	foreach my $f\
ilename (@filename_list) {\n		my $INFILE;\n		if($f\
ilename eq '-') { # STDIN.\n			open( $INFILE, '<-'\
 )\n			  or die 'Error: unable to STDIN (' . $! . \
')';\n		} else { # File.\n			open( $INFILE, '<', $\
filename )\n			  or die 'Error: unable to open fil\
e ' . $filename . ' (' . $! . ')';\n		}\n		while (\
<$INFILE>) {\n			my $seq = $_;\n			$seq =~ s/>$//;\
\n			if ( $seq =~ m/(\\S+)/ ) {\n				print STDERR \
\"Submitting job for: $1\\n\"\n				  if ( $outputL\
evel > 0 );\n				$seq = '>' . $seq;\n				&print_de\
bug_message( 'multi_submit_job', $seq, 11 );\n				\
&submit_job($seq);\n				$params{'outfile'} = undef\
 if ( $jobIdForFilename == 1 );\n			}\n		}\n		clos\
e $INFILE;\n	}\n	print_debug_message( 'multi_submi\
t_job', 'End', 1 );\n}\n\n=head2 list_file_submit_\
job()\n\nSubmit multiple jobs using a file contain\
ing a list of entry identifiers as \ninput.\n\n  &\
list_file_submit_job($list_filename)\n\n=cut\n\nsu\
b list_file_submit_job {\n	my $filename         = \
shift;\n	my $jobIdForFilename = 1;\n	$jobIdForFile\
name = 0 if ( defined( $params{'outfile'} ) );\n\n\
	# Iterate over identifiers, submitting each job\n\
	my $LISTFILE;\n	if($filename eq '-') { # STDIN.\n\
		open( $LISTFILE, '<-' )\n		  or die 'Error: unab\
le to STDIN (' . $! . ')';\n	} else { # File.\n		o\
pen( $LISTFILE, '<', $filename )\n		  or die 'Erro\
r: unable to open file ' . $filename . ' (' . $! .\
 ')';\n	}\n	while (<$LISTFILE>) {\n		my $line = $_\
;\n		chomp($line);\n		if ( $line ne '' ) {\n			&pr\
int_debug_message( 'list_file_submit_job', 'line: \
' . $line, 2 );\n			if ( $line =~ m/\\w:\\w/ ) {  \
  # Check this is an identifier\n				print STDERR \
\"Submitting job for: $line\\n\"\n				  if ( $outp\
utLevel > 0 );\n				&submit_job($line);\n			}\n			\
else {\n				print STDERR\n\"Warning: line \\\"$lin\
e\\\" is not recognised as an identifier\\n\";\n		\
	}\n		}\n		$params{'outfile'} = undef if ( $jobIdF\
orFilename == 1 );\n	}\n	close $LISTFILE;\n}\n\n=h\
ead2 load_data()\n\nLoad sequence data from file o\
r option specified on the command-line.\n\n  &load\
_data();\n\n=cut\n\nsub load_data {\n	print_debug_\
message( 'load_data', 'Begin', 1 );\n	my $retSeq;\\
n\n	# Query sequence\n	if ( defined( $ARGV[0] ) ) \
{    # Bare option\n		if ( -f $ARGV[0] || $ARGV[0]\
 eq '-' ) {    # File\n			$retSeq = &read_file( $A\
RGV[0] );\n		}\n		else {                          \
           # DB:ID or sequence\n			$retSeq = $ARGV\
[0];\n		}\n	}\n	if ( $params{'sequence'} ) {      \
             # Via --sequence\n		if ( -f $params{'\
sequence'} || $params{'sequence'} eq '-' ) {    # \
File\n			$retSeq = &read_file( $params{'sequence'}\
 );\n		}\n		else {    # DB:ID or sequence\n			$ret\
Seq = $params{'sequence'};\n		}\n	}\n	print_debug_\
message( 'load_data', 'End', 1 );\n	return $retSeq\
;\n}\n\n=head2 load_params()\n\nLoad job parameter\
s from command-line options.\n\n  &load_params();\\
n\n=cut\n\nsub load_params {\n	print_debug_message\
( 'load_params', 'Begin', 1 );\n\n	# Database(s) t\
o search\n	my (@dbList) = split /[ ,]/, $params{'d\
atabase'};\n	$tool_params{'database'} = \\@dbList;\
\n\n	# Match/missmatch\n	if ( $params{'match'} && \
$params{'missmatch'} ) {\n		$tool_params{'match_sc\
ores'} =\n		  $params{'match'} . ',' . $params{'mi\
ssmatch'};\n	}\n	\n	# Compatability options, old c\
ommand-line\n	if(!$tool_params{'alignments'} && $p\
arams{'numal'}) {\n		$tool_params{'alignments'} = \
$params{'numal'};\n	}\n	if(!$tool_params{'gapopen'\
} && $params{'opengap'}) {\n		$tool_params{'gapope\
n'} = $params{'opengap'};\n	}\n	if(!$tool_params{'\
gapext'} && $params{'extendgap'}) {\n		$tool_param\
s{'gapext'} = $params{'extendgap'};\n	}\n\n	print_\
debug_message( 'load_params', 'End', 1 );\n}\n\n=h\
ead2 client_poll()\n\nClient-side job polling.\n\n\
  &client_poll($job_id);\n\n=cut\n\nsub client_pol\
l {\n	print_debug_message( 'client_poll', 'Begin',\
 1 );\n	my $jobid  = shift;\n	my $status = 'PENDIN\
G';\n\n	my $errorCount = 0;\n	while ($status eq 'R\
UNNING'\n		|| $status eq 'PENDING'\n		|| ( $status\
 eq 'ERROR' && $errorCount < 2 ) )\n	{\n		$status \
= rest_get_status($jobid);\n		print STDERR \"$stat\
us\\n\" if ( $outputLevel > 0 );\n		if ( $status e\
q 'ERROR' ) {\n			$errorCount++;\n		}\n		elsif ( $\
errorCount > 0 ) {\n			$errorCount--;\n		}\n		if (\
   $status eq 'RUNNING'\n			|| $status eq 'PENDING\
'\n			|| $status eq 'ERROR' )\n		{\n\n			# Wait be\
fore polling again.\n			sleep $checkInterval;\n		}\
\n	}\n	print_debug_message( 'client_poll', 'End', \
1 );\n	return $status;\n}\n\n=head2 get_results()\\
n\nGet the results for a job identifier.\n\n  &get\
_results($job_id);\n\n=cut\n\nsub get_results {\n	\
print_debug_message( 'get_results', 'Begin', 1 );\\
n	my $jobid = shift;\n	print_debug_message( 'get_r\
esults', 'jobid: ' . $jobid, 1 );\n\n	# Verbose\n	\
if ( $outputLevel > 1 ) {\n		print 'Getting result\
s for job ', $jobid, \"\\n\";\n	}\n\n	# Check stat\
us, and wait if not finished\n	client_poll($jobid)\
;\n\n	# Use JobId if output file name is not defin\
ed\n	unless ( defined( $params{'outfile'} ) ) {\n	\
	$params{'outfile'} = $jobid;\n	}\n\n	# Get list o\
f data types\n	my (@resultTypes) = rest_get_result\
_types($jobid);\n\n	# Get the data and write it to\
 a file\n	if ( defined( $params{'outformat'} ) ) {\
    # Specified data type\n		my $selResultType;\n	\
	foreach my $resultType (@resultTypes) {\n			if ( \
$resultType->{'identifier'} eq $params{'outformat'\
} ) {\n				$selResultType = $resultType;\n			}\n		\
}\n		if ( defined($selResultType) ) {\n			my $resu\
lt =\n			  rest_get_result( $jobid, $selResultType\
->{'identifier'} );\n			if ( $params{'outfile'} eq\
 '-' ) {\n				write_file( $params{'outfile'}, $res\
ult );\n			}\n			else {\n				write_file(\n					$pa\
rams{'outfile'} . '.'\n					  . $selResultType->{'\
identifier'} . '.'\n					  . $selResultType->{'fil\
eSuffix'},\n					$result\n				);\n			}\n		}\n		els\
e {\n			die 'Error: unknown result format \"' . $p\
arams{'outformat'} . '\"';\n		}\n	}\n	else {    # \
Data types available\n		      # Write a file for e\
ach output type\n		for my $resultType (@resultType\
s) {\n			if ( $outputLevel > 1 ) {\n				print STDE\
RR 'Getting ', $resultType->{'identifier'}, \"\\n\\
";\n			}\n			my $result = rest_get_result( $jobid,\
 $resultType->{'identifier'} );\n			if ( $params{'\
outfile'} eq '-' ) {\n				write_file( $params{'out\
file'}, $result );\n			}\n			else {\n				write_fil\
e(\n					$params{'outfile'} . '.'\n					  . $resul\
tType->{'identifier'} . '.'\n					  . $resultType-\
>{'fileSuffix'},\n					$result\n				);\n			}\n		}\\
n	}\n	print_debug_message( 'get_results', 'End', 1\
 );\n}\n\n=head2 read_file()\n\nRead a file into a\
 scalar. The special filename '-' can be used to r\
ead from \nstandard input (STDIN).\n\n  my $data =\
 &read_file($filename);\n\n=cut\n\nsub read_file {\
\n	print_debug_message( 'read_file', 'Begin', 1 );\
\n	my $filename = shift;\n	print_debug_message( 'r\
ead_file', 'filename: ' . $filename, 2 );\n	my ( $\
content, $buffer );\n	if ( $filename eq '-' ) {\n	\
	while ( sysread( STDIN, $buffer, 1024 ) ) {\n			$\
content .= $buffer;\n		}\n	}\n	else {    # File\n	\
	open( my $FILE, '<', $filename )\n		  or die \"Er\
ror: unable to open input file $filename ($!)\";\n\
		while ( sysread( $FILE, $buffer, 1024 ) ) {\n			\
$content .= $buffer;\n		}\n		close($FILE);\n	}\n	p\
rint_debug_message( 'read_file', 'End', 1 );\n	ret\
urn $content;\n}\n\n=head2 write_file()\n\nWrite d\
ata to a file. The special filename '-' can be use\
d to write to \nstandard output (STDOUT).\n\n  &wr\
ite_file($filename, $data);\n\n=cut\n\nsub write_f\
ile {\n	print_debug_message( 'write_file', 'Begin'\
, 1 );\n	my ( $filename, $data ) = @_;\n	print_deb\
ug_message( 'write_file', 'filename: ' . $filename\
, 2 );\n	if ( $outputLevel > 0 ) {\n		print STDERR\
 'Creating result file: ' . $filename . \"\\n\";\n\
	}\n	if ( $filename eq '-' ) {\n		print STDOUT $da\
ta;\n	}\n	else {\n		open( my $FILE, '>', $filename\
 )\n		  or die \"Error: unable to open output file\
 $filename ($!)\";\n		syswrite( $FILE, $data );\n	\
	close($FILE);\n	}\n	print_debug_message( 'write_f\
ile', 'End', 1 );\n}\n\n=head2 usage()\n\nPrint pr\
ogram usage message.\n\n  &usage();\n\n=cut\n\nsub\
 usage {\n	print STDERR <<EOF\nNCBI BLAST\n=======\
===\n   \nRapid sequence database search programs \
utilizing the BLAST algorithm\n    \n[Required]\n\\
n  -p, --program      : str  : BLAST program to us\
e, see --paramDetail program\n  -D, --database    \
 : str  : database(s) to search, space separated. \
See\n                              --paramDetail d\
atabase\n      --stype        : str  : query seque\
nce type, see --paramDetail stype\n  seqFile      \
      : file : query sequence (\"-\" for STDIN, \\\
@filename for\n                              ident\
ifier list file)\n\n[Optional]\n\n  -m, --matrix  \
     : str  : scoring matrix, see --paramDetail ma\
trix\n  -e, --exp          : real : 0<E<= 1000. St\
atistical significance threshold \n               \
               for reporting database sequence mat\
ches.\n  -f, --filter       :      : filter the qu\
ery sequence for low complexity \n                \
              regions, see --paramDetail filter\n \
 -A, --align        : int  : pairwise alignment fo\
rmat, see --paramDetail align\n  -s, --scores     \
  : int  : number of scores to be reported\n  -n, \
--alignments   : int  : number of alignments to re\
port\n  -u, --match        : int  : Match score (B\
LASTN only)\n  -v, --mismatch     : int  : Mismatc\
h score (BLASTN only)\n  -o, --gapopen      : int \
 : Gap open penalty\n  -x, --gapext       : int  :\
 Gap extension penalty\n  -d, --dropoff      : int\
  : Drop-off\n  -g, --gapalign     :      : Optimi\
se gapped alignments\n      --seqrange     : str  \
: region within input to use as query\n      --mul\
tifasta   :      : treat input as a set of fasta f\
ormatted sequences\n\n[General]\n\n  -h, --help   \
      :      : prints this help text\n      --asyn\
c        :      : forces to make an asynchronous q\
uery\n      --email        : str  : e-mail address\
\n      --title        : str  : title for job\n   \
   --status       :      : get job status\n      -\
-resultTypes  :      : get available result types \
for job\n      --polljob      :      : poll for th\
e status of a job\n      --jobid        : str  : j\
obid that was returned when an asynchronous job \n\
                              was submitted.\n    \
  --outfile      : str  : file name for results (d\
efault is jobid;\n                              \"\
-\" for STDOUT)\n      --outformat    : str  : res\
ult format to retrieve\n      --params       :    \
  : list input parameters\n      --paramDetail  : \
str  : display details for input parameter\n      \
--quiet        :      : decrease output\n      --v\
erbose      :      : increase output\n   \nSynchro\
nous job:\n\n  The results/errors are returned as \
soon as the job is finished.\n  Usage: $scriptName\
 --email <your\\@email> [options...] seqFile\n  Re\
turns: results as an attachment\n\nAsynchronous jo\
b:\n\n  Use this if you want to retrieve the resul\
ts at a later time. The results \n  are stored for\
 up to 24 hours. 	\n  Usage: $scriptName --async -\
-email <your\\@email> [options...] seqFile\n  Retu\
rns: jobid\n\n  Use the jobid to query for the sta\
tus of the job. If the job is finished, \n  it als\
o returns the results/errors.\n  Usage: $scriptNam\
e --polljob --jobid <jobId> [--outfile string]\n  \
Returns: string indicating the status of the job a\
nd if applicable, results \n  as an attachment.\n\\
nFurther information:\n\n  http://www.ebi.ac.uk/To\
ols/webservices/services/sss/ncbi_blast_rest\n  ht\
tp://www.ebi.ac.uk/Tools/webservices/tutorials/per\
l\n\nSupport/Feedback:\n\n  http://www.ebi.ac.uk/s\
upport/\nEOF\n}\n\n=head1 FEEDBACK/SUPPORT\n\nPlea\
se contact us at L<http://www.ebi.ac.uk/support/> \
if you have any \nfeedback, suggestions or issues \
with the service or this client.\n\n=cut\n","\n=he\
ad1 NAME\n\nwublast_lwp.pl\n\n=head1 DESCRIPTION\n\
\nWU-BLAST (REST) web service Perl client using L<\
LWP>.\n\nTested with:\n\n=over\n\n=item *\nL<LWP> \
5.79, L<XML::Simple> 2.12 and Perl 5.8.3\n\n=item \
*\nL<LWP> 5.808, L<XML::Simple> 2.18 and Perl 5.8.\
8 (Ubuntu 8.04 LTS)\n\n=item *\nL<LWP> 5.834, L<XM\
L::Simple> 2.18 and Perl 5.10.1 (Ubuntu 10.04 LTS)\
\n\n=item *\nL<LWP> 6.03, L<XML::Simple> 2.18 and \
Perl 5.14.2 (Ubuntu 12.04 LTS)\n\n=back\n\nFor fur\
ther information see:\n\n=over\n\n=item *\nL<http:\
//www.ebi.ac.uk/Tools/webservices/services/sss/wu_\
blast_rest>\n\n=item *\nL<http://www.ebi.ac.uk/Too\
ls/webservices/tutorials/perl>\n\n=back\n\n=head1 \
LICENSE\n\nCopyright 2012-2013 EMBL - European Bio\
informatics Institute\n\nLicensed under the Apache\
 License, Version 2.0 (the \"License\");\nyou may \
not use this file except in compliance with the Li\
cense.\nYou may obtain a copy of the License at\n\\
n    http://www.apache.org/licenses/LICENSE-2.0\n\\
nUnless required by applicable law or agreed to in\
 writing, software\ndistributed under the License \
is distributed on an \"AS IS\" BASIS,\nWITHOUT WAR\
RANTIES OR CONDITIONS OF ANY KIND, either express \
or implied.\nSee the License for the specific lang\
uage governing permissions and\nlimitations under \
the License.\n\n=head1 VERSION\n\n$Id: wublast_lwp\
.pl 2560 2013-03-20 12:56:31Z hpm $\n\n=cut\n\nuse\
 strict;\nuse warnings;\n\nuse English;\nuse LWP;\\
nuse XML::Simple;\nuse Getopt::Long qw(:config no_\
ignore_case bundling);\nuse File::Basename;\nuse D\
ata::Dumper;\n\nmy $baseUrl = 'http://www.ebi.ac.u\
k/Tools/services/rest/wublast';\n\nmy $checkInterv\
al = 3;\n\nmy $outputLevel = 1;\n\nmy $numOpts = s\
calar(@ARGV);\nmy %params = ( 'debugLevel' => 0 );\
\n\nmy %tool_params = ();\nGetOptions(\n\n	# Tool \
specific options\n	'program|p=s'     => \\$tool_pa\
rams{'program'},      # BLAST program\n	'database|\
D=s'    => \\$params{'database'},     # Search dat\
abase\n	'matrix|m=s'      => \\$tool_params{'matri\
x'},       # Scoring matrix\n	'exp|E=f'         =>\
 \\$tool_params{'exp'},          # E-value thresho\
ld\n	'viewfilter|e'    => \\$tool_params{'viewfilt\
er'},   # Display filtered sequence\n	'filter|f=s'\
      => \\$tool_params{'filter'},       # Low com\
plexity filter name\n	'alignments|n=i'  => \\$tool\
_params{'alignments'},   # Number of alignments\n	\
'scores|s=i'      => \\$tool_params{'scores'},    \
   # Number of scores\n	'sensitivity|S=s' => \\$to\
ol_params{'sensitivity'},  # Search sensitivity\n	\
'sort|t=s'        => \\$tool_params{'sort'},      \
   # Sort hits by...\n	'stats|T=s'       => \\$too\
l_params{'stats'},        # Scoring statistic to u\
se\n	'strand|d=s'      => \\$tool_params{'strand'}\
,       # Strand to use\n	'topcombon|c=i'   => \\$\
tool_params{'topcombon'},    # Consistent sets of \
HSPs\n	'align|A=i'       => \\$tool_params{'align'\
},   # Pairwise alignment format\n	'stype=s' => \\\
$tool_params{'stype'},    # Sequence type 'protein\
' or 'dna'\n	'sequence=s' => \\$params{'sequence'}\
,         # Query sequence file or DB:ID\n	'multif\
asta' => \\$params{'multifasta'},       # Multiple\
 fasta input\n\n	# Compatability options, old comm\
and-line.\n	'echofilter|e'    => \\$params{'echofi\
lter'},   # Display filtered sequence\n	'b=i'  => \
\\$params{'numal'},        # Number of alignments\\
n	'appxml=s'        => \\$params{'appxml'},       \
# Application XML\n\n	# Generic options\n	'email=s\
'       => \\$params{'email'},          # User e-m\
ail address\n	'title=s'       => \\$params{'title'\
},          # Job title\n	'outfile=s'     => \\$pa\
rams{'outfile'},        # Output file name\n	'outf\
ormat=s'   => \\$params{'outformat'},      # Outpu\
t file type\n	'jobid=s'       => \\$params{'jobid'\
},          # JobId\n	'help|h'        => \\$params\
{'help'},           # Usage help\n	'async'        \
 => \\$params{'async'},          # Asynchronous su\
bmission\n	'polljob'       => \\$params{'polljob'}\
,        # Get results\n	'resultTypes'   => \\$par\
ams{'resultTypes'},    # Get result types\n	'statu\
s'        => \\$params{'status'},         # Get st\
atus\n	'params'        => \\$params{'params'},    \
     # List input parameters\n	'paramDetail=s' => \
\\$params{'paramDetail'},    # Get details for par\
ameter\n	'quiet'         => \\$params{'quiet'},   \
       # Decrease output level\n	'verbose'       =\
> \\$params{'verbose'},        # Increase output l\
evel\n	'debugLevel=i'  => \\$params{'debugLevel'},\
     # Debug output level\n	'baseUrl=s'     => \\$\
baseUrl,                  # Base URL for service.\\
n);\nif ( $params{'verbose'} ) { $outputLevel++ }\\
nif ( $params{'quiet'} )  { $outputLevel-- }\n\n&p\
rint_debug_message( 'MAIN', 'LWP::VERSION: ' . $LW\
P::VERSION,\n	1 );\n\n&print_debug_message( 'MAIN'\
, \"params:\\n\" . Dumper( \\%params ),           \
11 );\n&print_debug_message( 'MAIN', \"tool_params\
:\\n\" . Dumper( \\%tool_params ), 11 );\n\nmy $ua\
;\n\nmy $scriptName = basename( $0, () );\n\nif ( \
$params{'help'} || $numOpts == 0 ) {\n	&usage();\n\
	exit(0);\n}\n\n&print_debug_message( 'MAIN', 'bas\
eUrl: ' . $baseUrl, 1 );\n\nif (\n	!(\n		   $param\
s{'polljob'}\n		|| $params{'resultTypes'}\n		|| $p\
arams{'status'}\n		|| $params{'params'}\n		|| $par\
ams{'paramDetail'}\n	)\n	&& !( defined( $ARGV[0] )\
 || defined( $params{'sequence'} ) )\n  )\n{\n\n	#\
 Bad argument combination, so print error message \
and usage\n	print STDERR 'Error: bad option combin\
ation', \"\\n\";\n	&usage();\n	exit(1);\n}\n\nelsi\
f ( $params{'params'} ) {\n	&print_tool_params();\\
n}\n\nelsif ( $params{'paramDetail'} ) {\n	&print_\
param_details( $params{'paramDetail'} );\n}\n\nels\
if ( $params{'status'} && defined( $params{'jobid'\
} ) ) {\n	&print_job_status( $params{'jobid'} );\n\
}\n\nelsif ( $params{'resultTypes'} && defined( $p\
arams{'jobid'} ) ) {\n	&print_result_types( $param\
s{'jobid'} );\n}\n\nelsif ( $params{'polljob'} && \
defined( $params{'jobid'} ) ) {\n	&get_results( $p\
arams{'jobid'} );\n}\n\nelse {\n\n	# Multiple inpu\
t sequence mode, assume fasta format.\n	if ( $para\
ms{'multifasta'} ) {\n		&multi_submit_job();\n	}\n\
\n	# Entry identifier list file.\n	elsif (( define\
d( $params{'sequence'} ) && $params{'sequence'} =~\
 m/^\\@/ )\n		|| ( defined( $ARGV[0] ) && $ARGV[0]\
 =~ m/^\\@/ ) )\n	{\n		my $list_filename = $params\
{'sequence'} || $ARGV[0];\n		$list_filename =~ s/^\
\\@//;\n		&list_file_submit_job($list_filename);\n\
	}\n\n	# Default: single sequence/identifier.\n	el\
se {\n\n		# Load the sequence data and submit.\n		\
&submit_job( &load_data() );\n	}\n}\n\n=head1 FUNC\
TIONS\n\n=cut\n\n\n=head2 rest_user_agent()\n\nGet\
 a LWP UserAgent to use to perform REST requests.\\
n\n  my $ua = &rest_user_agent();\n\n=cut\n\nsub r\
est_user_agent() {\n	print_debug_message( 'rest_us\
er_agent', 'Begin', 21 );\n	# Create an LWP UserAg\
ent for making HTTP calls.\n	my $ua = LWP::UserAge\
nt->new();\n	# Set 'User-Agent' HTTP header to ide\
ntifiy the client.\n	'$Revision: 2560 $' =~ m/(\\d\
+)/;\n	$ua->agent(\"EBI-Sample-Client/$1 ($scriptN\
ame; $OSNAME) \" . $ua->agent());\n	# Configure HT\
TP proxy support from environment.\n	$ua->env_prox\
y;\n	print_debug_message( 'rest_user_agent', 'End'\
, 21 );\n	return $ua;\n}\n\n=head2 rest_error()\n\\
nCheck a REST response for an error condition. An \
error is mapped to a die.\n\n  &rest_error($respon\
se, $content_data);\n\n=cut\n\nsub rest_error() {\\
n	print_debug_message( 'rest_error', 'Begin', 21 )\
;\n	my $response = shift;\n	my $contentdata;\n	if(\
scalar(@_) > 0) {\n		$contentdata = shift;\n	}\n	i\
f(!defined($contentdata) || $contentdata eq '') {\\
n		$contentdata = $response->content();\n	}\n	# Ch\
eck for HTTP error codes\n	if ( $response->is_erro\
r ) {\n		my $error_message = '';\n		# HTML respons\
e.\n		if(	$contentdata =~ m/<h1>([^<]+)<\\/h1>/ ) \
{\n			$error_message = $1;\n		}\n		#  XML response\
.\n		elsif($contentdata =~ m/<description>([^<]+)<\
\\/description>/) {\n			$error_message = $1;\n		}\\
n		die 'http status: ' . $response->code . ' ' . $\
response->message . '  ' . $error_message;\n	}\n	p\
rint_debug_message( 'rest_error', 'End', 21 );\n}\\
n\n=head2 rest_request()\n\nPerform a REST request\
 (HTTP GET).\n\n  my $response_str = &rest_request\
($url);\n\n=cut\n\nsub rest_request {\n	print_debu\
g_message( 'rest_request', 'Begin', 11 );\n	my $re\
questUrl = shift;\n	print_debug_message( 'rest_req\
uest', 'URL: ' . $requestUrl, 11 );\n\n	# Get an L\
WP UserAgent.\n	$ua = &rest_user_agent() unless de\
fined($ua);\n	# Available HTTP compression methods\
.\n	my $can_accept;\n	eval {\n	    $can_accept = H\
TTP::Message::decodable();\n	};\n	$can_accept = ''\
 unless defined($can_accept);\n	# Perform the requ\
est\n	my $response = $ua->get($requestUrl,\n		'Acc\
ept-Encoding' => $can_accept, # HTTP compression.\\
n	);\n	print_debug_message( 'rest_request', 'HTTP \
status: ' . $response->code,\n		11 );\n	print_debu\
g_message( 'rest_request',\n		'response length: ' \
. length($response->content()), 11 );\n	print_debu\
g_message( 'rest_request',\n		'request:' .\"\\n\" \
. $response->request()->as_string(), 32 );\n	print\
_debug_message( 'rest_request',\n		'response: ' . \
\"\\n\" . $response->as_string(), 32 );\n	# Unpack\
 possibly compressed response.\n	my $retVal;\n	if \
( defined($can_accept) && $can_accept ne '') {\n	 \
   $retVal = $response->decoded_content();\n	}\n	#\
 If unable to decode use orginal content.\n	$retVa\
l = $response->content() unless defined($retVal);\\
n	# Check for an error.\n	&rest_error($response, $\
retVal);\n	print_debug_message( 'rest_request', 'r\
etVal: ' . $retVal, 12 );\n	print_debug_message( '\
rest_request', 'End', 11 );\n\n	# Return the respo\
nse data\n	return $retVal;\n}\n\n=head2 rest_get_p\
arameters()\n\nGet list of tool parameter names.\n\
\n  my (@param_list) = &rest_get_parameters();\n\n\
=cut\n\nsub rest_get_parameters {\n	print_debug_me\
ssage( 'rest_get_parameters', 'Begin', 1 );\n	my $\
url                = $baseUrl . '/parameters/';\n	\
my $param_list_xml_str = rest_request($url);\n	my \
$param_list_xml     = XMLin($param_list_xml_str);\\
n	my (@param_list)       = @{ $param_list_xml->{'i\
d'} };\n	print_debug_message( 'rest_get_parameters\
', 'End', 1 );\n	return (@param_list);\n}\n\n=head\
2 rest_get_parameter_details()\n\nGet details of a\
 tool parameter.\n\n  my $paramDetail = &rest_get_\
parameter_details($param_name);\n\n=cut\n\nsub res\
t_get_parameter_details {\n	print_debug_message( '\
rest_get_parameter_details', 'Begin', 1 );\n	my $p\
arameterId = shift;\n	print_debug_message( 'rest_g\
et_parameter_details',\n		'parameterId: ' . $param\
eterId, 1 );\n	my $url                  = $baseUrl\
 . '/parameterdetails/' . $parameterId;\n	my $para\
m_detail_xml_str = rest_request($url);\n	my $param\
_detail_xml     = XMLin($param_detail_xml_str);\n	\
print_debug_message( 'rest_get_parameter_details',\
 'End', 1 );\n	return ($param_detail_xml);\n}\n\n=\
head2 rest_run()\n\nSubmit a job.\n\n  my $job_id \
= &rest_run($email, $title, \\%params );\n\n=cut\n\
\nsub rest_run {\n	print_debug_message( 'rest_run'\
, 'Begin', 1 );\n	my $email  = shift;\n	my $title \
 = shift;\n	my $params = shift;\n	print_debug_mess\
age( 'rest_run', 'email: ' . $email, 1 );\n	if ( d\
efined($title) ) {\n		print_debug_message( 'rest_r\
un', 'title: ' . $title, 1 );\n	}\n	print_debug_me\
ssage( 'rest_run', 'params: ' . Dumper($params), 1\
 );\n\n	# Get an LWP UserAgent.\n	$ua = &rest_user\
_agent() unless defined($ua);\n\n	# Clean up param\
eters\n	my (%tmp_params) = %{$params};\n	$tmp_para\
ms{'email'} = $email;\n	$tmp_params{'title'} = $ti\
tle;\n	foreach my $param_name ( keys(%tmp_params) \
) {\n		if ( !defined( $tmp_params{$param_name} ) )\
 {\n			delete $tmp_params{$param_name};\n		}\n	}\n\
\n	# Submit the job as a POST\n	my $url = $baseUrl\
 . '/run';\n	my $response = $ua->post( $url, \\%tm\
p_params );\n	print_debug_message( 'rest_run', 'HT\
TP status: ' . $response->code, 11 );\n	print_debu\
g_message( 'rest_run',\n		'request:' .\"\\n\" . $r\
esponse->request()->as_string(), 11 );\n	print_deb\
ug_message( 'rest_run',\n		'response: ' . length($\
response->as_string()) . \"\\n\" . $response->as_s\
tring(), 11 );\n\n	# Check for an error.\n	&rest_e\
rror($response);\n\n	# The job id is returned\n	my\
 $job_id = $response->content();\n	print_debug_mes\
sage( 'rest_run', 'End', 1 );\n	return $job_id;\n}\
\n\n=head2 rest_get_status()\n\nCheck the status o\
f a job.\n\n  my $status = &rest_get_status($job_i\
d);\n\n=cut\n\nsub rest_get_status {\n	print_debug\
_message( 'rest_get_status', 'Begin', 1 );\n	my $j\
ob_id = shift;\n	print_debug_message( 'rest_get_st\
atus', 'jobid: ' . $job_id, 2 );\n	my $status_str \
= 'UNKNOWN';\n	my $url        = $baseUrl . '/statu\
s/' . $job_id;\n	$status_str = &rest_request($url)\
;\n	print_debug_message( 'rest_get_status', 'statu\
s_str: ' . $status_str, 2 );\n	print_debug_message\
( 'rest_get_status', 'End', 1 );\n	return $status_\
str;\n}\n\n=head2 rest_get_result_types()\n\nGet l\
ist of result types for finished job.\n\n  my (@re\
sult_types) = &rest_get_result_types($job_id);\n\n\
=cut\n\nsub rest_get_result_types {\n	print_debug_\
message( 'rest_get_result_types', 'Begin', 1 );\n	\
my $job_id = shift;\n	print_debug_message( 'rest_g\
et_result_types', 'jobid: ' . $job_id, 2 );\n	my (\
@resultTypes);\n	my $url                      = $b\
aseUrl . '/resulttypes/' . $job_id;\n	my $result_t\
ype_list_xml_str = &rest_request($url);\n	my $resu\
lt_type_list_xml     = XMLin($result_type_list_xml\
_str);\n	(@resultTypes) = @{ $result_type_list_xml\
->{'type'} };\n	print_debug_message( 'rest_get_res\
ult_types',\n		scalar(@resultTypes) . ' result typ\
es', 2 );\n	print_debug_message( 'rest_get_result_\
types', 'End', 1 );\n	return (@resultTypes);\n}\n\\
n=head2 rest_get_result()\n\nGet result data of a \
specified type for a finished job.\n\n  my $result\
 = rest_get_result($job_id, $result_type);\n\n=cut\
\n\nsub rest_get_result {\n	print_debug_message( '\
rest_get_result', 'Begin', 1 );\n	my $job_id = shi\
ft;\n	my $type   = shift;\n	print_debug_message( '\
rest_get_result', 'jobid: ' . $job_id, 1 );\n	prin\
t_debug_message( 'rest_get_result', 'type: ' . $ty\
pe,    1 );\n	my $url    = $baseUrl . '/result/' .\
 $job_id . '/' . $type;\n	my $result = &rest_reque\
st($url);\n	print_debug_message( 'rest_get_result'\
, length($result) . ' characters',\n		1 );\n	print\
_debug_message( 'rest_get_result', 'End', 1 );\n	r\
eturn $result;\n}\n\n\n=head2 print_debug_message(\
)\n\nPrint debug message at specified debug level.\
\n\n  &print_debug_message($method_name, $message,\
 $level);\n\n=cut\n\nsub print_debug_message {\n	m\
y $function_name = shift;\n	my $message       = sh\
ift;\n	my $level         = shift;\n	if ( $level <=\
 $params{'debugLevel'} ) {\n		print STDERR '[', $f\
unction_name, '()] ', $message, \"\\n\";\n	}\n}\n\\
n=head2 print_tool_params()\n\nPrint list of tool \
parameters.\n\n  &print_tool_params();\n\n=cut\n\n\
sub print_tool_params {\n	print_debug_message( 'pr\
int_tool_params', 'Begin', 1 );\n	my (@param_list)\
 = &rest_get_parameters();\n	foreach my $param ( s\
ort(@param_list) ) {\n		print $param, \"\\n\";\n	}\
\n	print_debug_message( 'print_tool_params', 'End'\
, 1 );\n}\n\n=head2 print_param_details()\n\nPrint\
 details of a tool parameter.\n\n  &print_param_de\
tails($param_name);\n\n=cut\n\nsub print_param_det\
ails {\n	print_debug_message( 'print_param_details\
', 'Begin', 1 );\n	my $paramName = shift;\n	print_\
debug_message( 'print_param_details', 'paramName: \
' . $paramName, 2 );\n	my $paramDetail = &rest_get\
_parameter_details($paramName);\n	print $paramDeta\
il->{'name'}, \"\\t\", $paramDetail->{'type'}, \"\\
\n\";\n	print $paramDetail->{'description'}, \"\\n\
\";\n	if(defined($paramDetail->{'values'}->{'value\
'})) {\n		if(ref($paramDetail->{'values'}->{'value\
'}) eq 'ARRAY') {\n			foreach my $value ( @{ $para\
mDetail->{'values'}->{'value'} } ) {\n				&print_p\
aram_value($value);\n			}\n		}\n		else {\n				&pri\
nt_param_value($paramDetail->{'values'}->{'value'}\
);\n		}\n	}\n	print_debug_message( 'print_param_de\
tails', 'End', 1 );\n}\n\n=head2 print_param_value\
()\n\nPrint details of a tool parameter value.\n\n\
  &print_param_details($param_value);\n\nUsed by p\
rint_param_details() to handle both singluar and a\
rray values.\n\n=cut\n\nsub print_param_value {\n	\
my $value = shift;\n	print $value->{'value'};\n	if\
 ( $value->{'defaultValue'} eq 'true' ) {\n		print\
 \"\\t\", 'default';\n	}\n	print \"\\n\";\n	print \
\"\\t\", $value->{'label'}, \"\\n\";\n	if ( define\
d( $value->{'properties'} ) ) {\n		foreach\n		  my\
 $key ( sort( keys( %{ $value->{'properties'}{'pro\
perty'} } ) ) )\n		{\n			if ( ref( $value->{'prope\
rties'}{'property'}{$key} ) eq 'HASH'\n				&& defi\
ned( $value->{'properties'}{'property'}{$key}{'val\
ue'} )\n			  )\n			{\n				print \"\\t\", $key, \"\\
\t\",\n				  $value->{'properties'}{'property'}{$k\
ey}{'value'}, \"\\n\";\n			}\n			else {\n				print\
 \"\\t\", $value->{'properties'}{'property'}{'key'\
},\n				  \"\\t\", $value->{'properties'}{'propert\
y'}{'value'}, \"\\n\";\n				last;\n			}\n		}\n	}\n\
}\n\n=head2 print_job_status()\n\nPrint status of \
a job.\n\n  &print_job_status($job_id);\n\n=cut\n\\
nsub print_job_status {\n	print_debug_message( 'pr\
int_job_status', 'Begin', 1 );\n	my $jobid = shift\
;\n	print_debug_message( 'print_job_status', 'jobi\
d: ' . $jobid, 1 );\n	if ( $outputLevel > 0 ) {\n	\
	print STDERR 'Getting status for job ', $jobid, \\
"\\n\";\n	}\n	my $result = &rest_get_status($jobid\
);\n	print \"$result\\n\";\n	if ( $result eq 'FINI\
SHED' && $outputLevel > 0 ) {\n		print STDERR \"To\
 get results: $scriptName --polljob --jobid \" . $\
jobid\n		  . \"\\n\";\n	}\n	print_debug_message( '\
print_job_status', 'End', 1 );\n}\n\n=head2 print_\
result_types()\n\nPrint available result types for\
 a job.\n\n  &print_result_types($job_id);\n\n=cut\
\n\nsub print_result_types {\n	print_debug_message\
( 'result_types', 'Begin', 1 );\n	my $jobid = shif\
t;\n	print_debug_message( 'result_types', 'jobid: \
' . $jobid, 1 );\n	if ( $outputLevel > 0 ) {\n		pr\
int STDERR 'Getting result types for job ', $jobid\
, \"\\n\";\n	}\n	my $status = &rest_get_status($jo\
bid);\n	if ( $status eq 'PENDING' || $status eq 'R\
UNNING' ) {\n		print STDERR 'Error: Job status is \
', $status,\n		  '. To get result types the job mu\
st be finished.', \"\\n\";\n	}\n	else {\n		my (@re\
sultTypes) = &rest_get_result_types($jobid);\n		if\
 ( $outputLevel > 0 ) {\n			print STDOUT 'Availabl\
e result types:', \"\\n\";\n		}\n		foreach my $res\
ultType (@resultTypes) {\n			print STDOUT $resultT\
ype->{'identifier'}, \"\\n\";\n			if ( defined( $r\
esultType->{'label'} ) ) {\n				print STDOUT \"\\t\
\", $resultType->{'label'}, \"\\n\";\n			}\n			if \
( defined( $resultType->{'description'} ) ) {\n			\
	print STDOUT \"\\t\", $resultType->{'description'\
}, \"\\n\";\n			}\n			if ( defined( $resultType->{\
'mediaType'} ) ) {\n				print STDOUT \"\\t\", $res\
ultType->{'mediaType'}, \"\\n\";\n			}\n			if ( de\
fined( $resultType->{'fileSuffix'} ) ) {\n				prin\
t STDOUT \"\\t\", $resultType->{'fileSuffix'}, \"\\
\n\";\n			}\n		}\n		if ( $status eq 'FINISHED' && \
$outputLevel > 0 ) {\n			print STDERR \"\\n\", 'To\
 get results:', \"\\n\",\n			  \"  $scriptName --p\
olljob --jobid \" . $params{'jobid'} . \"\\n\",\n	\
		  \"  $scriptName --polljob --outformat <type> -\
-jobid \"\n			  . $params{'jobid'} . \"\\n\";\n		}\
\n	}\n	print_debug_message( 'result_types', 'End',\
 1 );\n}\n\n=head2 submit_job()\n\nSubmit a job to\
 the service.\n\n  &submit_job($seq);\n\n=cut\n\ns\
ub submit_job {\n	print_debug_message( 'submit_job\
', 'Begin', 1 );\n\n	# Set input sequence\n	$tool_\
params{'sequence'} = shift;\n\n	# Load parameters\\
n	&load_params();\n\n	# Submit the job\n	my $jobid\
 = &rest_run( $params{'email'}, $params{'title'}, \
\\%tool_params );\n\n	# Simulate sync/async mode\n\
	if ( defined( $params{'async'} ) ) {\n		print STD\
OUT $jobid, \"\\n\";\n		if ( $outputLevel > 0 ) {\\
n			print STDERR\n			  \"To check status: $scriptN\
ame --status --jobid $jobid\\n\";\n		}\n	}\n	else \
{\n		if ( $outputLevel > 0 ) {\n			print STDERR \"\
JobId: $jobid\\n\";\n		}\n		sleep 1;\n		&get_resul\
ts($jobid);\n	}\n	print_debug_message( 'submit_job\
', 'End', 1 );\n}\n\n=head2 multi_submit_job()\n\n\
Submit multiple jobs assuming input is a collectio\
n of fasta formatted sequences.\n\n  &multi_submit\
_job();\n\n=cut\n\nsub multi_submit_job {\n	print_\
debug_message( 'multi_submit_job', 'Begin', 1 );\n\
	my $jobIdForFilename = 1;\n	$jobIdForFilename = 0\
 if ( defined( $params{'outfile'} ) );\n	my (@file\
name_list) = ();\n\n	# Query sequence\n	if ( defin\
ed( $ARGV[0] ) ) {    # Bare option\n		if ( -f $AR\
GV[0] || $ARGV[0] eq '-' ) {    # File\n			push( @\
filename_list, $ARGV[0] );\n		}\n		else {\n			warn\
 'Warning: Input file \"' . $ARGV[0] . '\" does no\
t exist'\n		}\n	}\n	if ( $params{'sequence'} ) {  \
                 # Via --sequence\n		if ( -f $para\
ms{'sequence'} || $params{'sequence'} eq '-' ) {  \
  # File\n			push( @filename_list, $params{'sequen\
ce'} );\n		}\n		else {\n			warn 'Warning: Input fi\
le \"' . $params{'sequence'} . '\" does not exist'\
\n		}\n	}\n\n	$/ = '>';\n	foreach my $filename (@f\
ilename_list) {\n		my $INFILE;\n		if($filename eq \
'-') { # STDIN.\n			open( $INFILE, '<-' )\n			  or\
 die 'Error: unable to STDIN (' . $! . ')';\n		} e\
lse { # File.\n			open( $INFILE, '<', $filename )\\
n			  or die 'Error: unable to open file ' . $file\
name . ' (' . $! . ')';\n		}\n		while (<$INFILE>) \
{\n			my $seq = $_;\n			$seq =~ s/>$//;\n			if ( $\
seq =~ m/(\\S+)/ ) {\n				print STDERR \"Submittin\
g job for: $1\\n\"\n				  if ( $outputLevel > 0 );\
\n				$seq = '>' . $seq;\n				&print_debug_message\
( 'multi_submit_job', $seq, 11 );\n				&submit_job\
($seq);\n				$params{'outfile'} = undef if ( $jobI\
dForFilename == 1 );\n			}\n		}\n		close $INFILE;\\
n	}\n	print_debug_message( 'multi_submit_job', 'En\
d', 1 );\n}\n\n=head2 list_file_submit_job()\n\nSu\
bmit multiple jobs using a file containing a list \
of entry identifiers as \ninput.\n\n  &list_file_s\
ubmit_job($list_filename)\n\n=cut\n\nsub list_file\
_submit_job {\n	print_debug_message( 'list_file_su\
bmit_job', 'Begin', 11 );\n	my $filename         =\
 shift;\n	my $jobIdForFilename = 1;\n	$jobIdForFil\
ename = 0 if ( defined( $params{'outfile'} ) );\n\\
n	# Iterate over identifiers, submitting each job\\
n	my $LISTFILE;\n	if($filename eq '-') { # STDIN.\\
n		open( $LISTFILE, '<-' )\n		  or die 'Error: una\
ble to STDIN (' . $! . ')';\n	} else { # File.\n		\
open( $LISTFILE, '<', $filename )\n		  or die 'Err\
or: unable to open file ' . $filename . ' (' . $! \
. ')';\n	}\n	while (<$LISTFILE>) {\n		my $line = $\
_;\n		chomp($line);\n		if ( $line ne '' ) {\n			&p\
rint_debug_message( 'list_file_submit_job', 'line:\
 ' . $line, 2 );\n			if ( $line =~ m/\\w:\\w/ ) { \
   # Check this is an identifier\n				print STDERR\
 \"Submitting job for: $line\\n\"\n				  if ( $out\
putLevel > 0 );\n				&submit_job($line);\n			}\n		\
	else {\n				print STDERR\n\"Warning: line \\\"$li\
ne\\\" is not recognised as an identifier\\n\";\n	\
		}\n		}\n		$params{'outfile'} = undef if ( $jobId\
ForFilename == 1 );\n	}\n	close $LISTFILE;\n	print\
_debug_message( 'list_file_submit_job', 'End', 11 \
);\n}\n\n=head2 load_data()\n\nLoad sequence data \
from file or option specified on the command-line.\
\n\n  &load_data();\n\n=cut\n\nsub load_data {\n	p\
rint_debug_message( 'load_data', 'Begin', 1 );\n	m\
y $retSeq;\n\n	# Query sequence\n	if ( defined( $A\
RGV[0] ) ) {    # Bare option\n		if ( -f $ARGV[0] \
|| $ARGV[0] eq '-' ) {    # File\n			$retSeq = &re\
ad_file( $ARGV[0] );\n		}\n		else {               \
                      # DB:ID or sequence\n			$ret\
Seq = $ARGV[0];\n		}\n	}\n	if ( $params{'sequence'\
} ) {                   # Via --sequence\n		if ( -\
f $params{'sequence'} || $params{'sequence'} eq '-\
' ) {    # File\n			$retSeq = &read_file( $params{\
'sequence'} );\n		}\n		else {    # DB:ID or sequen\
ce\n			$retSeq = $params{'sequence'};\n		}\n	}\n	p\
rint_debug_message( 'load_data', 'End', 1 );\n	ret\
urn $retSeq;\n}\n\n=head2 load_params()\n\nLoad jo\
b parameters from command-line options.\n\n  &load\
_params();\n\n=cut\n\nsub load_params {\n	print_de\
bug_message( 'load_params', 'Begin', 1 );\n\n	# Da\
tabase(s) to search\n	my (@dbList) = split /[ ,]/,\
 $params{'database'};\n	$tool_params{'database'} =\
 \\@dbList;\n\n	# Compatability options, old comma\
nd-line.\n	if(!$tool_params{'viewfilter'} && $para\
ms{'echofilter'}) {\n		$tool_params{'viewfilter'} \
= 'true';\n	}\n	if(!$tool_params{'alignments'} && \
$params{'numal'}) {\n		$tool_params{'alignments'} \
= $params{'numal'};\n	}\n	# TODO: set alignment fo\
rmat option to get NCBI BLAST XML.\n	if($params{'a\
ppxml'}) {\n		$tool_params{'align'} = '';\n	}\n\n	\
print_debug_message( 'load_params', 'End', 1 );\n}\
\n\n=head2 client_poll()\n\nClient-side job pollin\
g.\n\n  &client_poll($job_id);\n\n=cut\n\nsub clie\
nt_poll {\n	print_debug_message( 'client_poll', 'B\
egin', 1 );\n	my $jobid  = shift;\n	my $status = '\
PENDING';\n\n	my $errorCount = 0;\n	while ($status\
 eq 'RUNNING'\n		|| $status eq 'PENDING'\n		|| ( $\
status eq 'ERROR' && $errorCount < 2 ) )\n	{\n		$s\
tatus = rest_get_status($jobid);\n		print STDERR \\
"$status\\n\" if ( $outputLevel > 0 );\n		if ( $st\
atus eq 'ERROR' ) {\n			$errorCount++;\n		}\n		els\
if ( $errorCount > 0 ) {\n			$errorCount--;\n		}\n\
		if (   $status eq 'RUNNING'\n			|| $status eq 'P\
ENDING'\n			|| $status eq 'ERROR' )\n		{\n\n			# W\
ait before polling again.\n			sleep $checkInterval\
;\n		}\n	}\n	print_debug_message( 'client_poll', '\
End', 1 );\n	return $status;\n}\n\n=head2 get_resu\
lts()\n\nGet the results for a job identifier.\n\n\
  &get_results($job_id);\n\n=cut\n\nsub get_result\
s {\n	print_debug_message( 'get_results', 'Begin',\
 1 );\n	my $jobid = shift;\n	print_debug_message( \
'get_results', 'jobid: ' . $jobid, 1 );\n\n	# Verb\
ose\n	if ( $outputLevel > 1 ) {\n		print 'Getting \
results for job ', $jobid, \"\\n\";\n	}\n\n	# Chec\
k status, and wait if not finished\n	client_poll($\
jobid);\n\n	# Use JobId if output file name is not\
 defined\n	unless ( defined( $params{'outfile'} ) \
) {\n		$params{'outfile'} = $jobid;\n	}\n\n	# Get \
list of data types\n	my (@resultTypes) = rest_get_\
result_types($jobid);\n\n	# Get the data and write\
 it to a file\n	if ( defined( $params{'outformat'}\
 ) ) {    # Specified data type\n		my $selResultTy\
pe;\n		foreach my $resultType (@resultTypes) {\n		\
	if ( $resultType->{'identifier'} eq $params{'outf\
ormat'} ) {\n				$selResultType = $resultType;\n		\
	}\n		}\n		if ( defined($selResultType) ) {\n			my\
 $result =\n			  rest_get_result( $jobid, $selResu\
ltType->{'identifier'} );\n			if ( $params{'outfil\
e'} eq '-' ) {\n				write_file( $params{'outfile'}\
, $result );\n			}\n			else {\n				write_file(\n		\
			$params{'outfile'} . '.'\n					  . $selResultTy\
pe->{'identifier'} . '.'\n					  . $selResultType-\
>{'fileSuffix'},\n					$result\n				);\n			}\n		}\\
n		else {\n			die 'Error: unknown result format \"\
' . $params{'outformat'} . '\"';\n		}\n	}\n	else {\
    # Data types available\n		      # Write a file\
 for each output type\n		for my $resultType (@resu\
ltTypes) {\n			if ( $outputLevel > 1 ) {\n				prin\
t STDERR 'Getting ', $resultType->{'identifier'}, \
\"\\n\";\n			}\n			my $result = rest_get_result( $\
jobid, $resultType->{'identifier'} );\n			if ( $pa\
rams{'outfile'} eq '-' ) {\n				write_file( $param\
s{'outfile'}, $result );\n			}\n			else {\n				wri\
te_file(\n					$params{'outfile'} . '.'\n					  . \
$resultType->{'identifier'} . '.'\n					  . $resul\
tType->{'fileSuffix'},\n					$result\n				);\n			}\
\n		}\n	}\n	print_debug_message( 'get_results', 'E\
nd', 1 );\n}\n\n=head2 read_file()\n\nRead a file \
into a scalar. The special filename '-' can be use\
d to read from \nstandard input (STDIN).\n\n  my $\
data = &read_file($filename);\n\n=cut\n\nsub read_\
file {\n	print_debug_message( 'read_file', 'Begin'\
, 1 );\n	my $filename = shift;\n	print_debug_messa\
ge( 'read_file', 'filename: ' . $filename, 2 );\n	\
my ( $content, $buffer );\n	if ( $filename eq '-' \
) {\n		while ( sysread( STDIN, $buffer, 1024 ) ) {\
\n			$content .= $buffer;\n		}\n	}\n	else {    # F\
ile\n		open( my $FILE, '<', $filename )\n		  or di\
e \"Error: unable to open input file $filename ($!\
)\";\n		while ( sysread( $FILE, $buffer, 1024 ) ) \
{\n			$content .= $buffer;\n		}\n		close($FILE);\n\
	}\n	print_debug_message( 'read_file', 'End', 1 );\
\n	return $content;\n}\n\n=head2 write_file()\n\nW\
rite data to a file. The special filename '-' can \
be used to write to \nstandard output (STDOUT).\n\\
n  &write_file($filename, $data);\n\n=cut\n\nsub w\
rite_file {\n	print_debug_message( 'write_file', '\
Begin', 1 );\n	my ( $filename, $data ) = @_;\n	pri\
nt_debug_message( 'write_file', 'filename: ' . $fi\
lename, 2 );\n	if ( $outputLevel > 0 ) {\n		print \
STDERR 'Creating result file: ' . $filename . \"\\\
n\";\n	}\n	if ( $filename eq '-' ) {\n		print STDO\
UT $data;\n	}\n	else {\n		open( my $FILE, '>', $fi\
lename )\n		  or die \"Error: unable to open outpu\
t file $filename ($!)\";\n		syswrite( $FILE, $data\
 );\n		close($FILE);\n	}\n	print_debug_message( 'w\
rite_file', 'End', 1 );\n}\n\n=head2 usage()\n\nPr\
int program usage message.\n\n  &usage();\n\n=cut\\
n\nsub usage {\n	print STDERR <<EOF\nWU-BLAST\n===\
=====\n   \nRapid sequence database search program\
s utilizing the BLAST algorithm\n    \n[Required]\\
n\n  -p, --program      : str  : BLAST program to \
use, see --paramDetail program\n  -D, --database  \
   : str  : database(s) to search, space separated\
. See\n                              --paramDetail\
 database\n      --stype        : str  : query seq\
uence type, see --paramDetail stype\n  seqFile    \
        : file : query sequence (\"-\" for STDIN, \
\\@filename for\n                              ide\
ntifier list file)\n\n[Optional]\n\n  -m, --matrix\
       : str  : scoring matrix, see --paramDetail \
matrix\n  -e, --exp          : real : 0<E<= 1000. \
Statistical significance threshold \n             \
                 for reporting database sequence m\
atches.\n  -e, --viewfilter   :      : display the\
 filtered query sequence\n  -f, --filter       : s\
tr  : filter the query sequence for low complexity\
 \n                              regions, see --pa\
ramDetail filter\n  -A, --align        : int  : pa\
irwise alignment format, see --paramDetail align\n\
  -s, --scores       : int  : number of scores to \
be reported\n  -b, --alignments   : int  : number \
of alignments to report\n  -S, --sensitivity  : st\
r  : sensitivity of the search, \n                \
              see --paramDetail sensitivity\n  -t,\
 --sort	     : str  : sort order for hits, see --p\
aramDetail sort\n  -T, --stats        : str  : sta\
tistical model, see --paramDetail stats\n  -d, --s\
trand       : str  : DNA strand to search with,\n \
                             see --paramDetail str\
and\n  -c, --topcombon    : str  : consistent sets\
 of HSPs\n      --multifasta   :      : treat inpu\
t as a set of fasta formatted sequences\n\n[Genera\
l]\n\n  -h, --help         :      : prints this he\
lp text\n      --async        :      : forces to m\
ake an asynchronous query\n      --email        : \
str  : e-mail address\n      --title        : str \
 : title for job\n      --status       :      : ge\
t job status\n      --resultTypes  :      : get av\
ailable result types for job\n      --polljob     \
 :      : poll for the status of a job\n      --jo\
bid        : str  : jobid that was returned when a\
n asynchronous job \n                             \
 was submitted.\n      --outfile      : str  : fil\
e name for results (default is jobid;\n           \
                   \"-\" for STDOUT)\n      --outf\
ormat    : str  : result format to retrieve\n     \
 --params       :      : list input parameters\n  \
    --paramDetail  : str  : display details for in\
put parameter\n      --quiet        :      : decre\
ase output\n      --verbose      :      : increase\
 output\n   \nSynchronous job:\n\n  The results/er\
rors are returned as soon as the job is finished.\\
n  Usage: $scriptName --email <your\\@email> [opti\
ons...] seqFile\n  Returns: results as an attachme\
nt\n\nAsynchronous job:\n\n  Use this if you want \
to retrieve the results at a later time. The resul\
ts \n  are stored for up to 24 hours. 	\n  Usage: \
$scriptName --async --email <your\\@email> [option\
s...] seqFile\n  Returns: jobid\n\n  Use the jobid\
 to query for the status of the job. If the job is\
 finished, \n  it also returns the results/errors.\
\n  Usage: $scriptName --polljob --jobid <jobId> [\
--outfile string]\n  Returns: string indicating th\
e status of the job and if applicable, results \n \
 as an attachment.\n\nFurther information:\n\n  ht\
tp://www.ebi.ac.uk/Tools/webservices/services/sss/\
wu_blast_rest\n  http://www.ebi.ac.uk/Tools/webser\
vices/tutorials/perl\n\nSupport/Feedback:\n\n  htt\
p://www.ebi.ac.uk/support/\nEOF\n}\n\n=head1 FEEDB\
ACK/SUPPORT\n\nPlease contact us at L<http://www.e\
bi.ac.uk/support/> if you have any \nfeedback, sug\
gestions or issues with the service or this client\
.\n\n=cut\n","\n\n\nmy $PROBTRESH = 0.3;# base pai\
rs below this prob threshold will be ignored\nmy $\
WEIGHT = 100.0; # float!!\nmy $NUCALPH = \"ACGTUNR\
YMKSWHBVD\";\nuse vars qw($NUCALPH $WEIGHT);\n\nmy\
 $myname = basename($0);\n\nuse strict;\nuse warni\
ngs;\n\nuse File::Basename;\nuse Getopt::Long;\nus\
e File::Glob ':glob';\nuse File::Spec;\nuse File::\
Temp qw/ tempfile tempdir /;\n\n\n\n\nsub tcoffeel\
ib_header($;$)\n{\n    my ($nseq, $fd) = @_;\n    \
if (! defined($fd)) {\n        $fd = *STDOUT;\n   \
 }\n    printf $fd \"! TC_LIB_FORMAT_01\\n\";\n   \
 printf $fd \"%d\\n\", $nseq;\n}\n\n\nsub tcoffeel\
ib_header_addseq($$;$)\n{\n    my ($id, $seq, $fd)\
 = @_;\n    if (! defined($fd)) {\n        $fd = *\
STDOUT;\n    }\n    printf $fd \"%s %d %s\\n\", $i\
d, length($seq), $seq;\n}\n\n\nsub tcoffeelib_comm\
ent($;$)\n{\n    my ($comment, $fd) = @_;\n    if \
(! defined($fd)) {\n        $fd = *STDOUT;\n    }\\
n    printf $fd \"!\" . $comment . \"\\n\";\n}\n\n\
\nsub tcoffeelib_struct($$$;$)\n{\n    my ($nseq, \
$len, $bpm, $fd) = @_;\n\n    if (! defined($fd)) \
{\n        $fd = *STDOUT;\n    }\n\n    # output b\
asepair indices with fixed weight\n    printf $fd \
\"#%d %d\\n\", $nseq, $nseq;\n    # output basepai\
rs (only once) and with unit-offset\n    for (my $\
i=0; $i<$len; $i++) {\n        for (my $j=$i+1; $j\
<$len; $j++) {\n            if (! defined($bpm->[$\
i][$j])) {\n                print STDERR \"ERROR: \
\\$bpm->[$i][$j] undefined\\n\";\n            }\n \
           if ($bpm->[$i][$j]>0) {\n              \
  print $fd $i+1;\n                print $fd \" \"\
;\n                print $fd $j+1;\n              \
  print $fd \" \" . $bpm->[$i][$j] . \"\\n\";\n   \
         }\n        }\n    }\n}\n\n\nsub tcoffeeli\
b_footer(;$)\n{\n    my ($fd) = @_;\n    if (! def\
ined($fd)) {\n        $fd = *STDOUT;\n    }\n    p\
rint $fd \"! SEQ_1_TO_N\\n\";\n}\n\n\n    \nsub pl\
fold($$$)\n{    \n    my ($id, $seq, $probtresh) =\
 @_;\n    my (@struct);# return\n    my ($templ, $\
fhtmp, $fnametmp, $cmd, $ctr, $window_size);\n    \
our $ntemp++;\n    \n    $templ = $myname . \".pid\
-\" . $$ .$ntemp .\".XXXXXX\";\n    ($fhtmp, $fnam\
etmp) = tempfile($templ, UNLINK => 1); \n    print\
 $fhtmp \">$id\\n$seq\\n\";\n\n    # --- init base\
pair array\n    #\n    for (my $i=0; $i<length($se\
q); $i++) {\n        for (my $j=$i+1; $j<length($s\
eq); $j++) {\n            $struct[$i][$j]=0;\n    \
    }\n    }\n\n\n    # --- call rnaplfold and dro\
p a readme\n    #\n    $window_size=(length($seq)<\
70)?length($seq):70;\n    $cmd = \"RNAplfold -W $w\
indow_size < $fnametmp >/dev/null\";\n    system($\
cmd);\n    \n    if ($? != 0) {\n        printf ST\
DERR \"ERROR: RNAplfold ($cmd) exited with error s\
tatus %d\\n\", $? >> 8;\n        return;\n    }\n \
   #unlink($fnametmp);\n    my $fps = sprintf(\"%s\
_dp.ps\", $id); # check long name\n    \n    if (!\
 -s $fps) {\n      {\n\n	$fps = sprintf(\"%s_dp.ps\
\", substr($id,0,12)); # check short name\n 	if (!\
 -s $fps)\n	  {\n	    die(\"couldn't find expected\
 file $fps\\n\");\n	    return;\n	  }\n      }\n  \
  }\n\n    \n    # --- read base pairs from create\
d postscript\n    #\n    open(FH, $fps);\n    whil\
e (my $line = <FH>) {\n        my ($nti, $ntj, $pr\
ob);\n        chomp($line);        \n        # lin\
e: bp bp sqrt-prob ubox\n        my @match = ($lin\
e =~ m/^([0-9]+) +([0-9]+) +([0-9\\.]+) +ubox$/);\\
n        if (scalar(@match)) {\n            $nti=$\
1;\n            $ntj=$2;\n            $prob=$3*$3;\
# prob stored as square root\n\n            if ($p\
rob>$probtresh) {\n                #printf STDERR \
\"\\$struct[$nti][$ntj] sqrtprob=$3 prob=$prob > $\
probtresh\\n\";\n                $struct[$nti-1][$\
ntj-1] = $WEIGHT\n            }\n            # sto\
re with zero-offset\n        }\n    }\n    close(F\
H);\n\n    # remove or gzi postscript\n    #\n    \
unlink($fps);\n    #\n    # or gzip\n    #$cmd = \\
"gzip -qf $fps\";\n    #system($cmd);\n    #if ($?\
 != 0) {\n    #    printf STDERR \"ERROR: gzip ($c\
md) exited with error status %d\\n\", $? >> 8;\n  \
  #}\n\n    return \\@struct;\n}\n\n\n\n\n\nsub rn\
aseqfmt($)\n{\n    my ($seq) = @_;\n    # remove g\
aps\n    $seq =~ s/-//g;\n    # uppercase RNA\n   \
 $seq = uc($seq);\n    # T -> U\n    $seq =~ s/T/U\
/g;\n    # check for invalid charaters\n    $_ = $\
seq;\n    s/[^$NUCALPH]//g;\n    return $_;\n}\n\n\
\n\n\nsub usage(;$)\n{    \n    my ($errmsg) = @_;\
\n    if ($errmsg) {\n        print STDERR \"ERROR\
: $errmsg\\n\";\n    }\n    print STDERR << \"EOF\\
";\n$myname:\n Creates a T-Coffee RNA structure li\
brary from RNAplfold prediction.\n See FIXME:citat\
ion\nUsage:\n $myname -in seq_file -out tcoffee_li\
b\nEOF\n    exit(1);\n}\n\nsub read_fasta_seq \n  \
{\n    my $f=$_[0];\n    my %hseq;\n    my (@seq, \
@com, @name);\n    my ($a, $s,$nseq);\n\n    open \
(F, $f);\n    while (<F>)\n      {\n	$s.=$_;\n    \
  }\n    close (F);\n\n    \n    @name=($s=~/>(\\S\
*).*\\n[^>]*/g);\n    \n    @seq =($s=~/>.*.*\\n([\
^>]*)/g);\n    @com =($s=~/>(\\S*)(.*)\\n([^>]*)/g\
);\n\n\n    $nseq=$#name+1;\n  \n    for ($a=0; $a\
<$nseq; $a++)\n      {\n	my $n=$name[$a];\n	my $s;\
\n	$hseq{$n}{name}=$n;\n	$s=$seq[$a];$s=~s/\\s//g;\
\n	\n	$hseq{$n}{seq}=$s;\n	$hseq{$n}{com}=$com[$a]\
;\n      }\n    return %hseq;\n  }\n\n\n\n\n\n\n\n\
my $fmsq = \"\";\nmy $flib = \"\";\nmy %OPTS;\nmy \
%seq;\nmy ($id, $nseq, $i);\nmy @nl;\n\nGetOptions\
(\"in=s\" => \\$fmsq, \"out=s\" => \\$flib);\n\nif\
 (! -s $fmsq) {\n    usage(\"empty or non-existant\
 file \\\"$fmsq\\\"\")\n}\nif (length($flib)==0) {\
\n    usage(\"empty out-filename\")\n}\n\n\n\n\n\n\
\n%seq=read_fasta_seq($fmsq);\n\n\n@nl=keys(%seq);\
\n\n$nseq=$#nl+1;\nopen FD_LIB, \">$flib\" or die \
\"can't open $flib!\";\ntcoffeelib_header($nseq, *\
FD_LIB);\nforeach $id (keys (%seq))\n  {\n    my (\
$seq, $fmtseq);\n    \n    $seq = $seq{$id}{seq};\\
n    \n    $fmtseq = rnaseqfmt($seq);# check here,\
 formatting for folding important later\n    if (l\
ength($seq)!=length($fmtseq)) {\n        print STD\
ERR \"ERROR: invalid sequence $id is not an RNA se\
quence. read seq is: $seq\\n\";\n        exit\n   \
   }\n   \n    tcoffeelib_header_addseq($id, uc($s\
eq), *FD_LIB);\n  }\ntcoffeelib_comment(\"generate\
d by $myname on \" . localtime(), *FD_LIB);\n\n\n\\
n$i=0;\nforeach $id (keys (%seq))\n  {\n    my ($c\
leanid, $seq, $bpm);\n    $seq=$seq{$id}{seq};\n  \
  $cleanid = $id;\n    $cleanid =~ s,[/ ],_,g;# ne\
eded for rnaplfold\n    $seq = rnaseqfmt($seq);\n \
   \n    $bpm = plfold($cleanid, rnaseqfmt($seq), \
$PROBTRESH);       \n    \n    tcoffeelib_struct($\
i+1, length($seq), $bpm, *FD_LIB);\n    $i++;\n}\n\
\n\ntcoffeelib_footer(*FD_LIB);\nclose FD_LIB;\nex\
it (0);\n\n","\n\n\n\n\n$cmd=join ' ', @ARGV;\nif \
($cmd=~/-infile=(\\S+)/){ $seqfile=$1;}\nif ($cmd=\
~/-outfile=(\\S+)/){ $libfile=$1;}\n\n\n\n%s=read_\
fasta_seq ($seqfile);\n\nopen (F, \">$libfile\");\\
nforeach $name (keys (%s))\n  {\n    my $tclib=\"$\
name.RNAplfold_tclib\";\n    print (F \">$name _F_\
 $tclib\\n\");\n    seq2RNAplfold2tclib ($name, $s\
{$name}{seq}, $tclib);\n  }\nclose (F);\nexit (EXI\
T_SUCCESS);\n\nsub seq2RNAplfold2tclib\n  {\n    m\
y ($name, $seq, $tclib)=@_;\n    my ($tmp);\n    $\
n++;\n    $tmp=\"tmp4seq2RNAplfold_tclib.$$.$n.pep\
\";\n    open (RF, \">$tmp\");\n    print (RF \">$\
name\\n$seq\\n\");\n    close (RF);\n    \n    sys\
tem \"t_coffee -other_pg RNAplfold2tclib.pl -in=$t\
mp -out=$tclib\";\n    \n    unlink ($tmp);\n    r\
eturn $tclib;\n  }\n    \n    \nsub read_fasta_seq\
 \n  {\n    my $f=@_[0];\n    my %hseq;\n    my (@\
seq, @com, @name);\n    my ($a, $s,$nseq);\n\n    \
open (F, $f);\n    while (<F>)\n      {\n	$s.=$_;\\
n      }\n    close (F);\n\n    \n    @name=($s=~/\
>(\\S*).*\\n[^>]*/g);\n    \n    @seq =($s=~/>.*.*\
\\n([^>]*)/g);\n    @com =($s=~/>\\S*(.*)\\n([^>]*\
)/g);\n\n    \n    $nseq=$#name+1;\n    \n    for \
($a=0; $a<$nseq; $a++)\n      {\n	my $n=$name[$a];\
\n	$hseq{$n}{name}=$n;\n	$hseq{$n}{seq}=$seq[$a];\\
n	$hseq{$n}{com}=$com[$a];\n      }\n    return %h\
seq;\n  }\n","use Getopt::Long;\nuse File::Path;\n\
use Env;\nuse FileHandle;\nuse Cwd;\nuse Sys::Host\
name;\nour $PIDCHILD;\nour $ERROR_DONE;\nour @TMPF\
ILE_LIST;\nour $EXIT_FAILURE=1;\nour $EXIT_SUCCESS\
=0;\n\nour $REFDIR=getcwd;\nour $EXIT_SUCCESS=0;\n\
our $EXIT_FAILURE=1;\n\nour $PROGRAM=\"tc_generic_\
method.pl\";\nour $CL=$PROGRAM;\n\nour $CLEAN_EXIT\
_STARTED;\nour $debug_lock=$ENV{\"DEBUG_LOCK\"};\n\
our $LOCKDIR=$ENV{\"LOCKDIR_4_TCOFFEE\"};\nif (!$L\
OCKDIR){$LOCKDIR=getcwd();}\nour $ERRORDIR=$ENV{\"\
ERRORDIR_4_TCOFFEE\"};\nour $ERRORFILE=$ENV{\"ERRO\
RFILE_4_TCOFFEE\"};\n&set_lock ($$);\nif (isshellp\
id(getppid())){lock4tc(getppid(), \"LLOCK\", \"LSE\
T\", \"$$\\n\");}\n      \nour $print;\nmy ($fmsq1\
, $fmsq2, $output, $outfile, $arch, $psv, $hmmtop_\
home, $trim, $cov, $sample, $mode, $gor_home, $gor\
_seq, $gor_obs);\n\nGetOptions(\"-in=s\" => \\$fms\
q1,\"-output=s\" =>\\$output ,\"-out=s\" => \\$out\
file, \"-arch=s\" => \\$arch,\"-psv=s\" => \\$psv,\
 \"-hmmtop_home=s\", \\$hmmtop_home,\"-trim=s\" =>\
\\$trim ,\"-print=s\" =>\\$print,\"-cov=s\" =>\\$c\
ov , \"-sample=s\" =>\\$sample, \"-mode=s\" =>\\$m\
ode, \"-gor_home=s\"=>\\$gor_home, \"-gor_seq=s\"=\
>\\$gor_seq,\"-gor_obs=s\"=>\\$gor_obs);\n\n\nif (\
!$mode){$mode = \"hmmtop\"}\nelsif ($mode eq \"hmm\
top\"){;}\nelsif ($mode eq \"gor\"){;}\nelse {myex\
it(flush_error (\"-mode=$mode is unknown\"));}\n\n\
\nour $HOME=$ENV{\"HOME\"};\nour $MCOFFEE=($ENV{\"\
MCOFFEE_4_TCOFFEE\"})?$ENV{\"MCOFFEE_4_TCOFFEE\"}:\
\"$HOME/.t_coffee/mcoffee\";\n\nif ($mode eq \"hmm\
top\")\n  {\n    \n    check_configuration (\"hmmt\
op\");\n    if (-e $arch){$ENV{'HMMTOP_ARCH'}=$arc\
h;}\n    elsif (-e $ENV{HMMTOP_ARCH}){$arch=$ENV{H\
MMTOP_ARCH};}\n    elsif (-e \"$MCOFFEE/hmmtop.arc\
h\"){$arch=$ENV{'HMMTOP_ARCH'}=\"$MCOFFEE/hmmtop.a\
rch\";}\n    elsif (-e \"$hmmtop_home/hmmtop.arch\\
"){$arch=$ENV{'HMMTOP_ARCH'}=\"$hmmtop_home/hmmtop\
.arch\";}\n    else {myexit(flush_error ( \"Could \
not find ARCH file for hmmtop\"));}\n    \n    \n \
   if (-e $psv){$ENV{'HMMTOP_PSV'}=$psv;}\n    els\
if (-e $ENV{HMMTOP_PSV}){$psv=$ENV{HMMTOP_PSV};}\n\
    elsif (-e \"$MCOFFEE/hmmtop.psv\"){$psv=$ENV{'\
HMMTOP_PSV'}=\"$MCOFFEE/hmmtop.psv\";}\n    elsif \
(-e \"$hmmtop_home/hmmtop.psv\"){$psv=$ENV{'HMMTOP\
_PSV'}=\"$hmmtop_home/hmmtop.psv\";}\n    else {my\
exit(flush_error ( \"Could not find PSV file for h\
mmtop\"));}\n\n  }\nelsif ($mode eq \"gor\")\n  {\\
n    our $GOR_SEQ;\n    our $GOR_OBS;\n    \n    c\
heck_configuration (\"gorIV\");\n    if (-e $gor_s\
eq){$GOR_SEQ=$gor_seq;}\n    elsif (-e $ENV{GOR_SE\
Q}){$GOR_SEQ=$ENV{GOR_SEQ};}\n    elsif (-e \"$MCO\
FFEE/New_KS.267.seq\"){$GOR_SEQ=\"$MCOFFEE/New_KS.\
267.seq\";}\n    elsif (-e \"$gor_home/New_KS.267.\
seq\"){$GOR_SEQ=\"$gor_home/New_KS.267.seq\";}\n  \
  else {myexit(flush_error ( \"Could not find SEQ \
file for gor\"));}\n\n    if (-e $gor_obs){$GOR_OB\
S=$gor_obs;}\n    elsif (-e $ENV{GOR_OBS}){$GOR_OB\
S=$ENV{GOR_OBS};}\n    elsif (-e \"$MCOFFEE/New_KS\
.267.obs\"){$GOR_OBS=\"$MCOFFEE/New_KS.267.obs\";}\
\n    elsif (-e \"$gor_home/New_KS.267.obs\"){$GOR\
_OBS=\"$gor_home/New_KS.267.obs\";}\n    else {mye\
xit(flush_error ( \"Could not find OBS file for go\
r\"));}\n  }\n\n\nif ( ! -e $fmsq1){myexit(flush_e\
rror (\"Could Not Read Input file $fmsq1\"));}\n\n\
\nmy $fmsq2=vtmpnam();\nmy $fmsq3=vtmpnam();\nmy $\
tmpfile=vtmpnam();\nmy $predfile=vtmpnam();\n\nif \
($trim){$trim_action=\" +trim _aln_%%$trim\\_K1 \"\
;}\nif ($cov) {$cov_action= \" +sim_filter _aln_c$\
cov \";}\n&safe_system(\"t_coffee -other_pg seq_re\
format -in $fmsq1 -action +convert 'BOUJXZ-' $cov_\
action $trim_action -output fasta_aln -out $fmsq2\\
");\nmy (%pred, %seq, %predA);\n\n\n%seq=read_fast\
a_seq($fmsq2);\n%seq=fasta2sample(\\%seq, $sample)\
;\n\nif (1==2 &&$mode eq \"hmmtop\" && $output eq \
\"cons\")\n  {\n    fasta2hmmtop_cons($outfile,\\%\
seq);\n  }\nelse\n  {\n   \n    %pred=fasta2pred(\\
\%seq, $mode);\n    %predA=pred2aln (\\%pred, \\%s\
eq);\n    \n    \n    if (!$output || $output eq \\
"prediction\"){output_fasta_seq (\\%predA, $outfil\
e);}\n    elsif ($output eq \"color_html\"){pred2c\
olor (\\%pred,\\%seq, $outfile);}\n    elsif ($out\
put eq \"cons\"){pred2cons($outfile,\\%predA);}\n \
   else {flush_error (\"$output is an unknown outp\
ut mode\");}\n  }\n\nsub fasta2sample\n  {\n    my\
 $SR=shift;\n    my $it=shift;\n    my %S=%$SR;\n \
   \n    my $seq=index2seq_name (\\%S, 1);\n    my\
 $l=length($S{$seq}{seq});\n    my @sl=keys(%S);\n\
    my $nseq=$#sl+1;\n    my $index=$nseq;\n  \n  \
  if (!$sample) {return %S;}\n    for (my $a=0; $a\
<$it; $a++)\n      {\n	my $newseq=\"\";\n	my $nnam\
e=\"$seq\\_sampled_$index\";\n	for (my $p=0; $p<$l\
; $p++)\n	  {\n	    my $i=int(rand($nseq));\n	    \
\n	    my $name = $sl[$i];\n	    my $seq=$S{$name}\
{seq};\n	    my $r=substr ($seq, $p, 1);\n	    $ne\
wseq.=$r;\n	  }\n	$S{$nname}{name}=$nname;\n	$S{$n\
name}{seq}=$newseq;\n	$S{$nname}{com}=\"sampled\";\
\n	$S{$nname}{index}=++$index;\n      }\n    retur\
n %S;\n  }\n	      \nsub fasta2pred\n  {\n    my $\
s=shift;\n    my $mode=shift;\n\n    if ( $mode eq\
 \"hmmtop\"){return fasta2hmmtop_pred($s);}\n    e\
lsif ($mode eq \"gor\"){return fasta2gor_pred ($s)\
;}\n  }\nsub fasta2hmmtop_cons\n  {\n    my $outfi\
le=shift;\n    my $SR=shift;\n    \n    my $o = ne\
w FileHandle;\n    my $i = new FileHandle;\n    my\
 $tmp_in =vtmpnam();\n    my $tmp_out=vtmpnam();\n\
    my %seq=%$SR;\n    my %pred;\n    my $N=keys(%\
seq);\n    \n    output_fasta_seq (\\%seq,$tmp_in,\
 \"seq\");\n    `hmmtop -pi=mpred -if=$tmp_in -sf=\
FAS -pl 2>/dev/null >$tmp_out`;\n    open ($o, \">\
$outfile\");\n    open ($i, \"$tmp_out\");\n    wh\
ile (<$i>)\n      {\n	my $l=$_;\n	if (($l=~/>HP\\:\
\\s+(\\d+)\\s+(.*)/)){my $line=\">$2 NSEQ: $N\\n\"\
;print $o \"$line\";}\n	elsif ( ($l=~/.*pred(.*)/)\
)  {my $line=\"$1\\n\";print $o \"$line\";}\n     \
 }\n    close ($o);\n    close ($i);\n    return r\
ead_fasta_seq($tmp);\n  }\nsub fasta2hmmtop_pred\n\
  {\n    my $SR=shift;\n    my $o = new FileHandle\
;\n    my $i = new FileHandle;\n    my $tmp    =vt\
mpnam();\n    my $tmp_in =vtmpnam();\n    my $tmp_\
out=vtmpnam();\n    my %seq=%$SR;\n    my %pred;\n\
    \n\n    output_fasta_seq (\\%seq,$tmp_in, \"se\
q\");\n\n    \n    `hmmtop -if=$tmp_in -sf=FAS -pl\
 2>/dev/null >$tmp_out`;\n    \n\n    \n    \n    \
open ($o, \">$tmp\");\n    open ($i, \"$tmp_out\")\
;\n    while (<$i>)\n      {\n	my $l=$_;\n	if (($l\
=~/>HP\\:\\s+(\\d+)\\s+(.*)/)){my $line=\">$2\\n\"\
;print $o \"$line\";}\n	elsif ( ($l=~/.*pred(.*)/)\
)  {my $line=\"$1\\n\";print $o \"$line\";}\n     \
 }\n    close ($o);\n    close ($i);\n    return r\
ead_fasta_seq($tmp);\n  }\n    \n	\n	\n	    \n	\n	\
\n\n	\nsub fasta2gor_pred\n  {\n    my $SR=shift;\\
n    my $o = new FileHandle;\n    my $i = new File\
Handle;\n    my $tmp    =vtmpnam();\n    my $tmp_i\
n =vtmpnam();\n    my $tmp_out=vtmpnam();\n    my \
%seq=%$SR;\n    my %pred;\n    \n\n    output_fast\
a_seq (\\%seq,$tmp_in, \"seq\");\n    `gorIV -prd \
$tmp_in -seq $GOR_SEQ -obs $GOR_OBS >$tmp_out`;\n \
   open ($o, \">$tmp\");\n    open ($i, \"$tmp_out\
\");\n    while (<$i>)\n      {\n	my $l=$_;\n\n	\n\
	if ( $l=~/>/){print $o \"$l\";}\n	elsif ( $l=~/Pr\
edicted Sec. Struct./){$l=~s/Predicted Sec. Struct\
\\.//;print $o \"$l\";}\n      }\n    close ($o);\\
n    close ($i);\n    return read_fasta_seq($tmp);\
\n  }\n			\n			     \nsub index2seq_name\n  {\n   \
 \n    my $SR=shift;\n    my $index=shift;\n    \n\
    \n    my %S=%$SR;\n    \n    foreach my $s (%S\
)\n      {\n	if ( $S{$s}{index}==$index){return $s\
;}\n      }\n    return \"\";\n  }\n\nsub pred2con\
s\n  {\n    my $outfile=shift;\n    my $predR=shif\
t;\n    my $seq=shift;\n    my %P=%$predR;\n    my\
 %C;\n    my ($s,@r,$nseq);\n    my $f= new FileHa\
ndle;\n\n    open ($f, \">$outfile\");\n\n    if (\
!$seq){$seq=index2seq_name(\\%P,1);}\n    foreach \
my $s (keys(%P))\n      {\n	$nseq++;\n	$string= $P\
{$s}{seq};\n	$string = uc $string;\n	my @r=split (\
//,$string);\n	for (my $a=0; $a<=$#r; $a++)\n	  {\\
n	    if (($r[$a]=~/[OHICE]/)){$C{$a}{$r[$a]}++;}\\
n	  }\n      }\n    @l=keys(%C);\n    \n    \n    \
$s=$P{$seq}{seq};\n    print $f \">$seq pred based\
 on $nseq\\n\";\n    @r=split (//,$s);\n    \n    \
for (my $x=0; $x<=$#r; $x++)\n      {\n	if ($r[$x]\
 ne \"-\")\n	  {\n	    my $h=$C{$x}{H};\n	    my $\
i=$C{$x}{I};\n	    my $o=$C{$x}{O};\n	    my $c=$C\
{$x}{C};\n	    my $e=$C{$x}{E};\n	    my $l=$i+$o;\
\n	    \n	    if ($h>=$i && $h>=$o && $h>=$c && $h\
>=$e){$r[$x]='H';}\n	    elsif ($i>=$o && $i>=$c &\
& $i>=$e){$r[$x]='I';}\n	    elsif ($o>=$c && $o>=\
$e){$r[$x]='O';}\n	    elsif ($c>=$e){$r[$x]='C';}\
\n	    else {$r[$x]='E';}\n	  }\n      }\n    $j=j\
oin ('', @r);\n    print $f \"$j\\n\";\n    close \
($f);\n    return $j;\n  }\n\nsub pred2aln\n  {\n \
   my $PR=shift;\n    my $AR=shift;\n    \n    my \
$f=new FileHandle;\n    my %P=%$PR;\n    my %A=%$A\
R;\n    my %PA;\n    my $tmp=vtmpnam();\n    my $f\
= new FileHandle;\n    \n    open ($f, \">$tmp\");\
\n    foreach my $s (sort{$A{$a}{index}<=>$A{$b}{i\
ndex}}(keys (%A)))\n      {\n	my (@list, $seq, @pl\
ist, @pseq, $L, $PL, $c, $w);\n	my $seq;\n	my $seq\
=$A{$s}{seq};\n	my $pred=$P{$s}{seq};\n	$seq=pred2\
alnS($P{$s}{seq},$A{$s}{seq});\n	print $f \">$s\\n\
$seq\\n\";\n      }\n    close ($f);\n    return r\
ead_fasta_seq ($tmp);\n  }\nsub pred2alnS\n  {\n  \
  my $pred=shift;\n    my $aln= shift;\n    my ($j\
,$a,$b);\n    my @P=split (//, $pred);\n    my @A=\
split (//, $aln);\n    for ($a=$b=0;$a<=$#A; $a++)\
\n      {\n	if ($A[$a] ne \"-\"){$A[$a]=$P[$b++];}\
\n      }\n    if ($b!= ($#P+1)){add_warning (\"Co\
uld not thread sequence: $b $#P\");}\n    \n    $j\
= join ('', @A);\n    return $j;\n  }\nsub pred2co\
lor\n  {\n    my $predP=shift;\n    my $alnP=shift\
;\n    my $out=shift;\n    my $F=new FileHandle;\n\
    my $struc=vtmpnam();\n    my $aln=vtmpnam();\n\
    \n\n    output_fasta_seq ($alnP, $aln);\n    m\
y %p=%$predP;\n    \n    open ($F, \">$struc\");\n\
    \n    \n    foreach my $s (keys(%p))\n      {\\
n	\n	print $F \">$s\\n\";\n	my $s=uc($p{$s}{seq});\
\n	\n	$s=~s/[Oo]/0/g;\n	$s=~s/[Ee]/0/g;\n	\n	$s=~s\
/[Ii]/5/g;\n	$s=~s/[Cc]/5/g;\n	\n	$s=~s/[Hh]/9/g;\\
n	\n	print $F \"$s\\n\";\n      }\n    close ($F);\
\n    \n    \n    \n    safe_system ( \"t_coffee -\
other_pg seq_reformat -in $aln -struc_in $struc -s\
truc_in_f number_fasta -output color_html -out $ou\
t\");\n    return;\n  }\n	  \n    \nsub display_fa\
sta_seq\n  {\n    my $SR=shift;\n    my %S=%$SR;\n\
    \n    foreach my $s (sort{$S{$a}{index}<=>$S{$\
b}{index}}(keys (%S)))\n      {\n	print STDERR \">\
$s\\n$S{$s}{seq}\\n\";\n      }\n    close ($f);\n\
  }\nsub output_fasta_seq\n  {\n    my $SR=shift;\\
n    my $outfile=shift;\n    my $mode =shift;\n   \
 my $f= new FileHandle;\n    my %S=%$SR;\n    \n  \
  \n    open ($f, \">$outfile\");\n    foreach my \
$s (sort{$S{$a}{index}<=>$S{$b}{index}}(keys (%S))\
)\n      {\n	my $seq=$S{$s}{seq};\n	if ( $mode eq \
\"seq\"){$seq=~s/\\-//g;}\n	print $f \">$s\\n$seq\\
\n\";\n      }\n    close ($f);\n  }\n      \nsub \
read_fasta_seq \n  {\n    my $f=$_[0];\n    my %hs\
eq;\n    my (@seq, @com, @name);\n    my ($a, $s,$\
nseq);\n    my $index;\n    open (F, $f);\n    whi\
le (<F>)\n      {\n	$s.=$_;\n      }\n    close (F\
);\n\n    \n    @name=($s=~/>(\\S*).*\\n[^>]*/g);\\
n    \n    @seq =($s=~/>.*.*\\n([^>]*)/g);\n    @c\
om =($s=~/>.*(.*)\\n([^>]*)/g);\n\n\n    $nseq=$#n\
ame+1;\n    \n  \n    for ($a=0; $a<$nseq; $a++)\n\
      {\n	my $n=$name[$a];\n	my $s;\n	$hseq{$n}{na\
me}=$n;\n	$s=$seq[$a];$s=~s/\\s//g;\n	$hseq{$n}{in\
dex}=++$index;\n	$hseq{$n}{seq}=$s;\n	$hseq{$n}{co\
m}=$com[$a];\n      }\n    return %hseq;\n  }\n\n\\
nsub file2head\n      {\n	my $file = shift;\n	my $\
size = shift;\n	my $f= new FileHandle;\n	my $line;\
\n	open ($f,$file);\n	read ($f,$line, $size);\n	cl\
ose ($f);\n	return $line;\n      }\nsub file2tail\\
n      {\n	my $file = shift;\n	my $size = shift;\n\
	my $f= new FileHandle;\n	my $line;\n	\n	open ($f,\
$file);\n	seek ($f,$size*-1, 2);\n	read ($f,$line,\
 $size);\n	close ($f);\n	return $line;\n      }\n\\
n\nsub vtmpnam\n      {\n	my $r=rand(100000);\n	my\
 $f=\"file.$r.$$\";\n	while (-e $f)\n	  {\n	    $f\
=vtmpnam();\n	  }\n	push (@TMPFILE_LIST, $f);\n	re\
turn $f;\n      }\n\nsub myexit\n  {\n    my $code\
=@_[0];\n    if ($CLEAN_EXIT_STARTED==1){return;}\\
n    else {$CLEAN_EXIT_STARTED=1;}\n    ### ONLY B\
ARE EXIT\n    exit ($code);\n  }\nsub set_error_lo\
ck\n    {\n      my $name = shift;\n      my $pid=\
$$;\n\n      \n      &lock4tc ($$,\"LERROR\", \"LS\
ET\", \"$$ -- ERROR: $name $PROGRAM\\n\");\n      \
return;\n    }\nsub set_lock\n  {\n    my $pid=shi\
ft;\n    my $msg= shift;\n    my $p=getppid();\n  \
  &lock4tc ($pid,\"LLOCK\",\"LRESET\",\"$p$msg\\n\\
");\n  }\nsub unset_lock\n   {\n     \n    my $pid\
=shift;\n    &lock4tc ($pid,\"LLOCK\",\"LRELEASE\"\
,\"\");\n  }\nsub shift_lock\n  {\n    my $from=sh\
ift;\n    my $to=shift;\n    my $from_type=shift;\\
n    my $to_type=shift;\n    my $action=shift;\n  \
  my $msg;\n    \n    if (!&lock4tc($from, $from_t\
ype, \"LCHECK\", \"\")){return 0;}\n    $msg=&lock\
4tc ($from, $from_type, \"LREAD\", \"\");\n    &lo\
ck4tc ($from, $from_type,\"LRELEASE\", $msg);\n   \
 &lock4tc ($to, $to_type, $action, $msg);\n    ret\
urn;\n  }\nsub isshellpid\n  {\n    my $p=shift;\n\
    if (!lock4tc ($p, \"LLOCK\", \"LCHECK\")){retu\
rn 0;}\n    else\n      {\n	my $c=lock4tc($p, \"LL\
OCK\", \"LREAD\");\n	if ( $c=~/-SHELL-/){return 1;\
}\n      }\n    return 0;\n  }\nsub isrootpid\n  {\
\n    if(lock4tc (getppid(), \"LLOCK\", \"LCHECK\"\
)){return 0;}\n    else {return 1;}\n  }\nsub lock\
4tc\n	{\n	  my ($pid,$type,$action,$value)=@_;\n	 \
 my $fname;\n	  my $host=hostname;\n	  \n	  if ($t\
ype eq \"LLOCK\"){$fname=\"$LOCKDIR/.$pid.$host.lo\
ck4tcoffee\";}\n	  elsif ( $type eq \"LERROR\"){ $\
fname=\"$LOCKDIR/.$pid.$host.error4tcoffee\";}\n	 \
 elsif ( $type eq \"LWARNING\"){ $fname=\"$LOCKDIR\
/.$pid.$host.warning4tcoffee\";}\n	  \n	  if ($deb\
ug_lock)\n	    {\n	      print STDERR \"\\n\\t---l\
ock4tc(tcg): $action => $fname =>$value (RD: $LOCK\
DIR)\\n\";\n	    }\n\n	  if    ($action eq \"LCHEC\
K\") {return -e $fname;}\n	  elsif ($action eq \"L\
READ\"){return file2string($fname);}\n	  elsif ($a\
ction eq \"LSET\") {return string2file ($value, $f\
name, \">>\");}\n	  elsif ($action eq \"LRESET\") \
{return string2file ($value, $fname, \">\");}\n	  \
elsif ($action eq \"LRELEASE\") \n	    {\n	      i\
f ( $debug_lock)\n		{\n		  my $g=new FileHandle;\n\
		  open ($g, \">>$fname\");\n		  print $g \"\\nDe\
stroyed by $$\\n\";\n		  close ($g);\n		  safe_sys\
tem (\"mv $fname $fname.old\");\n		}\n	      else\\
n		{\n		  unlink ($fname);\n		}\n	    }\n	  return\
 \"\";\n	}\n	\nsub file2string\n	{\n	  my $file=@_\
[0];\n	  my $f=new FileHandle;\n	  my $r;\n	  open\
 ($f, \"$file\");\n	  while (<$f>){$r.=$_;}\n	  cl\
ose ($f);\n	  return $r;\n	}\nsub string2file \n  \
  {\n    my ($s,$file,$mode)=@_;\n    my $f=new Fi\
leHandle;\n    \n    open ($f, \"$mode$file\");\n \
   print $f  \"$s\";\n    close ($f);\n  }\n\nBEGI\
N\n    {\n      srand;\n    \n      $SIG{'SIGUP'}=\
'signal_cleanup';\n      $SIG{'SIGINT'}='signal_cl\
eanup';\n      $SIG{'SIGQUIT'}='signal_cleanup';\n\
      $SIG{'SIGILL'}='signal_cleanup';\n      $SIG\
{'SIGTRAP'}='signal_cleanup';\n      $SIG{'SIGABRT\
'}='signal_cleanup';\n      $SIG{'SIGEMT'}='signal\
_cleanup';\n      $SIG{'SIGFPE'}='signal_cleanup';\
\n      \n      $SIG{'SIGKILL'}='signal_cleanup';\\
n      $SIG{'SIGPIPE'}='signal_cleanup';\n      $S\
IG{'SIGSTOP'}='signal_cleanup';\n      $SIG{'SIGTT\
IN'}='signal_cleanup';\n      $SIG{'SIGXFSZ'}='sig\
nal_cleanup';\n      $SIG{'SIGINFO'}='signal_clean\
up';\n      \n      $SIG{'SIGBUS'}='signal_cleanup\
';\n      $SIG{'SIGALRM'}='signal_cleanup';\n     \
 $SIG{'SIGTSTP'}='signal_cleanup';\n      $SIG{'SI\
GTTOU'}='signal_cleanup';\n      $SIG{'SIGVTALRM'}\
='signal_cleanup';\n      $SIG{'SIGUSR1'}='signal_\
cleanup';\n\n\n      $SIG{'SIGSEGV'}='signal_clean\
up';\n      $SIG{'SIGTERM'}='signal_cleanup';\n   \
   $SIG{'SIGCONT'}='signal_cleanup';\n      $SIG{'\
SIGIO'}='signal_cleanup';\n      $SIG{'SIGPROF'}='\
signal_cleanup';\n      $SIG{'SIGUSR2'}='signal_cl\
eanup';\n\n      $SIG{'SIGSYS'}='signal_cleanup';\\
n      $SIG{'SIGURG'}='signal_cleanup';\n      $SI\
G{'SIGCHLD'}='signal_cleanup';\n      $SIG{'SIGXCP\
U'}='signal_cleanup';\n      $SIG{'SIGWINCH'}='sig\
nal_cleanup';\n      \n      $SIG{'INT'}='signal_c\
leanup';\n      $SIG{'TERM'}='signal_cleanup';\n  \
    $SIG{'KILL'}='signal_cleanup';\n      $SIG{'QU\
IT'}='signal_cleanup';\n      \n      our $debug_l\
ock=$ENV{\"DEBUG_LOCK\"};\n      \n      \n      \\
n      \n      foreach my $a (@ARGV){$CL.=\" $a\";\
}\n      if ( $debug_lock ){print STDERR \"\\n\\n\\
\n********** START PG: $PROGRAM *************\\n\"\
;}\n      if ( $debug_lock ){print STDERR \"\\n\\n\
\\n**********(tcg) LOCKDIR: $LOCKDIR $$ **********\
***\\n\";}\n      if ( $debug_lock ){print STDERR \
\"\\n --- $$ -- $CL\\n\";}\n      \n	     \n      \
\n      \n    }\nsub flush_error\n  {\n    my $msg\
=shift;\n    return add_error ($EXIT_FAILURE,$$, $\
$,getppid(), $msg, $CL);\n  }\nsub add_error \n  {\
\n    my $code=shift;\n    my $rpid=shift;\n    my\
 $pid=shift;\n    my $ppid=shift;\n    my $type=sh\
ift;\n    my $com=shift;\n    \n    $ERROR_DONE=1;\
\n    lock4tc ($rpid, \"LERROR\",\"LSET\",\"$pid -\
- ERROR: $type\\n\");\n    lock4tc ($$, \"LERROR\"\
,\"LSET\", \"$pid -- COM: $com\\n\");\n    lock4tc\
 ($$, \"LERROR\",\"LSET\", \"$pid -- STACK: $ppid \
-> $pid\\n\");\n   \n    return $code;\n  }\nsub a\
dd_warning \n  {\n    my $rpid=shift;\n    my $pid\
 =shift;\n    my $command=shift;\n    my $msg=\"$$\
 -- WARNING: $command\\n\";\n    print STDERR \"$m\
sg\";\n    lock4tc ($$, \"LWARNING\", \"LSET\", $m\
sg);\n  }\n\nsub signal_cleanup\n  {\n    print dt\
derr \"\\n**** $$ (tcg) was killed\\n\";\n    &cle\
anup;\n    exit ($EXIT_FAILURE);\n  }\nsub clean_d\
ir\n  {\n    my $dir=@_[0];\n    if ( !-d $dir){re\
turn ;}\n    elsif (!($dir=~/tmp/)){return ;}#safe\
ty check 1\n    elsif (($dir=~/\\*/)){return ;}#sa\
fety check 2\n    else\n      {\n	`rm -rf $dir`;\n\
      }\n    return;\n  }\nsub cleanup\n  {\n    #\
print stderr \"\\n----tc: $$ Kills $PIDCHILD\\n\";\
\n    #kill (SIGTERM,$PIDCHILD);\n    my $p=getppi\
d();\n    $CLEAN_EXIT_STARTED=1;\n    \n    \n    \
\n    if (&lock4tc($$,\"LERROR\", \"LCHECK\", \"\"\
))\n      {\n	my $ppid=getppid();\n	if (!$ERROR_DO\
NE) \n	  {\n	    &lock4tc($$,\"LERROR\", \"LSET\",\
 \"$$ -- STACK: $p -> $$\\n\");\n	    &lock4tc($$,\
\"LERROR\", \"LSET\", \"$$ -- COM: $CL\\n\");\n	  \
}\n      }\n    my $warning=&lock4tc($$, \"LWARNIN\
G\", \"LREAD\", \"\");\n    my $error=&lock4tc($$,\
  \"LERROR\", \"LREAD\", \"\");\n    #release erro\
r and warning lock if root\n    \n    if (isrootpi\
d() && ($warning || $error) )\n      {\n	\n	print \
STDERR \"**************** Summary *************\\n\
$error\\n$warning\\n\";\n\n	&lock4tc($$,\"LERROR\"\
,\"RELEASE\",\"\");\n	&lock4tc($$,\"LWARNING\",\"R\
ELEASE\",\"\");\n      } \n    \n    \n    foreach\
 my $f (@TMPFILE_LIST)\n      {\n	if (-e $f){unlin\
k ($f);} \n      }\n    foreach my $d (@TMPDIR_LIS\
T)\n      {\n	clean_dir ($d);\n      }\n    #No Mo\
re Lock Release\n    #&lock4tc($$,\"LLOCK\",\"LREL\
EASE\",\"\"); #release lock \n\n    if ( $debug_lo\
ck ){print STDERR \"\\n\\n\\n********** END PG: $P\
ROGRAM ($$) *************\\n\";}\n    if ( $debug_\
lock ){print STDERR \"\\n\\n\\n**********(tcg) LOC\
KDIR: $LOCKDIR $$ *************\\n\";}\n  }\nEND \\
n  {\n    \n    &cleanup();\n  }\n   \n\nsub safe_\
system \n{\n  my $com=shift;\n  my $ntry=shift;\n \
 my $ctry=shift;\n  my $pid;\n  my $status;\n  my \
$ppid=getppid();\n  if ($com eq \"\"){return 1;}\n\
  \n  \n\n  if (($pid = fork ()) < 0){return (-1);\
}\n  if ($pid == 0)\n    {\n      set_lock($$, \" \
-SHELL- $com (tcg)\");\n      exec ($com);\n    }\\
n  else\n    {\n      lock4tc ($$, \"LLOCK\", \"LS\
ET\", \"$pid\\n\");#update parent\n      $PIDCHILD\
=$pid;\n    }\n  if ($debug_lock){printf STDERR \"\
\\n\\t .... safe_system (fasta_seq2hmm)  p: $$ c: \
$pid COM: $com\\n\";}\n\n  waitpid ($pid,WTERMSIG)\
;\n\n  shift_lock ($pid,$$, \"LWARNING\",\"LWARNIN\
G\", \"LSET\");\n\n  if ($? == $EXIT_FAILURE || lo\
ck4tc($pid, \"LERROR\", \"LCHECK\", \"\"))\n    {\\
n      if ($ntry && $ctry <$ntry)\n	{\n	  add_warn\
ing ($$,$$,\"$com failed [retry: $ctry]\");\n	  lo\
ck4tc ($pid, \"LRELEASE\", \"LERROR\", \"\");\n	  \
return safe_system ($com, $ntry, ++$ctry);\n	}\n  \
    elsif ($ntry == -1)\n	{\n	  if (!shift_lock ($\
pid, $$, \"LERROR\", \"LWARNING\", \"LSET\"))\n	  \
  {\n	      add_warning ($$,$$,\"$com failed\");\n\
	    }\n	  else\n	    {\n	      lock4tc ($pid, \"L\
RELEASE\", \"LERROR\", \"\");\n	    }\n	  return $\
?;}\n      else\n	{\n	  if (!shift_lock ($pid,$$, \
\"LERROR\",\"LERROR\", \"LSET\"))\n	    {\n	      \
myexit(add_error ($EXIT_FAILURE,$$,$pid,getppid(),\
 \"UNSPECIFIED system\", $com));\n	    }\n	}\n    \
}\n  return $?;\n}\n\nsub check_configuration \n  \
  {\n      my @l=@_;\n      my $v;\n      foreach \
my $p (@l)\n	{\n	  \n	  if   ( $p eq \"EMAIL\")\n	\
    { \n	      if ( !($EMAIL=~/@/))\n		{\n		add_wa\
rning($$,$$,\"Could Not Use EMAIL\");\n		myexit(ad\
d_error ($EXIT_FAILURE,$$,$$,getppid(),\"EMAIL\",\\
"$CL\"));\n	      }\n	    }\n	  elsif( $p eq \"INT\
ERNET\")\n	    {\n	      if ( !&check_internet_con\
nection())\n		{\n		  myexit(add_error ($EXIT_FAILU\
RE,$$,$$,getppid(),\"INTERNET\",\"$CL\"));\n		}\n	\
    }\n	  elsif( $p eq \"wget\")\n	    {\n	      i\
f (!&pg_is_installed (\"wget\") && !&pg_is_install\
ed (\"curl\"))\n		{\n		  myexit(add_error ($EXIT_F\
AILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:wget\",\\
"$CL\"));\n		}\n	    }\n	  elsif( !(&pg_is_install\
ed ($p)))\n	    {\n	      myexit(add_error ($EXIT_\
FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:$p\",\"\
$CL\"));\n	    }\n	}\n      return 1;\n    }\nsub \
pg_is_installed\n  {\n    my @ml=@_;\n    my $r, $\
p, $m;\n    my $supported=0;\n    \n    my $p=shif\
t (@ml);\n    if ($p=~/::/)\n      {\n	if (safe_sy\
stem (\"perl -M$p -e 1\")==$EXIT_SUCCESS){return 1\
;}\n	else {return 0;}\n      }\n    else\n      {\\
n	$r=`which $p 2>/dev/null`;\n	if ($r eq \"\"){ret\
urn 0;}\n	else {return 1;}\n      }\n  }\n\n\n\nsu\
b check_internet_connection\n  {\n    my $internet\
;\n    my $tmp;\n    &check_configuration ( \"wget\
\"); \n    \n    $tmp=&vtmpnam ();\n    \n    if  \
   (&pg_is_installed    (\"wget\")){`wget www.goog\
le.com -O$tmp >/dev/null 2>/dev/null`;}\n    elsif\
  (&pg_is_installed    (\"curl\")){`curl www.googl\
e.com -o$tmp >/dev/null 2>/dev/null`;}\n    \n    \
if ( !-e $tmp || -s $tmp < 10){$internet=0;}\n    \
else {$internet=1;}\n    if (-e $tmp){unlink $tmp;\
}\n\n    return $internet;\n  }\nsub check_pg_is_i\
nstalled\n  {\n    my @ml=@_;\n    my $r=&pg_is_in\
stalled (@ml);\n    if (!$r && $p=~/::/)\n      {\\
n	print STDERR \"\\nYou Must Install the perl pack\
age $p on your system.\\nRUN:\\n\\tsudo perl -MCPA\
N -e 'install $pg'\\n\";\n      }\n    elsif (!$r)\
\n      {\n	myexit(flush_error(\"\\nProgram $p Sup\
ported but Not Installed on your system\"));\n    \
  }\n    else\n      {\n	return 1;\n      }\n  }\n\
\n\n\n","\n\n\n\n\nmy $FMODEL =\"\"; \nmy $TMPDIR \
= \"/tmp\";\n\n\n\n\nmy $NUCALPH = \"ACGTUNRYMKSWH\
BVD\";\nmy $PRIMNUCALPH = \"ACGTUN\";\nuse vars qw\
($NUCALPH $PRIMNUCALPH $TMPDIR);\n\n\nmy $errmsg;\\
nuse vars qw($errmsg);\n\n\n\nuse Getopt::Long;\nu\
se Cwd;\nuse File::Basename;\nuse File::Temp qw/ t\
empfile tempdir /;\nuse File::Copy;\nuse File::Pat\
h;\n\n\n\nsub usage(;$)\n{\n    my ($errmsg) = @_;\
\n    my $myname = basename($0);\n\n    if ($errms\
g) {\n        print STDERR \"ERROR: $errmsg\\n\";\\
n    }\n\n    print STDERR << \"EOF\";\n    \n$myn\
ame: align two sequences by means of consan\\'s sf\
old\nUsage:\n $myname -i file -o file -d path\nOpt\
ions:\n -i|--in : pairwise input sequence file\n -\
o|--out: output alignment\n -d|--directory contain\
ing data\n\nEOF\n}\n\nsub read_stk_aln \n  {\n    \
my $f=$_[0];\n    my ($seq, $id);\n    \n    my %h\
seq;\n\n    open (STK, \"$f\");\n    while (<STK>)\
\n      {\n	if ( /^#/ || /^\\/\\// || /^\\s*$/){;}\
\n	else\n	  {\n	    ($id,$seq)=/(\\S+)\\s+(\\S+)/;\
\n	    $hseq{$id}{'seq'}.=$seq;\n	  }\n      }\n  \
  close (STK);\n    return %hseq;\n  }\nsub read_f\
asta_seq \n  {\n    my $f=$_[0];\n    my %hseq;\n \
   my (@seq, @com, @name);\n    my ($a, $s,$nseq);\
\n\n    open (F, $f);\n    while (<F>)\n      {\n	\
$s.=$_;\n      }\n    close (F);\n\n    \n    @nam\
e=($s=~/>(.*).*\\n[^>]*/g);\n    \n    @seq =($s=~\
/>.*.*\\n([^>]*)/g);\n    @com =($s=~/>.*(.*)\\n([\
^>]*)/g);\n\n    \n    $nseq=$#name+1;\n    \n    \
for ($a=0; $a<$nseq; $a++)\n      {\n	my $n=$name[\
$a];\n	$hseq{$n}{name}=$n;\n	$hseq{$n}{seq}=$seq[$\
a];\n	$hseq{$n}{com}=$com[$a];\n      }\n    retur\
n %hseq;\n  }\n\n\n\nsub sfold_parseoutput($$)\n{\\
n    my ($frawout, $foutfa) = @_;\n    my %haln;\n\
    my ($fstk, $cmd, $id);\n    open FOUTFA, \">$f\
outfa\";\n    \n    $fstk = $frawout . \".stk\";\n\
    \n    # first line of raw out contains info\n \
   # remaining stuff is stockholm formatted\n    $\
cmd = \"sed -e '1d' $frawout\";\n    system(\"$cmd\
 > $fstk\");\n    if ($? != 0) {\n        $errmsg \
= \"command failed with exit status $?.\";\n      \
  $errmsg .=  \"Command was \\\"$cmd\\\"\";\n     \
   return -1;\n    }\n\n    # this gives an error \
message. just ignore it...\n    %haln=read_stk_aln\
 ( $fstk);\n    foreach $i (keys (%haln))\n      {\
\n	my $s;\n	$s=$haln{$i}{'seq'};\n	$s =~ s/\\./-/g\
;\n	print FOUTFA \">$i\\n$s\\n\";\n      }\n    cl\
ose FOUTFA;\n    return 0;\n}\n\n\n\n\nsub sfold_w\
rapper($$$$)\n{\n    \n    my ($fs1, $fs2, $fmodel\
, $foutfa) = @_;\n    \n\n    my ($cmd, $frawout, \
$ferrlog, $freadme, $ftimelog, $fstk);\n\n    # ad\
d  basename($fmsqin) (unknown here!)\n    $frawout\
 = \"sfold.log\";\n    $ferrlog = \"sfold.err\";\n\
    $ftimelog = \"sfold.time\";\n    $freadme =  \\
"sfold.README\";\n    $fstk = \"sfold.stk\";\n    \
\n    # prepare execution...\n    #\n    # ./tmp i\
s essential for dswpalign\n    # otherwise you'll \
get a segfault\n    mkdir \"./tmp\";\n    \n    $c\
md = \"sfold -m $fmodel $fs1 $fs2\";\n    open(FRE\
ADME,\">$freadme\");\n    print FREADME \"$cmd\\n\\
"; \n    close(FREADME);\n\n    # and go\n    #\n \
   system(\"/usr/bin/time -p -o $ftimelog $cmd >$f\
rawout 2>$ferrlog\");\n    if ($? != 0) {\n       \
 $errmsg = \"command failed with exit status $?\";\
\n        $errmsg .= \"command was \\\"$cmd\\\". S\
ee \" . getcwd . \"\\n\";\n        return -1;\n   \
 }\n\n    return sfold_parseoutput($frawout, $fout\
fa);\n}\n\n\n\n\n\n\n\nmy ($help, $fmsqin, $fmsaou\
t);\nGetOptions(\"help\"  => \\$help,\n           \
\"in=s\" => \\$fmsqin,\n           \"out=s\" => \\\
$fmsaout,\n	   \"data=s\" => \\$ref_dir);\n\n\n\ni\
f ($help) {\n    usage();\n    exit(0);\n}\nif (! \
defined($fmsqin)) {\n    usage('missing input file\
name');\n    exit(1);\n}\nif (! defined($fmsaout))\
 {\n    usage('missing output filename');\n    exi\
t(1);\n\n}\nif (scalar(@ARGV)) {\n    usage('Unkno\
wn remaining args');\n    exit(1);\n}\n\n$FMODEL =\
 \"$ref_dir/mix80.mod\";\nif (! -e \"$FMODEL\") {\\
n    die(\"couldn't find sfold grammar model file.\
 Expected $FMODEL\\n\");\n}\n\n\nmy %hseq=read_fas\
ta_seq ($fmsqin);\nmy $id;\n\nforeach $id (keys(%h\
seq))\n  {\n    push(@seq_array, $hseq{$id});\n  }\
\n\nif ( scalar(@seq_array) != 2 ) {\n    die(\"Ne\
ed *exactly* two sequences as input (pairwise alig\
nment!).\")\n}\n\n\n\nmy ($sec, $min, $hour, $mday\
, $mon, $year, $wday, $yday, $isdst) = localtime(t\
ime);\nmy $datei = sprintf(\"%4d-%02d-%02d\", $yea\
r+1900, $mon+1, $mday);\nmy $templ = basename($0) \
. \".\" . $datei . \".pid-\" . $$ . \".XXXXXX\";\n\
my $wd = tempdir ( $templ, DIR => $TMPDIR);\n\ncop\
y($fmsqin, \"$wd/\" . basename($fmsqin) . \".org\"\
); # for reproduction\ncopy($FMODEL, \"$wd\");\nmy\
 $fmodel = basename($FMODEL);\nmy $orgwd = getcwd;\
\nchdir $wd;\n\n\n\nmy @sepseqfiles;\nforeach $id \
(keys(%hseq)) {\n    my ($seq, $orgseq, $fname, $s\
out);\n    $seq=$hseq{$id}{'seq'};\n    \n    $fna\
me = basename($fmsqin) . \"_$id.fa\";\n    # repla\
ce funnies in file/id name (e.g. \"/\" \" \" etc)\\
n    $fname =~ s,[/ ],_,g;\n    open (PF, \">$fnam\
e\");\n    print (PF \">$id\\n$seq\\n\");\n    clo\
se (PF);\n\n    push(@sepseqfiles, $fname);\n}\n\n\
my ($f1, $f2, $fout);\n$f1 = $sepseqfiles[0];\n$f2\
 = $sepseqfiles[1];\n$fout = $wd . basename($fmsqi\
n) . \".out.fa\";\nif (sfold_wrapper($f1, $f2, $fm\
odel, \"$fout\") != 0) {\n    printf STDERR \"ERRO\
R: See logs in $wd\\n\";\n    exit(1);\n} else {\n\
    chdir $orgwd;\n    copy($fout, $fmsaout);\n   \
 rmtree($wd);\n   exit(0);\n}\n","\nuse Env qw(HOS\
T);\nuse Env qw(HOME);\nuse Env qw(USER);\n\n\n$tm\
p=clean_cr ($ARGV[0]);\nopen (F, $tmp);\n\nwhile (\
 <F>)\n  {\n    my $l=$_;\n    if ( $l=~/^# STOCKH\
OLM/){$stockholm=1;}\n    elsif ( $stockholm && $l\
=~/^#/)\n      {\n	$l=~/^#(\\S+)\\s+(\\S+)\\s+(\\S\
*)/g;\n	$l=\"_stockholmhasch_$1\\_stockholmspace_$\
2 $3\\n\";\n      }\n    $file.=$l;\n  }\nclose (F\
);\nunlink($tmp);\n$file1=$file;\n\n$file=~s/\\#/_\
hash_symbol_/g;\n$file=~s/\\@/_arobase_symbol_/g;\\
n\n\n$file=~s/\\n[\\.:*\\s]+\\n/\\n\\n/g;\n\n$file\
=~s/\\n[ \\t\\r\\f]+(\\b)/\\n\\1/g;\n\n\n$file=~s/\
(\\n\\S+)(\\s+)(\\S)/\\1_blank_\\3/g;\n\n$file=~s/\
[ ]//g;\n$file=~s/_blank_/ /g;\n\n\n\n$file =~s/\\\
n\\s*\\n/#/g;\n\n$file.=\"#\";\n$file =~s/\\n/@/g;\
\n\n\n\n\n@blocks=split /\\#/, $file;\nshift (@blo\
cks);\n@s=split /\\@/, $blocks[0];\n$nseq=$#s+1;\n\
\n\n\n$file=join '@', @blocks;\n@lines=split /\\@/\
,$file;\n\n$c=0;\n\nforeach $l (@lines)\n  {\n    \
if (!($l=~/\\S/)){next;}\n    elsif ($stockholm &&\
 ($l=~/^\\/\\// || $l=~/STOCKHOLM/)){next;}#get re\
ad of STOCHOLM Terminator\n   \n    $l=~/(\\S+)\\s\
+(\\S*)/g;\n    $n=$1; $s=$2;\n    \n    $seq[$c].\
=$s;\n    $name[$c]=$n;\n    $c++;\n    \n    if (\
 $c==$nseq){$c=0;}\n    \n  } \n\nif ( $c!=0)\n   \
   {\n	print STDERR \"ERROR: $ARGV[0] is NOT an MS\
A in Clustalw format: make sure there is no blank \
line within a block [ERROR]\\n\";\n	exit (EXIT_FAI\
LURE);\n      }\n\nfor ($a=0; $a< $nseq; $a++)\n  \
{\n    $name[$a]=cleanstring ($name[$a]);\n    $se\
q[$a]=cleanstring ($seq[$a]);\n    $seq[$a]=breaks\
tring($seq[$a], 60);\n    \n    $line=\">$name[$a]\
\\n$seq[$a]\\n\";\n    \n    print \"$line\";\n  }\
\nexit (EXIT_SUCCESS);\n\nsub cleanstring\n  {\n  \
  my $s=@_[0];\n    $s=~s/_hash_symbol_/\\#/g;\n  \
  $s=~s/_arobase_symbol_/\\@/g;\n    $s=~s/[ \\t]/\
/g;\n    return $s;\n  }\nsub breakstring\n  {\n  \
  my $s=@_[0];\n    my $size=@_[1];\n    my @list;\
\n    my $n,$ns, $symbol;\n    \n    @list=split /\
/,$s;\n    $n=0;$ns=\"\";\n    foreach $symbol (@l\
ist)\n      {\n	if ( $n==$size)\n	  {\n	    $ns.=\\
"\\n\";\n	    $n=0;\n	  }\n	$ns.=$symbol;\n	$n++;\\
n      }\n    return $ns;\n    }\n\nsub clean_cr\n\
  {\n    my $f=@_[0];\n    my $file;\n    \n    $t\
mp=\"f$.$$\";\n    \n    \n    open (IN, $f);\n   \
 open (OUT, \">$tmp\");\n    \n    while ( <IN>)\n\
      {\n	$file=$_;\n	$file=~s/\\r\\n/\\n/g;\n	$fi\
le=~s/\\n\\r/\\n/g;\n	$file=~s/\\r\\r/\\n/g;\n	$fi\
le=~s/\\r/\\n/g;\n	print OUT \"$file\";\n      }\n\
    \n    close (IN);\n    close (OUT);\n    retur\
n $tmp;\n  }\n","use strict;\nuse FileHandle;\nuse\
 Env qw(HOST);\nuse Env qw(HOME);\nuse Env qw(USER\
);\n\nmy $format=file2format ($ARGV[0]);\nif    ($\
format eq \"clustalw\"){clustalw2name_seq($ARGV[0]\
);}\nelsif ($format eq \"fasta\")   {fasta2name_se\
q($ARGV[0]);}\nelsif ($format eq \"msf\")   {msf2n\
ame_seq($ARGV[0]);}\nelsif ($format eq \"phylip\")\
   {phylip2name_seq($ARGV[0]);}\nelsif ($format eq\
 \"nameseq\") {display_file ($ARGV[0]);}\n \nexit \
(0);\n\nsub file2format\n  {\n    my $f=shift;\n  \
  \n    my $l=file2n_lines($f,2);\n    \n    if ( \
$l=~/^CLUSTAL/){return \"clustalw\";}\n    elsif (\
$l=~/^SAGA/){return \"clustalw\";}\n    elsif ($l=\
~/^>/){return \"fasta\";}\n    elsif ($l=~/^PileUp\
/){return \"msf\";}\n    elsif ($l=~/\\s+\\d+\\s+\\
\d+\\s/){return \"phylip\";}\n    elsif ($l=~/\\#N\
AMESEQ_01/){return \"nameseq\";}\n    else \n     \
 {\n	print STDERR \"ERROR: $f FILE is NOT a suppor\
ted format [ERROR]\\n\";\n	exit (1);\n      }\n  }\
\nsub display_file\n    {\n       my $file=shift;\\
n       my $F= new FileHandle;\n       open ($F, $\
file);\n       while (<$F>){print \"$_\";}\n      \
 close ($F);\n     }\nsub phylip2name_seq\n    {\n\
      my $file=shift;\n      my $F= new FileHandle\
;\n      my ($seq, $name,$seq);\n      my $query_s\
tart=-1;\n      my $query_end=-1;\n      my $in_al\
n=0;\n      my %list;\n      my ($first,$seq,$name\
, $cn, $nseq, $l,%len);\n      \n      open ($F, $\
file);\n      <$F>;\n      $l=$_;\n      $l=~/\\s*\
(\\d+)\\s*(\\d+)/;\n      $first=1;\n      $cn=0;\\
n      while (<$F>)\n	{\n	  my $l=$_;\n	  if (!($l\
=~/\\S/))\n	    {\n	      $cn=0;\n	      $first=0;\
\n	    }\n	  elsif ($first==1)\n	    {\n	      $l=\
~/\\s*(\\S+)(.*)/;\n	      my $name=$1;\n	      my\
 $seq=$2;\n	      chomp ($seq);\n	      $seq=~s/\\\
s//g;\n	      $list{$cn}{'name'}=$name;\n	      $l\
ist{$cn}{'seq'}.=$seq;\n	      $cn++;\n	      $nse\
q++;\n	    }\n	  else\n	    {\n	      chomp ($l);\\
n	      $l=~s/\\s//g;\n	      $list{$cn}{'seq'}.=$\
l;\n	      $cn++;\n	    }\n	}\n      close ($F);\n\
      \n      print \"#NAMESEQ_01\\n\";\n      pri\
nt \"# $nseq\\n\";\n      for (my $a=0; $a<$nseq; \
$a++)\n	{\n	  my $nl=length ($list{$a}{'name'});\n\
	  my $sl=length ($list{$a}{'seq'});\n	  print \">\
$nl $sl $list{$a}{'name'} $list{$a}{'seq'}\\n\";\n\
	}\n    }\n      \nsub msf2name_seq\n    {\n      \
my $file=shift;\n      my $F= new FileHandle;\n   \
   my ($seq, $name,$seq);\n      my $query_start=-\
1;\n      my $query_end=-1;\n      my $in_aln=0;\n\
      my %list;\n      my ($seq,$name, $n, $nseq, \
$l,%len);\n      \n      open ($F, $file);\n      \
while (<$F>)\n	{\n	  if ( /\\/\\//){$in_aln=1;}\n	\
  elsif ( $in_aln && /(\\S+)\\s+(.*)/)\n	    {\n	 \
     $name=$1;\n	      $seq=$2;\n	      $seq=~s/\\\
s//g;\n	      $seq=~s/\\~/\\-/g;\n	      $seq=~s/\\
\./\\-/g;\n	      if ( $list{$n}{'name'} && $list{\
$n}{'name'} ne $name)\n		{\n		  print \"$list{$n}{\
'name'} Vs $name\";\n		  \n		  exit (1);\n		}\n	  \
    else\n		{\n		  $list{$n}{'name'}= $name;\n		}\\
n	      \n	      $list{$n}{'seq'}=$list{$n}{'seq'}\
.$seq;\n	      \n	      $nseq=++$n;\n	      \n	   \
 }\n	  else\n	    {$n=0;}\n	}\n      close ($F);\n\
      print \"#NAMESEQ_01\\n\";\n      print \"# $\
nseq\\n\";\n      for (my $a=0; $a<$nseq; $a++)\n	\
{\n	  my $nl=length ($list{$a}{'name'});\n	  my $s\
l=length ($list{$a}{'seq'});\n	  print \">$nl $sl \
$list{$a}{'name'} $list{$a}{'seq'}\\n\";\n	}\n    \
}\n    \nsub fasta2name_seq\n    {\n      my $file\
=shift;\n      my $F= new FileHandle;\n      my ($\
seq, $name,$n,$l,%len);\n      \n      open ($F, $\
file);\n      while (<$F>)\n	{\n	  if ( /^>(\\S+)/\
){$n++;$seq=\"\";$name=$1;}\n	  else\n	    {\n	   \
   $l=$_;\n	      chomp ($l);\n	      $seq.=$l;\n	\
      $len{$name}=length($seq);\n	    }\n	}\n     \
 close ($F);\n      print \"#NAMESEQ_01\\n\";\n   \
   print \"# $n\";\n      \n      open ($F, $file)\
;\n      while (<$F>)\n	{\n	  if ( /^>(\\S+)(.*)\\\
n/)\n	    {\n	      my $name=$1;\n	      my $comme\
nt=$2;\n	      my $nl=length ($name);\n	      my $\
sl=$len{$name};\n	      if ($comment)\n		{\n		  $c\
omment=~s/^\\s+//g;\n		  my $cl=length ($comment);\
\n		  print \"\\n#$cl $comment\\n\";\n		}\n	      \
print \"\\n>$nl $sl $name \";\n	    }\n	  else\n	 \
   {\n	      $l=$_;\n	      chomp ($l);\n	      pr\
int \"$l\";\n	    }\n	}\n      print \"\\n\";\n   \
   close ($F);\n    }\nsub clustalw2name_seq\n  {\\
n    my $fname=shift;\n    my ($file1, $file);\n  \
  my (@blocks, @lines,@s, $n,$nseq, $c);\n    my (\
@name, @seq);\n    my $F= new FileHandle;\n    my \
$stockholm;\n    \n    open ($F, $fname);\n    \n \
   while ( <$F>)\n      {\n	my $l=$_;\n	$l=clean_c\
r($l);\n	if ( $l=~/^# STOCKHOLM/){$stockholm=1;}\n\
	elsif ( $stockholm && $l=~/^#/)\n	  {\n	    $l=~/\
^#(\\S+)\\s+(\\S+)\\s+(\\S*)/g;\n	    $l=\"_stockh\
olmhasch_$1\\_stockholmspace_$2 $3\\n\";\n	  }\n	$\
file.=$l;\n      }\n    close ($F);\n        \n   \
 #Protect # and @\n    $file=~s/\\#/_hash_symbol_/\
g;\n    $file=~s/\\@/_arobase_symbol_/g;\n    \n  \
  \n    #Remove annotation\n    $file=~s/\\n[\\.:*\
\\s]+\\n/\\n\\n/g;\n    \n    #Remove White spaces\
 before the sequence name\n    $file=~s/\\n[ \\t\\\
r\\f]+(\\b)/\\n\\1/g;\n    \n    \n    #Remove Int\
ernal Blanks\n    $file=~s/(\\n\\S+)(\\s+)(\\S)/\\\
1_blank_\\3/g;\n    \n    $file=~s/[ ]//g;\n    $f\
ile=~s/_blank_/ /g;\n    \n    \n    #Identify Dou\
ble Blank lines\n    \n    $file =~s/\\n\\s*\\n/#/\
g;\n    \n    $file.=\"#\";\n    $file =~s/\\n/@/g\
;\n    \n    \n    \n    \n    #count nseq\n    @b\
locks=split /\\#/, $file;\n    shift (@blocks);\n \
   @s=split /\\@/, $blocks[0];\n    $nseq=$#s+1;\n\
    \n    #Merge all the sequences and split every\
 Nseq\n    \n    \n    $file=join '@', @blocks;\n \
   @lines=split /\\@/,$file;\n    \n    $c=0;\n   \
 \n    foreach my $l (@lines)\n      {\n	my ($n, $\
s);\n	\n	if (!($l=~/\\S/)){next;}\n	elsif ($stockh\
olm && ($l=~/^\\/\\// || $l=~/STOCKHOLM/)){next;}#\
get read of STOCHOLM Terminator\n	\n	$l=~/(\\S+)\\\
s+(\\S*)/g;\n	$n=$1; $s=$2;\n	\n	$seq[$c].=$s;\n	$\
name[$c]=$n;\n	$c++;\n	\n	if ( $c==$nseq){$c=0;}\n\
	\n      } \n    \n    if ( $c!=0)\n      {\n	prin\
t STDERR \"ERROR: $fname is NOT an MSA in Clustalw\
 format: make sure there is no blank line within a\
 block [ERROR]\\n\";\n	exit (1);\n      }\n    pri\
nt \"#NAMESEQ_01\\n\";\n    print \"# $nseq\\n\";\\
n    for (my $a=0; $a< $nseq; $a++)\n      {\n	$na\
me[$a]=cleanstring ($name[$a]);\n	$seq[$a]=cleanst\
ring ($seq[$a]);\n	my $ln=length ($name[$a]);\n	my\
 $ls=length ($seq[$a]);\n	print \">$ln $ls $name[$\
a] $seq[$a]\\n\";\n      }\n  }\nsub cleanstring\n\
    {\n      my $s=@_[0];\n      $s=~s/_hash_symbo\
l_/\\#/g;\n      $s=~s/_arobase_symbol_/\\@/g;\n  \
    $s=~s/[ \\t]//g;\n      return $s;\n    }\n\ns\
ub clean_cr\n  {\n    my $f=shift;\n    $f=~s/\\r\\
\n/\\n/g;\n    $f=~s/\\n\\r/\\n/g;\n    $f=~s/\\r\\
\r/\\n/g;\n    $f=~s/\\r/\\n/g;\n    return $f;\n \
 }\n\nsub file2n_lines\n    {\n      my $file=shif\
t;\n      my $nl=shift;\n      my $ret;\n      my \
$F=new FileHandle;\n      my $n=0;\n      open ($F\
, $file);\n\n      while (<$F>)\n	{\n	  $ret.=$_;\\
n	  $n++;\n	  \n	  if ($n>=$n){close ($F); return \
$ret;}\n	}\n      close ($F);\n      return $ret;\\
n    }\n","use strict;\nuse FileHandle;\nuse Env q\
w(HOST);\nuse Env qw(HOME);\nuse Env qw(USER);\nmy\
 %name;\nmy $nseq;\nmy $fasta;\nif ($ARGV[2] eq \"\
-fasta\"){$fasta=1;}\nmy $F= new FileHandle;\n\nop\
en ($F, $ARGV[1]);\nwhile(<$F>)\n  {\n    my $l=$_\
;\n    if ($l=~/^#/){;}\n    elsif (($l=~/\\d+\\s+\
\\d+\\s+(\\S+)\\s+(\\S+)/))\n      {\n	my $n=$1;\n\
	$name{$1}++;\n      }\n  }\nclose ($F);\n\nopen (\
$F, $ARGV[0]);\nwhile(<$F>)\n  {\n    my $l=$_;\n \
   if ($l=~/^#/){;}\n    elsif ($l=~/\\d+\\s+\\d+\\
\s+(\\S+)\\s+(\\S+)/)\n      {\n	my $n=$1;\n	$name\
{$n}++;\n	if ($name{$n}==2){$nseq++;}\n      }\n  \
}\nclose ($F);\n\nif (!$fasta && $nseq>0)\n  {\n  \
  print \"#NAMESEQ_01\\n\";\n    print \"# $nseq\\\
n\";\n  }\nopen ($F, $ARGV[0]);\nwhile(<$F>)\n  {\\
n    my $l=$_;\n    if ($l=~/^#/){;}\n    elsif ($\
l=~/.\\d+\\s+\\d+\\s+(\\S+)\\s+(\\S+)/)\n      {\n\
	my $n=$1;\n	my $s=$2;\n	if ($name{$n}==2)\n	  {\n\
	    if ($fasta)\n	      {\n		print \">$n\\n$s\\n\\
";\n	      }\n	    else\n	      {\n		print \"$l\";\
\n	      }\n	  }\n      }\n  }\nclose ($F);\nexit \
(0);\n\n\n","use Env qw(HOST);\nuse Env qw(HOME);\\
nuse Env qw(USER);\n\n\n$query_start=-1;\n$query_e\
nd=-1;\n\nwhile (<>)\n  {\n    if ( /\\/\\//){$in_\
aln=1;}\n    elsif ( $in_aln && /(\\S+)\\s+(.*)/)\\
n      {\n\n\n	$name=$1;\n	\n\n	$seq=$2;\n	$seq=~s\
/\\s//g;\n        $seq=~s/\\~/\\-/g;\n	$seq=~s/\\.\
/\\-/g;\n	if ( $list{$n}{'name'} && $list{$n}{'nam\
e'} ne $name)\n	  {\n	    print \"$list{$n}{'name'\
} Vs $name\";\n	    \n	    exit (EXIT_FAILURE);\n	\
  }\n	else\n	  {\n	    $list{$n}{'name'}= $name;\n\
	  }\n\n	$list{$n}{'seq'}=$list{$n}{'seq'}.$seq;\n\
	\n	$nseq=++$n;\n	\n      }\n    else\n      {$n=0\
;}\n  }\n\n\nfor ($a=0; $a<$nseq; $a++)\n  {\n    \
print \">$list{$a}{'name'}\\n$list{$a}{'seq'}\\n\"\
;\n  }\n      \n","$run_anyway=2;\nmy $msaf=\"msa.\
in.tmp.$$\";\nmy $msaoutf=\"msa.out.tmp.$$\";\nmy \
$err=\"msa.out.err.$$\";\nopen  (F, $ARGV[0]);\nop\
en  (OUT, \">$msaf\");\n$nseq=0;\nwhile (<F>)\n  {\
\n    $l=$_;\n    if ( $l=~/^>(\\S+)/)\n      {\n	\
$s=$seqname{$nseq++}=$1;\n	print OUT \"$l\";\n	\n \
     }\n    else \n      {\n	$l=uc($l);\n	print OU\
T \"$l\";\n      }\n  }\n\nclose (F);\nclose(OUT);\
\n\nsystem (\"msa $msaf > $msaoutf 2>$err\");\nope\
n (F, \"$msaoutf\");\n$read=0;\n$cn=0;\nwhile (<F>\
)\n  {\n    $l=$_;\n    if ($read)\n      {\n	if (\
$l=~/End gaps not penalized/){$read=0;}\n	elsif (!\
($l=~/\\S/))\n	  {\n	    $cn=0;\n	  }\n	else\n	  {\
\n	    \n	    chomp ($l);\n	    $seqal{$cn++}.=$l;\
\n	    $tot++;\n	  }\n      }\n    elsif ($l=~/Opt\
imal Multiple Alignment/)\n      {\n	$read=1;\n   \
   }\n  }\nclose (F);\n\nif ($tot<1 && $run_anyway\
==1)\n  {\n    print STDERR \"\\nWarning: MSA retu\
rned a NULL file -- Use T-Coffee instead\\n\";\n  \
  open (F,$err);\n    while (<F>){print \"$_\";}\n\
      \n    system (\"t_coffee -seq $msaf -outfile\
 $ARGV[1]  -quiet\");\n  }\nelsif ($tot<1 && $run_\
anyway==2)\n  {\n    \n    \n    $nseq/=2;\n    $n\
seq=int ($nseq);\n    if ($nseq<2){$nseq=2;}\n    \
print \"RUN MSA with NSeq=$nseq\\n\";\n    #print \
(\"t_coffee -dpa -dpa_nseq $nseq -seq $ARGV[0] -dp\
a_tree codnd -outfile $ARGV[1] -dpa_method msa_msa\
\");\n    system (\"t_coffee -dpa -dpa_nseq $nseq \
-seq $ARGV[0] -dpa_tree codnd -outfile $ARGV[1] -d\
pa_method msa_msa>/dev/null\");\n\n  }\nelsif ($to\
t<1)\n  {\n    exit (EXIT_FAILURE);\n  }\nelse\n  \
{\n    open (OUT, \">$ARGV[1]\");\n    for ($a=0; \
$a<$nseq;$a++)\n      {\n	print OUT \">$seqname{$a\
}\\n$seqal{$a}\\n\";\n      }\n    close (OUT);\n \
 }\n\n\n\nunlink ($msaf);\nunlink ($msaoutf);\nunl\
ink ($err);\n","use strict;\nuse Cwd;\nuse File::B\
asename;\nmy $test=0;\n\nmy $tmpdir=\"/tmp/tco/ali\
gners/upp/\";\nmymkdir ($tmpdir);\n\n\n\nif ($ARGV\
[0] eq \"one\")\n  {\n    seq2msa ($ARGV[1], $ARGV\
[2]);\n  }\nelsif ($ARGV[0] eq \"all\")\n  {\n    \
listseq2listmsa ($ARGV[1]);\n  }\n\nsub listseq2li\
stmsa\n  {\n    my $list=shift;\n    my $cdir = ge\
tcwd;\n    my $dir=random_string(10);\n    $dir=\"\
$tmpdir/$dir/\";\n    my %h;\n    my $n;\n    mkdi\
r  ($dir);\n\n    open (F, $list);\n    while (<F>\
)\n      {\n        my $l=$_;\n\n        chomp($l)\
;\n        my @f=split (/\\s+/, $l);\n	if ( -e $f[\
0])\n          {\n            $h{$n}{in}=$f[0];\n \
           ($h{$n}{name},$h{$n}{path})=fileparse (\
$f[0]);\n            $h{$n}{NFin}= \"$dir/$h{$n}{n\
ame}.seq\";\n	    \n            $h{$n}{NFout}=\"$d\
ir/$h{$n}{name}.aln\";\n\n            $h{$n}{out}=\
$f[1];\n\n            fasta2fastaupp ($h{$n}{in}, \
$h{$n}{NFin});\n            $n++;\n          }\n  \
    }\n    close (F);\n    chdir ($dir);\n    \n  \
  if (!$test)\n      {\n	system (\"fbname=\\$(base\
name `ls *.seq` .seq); \\\n             run_upp.py\
 -s \\${fbname}.seq -m amino --cpu 1 -d outdir -o \
\\${fbname}.aln; \\\n             mv outdir/\\${fb\
name}.aln_alignment.fasta \\${fbname}.aln;\");\n  \
    }\n    \n    foreach my $n (keys (%h))\n      \
{\n	if ($test)\n	  {\n	    system (\"cp $h{$n}{NFi\
n} $h{$n}{NFout}\");\n	    print \"$h{$n}{NFin} $h\
{$n}{NFout} $h{$n}{out}\\n\";\n	  }\n        fasta\
upp2fasta ($h{$n}{NFout},$h{$n}{out});\n      }\n \
   chdir ($cdir);\n  }\n\nsub seq2msa\n    {\n    \
  my ($in, $out)=@_;\n      my $cdir=getcwd;\n    \
  \n      \n      if (!($in=~/\\//)){$in=$cdir.\"/\
\".$in;}\n      if (!($out=~/\\//)){$out=$cdir.\"/\
\".$out;}\n      \n      my $file=random_string(10\
);\n      $file=\"$tmpdir/$file\";\n      open (F,\
 \">$file\");\n      print F \"$in $out\\n\";\n   \
   close (F);\n      \n      return listseq2listms\
a ($file);\n    }\n	\nsub fasta2fastaupp\n  {\n   \
 my ($in, $out)=@_;\n    my ($name, $seq, $n);\n  \
  \n    if (!-e $in){return;}\n    \n    open (IN,\
 \"$in\");\n    open (OUT, \">$out\");\n    local \
$/ = \"\\n>\";  # read by FASTA record\n    \n    \
while (<IN>)\n      {\n	my $l=$_;\n	$l=~s/>//g;\n	\
$l=\">\".$l;\n	\n	$l=~/^>(.*)/;\n	$name=$1;\n	\n	$\
l=~s/^>*.+\\n//;\n	$l=~s/\\n//g;\n	$seq=$l;\n	\n	$\
seq=~s/u/x/g;\n	$seq=~s/U/X/g;\n	print OUT \">$nam\
e\\n$seq\\n\";\n	$n++;\n      }\n    if ($n==2)\n \
     {\n	print OUT \">fake_seq4upp\\n$seq\\n\";\n \
     }\n    close (IN);\n    close (OUT);\n    loc\
al $/=\"\\n\";\n  }\n\nsub fastaupp2fasta\n  {\n  \
  my ($in, $out)=@_;\n    my ($name, $seq, $n);\n \
   \n    if (!-e $in){return;}\n    \n    open (IN\
, \"$in\");\n    open (OUT, \">$out\");\n    local\
 $/ = \"\\n>\";  # read by FASTA record\n    \n   \
 while (<IN>)\n      {\n	my $l=$_;\n	$l=~s/>//g;\n\
	$l=\">\".$l;\n	\n	$l=~/^>(.*)/;\n	$name=$1;\n	\n	\
$l=~s/^>*.+\\n//;\n	$l=~s/\\n//g;\n	$seq=$l;\n	\n	\
$seq=~s/x/u/g;\n	$seq=~s/X/U/g;\n	\n	if (!($name=~\
/fake_seq4upp/))\n	  {\n	    print OUT \">$name\\n\
$seq\\n\";\n	  }\n      }\n    close (IN);\n    cl\
ose (OUT);\n    local $/=\"\\n\";\n  }\n\nsub rand\
om_string\n    {\n      my $len=shift;\n      my @\
chars = (\"A\"..\"Z\", \"a\"..\"z\");\n      my $s\
tring;\n      $string .= $chars[rand @chars] for 1\
..$len;\n      return $string;\n    }\n\nsub mymkd\
ir\n      {\n	my $d=shift;\n	my $cd='/';\n	\n	fore\
ach my $e (split (/\\//, $d))\n	  {\n	    $cd.=\"$\
e/\";\n	    if ( !-d $cd){mkdir ($cd);}\n	  }\n	re\
turn;\n      }\n      \n			  \n      \n","use stri\
ct;\nuse Cwd;\nuse File::Basename;\n\n\nmy $tmpdir\
=\"/tmp/tco/aligners/clustalo/\";\nmymkdir ($tmpdi\
r);\n\n\n\nif ($ARGV[0] eq \"one\")\n  {\n    seq2\
msa ($ARGV[1], $ARGV[2]);\n  }\nelsif ($ARGV[0] eq\
 \"all\")\n  {\n    listseq2listmsa ($ARGV[1]);\n \
 }\n\n\n\nsub listseq2listmsa\n  {\n    my $list=s\
hift;\n    my $cdir = getcwd;\n    my $dir=random_\
string(10);\n    $dir=\"$tmpdir/$dir/\";\n    my %\
h;\n    my $n;\n    mkdir  ($dir);\n    \n    open\
 (F, $list);\n    while (<F>)\n      {\n	my $l=$_;\
\n\n	chomp($l);\n	my @f=split (/\\s+/, $l);\n	#pri\
nt \"$l: 0:$f[0], 1:$f[1]\\n\";\n	if ( -e $f[0])\n\
	  {\n	    $h{$n}{in}=$f[0];\n	    ($h{$n}{name},$\
h{$n}{path})=fileparse ($f[0]);\n	    $h{$n}{NFin}\
= \"$dir/$h{$n}{name}.seq4nf\";\n	    $h{$n}{NFout\
}=\"$dir/$h{$n}{name}.aln\";\n	    \n	    $h{$n}{o\
ut}=$f[1];\n	    \n	    translate_fasta_seq (\"uU\\
", \"X\",$h{$n}{in}, $h{$n}{NFin});\n	    $n++;\n	\
  }\n      }\n    close (F);\n    \n    \n    chdi\
r ($dir);\n    dump_nf (\"nf\");\n    dump_config \
();\n   \n    #system (\"nextflow run nf  --name \\
\'*.seq4nf\\' >/dev/null 2>/dev/null\");\n    syst\
em (\"nextflow run nf  --name \\'*.seq4nf\\'\");\n\
    foreach my $n (keys (%h))\n      {\n	translate\
_fasta_seq (\"uU\", \"X\",$h{$n}{NFout},$h{$n}{out\
});\n      }\n    chdir ($cdir);\n  }\nsub seq2msa\
\n    {\n      my ($in, $out)=@_;\n      my $cdir=\
getcwd;\n      \n      \n      if (!($in=~/\\//)){\
$in=$cdir.\"/\".$in;}\n      if (!($out=~/\\//)){$\
out=$cdir.\"/\".$out;}\n      \n      my $file=ran\
dom_string(10);\n      $file=\"$tmpdir/$file\";\n \
     open (F, \">$file\");\n      print F \"$in $o\
ut\\n\";\n      close (F);\n      \n      return l\
istseq2listmsa ($file);\n    }\n	\nsub seq2msa_old\
\n  {\n    my ($in, $out)=@_;\n    my $cdir = getc\
wd;\n    my $dir=random_string(10);\n    $dir=\"/t\
mp/upp.nf4tcoffee/$dir\";\n    my $seq=random_stri\
ng(10);\n    $seq.=\".fa\";\n    my $aln=$seq;\n  \
  $aln.=\".aln\";\n    \n    mkdir ($dir);\n    tr\
anslate_fasta_seq (\"uU\", \"X\",$in, \"$dir/$seq\\
");\n    chdir ($dir);\n    \n    dump_nf (\"nf\")\
;\n    dump_config ();\n    print \"IN: $in OUT: $\
cdir/$out\\nDIR: $dir\\nnextflow run nf  --name \\\
'*.fa\\' \\n\";\n    system (\"nextflow run nf  --\
name \\'*.fa\\' \");\n    print \"$dir/$aln $cdir/\
$out\\n\";\n    translate_fasta_seq (\"xX\", \"U\"\
,$aln, \"$cdir/$out\");\n    chdir ($cdir);\n   } \
\nsub translate_fasta_seq\n  {\n    my ($from, $to\
, $in, $out)=@_;\n    my $n;\n    my $skip;\n    m\
y $l;\n    my $cseq;\n    if (!-e $in){return;}\n \
   \n    open (IN, \"$in\");\n    open (OUT, \">$o\
ut\");\n   \n    while (<IN>)\n      {\n	$l=$_;\n	\
if ($l=~\">\"){$n++;$cseq=\"\";}\n	else { $l=~s/[$\
from]/$to/;$cseq.=$l;}\n\n	if ($skip){$skip=0;}\n	\
elsif ($l=~/>fake_seq/){$skip=1;}\n	else\n	  {\n	 \
   print OUT \"$l\";\n	  }\n      }\n    if ($n==2\
 && $from eq \"uU\")\n      {\n	print OUT \">fake_\
seq\\n$cseq\\n\";\n      }\n    close (IN);\n    c\
lose (OUT);\n  }\n\nsub dump_config\n    {\n      \
open (F, \">nextflow.config\");\n\n      print F \\
"docker.enabled = true\\n\";\n      print F \"proc\
ess.container = \\'cbcrg/benchfam_large_scale\\'\\\
n\";\n      close (F);\n    }\n\nsub dump_nf\n  {\\
n    my $nff=shift;\n    open (F,\">$nff\");\n    \
print F \"#!/usr/bin/env nextflow\\n\";\n    print\
 F \"params.base_dir=\\\"./\\\"\\n\";\n    print F\
 \"params.out_dir=\\\"./\\\"\\n\";\n    print F \"\
Channel.fromPath(params.name)\\n\";\n    print F \\
"\\t.map{ tuple(it.baseName, it) }\\n\";\n    \n  \
  print F \"\\t.set{ file_names_1 }\\n\";\n    pri\
nt F \"process clustalo_align{\\n\";\n    print F \
\"\\tpublishDir params.out_dir, mode: \\\"copy\\\"\
\\n\";\n    print F \"tag \\\"\\${name}\\\"\";\n  \
  print F \"\\n\";\n    print F \"\\tinput:\\n\";\\
n    print F \"\\tset name, file(seq_file) from fi\
le_names_1\\n\";\n    print F \"\\toutput:\\n\";\n\
    print F \"\\tfile \\\"\\${name}.aln\\\"\\n\";\\
n    print F \"\\n\";\n    print F \" \\\"\\\"\\\"\
\\n\";\n    print F \" clustalo -i \\$seq_file -o \
\\${name}.aln\\n\";\n    print F \"\\\"\\\"\\\"\\n\
\\n\";\n    print F \"}\\n\";\n    close (F);\n  }\
\n\nsub random_string\n    {\n      my $len=shift;\
\n      my @chars = (\"A\"..\"Z\", \"a\"..\"z\");\\
n      my $string;\n      $string .= $chars[rand @\
chars] for 1..$len;\n      return $string;\n    }\\
n\nsub mymkdir\n      {\n	my $d=shift;\n	my $cd='/\
';\n	\n	foreach my $e (split (/\\//, $d))\n	  {\n	\
    $cd.=\"$e/\";\n	    if ( !-d $cd){mkdir ($cd);\
}\n	  }\n	return;\n      }\n      \n			  \n      \\
n","\nmy $msaf=\"msa.in.tmp.$$\";\nmy $msaoutf=\"m\
sa.out.tmp.$$\";\nmy $cost=\"blosum62.tmp.$$\";\n\\
nopen  (F, $ARGV[0]);\nopen  (OUT, \">$msaf\");\n$\
nseq=0;\nwhile (<F>)\n  {\n    $l=$_;\n    if ( $l\
=~/^>(\\S+)/)\n      {\n	my $simple=\"Seq$nseq\";\\
n	$s=$seqname{$nseq++}=$1;\n	$translate{$simple}=$\
s;\n	\n	print OUT \">$simple\\n\";\n	\n      }\n  \
  else\n      {\n	$l=uc($l);\n	print OUT \"$l\";\n\
      }\n  }\nclose (F);\nclose(OUT);\n\ndump_blos\
um ($cost);\nsystem (\"dca -c $cost -q $msaf> $msa\
outf 2>/dev/null\");\nopen (F, \"$msaoutf\");\nope\
n (OUT, \">$ARGV[1]\");\n\n$read=0;\nwhile (<F>)\n\
  {\n    $l=$_;\n    if ($l=~/^>(\\S+)/)\n      {\\
n	$read=1;\n	$name=$translate{$1};\n	print OUT \">\
$name\\n\";\n      }\n    elsif ($read && ($l=~/\\\
S/))\n      {\n	print OUT \"$l\";\n      }\n    el\
se\n      {\n	$read=0;\n      }\n  }\nclose (F);\n\
\nunlink ($cost);\nunlink ($msaf);\nunlink ($msaou\
tf);\n\nsub dump_blosum\n  {\n    my $f=shift;\n  \
  open (F, \">$f\");\n\n    print F \"6\\n\";\n   \
 print F \"- -   0\\n\";\n    print F \"W W   0\\n\
\";\n    print F \"Y Y   4\\n\";\n    print F \"F \
F   5\\n\";\n    print F \"V V   7\\n\";\n    prin\
t F \"L L   7\\n\";\n    print F \"I I   7\\n\";\n\
    print F \"M M   6\\n\";\n    print F \"K K   6\
\\n\";\nprint F \"R R   6\\n\";\n    print F \"H H\
   3\\n\";\n    print F \"Q Q   6\\n\";\n    print\
 F \"E E   6\\n\";\n    print F \"D D   5\\n\";\n \
   print F \"N N   5\\n\";\n    print F \"G G   5\\
\n\";\n    print F \"A A   7\\n\";\n    print F \"\
P P   4\\n\";\n    print F \"T T   6\\n\";\n    pr\
int F \"S S   7\\n\";\n    print F \"C C   2\\n\";\
\n    print F \"- C  10 \\n\";\n    print F \"- S \
 10\\n\";\n    print F \"- T  10 \\n\";\n    print\
 F \"- P  10\\n\";\n    print F \"- A  10 \\n\";\n\
    print F \"- G  10\\n\";\n    print F \"- N  10\
 \\n\";\n    print F \"- D  10\\n\";\n    print F \
\"- E  10 \\n\";\n    print F \"- Q  10\\n\";\npri\
nt F \"- H  10 \\n\";\n    print F \"- R  10\\n\";\
\n    print F \"- K  10 \\n\";\n    print F \"- M \
 10\\n\";\n    print F \"- I  10 \\n\";\n    print\
 F \"- L  10\\n\";\n    print F \"- V  10 \\n\";\n\
    print F \"- F  10\\n\";\n    print F \"- Y  10\
 \\n\";\n    print F \"- W  10\\n\";\n    print F \
\"W C  13 \\n\";\n    print F \"W S  14\\n\";\n   \
 print F \"W T  13 \\n\";\n    print F \"W P  15\\\
n\";\n    print F \"W A  14 \\n\";\n    print F \"\
W G  13\\n\";\n    print F \"W N  15 \\n\";\n    p\
rint F \"W D  15\\n\";\n    print F \"W E  14 \\n\\
";\n    print F \"W Q  13\\n\";\n    print F \"W H\
  13 \\n\";\n    print F \"W R  14\\n\";\n    prin\
t F \"W K  14 \\n\";\n    print F \"W M  12\\n\";\\
n    print F \"W I  14 \\n\";\n    print F \"W L  \
13\\n\";\n    print F \"W V  14 \\n\";\n    print \
F \"W F  10\\n\";\n    print F \"W Y   9 \\n\";\n \
   print F \"Y C  13\\n\";\n    print F \"Y S  13 \
\\n\";\n    print F \"Y T  13\\n\";\n    print F \\
"Y P  14 \\n\";\n    print F \"Y A  13\\n\";\n    \
print F \"Y G  14 \\n\";\n    print F \"Y N  13\\n\
\";\n    print F \"Y D  14 \\n\";\n    print F \"Y\
 E  13\\n\";\n    print F \"Y Q  12 \\n\";\n    pr\
int F \"Y H   9\\n\";\n    print F \"Y R  13 \\n\"\
;\n    print F \"Y K  13\\n\";\n    print F \"Y M \
 12 \\n\";\n    print F \"Y I  12\\n\";\n    print\
 F \"Y L  12 \\n\";\n    print F \"Y V  12\\n\";\n\
    print F \"Y F   8 \\n\";\n    print F \"F C  1\
3\\n\";\nprint F \"F S  13 \\n\";\n    print F \"F\
 T  13\\n\";\n    print F \"F P  15 \\n\";\n    pr\
int F \"F A  13\\n\";\n    print F \"F G  14 \\n\"\
;\n    print F \"F N  14\\n\";\n    print F \"F D \
 14 \\n\";\n    print F \"F E  14\\n\";\n    print\
 F \"F Q  14 \\n\";\n    print F \"F H  12\\n\";\n\
    print F \"F R  14 \\n\";\n    print F \"F K  1\
4\\n\";\n    print F \"F M  11 \\n\";\n    print F\
 \"F I  11\\n\";\n    print F \"F L  11 \\n\";\n  \
  print F \"F V  12\\n\";\n    print F \"V C  12 \\
\n\";\n    print F \"V S  13\\n\";\n    print F \"\
V T  11 \\n\";\n    print F \"V P  13\\n\";\n    p\
rint F \"V A  11 \\n\";\n    print F \"V G  14\\n\\
";\n    print F \"V N  14 \\n\";\n    print F \"V \
D  14\\n\";\nprint F \"V E  13 \\n\";\nprint F \"V\
 Q  13\\n\";\nprint F \"V H  14 \\n\";\nprint F \"\
V R  14\\n\";\nprint F \"V K  13 \\n\";\nprint F \\
"V M  10\\n\";\nprint F \"V I   8 \\n\";\nprint F \
\"V L  10\\n\";\nprint F \"L C  12 \\n\";\nprint F\
 \"L S  13\\n\";\nprint F \"L T  12 \\n\";\nprint \
F \"L P  14\\n\";\nprint F \"L A  12 \\n\";\nprint\
 F \"L G  15\\n\";\nprint F \"L N  14 \\n\";\nprin\
t F \"L D  15\\n\";\nprint F \"L E  14 \\n\";\npri\
nt F \"L Q  13\\n\";\nprint F \"L H  14 \\n\";\npr\
int F \"L R  13\\n\";\nprint F \"L K  13 \\n\";\np\
rint F \"L M   9\\n\";\nprint F \"L I   9 \\n\";\n\
print F \"I C  12\\n\";\nprint F \"I S  13 \\n\";\\
nprint F \"I T  12\\n\";\nprint F \"I P  14 \\n\";\
\nprint F \"I A  12\\n\";\nprint F \"I G  15 \\n\"\
;\nprint F \"I N  14\\n\";\nprint F \"I D  14 \\n\\
";\nprint F \"I E  14\\n\";\nprint F \"I Q  14 \\n\
\";\nprint F \"I H  14\\n\";\nprint F \"I R  14 \\\
n\";\nprint F \"I K  14\\n\";\nprint F \"I M  10 \\
\n\";\nprint F \"M C  12\\n\";\nprint F \"M S  12 \
\\n\";\nprint F \"M T  12\\n\";\nprint F \"M P  13\
 \\n\";\nprint F \"M A  12\\n\";\nprint F \"M G  1\
4 \\n\";\nprint F \"M N  13\\n\";\nprint F \"M D  \
14 \\n\";\nprint F \"M E  13\\n\";\nprint F \"M Q \
 11 \\n\";\nprint F \"M H  13\\n\";\nprint F \"M R\
  12 \\n\";\nprint F \"M K  12\\n\";\nprint F \"K \
C  14 \\n\";\nprint F \"K S  11\\n\";\nprint F \"K\
 T  12 \\n\";\nprint F \"K P  12\\n\";\nprint F \"\
K A  12 \\n\";\nprint F \"K G  13\\n\";\nprint F \\
"K N  11 \\n\";\nprint F \"K D  12\\n\";\nprint F \
\"K E  10 \\n\";\nprint F \"K Q  10\\n\";\nprint F\
 \"K H  12 \\n\";\nprint F \"K R   9\\n\";\nprint \
F \"R C  14 \\n\";\nprint F \"R S  12\\n\";\nprint\
 F \"R T  12 \\n\";\nprint F \"R P  13\\n\";\nprin\
t F \"R A  12 \\n\";\nprint F \"R G  13\\n\";\npri\
nt F \"R N  11 \\n\";\nprint F \"R D  13\\n\";\npr\
int F \"R E  11 \\n\";\nprint F \"R Q  10\\n\";\np\
rint F \"R H  11 \\n\";\nprint F \"H C  14\\n\";\n\
print F \"H S  12 \\n\";\nprint F \"H T  13\\n\";\\
nprint F \"H P  13 \\n\";\nprint F \"H A  13\\n\";\
\nprint F \"H G  13 \\n\";\nprint F \"H N  10\\n\"\
;\nprint F \"H D  12 \\n\";\nprint F \"H E  11\\n\\
";\nprint F \"H Q  11 \\n\";\nprint F \"Q C  14\\n\
\";\nprint F \"Q S  11 \\n\";\nprint F \"Q T  12\\\
n\";\nprint F \"Q P  12 \\n\";\nprint F \"Q A  12\\
\n\";\nprint F \"Q G  13 \\n\";\nprint F \"Q N  11\
\\n\";\nprint F \"Q D  11 \\n\";\nprint F \"Q E   \
9\\n\";\nprint F \"E C  15 \\n\";\nprint F \"E S  \
11\\n\";\nprint F \"E T  12 \\n\";\nprint F \"E P \
 12\\n\";\nprint F \"E A  12 \\n\";\nprint F \"E G\
  13\\n\";\nprint F \"E N  11 \\n\";\nprint F \"E \
D   9\\n\";\nprint F \"D C  14 \\n\";\nprint F \"D\
 S  11\\n\";\nprint F \"D T  12 \\n\";\nprint F \"\
D P  12\\n\";\nprint F \"D A  13 \\n\";\nprint F \\
"D G  12\\n\";\nprint F \"D N  10 \\n\";\nprint F \
\"N C  14\\n\";\nprint F \"N S  10 \\n\";\nprint F\
 \"N T  11\\n\";\nprint F \"N P  13 \\n\";\nprint \
F \"N A  13\\n\";\nprint F \"N G  11 \\n\";\nprint\
 F \"G C  14\\n\";\nprint F \"G S  11 \\n\";\nprin\
t F \"G T  13\\n\";\nprint F \"G P  13 \\n\";\npri\
nt F \"G A  11\\n\";\nprint F \"A C  11 \\n\";\npr\
int F \"A S  10\\n\";\nprint F \"A T  11 \\n\";\np\
rint F \"A P  12\\n\";\nprint F \"P C  14 \\n\";\n\
print F \"P S  12\\n\";\nprint F \"P T  12 \\n\";\\
nprint F \"T C  12\\n\";\nprint F \"T S  10 \\n\";\
\nprint F \"S C  12\\n\";\nclose (F);\n    return;\
\n  }\n    \n","\nuse Env qw(HOST);\nuse Env qw(HO\
ME);\nuse Env qw(USER);\n\n                       \
                                 \nuse strict;    \
                                         \nuse war\
nings;\nuse diagnostics;\n\nmy $in_hit_list, my $i\
n_aln=0, my(%name_list)=(),my (%list)=(),my $n_seq\
=0; my $test=0;\nmy($j)=0, my $n=0, my $nom, my $l\
g_query, my %vu=();\n\nopen (F, \">tmp\");\n\n$/=\\
"\\n\";\nwhile (<>)\n{\n    print F $_;\n    if($_\
 =~ /Query=\\s*(.+?)\\s/i) { $nom=$1;}\n\n    if (\
 /Sequences producing significant alignments/){$in\
_hit_list=1;}\n    \n    if ($_=~ /^pdb\\|/i) { $_\
=~ s/pdb\\|//g; }\n    if ($_=~ /^(1_\\d+)\\s+\\d+\
/) { $_=~ s/$1/QUERY/;}\n      \n    if ( /^(\\S+)\
.+?\\s+[\\d.]+\\s+([\\de.-]+)\\s+$/ && $in_hit_lis\
t)	\n    {\n	my($id)=$1; # \n	$id=~ s/\\|/_/g; #\n\
	if ($id =~ /.+_$/) { chop($id) }; #\n	$name_list{\
$n_seq++}=$id;\n	$name_list{$n_seq-1}=~ s/.*\\|//g\
;     \n    }\n  \n    if (/query/i) {$in_aln=1;}\\
n    if ( /^(\\S+)\\s+(\\d+)\\s+([a-zA-Z-]+)\\s+(\\
\d+)/ || /^(\\S+)(\\s+)(\\-+)(\\s+)/ && ($in_aln =\
= 1))\n    {\n	my $name=$1;\n	my $start=$2;\n	my $\
seq=$3;\n	my $end=$4;\n		\n	if ($name =~ /QUERY/i)\
 { $lg_query=length($seq); }\n\n	unless ($test > $\
n) #m\n	{\n	    my(@seqq)= split('',$seq);\n	    m\
y($gap_missing)= scalar(@seqq);\n	    \n	    while\
 ($gap_missing != $lg_query)  { unshift (@seqq,\"-\
\"); $gap_missing= scalar(@seqq); }\n	    $seq=joi\
n('',@seqq);  #m\n	}\n	\n	if ($name =~ /QUERY/i)\n\
	{\n	    $n=0; %vu=(); $j=0;\n	    $list{$n}{'real\
_name'}=\"$nom\";\n	}	\n	else\n	{\n	    unless (ex\
ists $vu{$name}) { ++$j;}	\n	    $list{$n}{'real_n\
ame'}=$name_list{$j-1};\n	}\n		\n	$list{$n}{'name'\
}=$name;\n\n	$seq=~tr/a-z/A-Z/;\n	$list{$n}{'seq'}\
=$list{$n}{'seq'};\n	$list{$n}{'seq'}.=$seq;\n\n	$\
n++;\n	$vu{$name}++;\n	$test++;\n   } \n    \n}\n\\
nmy @numero=();\n\nfor (my $a=0; $a<$n; $a++) #m\n\
{\n    my $long=length($list{0}{'seq'});  \n    my\
 $long1= length($list{$a}{'seq'});\n  \n    while \
($long1 ne $long)\n    {\n	$list{$a}{'seq'}.=\"-\"\
;\n	$long1= length ($list{$a}{'seq'});\n    } \n \\
n    push (@numero,\"$list{$a}{'name'} $list{$a}{'\
real_name'}\\n\");\n}\n\nmy %dejavu=();\n\n\nfor (\
my $i=0; $i<=$#numero; $i++)\n{\n    my $s=\">$lis\
t{$i}{'real_name'}\\n$list{$i}{'seq'}\\n\";\n    m\
y $k=0;\n    \n    if (exists $dejavu{$numero[$i]}\
) {next;}\n    else\n    {	\n	for ($j=0; $j<$n ; $\
j++)\n	{\n	    if (\"$numero[$i]\" eq \"$numero[$j\
]\" && $j != $i )\n	    {\n		++$k;\n		$s .=\">$lis\
t{$j}{'real_name'}\\n$list{$j}{'seq'}\\n\";\n	    \
}\n	}	\n    }\n    \n    if ($k>0) \n    {\n	my $c\
ons;\n	open (SOR,\">tempo_aln2cons\"); print SOR $\
s;  close SOR ;\n	open (COM,\"t_coffee -other_pg s\
eq_reformat -in tempo_aln2cons -action +aln2cons +\
upper |\") ; \n     	while (<COM>)\n	{	\n	    if (\
/^>/) { $cons =\">$list{$i}{'real_name'}\\n\"; nex\
t;}\n	    $_=~ s/\\n//g;\n	    $cons .=$_;\n	}\n	c\
lose COM; unlink (\"tempo_aln2cons\");\n	print $co\
ns,\"\\n\"; print F $cons,\"\\n\";\n    }	\n    el\
se  { print $s;  print F $s; }\n    \n    $dejavu{\
$numero[$i]}++;\n} #m\n\nexit;\n\n\n\n\n\n\n\n\n\n\
\n\n","use Env;\n\n\n$tmp_dir=\"\";\n$init_dir=\"\\
";\n$program=\"tc_generic_method.pl\";\n\n$blast=@\
ARGV[0];\n\n$name=\"query\";$seq=\"\";\n%p=blast_x\
ml2profile($name,$seq,100, 0, 0, $blast);\n&output\
_profile (%p);\n\n\nsub output_profile\n  {\n    m\
y (%profile)=(@_);\n    my ($a);\n    for ($a=0; $\
a<$profile{n}; $a++)\n      {\n	\n	print \">$profi\
le{$a}{name} $profile{$a}{comment}\\n$profile{$a}{\
seq}\\n\";\n      }\n    return;\n  }\nsub file_co\
ntains \n  {\n    my ($file, $tag, $max)=(@_);\n  \
  my ($n);\n    $n=0;\n    \n    if ( !-e $file &&\
 ($file =~/$tag/)) {return 1;}\n    elsif ( !-e $f\
ile){return 0;}\n    else \n      {\n	open (FC, \"\
$file\");\n	while ( <FC>)\n	  {\n	    if ( ($_=~/$\
tag/))\n	      {\n		close (FC);\n		return 1;\n	   \
   }\n	    elsif ($max && $n>$max)\n	      {\n		cl\
ose (FC);\n		return 0;\n	      }\n	    $n++;\n	  }\
\n      }\n    close (FC);\n    return 0;\n  }\n	 \
   \n	  \nsub file2string\n  {\n    my $f=@_[0];\n\
    my $string, $l;\n    open (F,\"$f\");\n    whi\
le (<F>)\n      {\n\n	$l=$_;\n	#chomp ($l);\n	$str\
ing.=$l;\n      }\n    close (F);\n    $string=~s/\
\\r\\n//g;\n    $string=~s/\\n//g;\n    return $st\
ring;\n  }\n\n\n\nsub tag2value \n  {\n    \n    m\
y $tag=(@_[0]);\n    my $word=(@_[1]);\n    my $re\
turn;\n    \n    $tag=~/$word=\"([^\"]+)\"/;\n    \
$return=$1;\n    return $return;\n  }\n      \nsub\
 hit_tag2pdbid\n  {\n    my $tag=(@_[0]);\n    my \
$pdbid;\n       \n    $tag=~/id=\"(\\S+)\"/;\n    \
$pdbid=$1;\n    $pdbid=~s/_//;\n    return $pdbid;\
\n  }\nsub id2pdbid \n  {\n    my $id=@_[0];\n  \n\
    if ($id =~/pdb/)\n      {\n	$id=~/pdb(.*)/;\n	\
$id=$1;\n      }\n    $id=~s/[|�_]//g;\n    return\
 $id;\n  }\nsub set_blast_type \n  {\n    my $file\
 =@_[0];\n    if (&file_contains ($file,\"EBIAppli\
cationResult\",100)){$BLAST_TYPE=\"EBI\";}\n    el\
sif (&file_contains ($file,\"NCBI_BlastOutput\",10\
0)) {$BLAST_TYPE=\"NCBI\";}\n    else\n      {\n	$\
BLAST_TYPE=\"\";\n      }\n    return $BLAST_TYPE;\
\n  }\nsub blast_xml2profile \n  {\n    my ($name,\
$seq,$maxid, $minid, $mincov, $file)=(@_);\n    my\
 (%p, $a, $string, $n);\n    \n\n\n    if ($BLAST_\
TYPE eq \"EBI\" || &file_contains ($file,\"EBIAppl\
icationResult\",100)){%p=ebi_blast_xml2profile(@_)\
;}\n    elsif ($BLAST_TYPE eq \"NCBI\" || &file_co\
ntains ($file,\"NCBI_BlastOutput\",100)){%p=ncbi_b\
last_xml2profile(@_);}\n    else \n      {\n	print\
 \"************ ERROR: Blast Returned an unknown X\
ML Format **********************\";\n	die;\n      \
}\n    for ($a=0; $a<$p{n}; $a++)\n      {\n	my $n\
ame=$p{$a}{name};\n	$p{$name}{seq}=$p{$a}{seq};\n \
     }\n    return %p;\n  }\nsub ncbi_blast_xml2pr\
ofile \n  {\n    my ($name,$seq,$maxid, $minid, $m\
incov, $string)=(@_);\n    my ($L,$l, $a,$b,$c,$d,\
$nhits,@identifyerL);\n    \n    \n    $seq=~s/[^a\
-zA-Z]//g;\n    $L=length ($seq);\n    \n    %hit=\
&xml2tag_list ($string, \"Hit\");\n    \n    \n   \
 for ($nhits=0,$a=0; $a<$hit{n}; $a++)\n      {\n	\
my ($ldb,$id, $identity, $expectation, $start, $en\
d, $coverage, $r);\n	my (%ID,%DE,%HSP);\n	\n	$ldb=\
\"\";\n\n	%ID=&xml2tag_list ($hit{$a}{body}, \"Hit\
_id\");\n	$identifyer=$ID{0}{body};\n	\n	%DE=&xml2\
tag_list ($hit{$a}{body}, \"Hit_def\");\n	$definit\
ion=$DE{0}{body};\n	\n	%HSP=&xml2tag_list ($hit{$a\
}{body}, \"Hsp\");\n	for ($b=0; $b<$HSP{n}; $b++)\\
n	  {\n	    my (%START,%END,%E,%I,%Q,%M);\n\n	 \n	\
    %START=&xml2tag_list ($HSP{$b}{body}, \"Hsp_qu\
ery-from\");\n	    %HSTART=&xml2tag_list ($HSP{$b}\
{body}, \"Hsp_hit-from\");\n	    \n	    %LEN=  &xm\
l2tag_list ($HSP{$b}{body}, \"Hsp_align-len\");\n	\
    %END=  &xml2tag_list ($HSP{$b}{body}, \"Hsp_qu\
ery-to\");\n	    %HEND=  &xml2tag_list ($HSP{$b}{b\
ody}, \"Hsp_hit-to\");\n	    %E=&xml2tag_list     \
($HSP{$b}{body}, \"Hsp_evalue\");\n	    %I=&xml2ta\
g_list     ($HSP{$b}{body}, \"Hsp_identity\");\n	 \
   %Q=&xml2tag_list     ($HSP{$b}{body}, \"Hsp_qse\
q\");\n	    %M=&xml2tag_list     ($HSP{$b}{body}, \
\"Hsp_hseq\");\n	    \n	    for ($e=0; $e<$Q{n}; $\
e++)\n\n	      {\n		$qs=$Q{$e}{body};\n		$ms=$M{$e\
}{body};\n		if ($seq eq\"\"){$seq=$qs;$L=length($s\
eq);}\n		\n		$expectation=$E{$e}{body};\n		$identi\
ty=($LEN{$e}{body}==0)?0:$I{$e}{body}/$LEN{$e}{bod\
y}*100;\n		$start=$START{$e}{body};\n		$end=$END{$\
e}{body};\n		$Hstart=$HSTART{$e}{body};\n		$Hend=$\
HEND{$e}{body};\n	\n		$coverage=(($end-$start)*100\
)/$L;\n\n	\n		if ($identity>$maxid || $identity<$m\
inid || $coverage<$mincov){next;}\n		@lr1=(split (\
//,$qs));\n		@lr2=(split (//,$ms));\n		$l=$#lr1+1;\
\n		for ($c=0;$c<$L;$c++){$p[$nhits][$c]=\"-\";}\n\
		for ($d=0,$c=0; $c<$l; $c++)\n		  {\n		    $r=$l\
r1[$c];\n		    if ( $r=~/[A-Za-z]/)\n		      {\n		\
	\n			$p[$nhits][$d + $start-1]=$lr2[$c];\n			$d++\
;\n		      }\n		  }\n		$Qseq[$nhits]=$qs;\n		$Hseq\
[$nhits]=$ms;\n		$QstartL[$nhits]=$start;\n		$Hsta\
rtL[$nhits]=$Hstart;\n		$identityL[$nhits]=$identi\
ty;\n		$endL[$nhits]=$end;\n		$definitionL[$nhits]\
=$definition;\n		$identifyerL[$nhits]=$identifyer;\
\n		$comment[$nhits]=\"$ldb|$identifyer [Eval=$exp\
ectation][id=$identity%][start=$Hstart end=$Hend]\\
";\n		$nhits++;\n	      }\n	  }\n      }\n    \n  \
  $profile{n}=0;\n    $profile{$profile{n}}{name}=\
$name;\n    $profile{$profile{n}}{seq}=$seq;\n    \
$profile {n}++;\n    \n    for ($a=0; $a<$nhits; $\
a++)\n      {\n	$n=$a+1;\n	\n	$profile{$n}{name}=\\
"$name\\_$a\";\n	$profile{$n}{seq}=\"\";\n	$profil\
e{$n}{Qseq}=$Qseq[$a];\n	$profile{$n}{Hseq}=$Hseq[\
$a];\n	$profile{$n}{Qstart}=$QstartL[$a];\n	$profi\
le{$n}{Hstart}=$HstartL[$a];\n	$profile{$n}{identi\
ty}=$identityL[$a];\n	$profile{$n}{definition}=$de\
finitionL[$a];\n	$profile{$n}{identifyer}=$identif\
yerL[$a];\n	$profile{$n}{comment}=$comment[$a];\n	\
for ($b=0; $b<$L; $b++)\n	  {\n	    if ($p[$a][$b]\
)\n	      {\n		$profile{$n}{seq}.=$p[$a][$b];\n	  \
    }\n	    else\n	      {\n		$profile{$n}{seq}.=\\
"-\";\n	      }\n	  }\n      }\n    \n    $profile\
{n}=$nhits+1;\n    return %profile;\n  }\nsub ebi_\
blast_xml2profile \n  {\n    my ($name,$seq,$maxid\
, $minid, $mincov, $string)=(@_);\n    my ($L,$l, \
$a,$b,$c,$d,$nhits,@identifyerL,$identifyer);\n   \
 \n\n    \n    $seq=~s/[^a-zA-Z]//g;\n    $L=lengt\
h ($seq);\n    %hit=&xml2tag_list ($string, \"hit\\
");\n    \n    for ($nhits=0,$a=0; $a<$hit{n}; $a+\
+)\n      {\n	my ($ldb,$id, $identity, $expectatio\
n, $start, $end, $coverage, $r);\n	my (%Q,%M,%E,%I\
);\n	\n	$ldb=&tag2value ($hit{$a}{open}, \"databas\
e\");\n	$identifyer=&tag2value ($hit{$a}{open}, \"\
id\");\n\n	$description=&tag2value ($hit{$a}{open}\
, \"description\");\n	\n	%Q=&xml2tag_list ($hit{$a\
}{body}, \"querySeq\");\n	%M=&xml2tag_list ($hit{$\
a}{body}, \"matchSeq\");\n	%E=&xml2tag_list ($hit{\
$a}{body}, \"expectation\");\n	%I=&xml2tag_list ($\
hit{$a}{body}, \"identity\");\n	\n\n	for ($b=0; $b\
<$Q{n}; $b++)\n	  {\n	    \n	    \n	    $qs=$Q{$b}\
{body};\n	    $ms=$M{$b}{body};\n	    if ($seq eq\\
"\"){$seq=$qs;$L=length($seq);}\n\n	    $expectati\
on=$E{$b}{body};\n	    $identity=$I{$b}{body};\n	 \
   \n	    	    \n	    $start=&tag2value ($Q{$b}{op\
en}, \"start\");\n	    $end=&tag2value ($Q{$b}{ope\
n}, \"end\");\n	    $startM=&tag2value ($M{$b}{ope\
n}, \"start\");\n	    $endM=&tag2value ($M{$b}{ope\
n}, \"end\");\n	    $coverage=(($end-$start)*100)/\
$L;\n	    \n	   # print \"$id: ID: $identity COV: \
$coverage [$start $end]\\n\";\n	    \n	    \n	    \
if ($identity>$maxid || $identity<$minid || $cover\
age<$mincov){next;}\n	    # print \"KEEP\\n\";\n\n\
	    \n	    @lr1=(split (//,$qs));\n	    @lr2=(spl\
it (//,$ms));\n	    $l=$#lr1+1;\n	    for ($c=0;$c\
<$L;$c++){$p[$nhits][$c]=\"-\";}\n	    for ($d=0,$\
c=0; $c<$l; $c++)\n	      {\n		$r=$lr1[$c];\n		if \
( $r=~/[A-Za-z]/)\n		  {\n		    \n		    $p[$nhits]\
[$d + $start-1]=$lr2[$c];\n		    $d++;\n		  }\n	  \
    }\n	  \n	    \n	    $identifyerL[$nhits]=$iden\
tifyer;\n	    $comment[$nhits]=\"$ldb|$identifyer \
[Eval=$expectation][id=$identity%][start=$startM e\
nd=$endM]\";\n	    $nhits++;\n	  }\n      }\n    \\
n    $profile{n}=0;\n    $profile{$profile{n}}{nam\
e}=$name;\n    $profile{$profile{n}}{seq}=$seq;\n \
   $profile {n}++;\n    \n    for ($a=0; $a<$nhits\
; $a++)\n      {\n	$n=$a+1;\n	$profile{$n}{name}=\\
"$name\\_$a\";\n	$profile{$n}{seq}=\"\";\n	$profil\
e{$n}{identifyer}=$identifyerL[$a];\n	\n	$profile{\
$n}{comment}=$comment[$a];\n	for ($b=0; $b<$L; $b+\
+)\n	  {\n	    if ($p[$a][$b])\n	      {\n		$profi\
le{$n}{seq}.=$p[$a][$b];\n	      }\n	    else\n	  \
    {\n		$profile{$n}{seq}.=\"-\";\n	      }\n	  }\
\n      }\n    $profile{n}=$nhits+1;\n    \n    re\
turn %profile;\n  }\n\nsub blast_xml2hit_list\n  {\
\n    my $string=(@_[0]);\n    return &xml2tag_lis\
t ($string, \"hit\");\n  }\nsub xml2tag_list  \n  \
{\n    my ($string_in,$tag)=@_;\n    my $tag_in, $\
tag_out;\n    my %tag;\n    \n    if (-e $string_i\
n)\n      {\n	$string=&file2string ($string_in);\n\
      }\n    else\n      {\n	$string=$string_in;\n\
      }\n    $tag_in1=\"<$tag \";\n    $tag_in2=\"\
<$tag>\";\n    $tag_out=\"/$tag>\";\n    $string=~\
s/>/>##1/g;\n    $string=~s/</##2</g;\n    $string\
=~s/##1/<#/g;\n    $string=~s/##2/#>/g;\n    @l=($\
string=~/(\\<[^>]+\\>)/g);\n    $tag{n}=0;\n    $i\
n=0;$n=-1;\n  \n \n\n    foreach $t (@l)\n      {\\
n\n	$t=~s/<#//;\n	$t=~s/#>//;\n	\n	if ( $t=~/$tag_\
in1/ || $t=~/$tag_in2/)\n	  {\n	 \n	    $in=1;\n	 \
   $tag{$tag{n}}{open}=$t;\n	    $n++;\n	    \n	  \
}\n	elsif ($t=~/$tag_out/)\n	  {\n	    \n\n	    $t\
ag{$tag{n}}{close}=$t;\n	    $tag{n}++;\n	    $in=\
0;\n	  }\n	elsif ($in)\n	  {\n	   \n	    $tag{$tag\
{n}}{body}.=$t;\n	  }\n      }\n  \n    return %ta\
g;\n  }\n\n\n\n\n","use Env qw(HOST);\nuse Env qw(\
HOME);\nuse Env qw(USER);\nwhile (<>)\n  {\n    if\
 ( /^>(\\S+)/)\n      {\n	if ($list{$1})\n	  {\n	 \
   print \">$1_$list{$1}\\n\";\n	    $list{$1}++;\\
n	  }\n	else\n	  {\n	    print $_;\n	    $list{$1}\
=1;\n	  }\n      }\n    else\n      {\n	print $_;\\
n      }\n  }\n      \n","\n\n\nuse Env qw(HOST);\\
nuse Env qw(HOME);\nuse Env qw(USER);\n\n\nopen (F\
,$ARGV[0]);\nwhile ( <>)\n  {\n    @x=/([^:,;\\)\\\
(\\s]+):[^:,;\\)\\(]*/g;\n    @list=(@list,@x);\n \
 }\n$n=$#list+1;\nforeach $n(@list){print \">$n\\n\
sequence\\n\";}\n\n\nclose (F);\n","\nopen (F, $AR\
GV[0]);\n\nwhile ( <F>)\n  {\n    @l=($_=~/(\\S+)/\
g);\n    \n    $name=shift @l;\n    \n    print ST\
DOUT \"\\n>$name\\n\";\n    foreach $e (@l){$e=($e\
 eq \"0\")?\"O\":\"I\";print \"$e\";}\n  }\nclose \
(F);\n\n		       \n    \n","use strict;\nuse FileH\
andle;\nuse Env qw(HOST);\nuse Env qw(HOME);\nuse \
Env qw(USER);\nmy %name;\nmy $nseq;\nmy $F= new Fi\
leHandle;\nopen ($F, $ARGV[0]);\nwhile(<$F>)\n  {\\
n    \n    my $l=$_;\n    if ($l=~/^#/){;}\n    el\
sif (($l=~/\\d+\\s+\\d+\\s+(\\S+)\\s+(\\S+)/))\n  \
    {\n	my $name=$1;\n	my $seq=$2;\n	print \">$nam\
e\\n$seq\\n\";\n      }\n  }\nclose ($F);\nexit (0\
);\n\n\n","use Env qw(HOST);\nuse Env qw(HOME);\nu\
se Env qw(USER);\n\n$tmp=\"$ARGV[0].$$\";\nopen (I\
N, $ARGV[0]);\nopen (OUT, \">$tmp\");\n\nwhile ( <\
IN>)\n  {\n    $file=$_;\n    $file=~s/\\r\\n/\\n/\
g;\n    $file=~s/\\n\\r/\\n/g;\n    $file=~s/\\r\\\
r/\\n/g;\n    $file=~s/\\r/\\n/g;\n    print OUT \\
"$file\";\n  }\nclose (IN);\nclose (OUT);\n\nopen \
(OUT, \">$ARGV[0]\");\nopen (IN, \"$tmp\");\n\nwhi\
le ( <IN>)\n{\n  print OUT \"$_\";\n}\nclose (IN);\
\nclose (OUT);\nunlink ($tmp);\n\n"};
