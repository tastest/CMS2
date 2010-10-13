#!/bin/env perl
# Autor: D.Kovalskyi (UCSB)
use warnings;
use strict;

die << 'EOF' unless (@ARGV==1);
Usage:
    pickEvents.pl <file>

The file should have the following format:
    * first line is dataset name
    * all other lines should be: RUN, LUMI, EVENT

EOF

my $template = << 'EOF';
import FWCore.ParameterSet.Config as cms
process = cms.Process("PE")
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(),
   eventsToProcess = cms.untracked.VEventRange( EVENTS )
)
process.source.fileNames.extend([ FILES ])
process.Out = cms.OutputModule("PoolOutputModule",
   fileName = cms.untracked.string(OUTPUT)
)
process.e = cms.EndPath(process.Out)
EOF

# Checking environment
die "CMS_PATH is not set. Abort\n" unless defined $ENV{CMS_PATH};
die "site-local-config.xml is not found\n" unless -e "$ENV{CMS_PATH}/SITECONF/local/JobConfig/site-local-config.xml";

# Check release consistency
die "CMSSW_VERSION is not set. Please set proper CMSSW environment. Abort\n" unless defined $ENV{CMSSW_VERSION};
my ($rv1,$rv2) = ($ENV{CMSSW_VERSION} =~ /CMSSW_(\d+)_(\d+)_/);

die "This CMSSW release is not supported. You have to use at least CMSSW_3_0_0 or later\n" unless ($rv1>2);
 
my @events = ();
my $dataset;
my %files = ();

open(IN,$ARGV[0]) || die "failed to open file $ARGV[0]\n$!\n";
print "querying DBS\n";
while (my $line = <IN>){
    $line =~ s/\n//;
    if ( !defined $dataset ){
	$dataset = $line;
	next;
    }
    my @entry = split(/\D+/,$line);
    die "Wrong file format. Expected 3 numbers got: $line\n" unless (@entry==3);
    push @events, \@entry;
    my @result = split(/\n/,`python \$DBSCMD_HOME/dbsCommandLine.py --noheader -c search --url= http://cmsdbsprod.cern.ch/cms_dbs_caf_analysis_01/servlet/DBSServlet --query="find file where run=$entry[0] and dataset=$dataset and lumi=$entry[1]"`);
    print "Warning no files are available at the local site for: $line\n" unless @result > 0;
    foreach my $file (@result){
	$file =~ s/\n//;
	$files{$file}++;
    }
    print ".";
}
print "\ndone\n";
print "Total number of files to process: ", scalar keys %files, "\n";

my $datasetRelease = `python \$DBSCMD_HOME/dbsCommandLine.py --noheader -c search --query="find dataset.release where dataset=$dataset"`;
$datasetRelease =~ s/\n//;
die "Failed to determine dataset release\n" unless ($datasetRelease =~ /CMSSW_\d+_\d+_/);
my ($dv1,$dv2) = ($datasetRelease =~ /CMSSW_(\d+)_(\d+)_/);

if ( $dv1>$rv1 || ( $dv1==$rv1 && $dv2>$rv2 )){
    die "ERROR: Working release is too old, make a project area with newer release.\n".
	"\tDataset release: $datasetRelease\n\tProject release: $ENV{CMSSW_VERSION}\n";
}

if ( $dv1<$rv1 || ( $dv1==$rv1 && $dv2<$rv2 )){
    print "WARNING: You may need to use an older release to match version.\n".
	"\tDataset release: $datasetRelease\n\tProject release: $ENV{CMSSW_VERSION}\n";
}

my $fileName = $dataset;
$fileName =~ s/^\/+//;
$fileName =~ s/\//_/g;

my @fileList = keys %files;

foreach my $i (0..$#fileList){
    $fileList[$i] =~ s/^[\'\"]*(.*?)[\'\"]*$/\'$1\'/;
}
my $vfiles = join(",",@fileList);

my $vevents;
foreach my $entry(@events){
    $vevents .= "," if ( defined $vevents );
    $vevents .= "\'$entry->[0]:$entry->[2]\'";
}

# prepare config file
$template =~ s/FILES/$vfiles/ms;
$template =~ s/EVENTS/$vevents/ms;
$template =~ s/OUTPUT/\'$fileName.root\'/ms;

open(OUT,">$fileName.py");
print OUT "$template\n";
close OUT;

print "Extracting events using cmsRun\n";
system("cmsRun $fileName.py");

exit
