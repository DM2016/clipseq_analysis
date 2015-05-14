package Clipseq;

use 5.010;
use strict;
use warnings;
use Carp; 
use File::Copy qw/ copy move /;
use feature qw/ state switch say /; #Modern perl features supported in 5.10 or newer.
require Exporter;


our @ISA = qw/ Exporter /;
our $VERSION = '0.01';
our @EXPORT_OK = qw/ is_valid_bed_line parse_bed_line read_bed_line write_bed_line discover_peaks build_genomic_coverage_strandwise widen_bed /; #Export nothing by default.

#Format structures, callable from other namespaces that use this module:
our @bed_fields   = qw/ chrom start end name score strand /; #The standard BED6 fields.
our @bed4_fields  = qw/ chrom start end strand /; #Some of pieces of bedtools pipe this out.
our @bed12_fields = qw/ chrom start end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts /; #Fields 7-12 named identically as described in the BED format specification on the UCSC FAQ, available at genome.ucsc.edu/faq at the time of this writing.

sub rewind {
  my $fh = shift;
  seek $fh, 0, 0;
}

sub parse_bed_line {
  my $line_ref = shift;
  chomp $$line_ref;
  my %line;
  @line{@bed_fields} = split /\t/ , $$line_ref;
  return \%line;
}

sub is_valid_bed_line {
  #Return true if BED line passed by scalar reference successfully parses.
  my $line_ref = shift;
  my $ignore_strands = shift; #Boolean to ignore strands.
  my $line = parse_bed_line( $line_ref );
  my $validity = 0; #Guilty until proven innocient.
  if ( ( $line->{chrom}  =~ /^chr/   ) &&
       ( $line->{start}  =~ /^\d+$/  ) &&
       ( $line->{end}    =~ /^\d+$/  ) &&
       ( $line->{end} > $line->{start} ) ) {
    $validity = 1; #Set to true: line looks good.
  }
  unless ( $ignore_strands ) {
    unless ( $line->{strand} =~ /^\+|-$/ ) {
      $validity = 0; #Invalid strand.
    }
  }
  $validity ? return $line : return $validity; #Returnary construct.
}

sub read_bed_line {
  #Given an opened filehandle reference, read the next BED line, parse it, and return a reference to a parsed hash of the BED line's contents.
  my $filehandle_ref = shift;
  my $ignore_strands = shift; #Optional Boolean to ignore strands.
  my $line = <$filehandle_ref>; #Scalar.
  return $line unless $line; #Return the unaltered undef return value from the readline operator <> if the EOF is reached.
  my $parsed_line = is_valid_bed_line( \$line , $ignore_strands ); #Will be a false value if the BED line is invalid.
  $parsed_line ? return $parsed_line : return read_bed_line( $filehandle_ref ); #Recurse if the line is invalid.
}

sub write_bed_line {
  #Given a reference to a parsed hash of a BED line's contents and an opened filehandle reference for writing, write the line to the BED file.
  my ( $filehandle_ref , $line ) = @_;
  print {$filehandle_ref} join "\t" , @{$line}{@bed_fields}; #Hash reference slice.
  print {$filehandle_ref} "\n";
}

sub count_bed_lines {
  #Given a BED file name as a scalar, read through the BED file and return the number of lines.
  my $bed_file = shift;
  my $lines = 0;
  open my $count, "<", $bed_file or croak "Failed to open $bed_file: $!\n";
  $lines++ while <$count>;
  close $count;
  return $lines;
}

sub extract_chromosomes {
  #Given a BED file, read through and return an array of chromosomes covered in the BED file.
  my $bed_file = shift;
  my %chroms;
  open my $chrom_collector , "<" , $bed_file or carp "Failed to open $bed_file: $!\n";
  while ( <$chrom_collector> ) {
    $chroms{ ( split /\t/ , $_ )[0] } = 1;
  }
  close $chrom_collector;
  return keys %chroms;
}

sub read_in_total_coverage {
  #Given a filehandle reference to a opened BED file, read in the coverage into a hierarchical data structurea nd return a reference.
  my $filehandle_ref = shift;
  my $ignore_strands = shift; #Optional second parameter to ignore strands in the coverage structure, useful for CHIP-seq peaks.
  my %coverage;
  while ( my $line = read_bed_line( $filehandle_ref , $ignore_strands ) ) {
    unless ( $ignore_strands ) {
      $coverage{ $line->{chrom} }{ $line->{strand} }{ $_ }++ for $line->{start} .. $line->{end};
    }
    else {
      $coverage{ $line->{chrom} }{ $_ }++ for $line->{start} .. $line->{end};
    }
  }
  return \%coverage;
}

sub discover_peaks {
  #Usage: discover_peaks( $parameters )

  #Discover peaks using a strandwise algorithm, comparing observed peak heights to iterative, transcript-wise randomizations.

  #Algorithm description:
  #1: Foreach chromosome, foreach strand, iterate through a transcriptome file (BED format) and extract transcripts on the current strand.
  #2: Extract all reads aligning to the current strand from the aligned reads file.
  #3: Foreach transcript on the current strand, iterate through reads on the strand.
  #4: If the read maps within that gene, push its length onto an array, and increment a coverage structure.
  #5: Compute p-value distribution for observed signal heights (number of signal intervals as tall or taller divided by total number of observed signal intervals).
  #6: For the specified number of iterations, randomly align intervals of identical length as the observed signal intervals to the current transcript.
  #7: Increment a random coverage array.
  #8: Collect a distribution of random signal interval heights and increment a height array.
  #9: Push a reference to that array onto an array of iterations.
  #10: Compute mean p-value and standard deviation for each signal height encountered across all randomization iterations.
  #11: Foreach observed signal interval, print it to output (score is max height) if it is at least as tall as the minimum height for the specified false discovery rate.

  my $parameters = shift; #Hash reference.

  #Assign local names to important parameters:
  my ( $bedfile , $gene_list , $total_iterations , $fdr , $outfile ) = ( $parameters->{aligned_reads_files_string} , $parameters->{transcriptome_file} , $parameters->{iterations} , $parameters->{fdr} , $parameters->{output_file} ); 

  $parameters->{verbose} and warn "Beginning peak discovery process.\n";

  #If $bedfile is a comma-separated list, concatenate the raw reads first.
  my ( $concatenated_bed_file , $did_concatenate , $bed );

  if ( $bedfile =~ /,/ ) {
    $parameters->{verbose} and warn "Detected replicate samples.\n";
    $did_concatenate = 1;
    my @bed_files = split /,/ , $bedfile;
    $concatenated_bed_file = $bed_files[0] . '_concatenated.' . rand(72585) . '.bed';
    if ( -f $concatenated_bed_file ) { #Overwrite file if it already exists (unlikely).
      `rm $concatenated_bed_file`;
    }
    open my $cat , ">>" , $concatenated_bed_file or croak "Failed to open file to write concatenated reads to: $concatenated_bed_file\nSystem says: $!\n";

    for my $bed_file ( @bed_files ) {
      open my $catenator , "<" , $bed_file or croak "Failed to open file to write concatenated reads from: $bed_file\nSystem says: $!\n";
      print {$cat} $_ while <$catenator>;
    }
    close $cat;

    $parameters->{verbose} and warn "Successfully concatenated replicates.\n";

    open $bed , "<" , $concatenated_bed_file or croak "Failed to open concatenated BED file named $concatenated_bed_file\nSystem says: $!\n";
  }
  else { #Only one replicate.
    open $bed , "<" , $bedfile or croak "Failed to open the reads BED file named $bedfile\nSystem says: $!\n";
  }

  open my $transcriptome, "<", $gene_list or croak "Failed to open $gene_list: $!\n";
  open my $out,   ">", $outfile   or croak "Failed to open $outfile: $!\n";

  #1: Foreach chromosome, foreach strand, iterate through a transcriptome file (BED format) and extract transcripts on the current strand.
  my @chromosomes = extract_chromosomes( $gene_list ); #Get chromosomes covered in the transcriptome.
  my $chromosome_count = scalar @chromosomes;
  $parameters->{verbose} and warn "Finished collecting chromosomes from the transcriptome.  Collected $chromosome_count chromosomes to search for peaks.\n";

  for my $current_chromosome ( @chromosomes ) {
    for my $strand ( '+' , '-' ) {
      #Get transcripts and aligned CLIP reads.
      my ( @transcripts, @reads, );
      rewind( $transcriptome );
      while ( my $transcript = read_bed_line($transcriptome) ) {
	push @transcripts, $transcript if $transcript->{chrom} eq $current_chromosome and $transcript->{strand} eq $strand;
      }

      #2: Extract all reads aligning to the current strand from the aligned reads file.
      rewind( $bed ); #To do: sort the BED file, then iterate through current strand of current chromosome only: no need to sweep over the entire BED file twice for each chromosome!
      while ( my $read = read_bed_line($bed) ) {
	push @reads, $read if $read->{chrom} eq $current_chromosome and $read->{strand} eq $strand;
      }

      #3: Iterate through the genes on that strand, for each of them,
      #   autovivify a coverage structure and save an array of the lengths of alignments that map to that gene.
      for my $current_transcript ( @transcripts ) {
	my ( @alignment_lengths, @observed_coverage, $hits, );

	for my $current_read ( @reads ) {
	  next unless $current_read->{start} >= $current_transcript->{start} and $current_read->{end} <= $current_transcript->{end};
	  my $aln_length = $current_read->{end} - $current_read->{start}; #Save this length for future randomization.
	  push @alignment_lengths, $aln_length;
	  my $cov_start = $current_read->{start} - $current_transcript->{start};
	  my $cov_end   = $current_read->{end}   - $current_transcript->{start};
	  $observed_coverage[$_]++ for $cov_start .. $cov_end;
	}

	#5: Compute p-value distribution for observed intervals of continuous CLIP signal:
	my ( %observed_peaks, $observed_peak_number, $total_observed_peaks );
	#Let "observed peaks" contain intervals of continuous CLIP signal bordered by signalless regions.
	#Discover interval starts:
	for my $pos ( 0 .. $#observed_coverage ) {
	  if ( defined $observed_coverage[$pos] and $pos == 0 ) { #First nt of gene has coverage (prevents $pos-1 becoming end of gene) 
	    $observed_peak_number++;
	    $observed_peaks{$observed_peak_number}{start} = $pos;
	    next;
	  }
	  elsif ( defined $observed_coverage[$pos] and ! defined $observed_coverage[$pos-1] ) {
	    $observed_peak_number++; #This is a nucleotide of positive coverage with nothing to the left.
	    $observed_peaks{$observed_peak_number}{start} = $pos;
	  }
	}

	#Discover interval ends, recording max peak height along the way:
	for my $peak ( sort { $a <=> $b } keys %observed_peaks ) {
	  $observed_peaks{$peak}{height} = 0;
	  my $peak_start = $observed_peaks{$peak}{start};
	  my $pos = $peak_start;
	  ++$pos until ! defined $observed_coverage[$pos+1]; #Walk to the next zero signal nucleotide.
	  $observed_peaks{$peak}{end} = $pos;
	  $total_observed_peaks++;
	  for ( $observed_peaks{$peak}{start} .. $observed_peaks{$peak}{end} ) {
	    $observed_coverage[$_] > $observed_peaks{$peak}{height} and $observed_peaks{$peak}{height} = $observed_coverage[$_];
	  }
	}

	#Collect distribution of observed CLIP signal interval heights:
	my @observed_heights;
	for my $peak ( sort { $a <=> $b } keys %observed_peaks ) {
	  my $height = $observed_peaks{$peak}{height};
	  $observed_heights[$height]++;
	}

	#Compute distribution of p-values for each height:
	my %observed_p_values;
	for my $height ( 1 .. $#observed_heights ) {
	  my $taller_peaks; #"As tall or taller", really.
	  for ( $height .. $#observed_heights ) {
	    if ( defined $observed_heights[$_] ) {
	      $taller_peaks += $observed_heights[$_];
	    }
	  }
	  $observed_p_values{$height} = $taller_peaks / $total_observed_peaks;
	}


	#6: Randomly align intervals of the observed alignment lengths for each of the specified iterations:
	my $i = 1;
	my %iteration;
	my $max_random_height = 0;

	until ( $i == ( $total_iterations + 1 ) ) {

	  #7: Randomly align reads to the current transcript, increment a random coverage array
	  my @random_coverage;
	  for my $aln_length ( @alignment_lengths ) {
	    my $last_mappable_pos = ( $current_transcript->{end} - $current_transcript->{start} ) - $aln_length;
	    my $random_start  = int( rand( $last_mappable_pos + 1 ) );
	    my $random_end    = $random_start + $aln_length;
	    for my $nt ( $random_start ..  $random_end ) {
	      $random_coverage[$nt]++;
	    }
	  }

	  #8: Collect distribution of random peak heights
	  my ( %random_peaks, $random_peak_number, @random_heights, $total_random_peaks );
	  #Collect randomized interval starts:
	  for my $pos ( 0 .. $#random_coverage ) {
	    if ( defined $random_coverage[$pos] and ! defined $random_coverage[$pos-1] ) {
	      ++$random_peak_number;
	      $random_peaks{$random_peak_number}{start} = $pos;
	    }
	  }

	  #Find each random interval's end and record maximum height along the way:
	  for my $random_peak ( sort { $a <=> $b } keys %random_peaks ) {
	    $random_peaks{$random_peak}{height} = 0;
	    my $random_peak_start = $random_peaks{$random_peak}{start};
	    my $pos = $random_peak_start;
	    ++$pos until ! defined $random_coverage[$pos+1];
	    $random_peaks{$random_peak}{end} = $pos;
	    ++$total_random_peaks;

	    for ( $random_peaks{$random_peak}{start} .. $random_peaks{$random_peak}{end} ) {
	      $random_coverage[$_] > $random_peaks{$random_peak}{height} and $random_peaks{$random_peak}{height} = $random_coverage[$_];
	    }
	  }

	  #Collect distribution of randomized intervals' heights:
	  for my $random_peak ( sort { $a <=> $b } keys %random_peaks ) {
	    my $random_peak_height = $random_peaks{$random_peak}{height};
	    $random_heights[$random_peak_height]++;
	  }

	  #Translate distribution of random peak heights to distribution of p-values
	  my %p_values;
	  for my $random_height ( 1 .. $#random_heights ) {
	    my $taller_random_peaks; #as tall or taller random peaks, really
	    for ( $random_height .. $#random_heights ) {
	      if ( defined $random_heights[$_] ) {
		$taller_random_peaks += $random_heights[$_];
	      }
	    }

	    $p_values{$random_height} = $taller_random_peaks / $total_random_peaks;
	  }

	  for my $height ( sort { $a <=> $b } keys %p_values ) {
	    $height > $max_random_height and $max_random_height = $height;
	  }

	  $iteration{$i} = \%p_values;
	  $i++;
	}

	#10. Compute average p-val for each random height:
	my ( %total_p_values, %mean_p_values, %stdevs, $total_squared_displacements );
	for my $height ( 1 .. $max_random_height ) {
	  #First, compute the sum of p-values for each random height
	  for my $j ( sort { $a <=> $b }  keys %iteration ) {
	    if ( defined $iteration{$j}{$height} ) {
	      $total_p_values{$height} += $iteration{$j}{$height};
	    }
	  }

	  #Then divide by the number of iterations to get average (mean) p-val for this height
	  $mean_p_values{$height} = $total_p_values{$height} / $total_iterations;

	  #Then compute the value of one standard deviation for this height:
	  my $total_squared_displacement;
	  for my $j ( sort { $a <=> $b } keys %iteration ) {
	    if ( defined $iteration{$j}{$height} ) {
	      my $displacement = $iteration{$j}{$height} - $mean_p_values{$height};
	      my $squared_displacement = $displacement ** 2;
	      $total_squared_displacements += $squared_displacement;
	    }
	  }

	  $stdevs{$height} = ( $total_squared_displacements / $total_iterations ) ** 0.5;
	}

	#10b: Compute FDR distribution:
	my ( %fdrs, $min_peak_height );
	for my $height ( 1 .. $max_random_height ) {
	  if ( defined $observed_p_values{$height} ) {
	    $fdrs{$height} = ( $mean_p_values{$height} + $stdevs{$height} ) / $observed_p_values{$height};
	  }

	  unless ( defined $min_peak_height ) {
	    if ( $fdrs{$height} ) {
	      if ( $fdrs{$height} < $fdr ) {
		$min_peak_height = $fdrs{$height};
	      }
	    }
	  }
	}

	$min_peak_height ||= $max_random_height;

	my $significant_peaks;
	for my $peak ( keys %observed_peaks ) {
	  if ( $observed_peaks{$peak}{height} > $min_peak_height ) {
	    $observed_peaks{$peak}{start} += $current_transcript->{start};
	    $observed_peaks{$peak}{end}   += $current_transcript->{start};
	    print {$out} "$current_chromosome\t$observed_peaks{$peak}{start}\t$observed_peaks{$peak}{end}\t$current_transcript->{name}\t$observed_peaks{$peak}{height}\t$strand\n";
	    $significant_peaks++;
	  }
	}
      } #End transcript.
    } #End strand.
  } #End chromosome.

  close $transcriptome;
  close $out;
  close $bed;

  my $peak_count = count_bed_lines( $outfile );
  $parameters->{verbose} and warn "Finished collecting $peak_count raw peaks on potentially overlapping transcripts. Merging overlapping peaks into unique genomic intervals...\n";

  #Post-process peaks: first, combine peaks that were discovered on more than one transcript and merge overlapping peaks.
  my ( $reads_in_peaks_file , $flattened_reads_file ) = ( $outfile . '.reads.bed' , $outfile . '.flattened_reads.bed' );
  if ( $did_concatenate ) {
    bedtools_intersect_u( $concatenated_bed_file , $outfile , $reads_in_peaks_file ); #Grab the reads that are in peaks from the concatenated reads file.
  }
  else {
    bedtools_intersect_u( $bedfile , $outfile , $reads_in_peaks_file ); #Grab the reads that are in peaks from the input aligned reads file.
  }

  $parameters->{verbose} and warn "Encountered " , count_bed_lines( $reads_in_peaks_file ) , " reads in peaks.\n";

  slice_genomic_coverage( $reads_in_peaks_file , 1 , $flattened_reads_file ); #Flatten the reads into continuous intervals.
  unlink $reads_in_peaks_file; #Remove temporary file containing reads that map to peaks.
  move $flattened_reads_file , $outfile; #Replace the peaks file with the post-processed peaks file.

  $parameters->{verbose} and warn "Finished merging overlapping peaks. Total peaks: " , count_bed_lines( $outfile ) , "\n";

  #Remove the concatenated BED file if it exists:
  if ( $did_concatenate ) {
    unink $concatenated_bed_file;
  }
} #Function return.

sub slice_genomic_coverage {
  #Print the coordinates of intervals of continous coverage with a maximum overlap greater than or equal to the specified height.
  #Usage: build_genomic_coverage_strandwise( $bedfile , $min_peak_height , $outfile );

  #1) Autovivify genomic coverage from a BED file
  #2) Interval from start to end of coverage is printed if max height is above the $min_peak_height threshold.
  #2) Iterate through each strand of each chromosome, autovivifying a coverage hash from just that strand

  my ( $bed_file , $min_peak_height , $out_file ) = @_;
  open my $bed, "<", $bed_file or croak "$bed_file: $!\n";
  my $out;
  open $out, ">", $out_file or croak "$out_file: $!\n"; #Clobber.

  my $peak_number;
  #Iterate by strand:
  my @chromosomes = extract_chromosomes( $bed_file ); #Hash reference.
  for my $chrom ( @chromosomes ) {
    for my $strand ( '+' , '-' ) {
      #Read through the BED file, autovivify genomic structure, collect coverage values:
      my @coverage;

      rewind $bed;
      while ( my $line = read_bed_line($bed) ) {
	next unless $line->{chrom} eq $chrom and $line->{strand} eq $strand;

	for my $nt ( $line->{start} .. $line->{end} ) {
	  $coverage[$nt]++;
	}
      }

      my %peaks; #Intervals of continuous nonzero coverage.

      if ( $coverage[0] ) { #Check to see if the first nucleotide has coverage.
	$peak_number++;
	$peaks{$peak_number}{start} = 0;
      }
      for my $position ( 1 .. $#coverage ) {
	my $peak_max_height = 0;
	if ( $coverage[$position] ) {
	  my $minus_one = $position - 1;
	  my $plus_one  = $position + 1;
	  unless ( $coverage[$minus_one] ) { #This is a peak start, the previous nucleotide wasn't covered.
	    $peak_number++;
	    $peaks{$peak_number}{start} = $position;
	    next;
	  }
	  elsif ( ! defined $coverage[$plus_one] ) {
	    #This is a peak end.
	    $peaks{$peak_number}{end} = $position;

	    #Find peak max height, save it as the score:
	    for my $nt ( $peaks{$peak_number}{start} .. $peaks{$peak_number}{end} ) {
	      $coverage[$nt] > $peak_max_height and $peak_max_height = $coverage[$nt];
	    }
	    if ( $peak_max_height >= $min_peak_height ) {
	      my $peak_name = 'peak_' . $peak_number;
	      print {$out} "$chrom\t$peaks{$peak_number}{start}\t$peaks{$peak_number}{end}\t$peak_name\t$peak_max_height\t$strand\n";
	    }
	  }
	}
      }
    }
  }
}

sub widen_bed {
  #Symmetrically widen (or shorten) a BED file by the specified number of nucleotides.
  #Usage: widen_bed $bed_file $nts $output_file
  #Warning: When shortening BED intervals, there is currently no protection from generating intervals with end coordinates less than start coordinates!
  my ( $input_file , $nts , $output_file ) = @_;
  open my $in  , "<" , $input_file  or carp "Failed to open $input_file:  $!\n";
  open my $out , ">" , $output_file or carp "Failed to open $output_file: $!\n";
  while ( my $line = read_bed_line( $in ) ) {
    $line->{start} -= $nts;
    $line->{end}   += $nts;
    write_bed_line( $line , $output_file );
  }
}

sub bedtools_intersect_u { 
  #Run bedtools intersect.
  #Usage: bedtools_intersect_u( $a_file , $b_file , $output_file );
  my ( $a_file , $b_file , $output_file ) = @_;
  `bedtools intersect -u -s -a $a_file -b $b_file > $output_file`;
}

sub bedtools_intersect_v {
  #Run bedtools interect negatively.
  #Usage: bedtools_intersect_v( $a_file , $b_file , $output_file );
  my ( $a_file , $b_file , $output_file ) = @_;
  `bedtools intersect -v -s -a $a_file -b $b_file > $output_file`;
}

1;

__END__

=head1 NAME

Clipseq.pm - Perl module for the analysis of CLIP-seq data.

=head1 SYNOPSIS

  This module provides an interface to low-level manipulations of CLIP-seq data and is meant to be called from a higher-level program, such as the discover_peaks script provided with the clipseq_analysis distribution.  See the embedded documentation for that script with `perldoc`.
  The functions in this module handle file I/O, format parsing, building genomic coverage data structures, and iterative randomization for CLIP-seq signal peak detection.

  Some functions in this module may use reimplementations of genome arithmetic that could be performed faster with bedtools.  The ability to integrate bedtools calls is provided, but this module does not require bedtools to function.

=head2 EXPORT

  None by default. Exportable functions are documented here:

=head1 AUTHOR

Brian Sebastian Cole, E<lt>colebr@mail.med.upenn.eduE<gt>
The empirical algorithm used to discover significant CLIP-seq sites ("peaks") was first described by Gene Yeo et al in the methods section of scientific manuscripts.  The code implementation contained in this distribution is entirely the author's work.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012,2013, and 2014 by Brian Sebastian Cole

This library is free software; you can redistribute it and/or modify it under the same terms as Perl itself, either Perl version 5.14.2 or, at your option, any later version of Perl 5 you have available.

=cut
