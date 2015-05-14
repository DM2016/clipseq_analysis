package Z_score;

use Carp;

#Subroutines for the computation of z_scores for each kmer in a given set of RNA binding protein sites relative to the bound transcripts.

#Contents:

# read_in_transcripts_on_strand() Reads through a reference to a filehandle from a file that contains transcripts in BED format, and returns a reference to an array of BED hashrefs.

# read_in_binding_sites_on_strand() Reads through a reference to a filehandle from a file that contains binding sites in BED format, and returns a reference to an array of BED hashrefs.

# collect_bound_transcripts_on_strand() Collects the transcripts specified in a reference to an array of BED hashrefs bound by the binding sites specified in a reference to an array of BED hashrefs and returns a reference to an array of BED hashrefs containing the transcripts which were bound.

# convert_bed_hashes_to_sequences() Given a reference to an array of BED hashrefs, return a reference to a hash of sequences keyed by the BED name

# randomize_binding_sites_within_transcripts() Given a reference to an array of BED hashrefs containing transcripts, randomize the binding sites specified in a refernce to an array of BED hashrefs without overlapping the binding sites, return a reference to an array of BED hashrefs representing the positions of the randomized binding sites.

# convert_bed_hash_to_sequence() Given a reference to a BED hash and a path to a directory containing the chromosome FASTA files, return a reference to a hash of the sequence keyed by the BED name

# randomize_binding_sites_within_transcript() Given a reference to a an array of BED hashrefs containing intervals within a transcript and a reference to a BED hash of the transcript, return a reference to an array of BED hashrefs containing randomized intervals, without any overlap of the original intervals

my @bed_fields = qw/ chrom start end name score strand /; #The six BED fields, an array used for slicing the fields out of BED files into hashes.

sub bed_to_fa {
    #use bedtools getfasta -s to convert BED intervals to FASTA, in hg19 (hardcoded paths, this is my shit)
    my ( $bed_file , $fasta_file ) = @_;
    `bedtools getfasta -s -fi /home/colebr/data/hg19/hg19.fa -bed $bed_file -fo $fasta_file`;
}

sub count_kmers {
    #counts kmers in a fasta file, returns a reference to an unsorted hash of their frequency
    my ( $fasta_file , $k ) = @_;

    open my $in, "<", $fasta_file or croak "$fasta_file: $!\n";
    
    my ( %total_kmer_hits , $total_kmers , %kmer_frequency );

    while ( <$in> ) {
	chomp;
	next if /^>/; #skip FASTA name lines
	my $seq = $_;
	my $length = length $seq;
	next unless $length >= $k; #skip sequences shorter than k
	my $last_start_index = $length - $k;
	
	for my $pos ( 0 .. $last_start_index ) {
	    my $sequence = uc( substr $seq , $pos , $k );
	    $total_kmers++;
	    $total_kmer_hits{ $sequence }++;
	}
    }
    for my $kmer ( keys %total_kmer_hits ) {
	$kmer_frequency{ $kmer } = $total_kmer_hits{ $kmer } / $total_kmers;
    }
    close $in;
    return \%kmer_frequency;
}
	    
	    


sub read_in_transcripts_on_strand {
    #usage: $refence_to_array_of_BED_hashrefs_containing_transcripts = read_in_transcripts_on_strand( $transcripts_BED_file , $chromosome , $strand );

    my @array_of_BED_hashrefs; #A container in which to push references to hashes for the transcripts on the strand of interest.

    my ( $transcripts_BED_file , $chromosome , $strand ) = @_;
    
    #1: Open transcripts BED file
    open my $transcripts, "<", $transcripts_BED_file or die "Failed to open $transcripts_bed_file to read in transcripts.  System says: $!\n";
    
    #2: Read through the file, parsing the lines, skipping the lines that aren't on the strand of interest

    while ( <$transcripts> ) {
	next unless /^chr/; #Skip non-BED lines.
	chomp;
	my %line;
	@line{@bed_fields} = split /\t/;
	next unless $line{strand} eq $strand and $line{chrom} eq $chromosome;
	push @array_of_BED_hashrefs, \%line;
    }

    close $transcripts;
    #3: Return a reference to the array of BED hashrefs containing the transcripts on the strand of interest
    return \@array_of_BED_hashrefs;
}

sub read_in_binding_sites_on_strand {
    my @array;
    my ( $binding_sites_file , $chromosome , $strand ) = @_;

    open my $binding_sites, "<", $binding_sites_file or die "Failed to open $binding_sites_file to read in binding sites.  System says: $!\n";

    while ( <$binding_sites> ) {
	next unless /^chr/;
	chomp;
	my %line;
	@line{@bed_fields} = split /\t/;

	next unless $line{strand} eq $strand;
	next unless $line{chrom} eq $chromosome;
	push @array, \%line;
    }

    close $binding_sites;
    return \@array;

}

sub collect_bound_transcripts_on_strand {
    #usage: $reference_to_array_of_bound_transcipts_BED_hashrefs = collect_bound_transcripts_on_strand( $reference_to_array_of_transcripts_BED_hashrefs , $reference_to_array_of_binding_site_hashrefs , $chrom , $strand );

    my ( $reference_to_array_of_transcripts_BED_hashrefs , $reference_to_array_of_binding_site_hashrefs ) = @_;
    
    my @array_of_bound_transcripts_BED_hashrefs;

    for my $transcript_hashref ( @$reference_to_array_of_transcripts_BED_hashrefs ) {
	my $binding_sites_within_transcript; #Scalar containing the total nubmer of binding sites within that transcript.
	for my $binding_site_hashref ( @$reference_to_array_of_binding_site_hashrefs ) {
	    $binding_sites_within_transcript++ if $binding_site_hashref->{start} >= $transcript_BED_hashref->{start} and $binding_site_hashref->{end} <= $transcript_BED_hashref->{end};
	}
	
	$binding_sites_within_transcript and push @array_of_bound_transcripts_BED_hashrefs, $transcript_hashref;
    }
    
    return \@array_of_bound_transcripts_BED_hashrefs;

}

sub randomize_binding_sites_within_transcript {
    #usage: $reference_to_array_of_randomized_binding_sites = randomize_binding_sites_within_transcript( $reference_transcript_BED_hash , $reference_to_array_of_binding_sites_BED_hashrefs );

    my ( $reference_to_transcript_BED_hash , $reference_to_array_of_binding_sites_BED_hashrefs ) = @_;
    
    my @binding_site_coverage; #An array of nucleotides indexed by the offset within the transcript, whose indexed values are the coverage at that position seen in the binding sites, to be used to check if randomized binding sites overlap the binding sites, and if they do, rerandomize them.
    
    my $transcript_length = $reference_to_transcript_BED_hash->{end} - $reference_to_transcript_BED_hash->{start};

    my @array_of_randomized_binding_sites; #An array onto which to push references to BED hashes of randomized binding sites that don't overlap the observed binding site.

    for my $binding_site_BED_hashref ( @$reference_to_array_of_binding_sites_BED_hashrefs ) {

	my $binding_site_length = $binding_site_BED_hashref->{end} - $binding_site_BED_hashref->{start};
	my $random_start_offset   = int( rand( $transcript_length ) );
	my $random_end_offset     = $random_start_offset + $binding_site_length;

	#Check if the randomized alignment position overlaps observed coverage in the binding sites, and if it does, rerandomize
	my $overlaps_observed_coverage; #Will be true if the randomized site overlaps observed coverage.
	for my $pos ( $random_start_offset .. $random_end_offset ) {
	    $binding_site_coverage[ $pos ] and $overlaps_observed_coverage += $binding_site_coverage[ $pos ];
	}
	
	redo if $overlaps_observed_coverage;
	
	my $random_start_position = $reference_to_transcript_BED_hash->{start} + $random_start_offset;
	my $random_end_position   = $random_start_position + $binding_site_length;
	
	my %randomized_binding_site;
	@randomized_binding_site{@bed_fields} = ( $refence_to_transcript_BED_hash->{chrom} , $random_start_position , $random_end_position , 1 , 1, $reference_to_transcript_BED_hash->{strand} );
	push @array_of_randomized_binding_sites, \%randomized_binding_site;
    }
    return \@array_of_randomized_binding_sites;
}

sub randomize_binding_sites_within_transcripts {
    #usage: $refence_to_array_of_randomized_binding_sites_BED_hashrefs = randomize_binding_sites_within_transcirpts( $refence_to_array_of_transcript_BED_hashrefs , $refence_to_array_of_binding_sites_BED_hahrefs );

    my ( $reference_to_array_of_transcript_BED_hashrefs , $refence_to_array_of_binding_site_hashrefs ) = @_;

    my @array_of_randomized_binding_sites; #An array onto which to push references to BED hashrefs of randomized binding sites.
    
    #For each transcript, call randomize_binding_sites_within_transcript() on the binding sites within that transcript.

    for my $transcript_BED_hashref ( @$refence_to_array_of_transcript_BED_hashrefs ) {
	my @binding_sites_within_transcript; #An array onto which to push references to BED hashrefs for the binding sites within the transcript of interest.
	for my $binding_site_BED_hashref ( @$reference_to_array_of_binding_site_hashrefs ) {
	    next unless $binding_site_BED_hashref->{start} >= $transcript_BED_hashref->{start} and $binding_site_BED_hashref->{end} <= $transcript_BED_hashref->{end};
	    push @binding_sites_within_transcript, $binding_site_BED_hashref;
	}
	
	push @array_of_randomized_binding_sites, \@binding_sites_within_transcript;
    }

    return \@array_of_randomized_binding_sites;
}

sub convert_BED_hash_to_sequence {
    #usage: $reference_to_sequence = convert_BED_hash_to_sequence( $reference_to_BED_hash , $directory_containing_chromosome_fastas );

    my ( $reference_to_BED_hash , $chromosome_dir ) = @_;
    
    my $chromosome_fasta_file = $chromosome_dir . '/' . $reference_to_BED_hash->{chrom} . '.fa';

    open my $chrom, "<", $chromosome_fasta_file or croak "convert_BED_hash_to_sequence failed to open $chromosome_fasta_file: $!\n";

    my $chromosome_sequence_string; #A single sequence string containing all the nucleotides of the chromosome.

    #Add the sequences from the FASTA file onto this string
    while ( <$chrom> ) {
	next unless /^[ACTGNactgn]+$/; #skip fasta headers
	chomp;
	$chromosome_sequence_string .= $_;
    }

    croak "convert_BED_hash_to_sequence failed to build a chromosome sequence string from $chromosome_fasta_file.\nFatality!\n" unless defined $chromosome_sequence_string;

    my $sequence;
    if ( $reference_to_BED_hash->{strand} eq "+" ) {
	$sequence = substr $chromosome_sequence_string, $reference_to_BED_hash->{start}, $reference_to_BED_hash->{end};
    }
    elsif ( $reference_to_BED_hash->{strand} eq "-" ) {
	my $reverse_compliment_of_sequence = substr $chromosome_sequence_string, $reference_to_BED_hash->{start}, $reference_to_BED_hash->{end};
	my $reverse_of_sequence = reverse( $reverse_complement_of_sequence );
	( $sequence = $reverse_of_sequence ) =~ tr/ACTGNactgn/TGACNtgacn/;
    }

    croak "convert_BED_hash_to_sequence failed to build a sequence string through substr.\nFatality!\n" unless defined $sequence;
    my $upper_case_sequence = uc( $sequence );
    return \$upper_case_sequence;
}

sub convert_BED_hashes_to_sequences {
  #usage: $reference_to_array_of_sequences = convert_BED_hashes_to_sequences( $reference_to_array_of_BED_hashes , $directory_containing_chromosome_fastas );

    my ( $reference_to_array_of_BED_hashes , $chromosome_dir ) = @_;

    my @array_of_references_to_sequences;
    
    for my $reference_to_BED_hash ( @$reference_to_array_of_BED_hashes ) {
	my $reference_to_sequence = convert_BED_hash_to_sequence( $reference_to_BED_hash , $chromosome_dir );
	push @array_of_references_to_sequences , $reference_to_sequence;
    }
    
    return \@array_of_references_to_sequences;
}


1;
