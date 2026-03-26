=head1 LICENSE

 Copyright 2013 EMBL - European Bioinformatics Institute

 Licensed under the Apache License, Version 2.0 (the
 "License"); you may not use this file except in
 compliance with the License.  You may obtain a copy of
 the License at

     http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, 
 software distributed under the License is distributed on 
 an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY 
 KIND, either express or implied. See the License for the 
 specific language governing permissions and limitations 
 under the License.
                                                                                                                     
=head1 CONTACT                                                                                                       

 Graham Ritchie <grsr@ebi.ac.uk>
    
=cut

=head1 NAME

 Gwava

=head1 SYNOPSIS

 mv Gwava.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin Gwava,tss,/path/to/gwava_scores.bed.gz

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that 
 retrieves precomputed Genome Wide Annotation of VAriants (GWAVA) scores
 for any  variant that overlaps a known variant from the Ensembl variation 
 database. It adds one new entry class to the VEP's Extra column, GWAVA 
 which is the GWAVA score. Two arguments are required: the classifier model 
 to use, one of "region", "tss", "unmatched", and the path to the score
 file (which is available at the same location as this plugin).

=cut

package Gwava;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub feature_types {
    return ['Feature', 'Intergenic'];
}

sub get_header_info {
    my $self = shift;
    return {
        GWAVA => "Genome Wide Annotation of VAriants score (".$self->{model}." model)",
    };
}

sub new {
    my $class = shift;
    
    my $self = $class->SUPER::new(@_);

    $self->{model} = $self->params->[0];
    $self->{file} = $self->params->[1];
   
    unless ($self->{model} =~ /^region|tss|unmatched$/) {
        die "Unrecognised model option, use one of 'region', 'tss', 'unmatched'\n";
    }

    unless (-e $self->{file}) {
        die "Can't find the score file: ".$self->{file}."\n";
    }

    unless (`which tabix`) {
        die "The tabix program is required to run this plugin\n";
    }

    return $self;
}

sub run {
    my ($self, $bvfoa) = @_;
 
    my $bvf = $bvfoa->base_variation_feature;

    my $chr = "chr".$bvf->seq_region_name;
    my $start = $bvf->seq_region_start;
    my $end = $bvf->seq_region_end;

    my $key = join("_", $chr, $start, $end);

    unless (exists $self->{cache}->{$key}) {
        my $cmd = "tabix ".$self->{file}." ${chr}:${start}-${end}";
        my @res = `$cmd`;
        if (@res > 0) {
            # just take the first result
            my ($c, $s, $e, $id, $reg, $tss, $unm) = split "\t", $res[0];
            my $score = $self->{model} eq 'region' ? $reg : $self->{model} eq 'tss' ? $tss : $unm; 
            $self->{cache}->{$key} = $score;
        }
        else {
            $self->{cache}->{$key} = undef;
        }
    }
   
    if (my $score = $self->{cache}->{$key}) {
        return { GWAVA => $score };
    }

    return {};
}

1;

