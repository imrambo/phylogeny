#!/usr/bin/python3
import argparse
import re
import logging
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
"""
Ian Rambo 2019
Purpose: preprocessing of multiple sequence alignment(s).
The user can remove columns or sequences with ambigous characters,
keep only unique sequences and headers, and change the MSA format.

Thirteen... that's a mighty unlucky number... for somebody!
"""
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument('--input', type = str, dest = 'input_file',
action = 'store', help = 'path to input MSA file.')
parser.add_argument('--output', type = str, dest = 'output_file',
action = 'store', help = 'path to output MSA file.')
parser.add_argument('--rm_ambigs', default = False, action = 'store_true',
dest = 'rm_ambig', help = 'Boolean. Optional. Remove sequences with ambiguous characters.')
parser.add_argument('--ambig_strategy', action = 'store', dest = 'ambig_strategy',
nargs = '?', type = str, help = 'Use along with rm_ambigs. Remove columns or sequences. Choose "column" or "sequence".')
parser.add_argument('--unique', default=False, action = 'store_true',
dest = 'unique', help = 'Boolean. Optional. Only store unique sequences and identifiers.')
parser.add_argument('--seq_type', type = str, action = 'store', dest = 'seq_type',
nargs = '?', help='Sequence type, specify "dna", "rna", or "amino". Optional, used for ambiguous character removal.')
parser.add_argument('--input_format', type = str, dest = 'infmt', action = 'store',
help="""input MSA format. Required. Valid formats: clustal, emboss, fasta, fasta-m10, ig,
maf, mauve, nexus, phylip, phylip-sequential, phylip-relaxed, stockholm""")
parser.add_argument('--output_format', type = str, dest = 'outfmt', action = 'store',
help="""output MSA format. Required. Valid formats: clustal, emboss, fasta, fasta-m10, ig,
maf, mauve, nexus, phylip, phylip-sequential, phylip-relaxed, stockholm""")
parser.add_argument('--logfile', type = str, dest = 'logfile', action = 'store',
nargs = '?', help = 'path to optional logfile')
opts = parser.parse_args()

#==============================================================================
#Valid and ambiguous amino acids
valid_aa = ['A','C','D','E','F','G','H','I','K','L','M','N',
'P','Q','R','S','T','V','W','Y']
ambig_aa = ['B','J','O','U','X','Z']
#Ambiguous nucleotides
ambig_nuc = ['R','Y','W','S','K','M','D','V','H','B']
#Characters illegal in header
illegal_chars = re.compile(r'\:|\,|\)|\(|\;|\,|\]|\[|\,|\'')
#MSA formats compatible with Bio.AlignIO
valid_formats = ['clustal', 'emboss', 'fasta', 'fasta-m10', 'ig', 'maf',
'mauve', 'nexus', 'phylip', 'phylip-sequential', 'phylip-relaxed', 'stockholm']
valid_seqtype = ['amino', 'nuc']
#==============================================================================
# def msa_convert(msain, infmt, msaout, outfmt):
#     '''
#     Change the format of an msa file.
#     '''
#     outfmt = outfmt.lower()
#     if not outfmt in valid_formats:
#         logging.warning('Specified MSA format is not valid! Choose from: %s' % ','.join(valid_formats))
#         return
#     AlignIO.convert(msain, infmt, msaout, outfmt)
#     print('converted %s from %s to %s' % (msain, infmt, outfmt))
#     return None
# #------------------------------------------------------------------------------
# def msa_format(alignment, format):
#     '''
#     Change the format of an alignment object.
#     '''
#     if not format in valid_formats:
#         logging.warning('Specified MSA format is not valid! Choose from: %s' % ','.join(valid_formats))
#         return
#     fmt_alignment = alignment.format(format)
#     print("Alignment changed from '%s' to '%s' format" % (opts.infmt, opts.outfmt))
#     return fmt_alignment
#------------------------------------------------------------------------------
def rm_ambig(alignment, strategy, ambig_list):
    """
    Remove either columns or entire sequences in an alignment
    object that contain ambiguous characters specified in a list.
    Returns a MultipleSequenceAlignment object.
    """
    if strategy == 'column':
        print('Strategy: remove columns with ambiguous characters')
        keep_cols = []
        rm_cols = []
        for i in range(alignment.get_alignment_length()):
            #Check if any ambiguous characters in the columns.
            if any([ambig in alignment[:, i] for ambig in ambig_list]):
                rm_cols.append(i)
            else:
                keep_cols.append(i)
        if rm_cols:
            print('ambiguous symbols found in column(s): %s' % ','.join([str(x) for x in rm_cols]))
        else:
            pass
        ali_records = []
        for record in alignment:
            sequence = ''.join([record.seq[i] for i in keep_cols])
            seqrec = SeqRecord(Seq(sequence), id = record.id, description = '')
            ali_records.append(seqrec)
        alignment_edit = MultipleSeqAlignment(ali_records)
        return alignment_edit
    elif strategy == 'sequence':
        print('Strategy: remove sequences with ambiguous characters')
        #Get the alignment records whose sequences do not contain ambiguous characters
        ali_records = [record for record in alignment if all([not char in record.seq for char in ambig_list])]
        removed = len(alignment) - len(ali_records)
        alignment_edit = MultipleSeqAlignment(ali_records)
        print('%d sequences containing ambiguous characters removed from alignment' % removed)
        return alignment_edit
    else:
        logging.warning('invalid ambig removal strategy "%s" specified. Choose "column" or "sequence"' % strategy)
        return
#------------------------------------------------------------------------------
def unique_records(alignment, illegals):
    """
    Create a dictionary from an alignment object containing unique sequences
    and identifiers. Characters in the identifier that often break alignment
    or tree software are substituted. Characters specified in a string literal
    are replaced in the identifiers. Returns a MultipleSequenceAlignment object.
    """
    sequence_dict = {}
    id_counts = {}
    identicals = {}
    omit = 0
    illegal_chars = re.compile(illegals)
    #Create the sequence_dict, which contains unique sequences and headers
    for record in alignment:
        sequence = str(record.seq)
        id_legal = illegal_chars.sub('_', record.id)
        #Store unique sequences in the dictionary
        if not sequence in sequence_dict.keys():
            #Make sure the identifier is unique
            if id_legal in sequence_dict and id_counts[id_legal] > 1:
                id_counts[id_legal] += 1
                id_legal = id_legal + '_%d' % id_counts[id_legal]
            else:
                pass
            sequence_dict[sequence] = id_legal
        else:
            #Non-unique sequence found. Store in a dictionary of identicals
            #to display a message about duplicate sequences
            if not sequence in identicals.keys():
                identicals[sequence] = [id_legal, sequence_dict[sequence]]
            else:
                identicals[sequence].append(id_legal)
            omit += 1
    #Generate a message detailing identical sequences that were removed,
    #and which of these were kept in the alignment
    identical_sequences = [identicals[key] for key in identicals.keys() if identicals[key]]
    if identical_sequences:
        for i in identical_sequences:
            if not any([v for v in i if v in sequence_dict.values()]):
                logging.warning('Of identical sequences [%s] none were kept!!' % ','.join(i))
            else:
                for v in i:
                    if v in sequence_dict.values():
                        print('####\nsequences:\n%s\nare identical. %s is kept in the alignment\n####' % (','.join(i), v))
                    else:
                        pass

    print('total non-unique sequences omitted: %d' % omit)
    #Convert the dictionary into a MultipleSequenceAlignment object
    ali_records = []
    for sequence in sequence_dict.keys():
        sequence = str(sequence)
        seqrec = SeqRecord(Seq(sequence), id = sequence_dict[sequence], description = '')
        ali_records.append(seqrec)

    alignment_uniq = MultipleSeqAlignment(ali_records)

    return alignment_uniq
#==============================================================================
####
#MAIN
####

#Read in the alignment as a MultipleSeqAlignment
alignment = AlignIO.read(opts.input_file, opts.infmt)

if opts.infmt != opts.outfmt:
    if opts.outfmt in valid_formats:
        print("Alignment changed from '%s' to '%s' format." % (opts.infmt, opts.outfmt))
    else:
        logging.warning('Specified MSA format is not valid! Choose from: %s' % ','.join(valid_formats))

if not opts.rm_ambig and not opts.unique:
    print("No preprocessing steps specified. Writing alignment in '%s' format to '%s' format." % (opts.infmt, opts.outfmt))
    AlignIO.write(alignment, opts.output_file, opts.outfmt)
else:
    #Remove ambiguous characters if specified
    if opts.rm_ambig:
        ambig_list = list()
        if opts.seq_type == 'amino':
            ambig_list = ambig_aa
        elif opts.seq_type == 'dna' or opts.seq_type == 'rna':
            ambig_list = ambig_nuc
        else:
            logging.warning('choose a valid --seq_type: "amino", "dna", or "rna"')
            pass
        if ambig_list:
            rm_ambig(alignment = alignment, strategy = opts.ambig_strategy, ambig_list = ambig_list)

    #Get unique sequences and identifiers
    if opts.unique:
        alignment = unique_records(alignment = alignment, illegals = illegal_chars)

    AlignIO.write(alignment, opts.output_file, opts.outfmt)
