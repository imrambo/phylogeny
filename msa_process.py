#!/usr/bin/python3
import argparse
import re
import logging
import os
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
parser.add_argument('--rm_ambigs', default = 'column', action = 'store',
                    dest = 'rm_ambig', nargs = '?', type = str,
                    help = '''Remove columns or sequences with ambiguous characters.
                    Choose "column" or "sequence". Default is "column".''')
# parser.add_argument('--ambig_strategy', action = 'store', dest = 'ambig_strategy',
# nargs = '?', type = str, help = 'Use along with rm_ambigs. Remove columns or sequences. Choose "column" or "sequence".')
parser.add_argument('--unique', default = False, action = 'store_true',
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
nargs = '?', help = 'path to logfile')
# parser.add_argument('--quiet', default = False, action = 'store_true',
# dest = 'quiet', help = 'suppress status messages' )
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
def rm_ambig(alignment, strategy, ambig_list):
    """
    Remove either columns or entire sequences in an alignment
    object that contain ambiguous characters specified in a list.
    Returns a MultipleSequenceAlignment object.
    """
    if strategy == 'column':
        strat_msg = 'Strategy: remove columns with ambiguous characters'
        #if not quiet:
        logging.info(strat_msg)
        keep_cols = []
        rm_cols = []
        for i in range(alignment.get_alignment_length()):
            #Check if any ambiguous characters in the columns.
            ambig_pos = [(alignment[pos].id, char, i) for pos, char in enumerate(str(alignment[:, i])) if char in ambig_list]
            #if any([ambig in alignment[:, i] for ambig in ambig_list]):
            if ambig_pos:
                ambig_msg = ['Sequence %s: ambiguous character "%s" found at position %d' % (a, b, c) for a, b, c in ambig_pos]
                logging.info(ambig_msg[0])
                #the position of the ambigs in the columns. which indices?
                #ambig_pos = [(alignment[pos].id, char, i) for pos, char in enumerate(str(alignment[:, i])) if char in ambig_list]
                rm_cols.append(i)
            else:
                keep_cols.append(i)
        if rm_cols:
            rmcol_msg = 'ambiguous symbols found in %d column(s): %s' % (len(rm_cols), ','.join([str(x) for x in rm_cols]))
            logging.info(rmcol_msg)
            # if not quiet:
            #     logging.info(rmcol_msg)
            # if logfile:
            #     logging.info(rmcol_msg)
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
        logging.info('Strategy: remove sequences with ambiguous characters')
        #Get the alignment records whose sequences do not contain ambiguous characters
        removed_records = []
        #ali_records = [alignment[i,:] else removed_records.append(alignment[i,:]) for i in range(len(alignment)) if all([not char in alignment[i,:].seq for char in ambig_list])]
        #ali_records = [record for record in alignment if all([not char in record.seq for char in ambig_list])]
        #ali_records = [record if all([not char in record.seq for char in ambig_list]) else removed_records.append(record) for record in alignment]
        ali_records = []
        for record in alignment:
            if all([not char in record.seq for char in ambig_list]):
                ali_records.append(record)
            else:
                removed_records.append(record)


        removed = len(alignment) - len(ali_records)
        alignment_edit = MultipleSeqAlignment(ali_records)
        logging.info('%d sequences containing ambiguous characters removed from alignment' % removed)
        logging.info('removed sequences %s' % ','.join([record.id for record in removed_records]))

        return alignment_edit
    else:
        inval_strat_msg = 'invalid ambig removal strategy "%s" specified. Choose "column" or "sequence"' % strategy
        logging.warning(inval_strat_msg)
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
                        #print('####\nsequences:\n%s\nare identical. %s is kept in the alignment\n####' % (','.join(i), v))
                        logging.info('sequences:\n%s\nare identical. %s is kept in the alignment' % (','.join(i), v))

                    else:
                        pass

    #print('total non-unique sequences omitted: %d' % omit)
    logging.info('total non-unique sequences omitted: %d' % omit)
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

#Open the logfile
if opts.logfile:
    #open(opts.logfile, 'w')
    logging.basicConfig(filename=opts.logfile,
                        filemode='w',
                        format='%(asctime)s,%(msecs)d %(name)s - %(levelname)s - %(message)s',
                        datefmt='%H:%M:%S',
                        level=logging.DEBUG)
else:
    pass

#Read in the alignment as a MultipleSeqAlignment
alignment = AlignIO.read(os.path.abspath(opts.input_file), opts.infmt)
logging.info('input alignment %s, nrow = %d, ncol = %d' % (os.path.abspath(opts.input_file), len(alignment[:]), alignment.get_alignment_length()))

#Specify that the alignment format is to be changed
if opts.infmt != opts.outfmt:
    if opts.outfmt in valid_formats:
        logging.info("Alignment changed from '%s' to '%s' format." % (opts.infmt, opts.outfmt))
    else:
        logging.error('Specified MSA format is not valid! Choose from: %s' % ','.join(valid_formats))

#If not removing columns/rows with ambiguous characters or non-unique sequences
if opts.infmt and opts.outfmt and not opts.rm_ambig and not opts.unique:
    logging.info("No preprocessing steps specified. Converting alignment in '%s' format to '%s' format." % (opts.infmt, opts.outfmt))
    logging.info('output alignment %s, nrow = %d, ncol = %d' % (os.path.abspath(opts.input_file), len(alignment[:]), alignment.get_alignment_length()))
    AlignIO.write(alignment, os.path.abspath(opts.output_file), opts.outfmt)
else:
    #Get unique sequences and identifiers
    if opts.unique:
        alignment = unique_records(alignment = alignment, illegals = illegal_chars)
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
            alignment = rm_ambig(alignment = alignment, strategy = opts.rm_ambig, ambig_list = ambig_list)
#------------------------------------------------------------------------------
    #Write the final alignment object
    final_msg = 'output alignment %s, nrow = %d, ncol = %d' % (os.path.abspath(opts.input_file), len(alignment[:]), alignment.get_alignment_length())
    logging.info(final_msg)
    AlignIO.write(alignment, os.path.abspath(opts.output_file), opts.outfmt)
