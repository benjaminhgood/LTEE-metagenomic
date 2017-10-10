import numpy
import sys
from math import fabs
import glob
import os


data_directory='data_files/'
figure_directory='manuscript/figures/'

default_min_depth = 5

ancestral_araA_mutation_location = 70867 #, araA, T->C @ 70867
ancestral_recD_mutation_location = 2847052 #, recD, A->G @ 2847052
        
#NUM_SYNONYMOUS_SITES = 941000.0

base_table = {'A':'T','T':'A','G':'C','C':'G'}

codon_table = { 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'CGT': 'R', 'CGC': 'R', 'CGA':'R',
'CGG':'R', 'AGA':'R', 'AGG':'R', 'AAT':'N', 'AAC':'N', 'GAT':'D', 'GAC':'D', 'TGT':'C', 'TGC':'D', 'CAA':'Q', 'CAG':'Q', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'CAT':'H', 'CAC':'H', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'AAA':'K', 'AAG':'K', 'ATG':'M', 'TTT':'F', 'TTC':'F', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'TGG':'W', 'TAT':'Y', 'TAC':'Y', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'TAA':'!', 'TGA':'!', 'TAG':'!' }

# calculate number of synonymous opportunities for each codon
codon_synonymous_opportunity_table = {}
for codon in codon_table.keys():
    codon_synonymous_opportunity_table[codon] = {}
    for i in xrange(0,3):
        codon_synonymous_opportunity_table[codon][i] = -1 # since G->G is by definition synonymous, but we don't want to count it
        codon_list = list(codon)
        for base in ['A','C','T','G']:
            codon_list[i]=base
            new_codon = "".join(codon_list)
            if codon_table[codon]==codon_table[new_codon]:
                # synonymous!
                codon_synonymous_opportunity_table[codon][i]+=1

bases = set(['A','C','T','G'])
substitutions = []
for b1 in bases:
    for b2 in bases:
        if b2==b1:
            continue
        
        substitutions.append( '%s->%s' % (b1,b2) )
        
codon_synonymous_substitution_table = {}
codon_nonsynonymous_substitution_table = {}
for codon in codon_table.keys():
    codon_synonymous_substitution_table[codon] = [[],[],[]]
    codon_nonsynonymous_substitution_table[codon] = [[],[],[]]
    
    for i in xrange(0,3):
        reference_base = codon[i]
        
        codon_list = list(codon)
        for derived_base in ['A','C','T','G']:
            if derived_base==reference_base:
                continue
            substitution = '%s->%s' % (reference_base, derived_base)
            codon_list[i]=derived_base
            new_codon = "".join(codon_list)
            if codon_table[codon]==codon_table[new_codon]:
                # synonymous!
                codon_synonymous_substitution_table[codon][i].append(substitution)
            else:
                codon_nonsynonymous_substitution_table[codon][i].append(substitution)


def parse_mask_list(filename="additional_data/REL606.L20.G15.P0.M35.mask.gd"):
    
    # Add masks calculated in Tenaillon et al (Nature, 2016)
    # Downloaded from barricklab/LTEE-Ecoli/reference/ github repository
    
    mask_start_positions = []
    mask_end_positions = []
    
    file = open(filename,"r")
    file.readline() # header
    for line in file:
        items = line.split()
        start = long(items[4])
        length = long(items[5])
        mask_start_positions.append(start)
        mask_end_positions.append(start+length-1)
    
    # Add masking of prophage elements (Methods section of Tenaillon et al (Nature, 2016))
    mask_start_positions.append(880528)
    mask_end_positions.append(904682)
    
    return numpy.array(mask_start_positions), numpy.array(mask_end_positions)

def is_repeat_masked(position, position_gene_map=None):

    if position_gene_map==None:
        position_gene_map,effective_gene_lengths, substitution_specific_synonymous_fraction = create_annotation_map()
        
    if position in position_gene_map:
        if position_gene_map[position]=='repeat':
            return True
            
    return False

def create_annotation_map(gene_data=None, repeat_data=None, mask_data=None):

    if gene_data==None:
        gene_data = parse_gene_list()
        repeat_data = parse_repeat_list()
        mask_data = parse_mask_list()
    

    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = gene_data
    repeat_names, repeat_start_positions, repeat_end_positions, repeat_complements = repeat_data
    mask_start_positions, mask_end_positions = mask_data
    
    position_gene_map = {}
    gene_position_map = {}

    num_masked_sites = 0

    # first mark things that are repeats
    # this takes precedence over all other annotations
    for start,end in zip(repeat_start_positions,repeat_end_positions):
        for position in xrange(start,end+1):
            if position not in position_gene_map:
                position_gene_map[position]='repeat'
                num_masked_sites+=1
    
    # then mark masked things
    for start,end in zip(mask_start_positions, mask_end_positions):
        for position in xrange(start,end+1):
            if position not in position_gene_map:
                position_gene_map[position]='repeat'
                num_masked_sites+=1
    
    
    # then greedily annotate genes at remaining sites
    for gene_name,start,end in zip(gene_names,gene_start_positions,gene_end_positions):
        for position in xrange(start,end+1):
            if position not in position_gene_map:
                position_gene_map[position] = gene_name
                if gene_name not in gene_position_map:
                    gene_position_map[gene_name]=[]
                gene_position_map[gene_name].append(position)
    
    # remove 'partial' genes that have < 10bp unmasked sites
    for gene_name in list(sorted(gene_position_map.keys())):
        if len(gene_position_map[gene_name]) < 10:
            for position in gene_position_map[gene_name]:
                position_gene_map[position] = 'repeat'
            del gene_position_map[gene_name]
    
    # count up number of synonymous opportunities
    effective_gene_synonymous_sites = {}
    effective_gene_nonsynonymous_sites = {}
    
    substitution_specific_synonymous_sites = {substitution: 0 for substitution in substitutions}
    substitution_specific_nonsynonymous_sites = {substitution: 0 for substitution in substitutions}
    
    for gene_name,start,end,gene_sequence,strand in zip(gene_names, gene_start_positions, gene_end_positions, gene_sequences, strands):
        
        if gene_name not in gene_position_map:
            continue
        
        if strand=='forward':
            oriented_gene_sequence = gene_sequence
        else:
            oriented_gene_sequence = calculate_reverse_complement_sequence(gene_sequence)
        
        for position in gene_position_map[gene_name]:
        
            if gene_name not in effective_gene_synonymous_sites:
                effective_gene_synonymous_sites[gene_name]=0
                effective_gene_nonsynonymous_sites[gene_name]=0
        
            if gene_name.startswith('tRNA') or gene_name.startswith('rRNA'):
                pass
            
            else:
        
                # calculate position in gene
                if strand=='forward':
                    position_in_gene = position-start
                else:
                    position_in_gene = end-position
        
                # calculate codon start
                codon_start = long(position_in_gene/3)*3
                codon = gene_sequence[codon_start:codon_start+3] 
                position_in_codon = position_in_gene%3       
            
                #print gene_name, start, end, position, codon,position_in_codon
            
                effective_gene_synonymous_sites[gene_name] += codon_synonymous_opportunity_table[codon][position_in_codon]/3.0
                effective_gene_nonsynonymous_sites[gene_name] += 1-codon_synonymous_opportunity_table[codon][position_in_codon]/3.0 
                
                for substitution in codon_synonymous_substitution_table[codon][position_in_codon]:
                    substitution_specific_synonymous_sites[substitution] += 1
                    
                for substitution in codon_nonsynonymous_substitution_table[codon][position_in_codon]:
                    substitution_specific_nonsynonymous_sites[substitution] += 1
    
    substitution_specific_synonymous_fraction = {substitution: substitution_specific_synonymous_sites[substitution]*1.0/(substitution_specific_synonymous_sites[substitution]+substitution_specific_nonsynonymous_sites[substitution]) for substitution in substitution_specific_synonymous_sites.keys()}
                    
    # then annotate promoter regions at remaining sites
    for gene_name,start,end in zip(gene_names,promoter_start_positions,promoter_end_positions):
        for position in xrange(start,end+1):
            if position not in position_gene_map:
                # position hasn't been annotated yet
                
                if gene_name not in gene_position_map:
                    # the gene itself has not been annotated
                    # so don't annotate the promoter
                    pass
                else:
                    position_gene_map[position] = gene_name
                    gene_position_map[gene_name].append(position)
    
    # calculate effective gene lengths
    effective_gene_lengths = {gene_name: len(gene_position_map[gene_name])-effective_gene_synonymous_sites[gene_name] for gene_name in gene_position_map.keys()}
    effective_gene_lengths['synonymous'] = sum([effective_gene_synonymous_sites[gene_name] for gene_name in gene_position_map.keys()])
    effective_gene_lengths['nonsynonymous'] = sum([effective_gene_nonsynonymous_sites[gene_name] for gene_name in gene_position_map.keys()])
    effective_gene_lengths['masked'] = num_masked_sites
    effective_gene_lengths['noncoding'] = calculate_genome_length()-effective_gene_lengths['synonymous']-effective_gene_lengths['nonsynonymous']-effective_gene_lengths['masked']
    
    return position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction      

def calculate_synonymous_nonsynonymous_target_sizes():
    position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction  = create_annotation_map()
    return effective_gene_lengths['synonymous'], effective_gene_lengths['nonsynonymous'], substitution_specific_synonymous_fraction   

def calculate_reverse_complement_sequence(dna_sequence):
    return "".join(base_table[base] for base in dna_sequence[::-1])
 
def calculate_codon_sequence(dna_sequence):
    return "".join(codon_table[dna_sequence[3*i:3*i+3]] for i in xrange(0,len(dna_sequence)/3))

def get_closest_gene_name(location, gene_data):
    gene_names, start_positions, end_positions, gene_sequences, strands = gene_data

    closest_start_idx = numpy.fabs(start_positions-location).argmin()
    closest_end_idx = numpy.fabs(end_positions-location).argmin()
    if fabs(start_positions[closest_start_idx]-location) < fabs(end_positions[closest_end_idx]-location):
        return gene_names[closest_start_idx]
    else:
        return gene_names[closest_end_idx]

var_types = ['synonymous','missense','nonsense','noncoding','indel','sv']

def annotate_gene(position, position_gene_map):
    
    if position in position_gene_map:
        gene_name = position_gene_map[position]
    else:
        gene_name = 'intergenic'
        
    return gene_name

def annotate_variant(position, allele, gene_data, position_gene_map):
    
    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = gene_data
    
    # get gene
    gene_name = annotate_gene(position, position_gene_map)
    
    if allele.startswith('Depth'):
        var_type = 'unknown'
    elif allele.startswith('MOB') or allele.startswith('junction'): 
        var_type = 'sv'
    elif allele.startswith('indel'):
            var_type = 'indel'
    elif allele[1:3]=='->':
        # a SNP, so annotate it
        if gene_name=='intergenic':
            var_type = 'noncoding'
        elif gene_name=='repeat':
            var_type = 'repeat'
        else:
            # must be in a real gene
            # so get it
            i = gene_names.index(gene_name)
            
            gene_start_position = gene_start_positions[i]
            gene_end_position = gene_end_positions[i]
            promoter_start_position = promoter_start_positions[i]
            promoter_end_position = promoter_end_positions[i]
            gene_sequence = gene_sequences[i]
            strand = strands[i]
            
            if position<gene_start_position or position>gene_end_position:
                #var_type='promoter'
                var_type='noncoding' # (promoter)
            else:
            
                if gene_name.startswith('tRNA') or gene_name.startswith('rRNA'):
                    var_type='noncoding'
                else:
                
                    # calculate position in gene
                    if strand=='forward':
                        position_in_gene = position-gene_start_position
                        oriented_gene_sequence = gene_sequence
                        new_base = allele[3]
                    else:
                        position_in_gene = gene_end_position-position
                        oriented_gene_sequence = calculate_reverse_complement_sequence(gene_sequence)
                        new_base = base_table[allele[3]]
            
                
                    # calculate codon start
                    codon_start = long(position_in_gene/3)*3
                    codon = oriented_gene_sequence[codon_start:codon_start+3]
                    codon_list = list(codon) 
                    position_in_codon = position_in_gene%3 
                    codon_list[position_in_codon]=new_base
                    new_codon="".join(codon_list)
                    if codon_table[codon]==codon_table[new_codon]:
                        var_type='synonymous'
                    else:
            
                        if codon_table[new_codon]=='!':
                            var_type='nonsense'
                        else:
                            var_type='missense'
    else:
        sys.stderr.write("Unknown: %s\n" % allele)
        var_type='unknown'
    
    return gene_name, var_type
        
def old_annotate_variant(location, allele, gene_data, repeat_data):
    
    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = gene_data
    repeat_names, repeat_start_positions, repeat_end_positions, repeat_complements = repeat_data
    
    if allele[0]=='+' or allele[0]=='-':
        location+=1
        # where the actual base change resides

    # get gene
    gene_idxs = numpy.nonzero((location >= gene_start_positions) * (location <= gene_end_positions))[0]
    repeat_idxs = numpy.nonzero((location >= repeat_start_positions) * (location <= repeat_end_positions))[0]
    promoter_idxs = numpy.nonzero((location >= promoter_start_positions) * (location <= promoter_end_positions))[0]
    
    #sys.stderr.write('%s\n' % str(promoter_idxs))
    #sys.stderr.write('%d %d\n' % (gene_start_positions[0], promoter_start_positions[0]))
    #sys.stderr.write('%d %d\n' % (gene_end_positions[0], promoter_end_positions[0]))
    
    if len(gene_idxs) > 0:
    
        # if it is inside a gene
    
        gene_idx = gene_idxs[0]
        gene_name = gene_names[gene_idx]
        gene_start = gene_start_positions[gene_idx]
        gene_end = gene_end_positions[gene_idx]
        gene_sequence = gene_sequences[gene_idx]
        gene_position = location-gene_start
        
        #print gene_sequence
        if allele.startswith('MOB') or allele.startswith('JC'):
            var_type = 'sv'
        
        if allele[1:3]=='->':
            # a SNP, so annotate it
        
            if gene_sequence=="":
                var_type = 'missense'
            else:    
            
                new_gene_sequence = list(gene_sequence)
                new_gene_sequence[gene_position] = allele[-1]
                new_gene_sequence = "".join(new_gene_sequence)
            
                if strands[gene_idx] == 'reverse':
                   gene_sequence = calculate_reverse_complement_sequence(gene_sequence)
                   new_gene_sequence = calculate_reverse_complement_sequence(new_gene_sequence)
                
                codon_sequence = calculate_codon_sequence(gene_sequence)
                new_codon_sequence = calculate_codon_sequence(new_gene_sequence)
                if codon_sequence == new_codon_sequence:
                    var_type = 'synonymous'
                else:
                    if new_codon_sequence.find('!') < len(new_codon_sequence)-1:
                        var_type = 'nonsense'
                        #sys.stderr.write("nonsense: %d\n" % (len(new_codon_sequence)-new_codon_sequence.find('!')-1))
                    else:
                        var_type = 'missense'
                
            
        elif allele.startswith('MOB') or allele.startswith('junction'):
            var_type = 'sv'
        else:
            var_type = 'indel'
            
        if len(repeat_idxs) > 0:
            # a repeat sequence 
            repeat_idx = repeat_idxs[0]
            #sys.stderr.write("Also a repeat sequence: %s %s\n" % (gene_name, repeat_names[repeat_idx]))
        
                
    elif len(promoter_idxs) > 0:
        # if it is inside a promoter
        promoter_idx = promoter_idxs[0]
        gene_name = gene_names[promoter_idx]
        if allele[1:3] == '->':
            var_type = 'promoter'
        elif allele.startswith('MOB') or allele.startswith('JC'):
            var_type = 'sv'
        else:
            var_type = 'indel'
        
    elif len(repeat_idxs) > 0: 
        # a repeat sequence 
        repeat_idx = repeat_idxs[0]
        gene_name = repeat_names[repeat_idx] 
        
        if allele[1:3] == '->':
            var_type = 'repeat'  
        elif allele.startswith('MOB') or allele.startswith('JC'):
            var_type = 'sv'
        else:
            var_type = 'indel' 
            
    else:
        # an intergenic mutation
            
        gene_name = 'intergenic'
        if allele[1:3]=='->':
            # SNP
            var_type = 'intergenic'
        elif allele.startswith('MOB') or allele.startswith('junction'):
            var_type = 'sv'
        else:
            var_type = 'indel'
            
    return gene_name, var_type

def get_closest_repeat_idx(location, repeat_data):
    repeat_names, start_positions, end_positions, complements = repeat_data

    closest_start_idx = numpy.fabs(start_positions-location).argmin()
    closest_end_idx = numpy.fabs(end_positions-location).argmin()
    if fabs(start_positions[closest_start_idx]-location) < fabs(end_positions[closest_end_idx]-location):
        return closest_start_idx
    else:
        return closest_end_idx

def get_repeat_idx(location, repeat_data):
    repeat_names, start_positions, end_positions, complements = repeat_data
    
    repeat_idxs = numpy.nonzero((location >= start_positions-100) * (location <= end_positions+100))[0]
    if len(repeat_idxs) > 0:
        return repeat_idxs[0]
    else:
        return -1
        

def in_repeat_region(location, repeat_data):
    
    repeat_names, start_positions, end_positions, complements = repeat_data
    
    repeat_idxs = numpy.nonzero((location >= start_positions) * (location <= end_positions))[0]
    if len(repeat_idxs) > 0:
        return True
    else:
        return False
        
    
def parse_reference_genome(filename="additional_data/REL606.6.gbk"):
    reference_sequences = []
    
    # GBK file
    if filename[-3:] == 'gbk':
        file = open(filename,"r")
        origin_reached = False
        for line in file:
            if line.startswith("ORIGIN"):
                origin_reached=True
            if origin_reached:
                items = line.split()
                if items[0].isdigit():
                    reference_sequences.extend(items[1:])    
        file.close()
    
    # FASTA file
    else:
        file = open(filename,"r")
        file.readline() # header
        for line in file:
            reference_sequences.append(line.strip())
        file.close()
    
    reference_sequence = "".join(reference_sequences).upper()
    return reference_sequence

def calculate_genome_length(reference_sequence=None):
    if reference_sequence==None:
        reference_sequence=parse_reference_genome()
    return len(reference_sequence)
    
def print_reference_fasta(reference_sequence):
    print ">chr1"
    for i in xrange(0,len(reference_sequence),70):
        print reference_sequence[i:min([len(reference_sequence),i+70])]
    
def annotate_operon(gene_name,operon_data):
    
    operon_gene_map, gene_operon_map, operon_size_map  = operon_data 
    
    if gene_name in gene_operon_map:
        return gene_operon_map[gene_name]
    else:
        return None    
    
def parse_operon_list(filename="additional_data/door_operon_list.txt"): #, gene_data=None, repeat_data=None, position_gene_map=None):

    gene_data = parse_gene_list()
    repeat_data = parse_repeat_list()
    mask_data = parse_mask_list()  
    position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction = create_annotation_map(gene_data, repeat_data, mask_data)
    gene_size_map = create_gene_size_map(effective_gene_lengths)

    
    # returns dictionary: operon name

    operons = {}
    raw_operon_lengths = {}
    
    file = open(filename,"r")
    file.readline() # header
    for line in file:
     
        if line.strip()=="":
            break
    
        items = line.split()
        operon_id = items[0].strip()
        
        start_position = long(items[3])
        end_position = long(items[4])
        
        for position in xrange(start_position,end_position+1):
            if position in position_gene_map:
                gene_name = position_gene_map[position]
                if gene_name!='repeat':
                    # add it to the list
                    if operon_id not in operons:
                        operons[operon_id] = set()
                        raw_operon_lengths[operon_id] = end_position-start_position
                    operons[operon_id].add(gene_name)
        
    # create (possibly nonunique) reverse map
    gene_operon_multimap = {}
    for operon_id in operons.keys():
        for gene_name in operons[operon_id]:
            if gene_name not in gene_operon_multimap:
                gene_operon_multimap[gene_name] = []
            gene_operon_multimap[gene_name].append(operon_id)
            
    # now choose unique operon for each gene
    for gene_name in gene_operon_multimap:
        if len(gene_operon_multimap[gene_name]) > 1:
            operon_lengths = numpy.array([raw_operon_lengths[operon_id] for operon_id in gene_operon_multimap[gene_name]])
            
            desired_idx = operon_lengths.argmax()
            
            for i in xrange(0,len(gene_operon_multimap[gene_name])):
                
                if i!=desired_idx:
                    operons[ gene_operon_multimap[gene_name][i] ].remove(gene_name)
            
            
    # everything should be unique now        
        
    # create reverse map
    gene_operon_map = {}
    operon_gene_map = {}
    for operon_id in operons.keys():
        operon_name = ";".join(operons[operon_id])
        operon_gene_map[operon_name] = operons[operon_id]
        for gene_name in operons[operon_id]:
            gene_operon_map[gene_name] = operon_name
    
    # add in genes that were left out before
    for gene_name in gene_size_map:
        if gene_name not in gene_operon_map:
            gene_operon_map[gene_name] = gene_name
            operon_gene_map[gene_name] = set([gene_name])
            
    # calculate operon sizes
    
    operon_size_map = {}
    
    for gene_name in gene_operon_map.keys():
        operon_name = gene_operon_map[gene_name]
        if operon_name not in operon_size_map:
            operon_size_map[operon_name]=0
            
        operon_size_map[operon_name] += gene_size_map[gene_name]
            
    return operon_gene_map, gene_operon_map, operon_size_map


def create_gene_size_map(effective_gene_lengths=None):

    if effective_gene_lengths==None:
        reference_sequence = parse_reference_genome()
        gene_data = parse_gene_list()
        repeat_data = parse_repeat_list()
        mask_data = parse_mask_list()
        gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = gene_data
    
        position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction = create_annotation_map(gene_data, repeat_data, mask_data)

    
    excluded_genes=set(['synonymous','nonsynonymous','noncoding','masked'])
    
    gene_size_map = {}
    for gene_name in effective_gene_lengths.keys():
        
        #if gene_name.startswith('tRNA'):
        #    print gene_name 
            
        if gene_name in excluded_genes:
            continue
            
        gene_size_map[gene_name] = effective_gene_lengths[gene_name]
        
    return gene_size_map
    

def old_create_gene_size_map(gene_data):
    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = gene_data
    
    gene_size_map = {}
    
    for i in xrange(0,len(gene_names)):
        gene_name = gene_names[i]
        if not gene_name.startswith('ins'): 
        
            if gene_name.startswith('tRNA') or gene_name.startswith('rRNA'):
                nonsynonymous_multiple = 1.0
            else:
                nonsynonymous_multiple = 3.0/4
         
            gene_size_map[gene_name] = (fabs(gene_start_positions[i]-gene_end_positions[i])*nonsynonymous_multiple + fabs(promoter_start_positions[i]-promoter_end_positions[i]))

    return gene_size_map
    

#####################################################################
#
# Reads through the Genbank file for the reference and 
# compiles a list of genes, tRNAs, etc.
#
#####################################################################        
def parse_gene_list(reference_sequence=None, filename="additional_data/REL606.6.gbk"):

    features = set(['CDS','gene','tRNA','rRNA','repeat_region'])

    if reference_sequence==None:
        reference_sequence=parse_reference_genome()

    observed_gene_names = set()

    gene_names = []
    start_positions = []
    end_positions = []
    promoter_start_positions = []
    promoter_end_positions = []
    gene_sequences = []
    strands = []

    file = open(filename,"r")
    line = file.readline()
    while line!="":
        
        items = line.split()
        feature = items[0]
        
        gene_name = ""
        feature_location = ""
        
        if feature=='CDS':
            feature_location = items[1]
        
            line = file.readline().strip()
            
            gene_name=""
            locus_name=""
            is_pseudo=False
            
            while line.split()[0] not in features:
                
                if line.startswith('/gene'):
                    gene_name = line.split('=')[1].strip()[1:-1]
                if line.startswith('/locus_tag'):
                    locus_name = line.split('=')[1].strip()[1:-1] 
                if line.startswith('/pseudo'):
                    is_pseudo=True   
        
                line = file.readline().strip()
                
            if gene_name=="":
                gene_name = locus_name
                
            if is_pseudo:
                gene_name = ""
                     
            # done here
          
        elif feature=='tRNA' or feature=='rRNA':
        
            feature_location = items[1]
            #print feature_location
        
            while not line.strip().startswith('/gene'):
                line = file.readline().strip()
            gene_name = line.split('=')[1].strip()[1:-1]
            gene_name = '%s:%s' % (feature, gene_name)
            
            
        else:
            # nothing to see here
            line = file.readline().strip()
        
        # If the element has a feature location string and a name
        # it should either be a gene, tRNA, or rRNA, so let's get details
        if feature_location!="" and gene_name!="":
        
            location_str = feature_location.lstrip("complement(").lstrip("join(").rstrip(")")
            location_strs = location_str.split(",")
            
            for location_str in location_strs:
            
                locations = [long(subitem) for subitem in location_str.split("..")]
            
                gene_start = locations[0]
                gene_end = locations[1]
                
                if feature=="CDS":
                    gene_sequence = reference_sequence[gene_start-1:gene_end]
                else:
                    gene_sequence = ""
                
                strand = 'forward'
                promoter_start = gene_start - 100 # by arbitrary definition, we treat the 100bp upstream as promoters
                promoter_end = gene_start - 1
                        
                
                if gene_sequence!="" and (not len(gene_sequence)%3==0):
                    pass
                    print gene_name, gene_start, "Not a multiple of 3"
                    
                if feature_location.startswith('complement'):
                    strand='reverse'
                    promoter_start = gene_end+1
                    promoter_end = gene_end+100
                    
                    if gene_sequence=="":
                        promoter_end = promoter_start
                
                
                # record information
                
                # first make sure gene name is unique
                i = 1
                old_gene_name = gene_name
                while gene_name in observed_gene_names:
                    i+=1
                    gene_name = "%s_%d" % (old_gene_name,i)
                
                start_positions.append(gene_start)
                end_positions.append(gene_end)
                promoter_start_positions.append(promoter_start)
                promoter_end_positions.append(promoter_end)
                gene_names.append(gene_name)
                gene_sequences.append(gene_sequence)
                strands.append(strand)
                observed_gene_names.add(gene_name)
        
    file.close()
    
    # sort genes based on start position
    
    gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = (list(x) for x in zip(*sorted(zip(gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands), key=lambda pair: pair[1])))
    
    return gene_names, numpy.array(start_positions), numpy.array(end_positions), numpy.array(promoter_start_positions), numpy.array(promoter_end_positions), gene_sequences, strands
    
def parse_repeat_list(filename="additional_data/REL606.6.gbk"):

    repeat_names = []
    start_positions = []
    end_positions = []
    complements = []
 
    file = open(filename,"r")
    line = file.readline()
    while line!="":
        items = line.split()
        feature = items[0]
        
        if feature == 'repeat_region':
            feature_location = items[1]
        
            # Get name of mobile element
            repeat_name = 'unknown'
            
            line = file.readline()
            while line.strip()[0]=='/':
                if line.strip().startswith('/mobile_element'):
                    repeat_name = line.split('=')[1].strip()[1:-1]
                line = file.readline()
            
            # Finished at next non '/' entry, presumably next feature

            if feature_location.startswith('complement'):
                complement = True
            else:
                complement = False            
            
            location_str = feature_location.lstrip("complement(").lstrip("join(").rstrip(")")
            location_strs = location_str.split(",")
            for location_str in location_strs:
            
                locations = [long(subitem) for subitem in location_str.split("..")]
                start_positions.append(locations[0])
                end_positions.append(locations[1])
                repeat_names.append(repeat_name)
                complements.append(complement)
        
        else:
                
            line = file.readline()
    file.close()
    
    return repeat_names, numpy.array(start_positions), numpy.array(end_positions), complements
    
def parse_timecourse(filename):

    mutations = []

    file = open(filename,"r")
    header_line = file.readline()
    items = header_line.strip().split(",")
    
    times = []
    for i in xrange(10,len(items),2):
        times.append(long(items[i].split(":")[1]))
    times = numpy.array(times) 
    
    for line in file:
        items = line.strip().split(",")
        location = long(items[0])
        gene_name = items[1].strip()
        ref_allele = items[2].strip()
        alt_allele = items[3].strip()
        var_type = items[4].strip()
        test_statistic = float(items[5])
        pvalue = float(items[6])
        cutoff_idx = long(items[7])
        depth_fold_change = float(items[8])
        depth_change_pvalue = float(items[9])
        
        alts = []
        depths = []
        for i in xrange(10,len(items),2):
            alts.append(long(float(items[i])))
            depths.append(long(float(items[i+1])))
            
        alts = numpy.array(alts)
        depths = numpy.array(depths)
        
        pop_times = times[times<1000000]
        pop_alts = alts[times<1000000]
        pop_depths = depths[times<1000000]
        
        clone_times = times[times>1000000]-1000000
        clone_alts = alts[times>1000000]
        clone_depths = depths[times>1000000]
        
        mutations.append((location, gene_name, ref_allele, alt_allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, pop_times, pop_alts, pop_depths, clone_times, clone_alts, clone_depths)) 
          
    file.close()
    
    # sort by position
    keys = [mutation[0] for mutation in mutations]
    keys, mutations = (list(t) for t in zip(*sorted(zip(keys, mutations))))
    return mutations    

def parse_annotated_timecourse(population, only_passed=True, min_coverage=5):

    mutations = []

    timecourse_filename = data_directory+("%s_annotated_timecourse.txt" % population)

    file = open(timecourse_filename, "r")
    
    header_line = file.readline()
    items = header_line.strip().split(",")
    
    times = []
    for i in xrange(13,len(items),2):
        times.append(long(items[i].split(":")[1]))
    times = numpy.array(times) 
    
    # depth line
    depth_line = file.readline()
    items = depth_line.strip().split(",")
    avg_depths = []
    for i in xrange(13,len(items),2):
        avg_depths.append(float(items[i+1]))
    avg_depths = numpy.array(avg_depths) 
    
    population_avg_depth_times = times[times<1000000]
    population_avg_depths = avg_depths[times<1000000]
    clone_avg_depth_times = times[times>1000000]-1000000
    clone_avg_depths = avg_depths[times>1000000]
    
    for line in file:
        items = line.strip().split(",")
        location = long(items[0])
        gene_name = items[1].strip()
        allele = items[2].strip()
        var_type = items[3].strip()
        test_statistic = float(items[4])
        pvalue = float(items[5])
        cutoff_idx = long(items[6])
        depth_fold_change = float(items[7])
        depth_change_pvalue = float(items[8])
        
        duplication_idx = long(items[9])
        fold_increase = float(items[10])
        duplication_pvalue = float(items[11])
        
        passed_str = items[12]
        if passed_str.strip()=='PASS':
            passed = True
        else:
            passed = False
        
        alts = []
        depths = []
        
        for i in xrange(13,len(items),2):
            alts.append(long(float(items[i])))
            depths.append(long(float(items[i+1])))
            
        alts = numpy.array(alts)
        depths = numpy.array(depths)
        
        # zero out timepoints with individual coverage lower than some threshold
        alts *= (depths>=min_coverage)*(avg_depths>=min_coverage)
        depths *= (depths>=min_coverage)*(avg_depths>=min_coverage)
        
        pop_times = times[(times<1000000)]
        pop_alts = alts[(times<1000000)]
        pop_depths = depths[(times<1000000)]
        
        clone_times = times[(times>1000000)]-1000000
        clone_alts = alts[(times>1000000)]
        clone_depths = depths[(times>1000000)]
        
        if passed or (not only_passed):
            mutations.append((location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, pop_times, pop_alts, pop_depths, clone_times, clone_alts, clone_depths)) 
          
    file.close()
    
    # sort by position
    keys = [mutation[0] for mutation in mutations]
    keys, mutations = (list(t) for t in zip(*sorted(zip(keys, mutations))))
    return mutations, (population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths)  

clade_hmm_states = {'A':0,'E':1,'FB':2,'FM':3, 'Fm':4,'PB':5,'PM':6,'Pm':7,'PB*':8}
well_mixed_hmm_states = {'A':0,'E':1,'F':2,'P':3}    

UNBORN = clade_hmm_states['A']
EXTINCT= clade_hmm_states['E']
ANCESTRAL_FIXED = clade_hmm_states['FB']
MINOR_FIXED=clade_hmm_states['Fm']
MAJOR_FIXED=clade_hmm_states['FM']
ANCESTRAL_POLYMORPHIC=clade_hmm_states['PB']
MINOR_POLYMORPHIC=clade_hmm_states['Pm']
MAJOR_POLYMORPHIC=clade_hmm_states['PM']

clade_extinct_states = set([clade_hmm_states['A'],clade_hmm_states['E']])
clade_fixed_states = set([clade_hmm_states['FB'], clade_hmm_states['FM'], clade_hmm_states['Fm'], clade_hmm_states['PB*']])

clade_majority_states = set([clade_hmm_states['FB'], clade_hmm_states['FM'], clade_hmm_states['PB'],clade_hmm_states['PM'],clade_hmm_states['PB*']])

clade_polymorphic_states = set([clade_hmm_states['PB'], clade_hmm_states['PM'], clade_hmm_states['Pm']])

well_mixed_extinct_states = set([well_mixed_hmm_states['A'], well_mixed_hmm_states['E']])
well_mixed_polymorphic_states = set([well_mixed_hmm_states['P']])
well_mixed_fixed_states = set([well_mixed_hmm_states['F']])

FIXED = well_mixed_hmm_states['F']
POLYMORPHIC = well_mixed_hmm_states['P']


def parse_haplotype_timecourse(population):

    haplotype_filename = data_directory+('%s_haplotype_timecourse.txt' % population)
 
    file = open(haplotype_filename,"r")
     
    times = numpy.array([float(item) for item in file.readline().split(",")])
    fmajors = numpy.array([float(item) for item in file.readline().split(",")])
    fminors = numpy.array([float(item) for item in file.readline().split(",")])
    file.readline()
    file.readline()
    file.readline()
    file.readline()
    file.readline()
    file.readline()
    haplotypes = []
    for line in file:
        Ls = numpy.array([float(item) for item in line.split(",")])
        haplotypes.append(Ls)
    file.close()    
    return times,fmajors,fminors,haplotypes
    
def parse_well_mixed_state_timecourse(population):
 
    haplotype_filename = data_directory+('%s_well_mixed_state_timecourse.txt' % population)
 
    file = open(haplotype_filename,"r")
     
    times = numpy.array([float(item) for item in file.readline().split(",")])
    num_unborn = numpy.array([float(item) for item in file.readline().split(",")])
    num_extinct = numpy.array([float(item) for item in file.readline().split(",")])
    num_fixed = numpy.array([float(item) for item in file.readline().split(",")])
    num_polymorphic = numpy.array([float(item) for item in file.readline().split(",")])
    
    states = []
    for line in file:
        Ls = numpy.array([float(item) for item in line.split(",")])
        states.append(Ls)
    file.close()    
    return times, states


def print_timecourse(mutations):

    print_strs = ['Position', 'Gene', 'Ref Allele', 'Alt Allele', 'Annotation', 'Test statistic', 'P-value', 'Cutoff index', 'Fold change', 'P-value']
    times = mutations[0][10]
    for t in zip(times):
            print_strs.append('AC:%d' % t)
            print_strs.append('DP:%d' % t)
            
    print ", ".join(print_strs)
    
    for location, gene_name, ref_allele, alt_allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths in mutations:
    
        print_strs = [str(location),gene_name, ref_allele, alt_allele, var_type, str(test_statistic), str(pvalue), str(cutoff_idx), str(depth_fold_change), str(depth_change_pvalue)]
        for alt,depth in zip(alts,depths):
            print_strs.append(str(alt))
            print_strs.append(str(depth))
             
        print ", ".join(print_strs) 

def parse_coverage(coverage_filename):
    coverage_file = open(coverage_filename)
    ts = numpy.array([long(item) for item in coverage_file.readline().split()])
    ds = numpy.array([float(item) for item in coverage_file.readline().split()])
    return ts,ds
    
import numpy
from scipy.optimize import newton
from math import log,exp,fabs
from bz2 import BZ2File
import sys


all_lines = ['m5','p2','p4','p1','m6','p5','m1','m2','m3','m4','p3','p6']

all_old_lines = ['%s_old' % population for population in all_lines]
 
all_line_colors = ['#5DA5DA', '#FAA43A', '#60BD68', '#B276B2', '#F15854', '#4D4D4D']*2
 
 
#nonmutator_line_colors = ['#d0d1e6', '#a6bddb', '#67a9cf', '#3690c0', '#02818a', '#016450']
# or reversed (blue to light)
#nonmutator_line_colors = ['#016450', '#02818a', '#3690c0', '#67a9cf', '#a6bddb', '#d0d1e6']
#mutator_line_colors = ['#fdd49e', '#fdbb84', '#fc8d59', '#ef6548', '#d7301f', '#990000']
 
nonmutator_line_colors = ['#4A1486','#807DBA','#084594','#4292C6','#005A32','#41AB5D']
 
mutator_line_colors = ['#8C2D04','#CC4C02','#B10026','#E31A16','#FC4E2A','#FD8D3C']

line_color_map = {'m5': '#4A1486', 'p2': '#807DBA', 'p4': '#084594','p1': '#4292C6', 'p5': '#005A32', 'm6': '#41AB5D', 'm1': '#8C2D04', 'm2': '#CC4C02', 'm3': '#B10026', 'm4': '#E31A16', 'p3': '#FC4E2A', 'p6': '#FD8D3C'}

def get_line_color(population):
    return line_color_map[population]

nonmutator_group_color = get_line_color('p4')
mutator_group_color = get_line_color('m3')

    
 
#4D4D4D (gray)
#5DA5DA (blue)
#FAA43A (orange)
#60BD68 (green)
#F17CB0 (pink)
#B2912F (brown)
#B276B2 (purple)
#DECF3F (yellow)
#F15854 (red)
 
complete_lines = ['m1', 'm4', 'm5', 'm6', 'p1', 'p2', 'p3', 'p4', 'p5']    
complete_nonmutator_lines = ['m5','m6','p1','p2','p4','p5']
complete_mutator_lines = ['m1','m4','p3']
complete_early_mutator_lines = ['m4','p3']
early_mutator_lines = ['m4','p3']
mutator_lines = ['m1','m2','m3','m4','p3','p6']

idx_10k = 15
idx_40k = 35
idx_50k = 40

mut_idx_10k = 3
mut_idx_40k = 7

def get_pretty_name(line):
    if line[0]=='m':
        return 'Ara-%s' % line[1]
    else:
        return 'Ara+%s' % line[1]


def calculate_W(Ne0,Nr0,Nef,Nrf,cycles=1):
    return log((Nef)*(100.0**(0))/(Ne0))/log((Nrf)*(100.0**(0))/(Nr0+1))

def calculate_X(Ne0,Nr0,Nef,Nrf,cycles=1):
    return log( ((Nef)/(Nrf)) / ((Ne0)/(Nr0)) ) / (6.64*cycles)

def w2x(w):
    x = newton(lambda x: x-(w-1)*log(100/(1+0.5*(exp(x*6.64)-1)))/6.64, (w-1)*log(100)/6.64)
    return x

def parse_ancestor_fitnesses(filename="Concatenated.LTEE.data.all.csv"):
    line_data = {line: [] for line in all_lines}
    file = open(filename,"r")
    file.readline() # skip headers
    for line in file:
        #print line
        items = line.split(",")
        t = long(items[0])
        
        subitems = items[3].split()
        idx=long(subitems[2])
        if subitems[1]=='-':
            line_idx='m%d' % idx
        else:
            line_idx='p%d' % idx
        
        if items[5] == '' or items[8]=='':
            continue
        
        
        R0c = float(items[5])
        W0c = float(items[6])
        D0 = float(items[7])
        Rfc = float(items[8])
        Wfc = float(items[9])
        Df = float(items[10])
        theirW = float(items[11])
        
        R0 = 1.0*R0c*D0
        W0 = 1.0*W0c*D0
        Rf = 1.0*Rfc*Df
        Wf = 1.0*Wfc*Df
        
        if items[1][0]=='A':
            # evolved clone is red pop
            Ne0 = R0
            Nef = Rf
            Nr0 = W0
            Nrf = Wf
        else:
            Ne0 = W0
            Nef = Wf
            Nr0 = R0
            Nrf = Rf
        
        W = calculate_W(Ne0,Nr0,Nef,Nrf)
            
        if fabs(W-theirW) > 1e-02:
            # should never happen!
            print "Fitnesses don't match up!"
            
        line_data[line_idx].append((t,Ne0,Nr0,Nef,Nrf))        
            
    file.close()
    
    # process data
    
    new_line_data = {line: [{},{}] for line in all_lines}
    for line in all_lines:
        for t,Ne0,Nr0,Nef,Nrf in line_data[line]:
            W = calculate_W(Ne0,Nr0,Nef,Nrf)
            X = calculate_X(Ne0,Nr0,Nef,Nrf)
            if not t in new_line_data[line][0]:
                new_line_data[line][0][t] = []
                new_line_data[line][1][t] = []
            new_line_data[line][0][t].append(X)
            new_line_data[line][1][t].append(W)
            
    trajectories = {line: None for line in all_lines}
    for line in all_lines:
        ts = new_line_data[line][0].keys()
        ts.sort()
        xs = numpy.zeros_like(ts)*1.0
        error_xs = numpy.zeros_like(ts)*1.0
        ws = numpy.zeros_like(ts)*1.0
        for i in xrange(0,len(ts)):
            xs[i] = numpy.array(new_line_data[line][0][ts[i]]).mean()
            error_xs[i] = new_line_data[line][0][ts[i]][0] - xs[i]
            ws[i] = numpy.array(new_line_data[line][1][ts[i]]).mean()
        trajectories[line] = [ts,xs,ws,error_xs]
        
    return trajectories, line_data

def parse_40k_50k_fitnesses(filename="additional_data/All.Hands.Concatenated.Data.csv"):
    line_data = {line: [] for line in complete_lines}
    
    file = open(filename,"r")
    file.readline() # skip headers
    
    for line in file:
        items = line.split(",")
        flask = long(items[0])
        red_pop_str = items[1].strip()
        white_pop_str = items[2].strip()
        t = long(items[3])
        dt = 3
        
        if items[5]=='' or items[6]=='' or items[7]=='' or items[8]=='':
            continue
            
        R0 = float(items[5])
        W0 = float(items[6])
        Rf = float(items[7])
        Wf = float(items[8])

        if red_pop_str.startswith('Ara'):
            # evolved clone is red pop
            Ne0 = R0
            Nef = Rf
            Nr0 = W0
            Nrf = Wf
            pop_str = red_pop_str
        else:
            Ne0 = W0
            Nef = Wf
            Nr0 = R0
            Nrf = Rf
            pop_str = white_pop_str
        
        if Ne0 < 10 or Nr0 < 10:
            continue
        
        W = calculate_W(Ne0,Nr0,Nef*(100**3),Nrf*(100**3))
        
        subitems = pop_str.split()
        idx = long(subitems[2])
        if subitems[1]=='-':
            line_idx='m%d' % idx
        else:
            line_idx='p%d' % idx
        
        line_data[line_idx].append((t,Ne0,Nr0,Nef,Nrf))
    
    file.close()
    
    # process data
    
    new_line_data = {line: [{},{}] for line in complete_lines}
    for line in complete_lines:
        for t,Ne0,Nr0,Nef,Nrf in line_data[line]:
            W = calculate_W(Ne0,Nr0,Nef*(100**3),Nrf*(100**3))
            X = calculate_X(Ne0,Nr0,Nef,Nrf,cycles=3)
            if not t in new_line_data[line][0]:
                new_line_data[line][0][t] = []
                new_line_data[line][1][t] = []
            new_line_data[line][0][t].append(X)
            new_line_data[line][1][t].append(W)

    trajectories = {line: None for line in complete_lines}
    residuals = []
    for line in complete_lines:
        ts = new_line_data[line][0].keys()
        ts.sort()
        xs = numpy.zeros_like(ts)*1.0
        std_xs = numpy.zeros_like(ts)*1.0
        ws = numpy.zeros_like(ts)*1.0
        for i in xrange(0,len(ts)):
            x_array = numpy.array(new_line_data[line][0][ts[i]]) 
            xs[i] = x_array.mean()
            std_xs[i] = x_array.std()/(len(x_array)**0.5)
            #print len(x_array)
        trajectories[line] = [ts,xs,std_xs]
        
    return trajectories

      

def parse_40k_fitnesses(filename="Concatenated.40k.v.50k.csv"):
    line_data = {line: [] for line in complete_lines}
    file = open(filename,"r")
    file.readline() # skip headers
    for line in file:
        #print line
        items = line.split(",")
        t = long(items[4])
        
        subitems = items[2].split()
        idx=long(subitems[2])
        if subitems[1]=='-':
            line_idx='m%d' % idx
        else:
            line_idx='p%d' % idx
        
        if items[5] == '' or items[8]=='':
            continue
        
        
        R0c = float(items[6])
        W0c = float(items[7])
        D0 = float(items[8])
        Rfc = float(items[9])
        Wfc = float(items[10])
        Df = float(items[11])
        theirW = float(items[12])
        
        R0 = 1.0*R0c*D0
        W0 = 1.0*W0c*D0
        Rf = 1.0*Rfc*Df
        Wf = 1.0*Wfc*Df
        
        if items[0][0]=='A':
            # evolved clone is red pop
            Ne0 = R0
            Nef = Rf
            Nr0 = W0
            Nrf = Wf
        else:
            Ne0 = W0
            Nef = Wf
            Nr0 = R0
            Nrf = Rf
        
        W = calculate_W(Ne0,Nr0,Nef,Nrf,cycles=3)
            
        if fabs(W-theirW) > 1e-02:
            # should never happen!
            print "Fitnesses don't match up!"
            
        line_data[line_idx].append((t,Ne0,Nr0,Nef,Nrf))        
            
    file.close()
    
    # process data
    
    new_line_data = {line: [{},{}] for line in complete_lines}
    for line in complete_lines:
        for t,Ne0,Nr0,Nef,Nrf in line_data[line]:
            W = calculate_W(Ne0,Nr0,Nef,Nrf,cycles=3)
            X = calculate_X(Ne0,Nr0,Nef,Nrf,cycles=3)
            if not t in new_line_data[line][0]:
                new_line_data[line][0][t] = []
                new_line_data[line][1][t] = []
            new_line_data[line][0][t].append(X)
            new_line_data[line][1][t].append(W)
            
    trajectories = {line: None for line in complete_lines}
    residuals = []
    for line in complete_lines:
        ts = new_line_data[line][0].keys()
        ts.sort()
        xs = numpy.zeros_like(ts)*1.0
        std_xs = numpy.zeros_like(ts)*1.0
        ws = numpy.zeros_like(ts)*1.0
        for i in xrange(0,len(ts)):
            x_array = numpy.array(new_line_data[line][0][ts[i]]) 
            xs[i] = x_array.mean()
            std_xs[i] = x_array.std()/(len(x_array)**0.5)
        trajectories[line] = [ts,xs,std_xs]
        
    return trajectories

def parse_simulation_trajectory(line):
    items = line.split(":")
    #sys.stderr.write(line)
    params = [float(subitem) for subitem in items[0].split(",")]
    xs = numpy.array([float(subitem) for subitem in items[1].split(",")])
    ps = numpy.array([float(subitem) for subitem in items[2].split(",")])
    
    if items[3].count(",") > 0:
        dxs = numpy.array([float(subitem) for subitem in items[3].split(",")])
    else:
        dxs = [0]
    #sys.stderr.write("%s \n" % items[4])
    if items[4].count(";") > 0:
        s_trajectories = []
        for subitem in items[4].split(";")[:-1]:
            s_trajectory = []
            if subitem.count(",") > 0:
                s_trajectory = [float(subsubitem) for subsubitem in subitem.split(",")]
            else:
                s_trajectory = []
            if len(s_trajectory) < 5:
                s_trajectory.extend([0]*(5-len(s_trajectory)))
            s_trajectories.append(numpy.array(s_trajectory)) 
    else:
        s_trajectories = [numpy.array([0,0,0,0,0])]
    
    if items[5].count(";") > 0:
        x_trajectories = numpy.array([[float(subsubitem) for subsubitem in subitem.split(",")] for subitem in items[5].split(";")[:-1]]) 
    else:
        print "No x trajectories!"
        x_trajectories = numpy.array([xs])
    
    if items[6].count(";") > 0:
        p_trajectories = numpy.array([[float(subsubitem) for subsubitem in subitem.split(",")] for subitem in items[6].split(";")[:-1]])
    else:
        print "No p trajectories!"
        p_trajectories = [ps]
    
    trajectory_data = [params,xs,ps,dxs,s_trajectories,x_trajectories,p_trajectories]
    for i in xrange(7,len(items)):
        trajectory_data.append([float(subitem) for subitem in items[i].split(",")])
    return trajectory_data
    
def parse_simulation_trajectories_bzip(filename):
    trajectories = []
    file = BZ2File(filename, "r")
    for line in file:
        #print line
        trajectories.append(parse_simulation_trajectory(line))
    return trajectories

def parse_simulation_trajectories(filename):
    trajectories = []
    file = open(filename, "r")
    for line in file:
        try:
            items = line.split(":")
            params = [float(subitem) for subitem in items[0].split(",")]
            xs = numpy.array([float(subitem) for subitem in items[1].split(",")])
            ps = numpy.array([float(subitem) for subitem in items[2].split(",")])
            
            if items[3].count(",") > 0:
                dxs = numpy.array([float(subitem) for subitem in items[3].split(",")])
            else:
                dxs = [0]
            if items[4].count(";") > 0:
                s_trajectories = [numpy.array([float(subsubitem) for subsubitem in subitem.split(",")]) for subitem in items[4].split(";")[:-1]] 
            else:
                s_trajectories = [numpy.array([0,0,0,0,0])]
            if items[5].count(";") > 0:
                x_trajectories = numpy.array([[float(subsubitem) for subsubitem in subitem.split(",")] for subitem in items[5].split(";")[:-1]]) 
            else:
                print "No x trajectories!"
                x_trajectories = numpy.array([xs])
            trajectories.append([params,xs,ps,dxs,s_trajectories,x_trajectories,times])
        except ValueError, e:
            pass
    file.close()
    return trajectories

def parse_simulation_timecourse_trajectories(filename):
    trajectories = []
    file = open(filename, "r")
    params = [float(subitem) for subitem in file.readline().split()]
    s = float(params[0])
    Ub = float(params[1])
    times = numpy.array([float(item) for item in file.readline().split()])
    avg_xs = numpy.zeros_like(times)
    avg_ps = numpy.zeros_like(times)
    x_trajectories = []
    p_trajectories = []
    s_trajectories = []
    for line in file:
        items = line.split(":")
        xs = numpy.array([float(subitem) for subitem in items[0].split()])
        ps = numpy.array([float(subitem) for subitem in items[1].split()])
        ss = numpy.array([float(subitem.split(",")[0]) for subitem in items[2].split()])
        avg_xs += xs
        avg_ps += ps
        x_trajectories.append(xs)
        p_trajectories.append(ps)
        s_trajectories.append(ss)
        
    file.close()
    return [[params,xs,ps,['0'],s_trajectories,numpy.array(x_trajectories),times]]

def get_time_indices(smaller_times, bigger_times):
    indices = []
    for t in smaller_times:
        i = numpy.fabs(bigger_times-t).argmin()
        indices.append(i)
    return numpy.array(indices)
    
def parse_barrick_data():
    mutations = {}
    barrick_mutation_times = [0, 2000,5000,10000,15000,20000]
    barrick_mutation_counts = [0, 6, 15, 28, 36, 45]
    for t,p in zip(barrick_mutation_times,barrick_mutation_counts):
        mutations[('m1',t)] = [[str(i) for i in xrange(0,p)]]
    return mutations

def parse_wielgoss_data(directory='wielgoss_genomes'):

    mutation_types = {'SNP': True, 'SUB': True, 'DEL': True, 'INS': True, 'MOB': True, 'AMP': True, 'CON': True, 'INV': True, 'RA': False, 'MC': False, 'JC': False, 'UN': False, 'TSEQ': False, 'PFLP': False, 'RFLP': False, 'PFGE': False, 'PHYL': False, 'CURA': False}

    old_filenames = ['A+1_40K_REL11008.gd', 'A+4_40K_REL10957.gd', 'A-1_20K_REL8593C.gd',	'A-5_40K_REL10948.gd', 'A+1_40K_REL11009.gd', 'A+5_40K_REL10982.gd', 'A-3_30K_ZDB16.gd', 'A-5_40K_REL10949.gd', 'A+2_40K_REL10950.gd', 'A+5_40K_REL10983.gd', 'A-3_30K_ZDB357.gd',	'A-6_40K_REL11005.gd', 'A+2_40K_REL10951.gd', 'A-1_20K_REL8593A.gd', 'A-3_40K_REL10988.gd', 'A-6_40K_REL11006.gd', 'A+4_40K_REL10956.gd', 'A-1_20K_REL8593B.gd', 'A-5_40K_REL10947.gd'] #has a citrate guy
    
    filenames = ['A-1_20K_REL8593C.gd', 'A-1_20K_REL8593A.gd', 'A-1_20K_REL8593B.gd', 'A-3_30K_ZDB16.gd', 'A-3_30K_ZDB357.gd', 'A-5_40K_REL10948.gd', 'A-5_40K_REL10949.gd',  'A-5_40K_REL10947.gd', 'A-6_40K_REL11006.gd', 'A-6_40K_REL11005.gd', 'A+1_40K_REL11008.gd', 'A+1_40K_REL11009.gd', 'A+2_40K_REL10950.gd', 'A+2_40K_REL10951.gd', 'A+4_40K_REL10957.gd', 'A+4_40K_REL10956.gd', 'A+5_40K_REL10982.gd', 'A+5_40K_REL10983.gd'] # w/o citrate sweep
    
    mutations = {}
    
    for filename in filenames:  
        
        if filename[1]=='-':
            population = 'm%s' % filename[2]
        else:
            population = 'p%s' % filename[2]
        
        if True:
            generation = long(filename[4:6])*1000
            
            if not (population,generation) in mutations:    
                mutations[(population, generation)] = []
                #print population, generation, filename
                
                
            clone_mutations = []

            file = open('%s/%s' % (directory, filename),"r")
            file.readline()
            for line in file:
                items = line.split()
                mutation_type = items[0]
                                  
                if mutation_types[mutation_type]:
                    for i in xrange(1,len(items)): 
                        if items[i] == 'REL606':
                            mutation_str = "%s %s" % (mutation_type, " ".join(items[i+1:]))
                            #print mutation_str
                    clone_mutations.append(mutation_str)
                else:
                    break
                    
            file.close()
            mutations[(population,generation)].append(clone_mutations)
    
    return mutations

def parse_tenaillon_data(directory='additional_data/tenaillon_genomes/LTEE-clone-curated'):

    mutation_types = {'SNP': True, 'SUB': True, 'DEL': True, 'INS': True, 'MOB': True, 'AMP': True, 'CON': True, 'INV': True, 'RA': False, 'MC': False, 'JC': False, 'UN': False, 'TSEQ': False, 'PFLP': False, 'RFLP': False, 'PFGE': False, 'PHYL': False, 'CURA': False}

    filenames = next(os.walk(directory))[2]
    mutations = {}
    
    for filename in filenames:  
        
        items = filename[:-3].split("_")
        population_str = items[0]
        generation = long(items[1][:-3])
        clone_str = items[2]
        
        
        if population_str.startswith('Anc'):
            if population_str[3]=='+':
                populations = ['p%d' % i for i in xrange(1,7)]
            else:
                populations = ['m%d' % i for i in xrange(1,7)]    
        else:
            if population_str[3]=='-':
                populations = ['m%s' % population_str[4]]
            else:
                populations = ['p%s' % population_str[4]]
        
        
        
        for population in populations:
            
            if not (population,generation) in mutations:    
                mutations[(population, generation)] = []
                #print population, generation, filename
                
                
            clone_mutations = []

            file = open('%s/%s' % (directory, filename),"r")
            line = file.readline()
            while line.startswith('#'):
                line=file.readline()
            while line!="":
                items = line.split()
                mutation_type = items[0]
                                  
                if mutation_types[mutation_type]:
                    for i in xrange(1,len(items)): 
                        if items[i] == 'REL606':
                            mutation_str = "%s %s" % (mutation_type, " ".join(items[i+1:]))
                            #print mutation_str
                    clone_mutations.append(mutation_str)
                else:
                    break
                
                line=file.readline()    
            file.close()
            mutations[(population,generation)].append(clone_mutations)
    
    return mutations
    
def count_differences(mutation_list_1, mutation_list_2):
    unique_mutations = set()
    unique_mutations.update(mutation_list_1)
    unique_mutations.update(mutation_list_2)
    #print unique_mutations
    #print len(unique_mutations),len(mutation_list_1),len(mutation_list_2)
    return 2*len(unique_mutations)-len(mutation_list_1)-len(mutation_list_2) 

def get_tenaillon_mutation_trajectories():
    tenaillon_data = parse_tenaillon_data()

    mutation_trajectories = {}

    final_mutation_trajectories = {}

    for population,t in tenaillon_data.keys():
        if population not in mutation_trajectories:
            mutation_trajectories[population] = {}
        if not t in mutation_trajectories[population]:
            mutation_trajectories[population][t] = []
        for i in xrange(0,len(tenaillon_data[(population,t)])):
            clone_mutations = tenaillon_data[(population,t)][i]
            mutation_trajectories[population][t].append( len(clone_mutations)*1.0 )
            #print len(clone_mutations)*1.0
    
    for population in mutation_trajectories.keys():
    
        
        ts = numpy.array(sorted(mutation_trajectories[population].keys()))
        avg_muts = numpy.array([numpy.array(mutation_trajectories[population][t]).mean() for t in ts])
        avg_muts -= avg_muts[0]    
        stderr_muts = numpy.array([numpy.array(mutation_trajectories[population][t]).std()/(1.0*len(mutation_trajectories[population][t]))**0.5 for t in ts])
        final_mutation_trajectories[population] = (ts,avg_muts,stderr_muts)
    
        #print population, ts
    
    return final_mutation_trajectories

def get_tenaillon_mutation_trajectory(masked_populations=[]):

    final_mutation_trajectories = get_tenaillon_mutation_trajectories()

    masked_populations = set(masked_populations)
    
    avg_ts = final_mutation_trajectories[final_mutation_trajectories.keys()[0]][0]
    avg_muts = numpy.zeros_like(avg_ts)*1.0
    avg_mut2s = numpy.zeros_like(avg_ts)*1.0
    num_samples = 0
    for population in complete_nonmutator_lines:
        if population not in masked_populations:
            avg_muts+=final_mutation_trajectories[population][1]
            avg_mut2s+=numpy.square(final_mutation_trajectories[population][1])
            num_samples += 1
            
    avg_muts /= num_samples
    avg_mut2s /= num_samples
    
    stderr_muts = numpy.sqrt((avg_mut2s-avg_muts*avg_muts)/num_samples)
    
    return avg_ts, avg_muts, stderr_muts

def get_tenaillon_fixed_mutation_trajectories():
    tenaillon_data = parse_tenaillon_data()

    mutation_trajectories = {}

    final_mutation_trajectories = {}

    for population,t in tenaillon_data.keys():
        if population not in mutation_trajectories:
            mutation_trajectories[population] = {}
        if not t in mutation_trajectories[population]:
            mutation_trajectories[population][t] = []
            
        fixed_mutations = None
        for i in xrange(0,len(tenaillon_data[(population,t)])):
            clone_mutations = set(tenaillon_data[(population,t)][i])
            if fixed_mutations == None:
                fixed_mutations = clone_mutations
            
            fixed_mutations &= clone_mutations
            
            
        mutation_trajectories[population][t].append( len(fixed_mutations)*1.0 )
            #print len(clone_mutations)*1.0
    
    for population in mutation_trajectories.keys():
    
        
        ts = numpy.array(sorted(mutation_trajectories[population].keys()))
        avg_muts = numpy.array([numpy.array(mutation_trajectories[population][t]).mean() for t in ts])
        avg_muts -= avg_muts[0]    
        stderr_muts = numpy.array([numpy.array(mutation_trajectories[population][t]).std()/(1.0*len(mutation_trajectories[population][t]))**0.5 for t in ts])
        final_mutation_trajectories[population] = (ts,avg_muts,stderr_muts)
    
        #print population, ts
    
    return final_mutation_trajectories

        
def get_mutation_trajectory(allowed_lines = all_lines): 
    muts =  {}
    wielgoss_data = parse_wielgoss_data()
    barrick_data = parse_barrick_data()

    for population,t in wielgoss_data.keys():
        if population in allowed_lines:
            if not t in muts:
                muts[t] = []
            for i in xrange(0,len(wielgoss_data[(population,t)])):
                clone_mutations = wielgoss_data[(population,t)][i]
                muts[t].append(len(clone_mutations))

    for population,t in barrick_data.keys():
        if population in allowed_lines:
            if not t in muts:
                muts[t] = []
            for i in xrange(0,len(barrick_data[(population,t)])):
                clone_mutations = barrick_data[(population,t)][i]
                muts[t].append(len(clone_mutations))

    mut_ts = numpy.array(sorted(muts.keys()))
    avg_muts = numpy.array([numpy.array(muts[t]).mean() for t in mut_ts])
    stderr_muts = numpy.array([numpy.array(muts[t]).std()/(1.0*len(muts[t]))**0.5 for t in mut_ts])
    return mut_ts, avg_muts, stderr_muts        

def get_fitness_trajectory():
    trajectories, line_data = parse_ancestor_fitnesses()

    n = len(complete_nonmutator_lines)

    avg_ts = numpy.zeros(41)
    avg_xs = numpy.zeros(41)
    avg_x2s = numpy.zeros(41)

    for line in complete_nonmutator_lines:
        ts = numpy.array(trajectories[line][0])*1.0
        xs = numpy.array(trajectories[line][1])
        avg_ts += ts
        avg_xs += xs
        avg_x2s += xs*xs
    
    avg_ts /= n
    avg_xs /= n
    avg_x2s /= n
    std_xs = numpy.sqrt((avg_x2s-avg_xs*avg_xs))*n/(n-1)
    return avg_ts, avg_xs, std_xs

def get_40k_gains():
    trajectories = parse_40k_fitnesses()

    ss = []
    ds2s = []
    for line,i in zip(complete_nonmutator_lines,range(1,1+len(complete_nonmutator_lines))):
        x40 = trajectories[line][1][0]
        dx40 = trajectories[line][2][0]
        x50 = trajectories[line][1][1]
        s = (x50-x40)
        ss.append(s)
        dx50 = trajectories[line][2][1]
        ds2 = (dx50**2+dx40**2)
        ds2s.append(ds2)
        ds = ds2**(0.5)
        
    
    avg_s = sum(ss)/len(ss)
    std_s = (sum(ds2s)/len(ds2s)**2)**(0.5) 
    return 40000,50000, avg_s, std_s

def get_individual_40k_gains():
    trajectories = parse_40k_fitnesses()

    individual_dxs = []
    individual_std_dxs = []
    
    for line,i in zip(complete_nonmutator_lines,range(1,1+len(complete_nonmutator_lines))):
        x40 = trajectories[line][1][0]
        dx40 = trajectories[line][2][0]
        x50 = trajectories[line][1][1]
        dx = (x50-x40)
        individual_dxs.append(dx)
        dx50 = trajectories[line][2][1]
        std_dx = (dx50**2+dx40**2)**0.5
        individual_std_dxs.append(std_dx)
        
        
    individual_dxs = numpy.array(individual_dxs)
    individual_std_dxs = numpy.array(individual_std_dxs)
    
    avg_dx = individual_dxs.mean()
    std_dx = individual_dxs.std()/(1.0*len(individual_dxs))**(0.5)
    std_dx_measurement = (numpy.square(individual_std_dxs).mean()/len(individual_std_dxs))**0.5
    
    return individual_dxs, individual_std_dxs, avg_dx, std_dx, std_dx_measurement



def get_all_40k_gains():
    trajectories = parse_40k_fitnesses()

    individual_dxs = []
    individual_std_dxs = []
    
    for line,i in zip(complete_lines,range(1,1+len(complete_lines))):
        x40 = trajectories[line][1][0]
        dx40 = trajectories[line][2][0]
        x50 = trajectories[line][1][1]
        dx = (x50-x40)
        individual_dxs.append(dx)
        dx50 = trajectories[line][2][1]
        std_dx = (dx50**2+dx40**2)**0.5
        individual_std_dxs.append(std_dx)
        
        
    individual_dxs = numpy.array(individual_dxs)
    individual_std_dxs = numpy.array(individual_std_dxs)
    
    avg_dx = individual_dxs.mean()
    std_dx = individual_dxs.std()/(1.0*len(individual_dxs))**(0.5)
    std_dx_measurement = (numpy.square(individual_std_dxs).mean()/len(individual_std_dxs))**0.5
    
    return individual_dxs, individual_std_dxs, avg_dx, std_dx, std_dx_measurement
    
def get_all_50k_fitnesses():
    trajectories = parse_40k_fitnesses()

    individual_xs = []
    individual_std_xs = []
    
    for line,i in zip(complete_lines,range(1,1+len(complete_lines))):
        x50 = trajectories[line][1][1]
        individual_xs.append(x50)
        dx50 = trajectories[line][2][1]
        std_x = dx50
        individual_std_xs.append(std_x)
        
        
    individual_xs = numpy.array(individual_xs)
    individual_std_xs = numpy.array(individual_std_xs)
    
    avg_x = individual_xs.mean()
    std_x = individual_xs.std()/(1.0*len(individual_xs))**(0.5)
    std_x_measurement = (numpy.square(individual_std_xs).mean()/len(individual_std_xs))**0.5
    
    return individual_xs, individual_std_xs, avg_x, std_x, std_x_measurement
 

def get_all_40k_50k_gains():
    trajectories = parse_40k_50k_fitnesses()

    individual_40k_dxs = []
    individual_40k_std_dxs = []
    individual_50k_dxs = []
    individual_50k_std_dxs = []
    
    
    for line,i in zip(complete_lines,range(1,1+len(complete_lines))):
        x40 = trajectories[line][1][0]
        std_x40 = trajectories[line][2][0]
        x50 = trajectories[line][1][1]
        std_x50 = trajectories[line][2][1]
        x60 = trajectories[line][1][2]
        std_x60 = trajectories[line][2][2]
        
        dx_40k = x50-x40
        dx_50k = x60-x50
        
        individual_40k_dxs.append(dx_40k)
        individual_50k_dxs.append(dx_50k)
        
        std_40k_dx = (std_x50**2+std_x40**2)**0.5
        std_50k_dx = (std_x60**2+std_x50**2)**0.5
        
        individual_40k_std_dxs.append(std_40k_dx)
        individual_50k_std_dxs.append(std_50k_dx)
        
    individual_40k_dxs = numpy.array(individual_40k_dxs)
    individual_40k_std_dxs = numpy.array(individual_40k_std_dxs)
    individual_50k_dxs = numpy.array(individual_50k_dxs)
    individual_50k_std_dxs = numpy.array(individual_50k_std_dxs)
    
    avg_40k_dx = individual_40k_dxs.mean()
    std_40k_dx = individual_40k_dxs.std()/(1.0*len(individual_40k_dxs))**(0.5)
    std_40k_dx_measurement = (numpy.square(individual_40k_std_dxs).mean()/len(individual_40k_std_dxs))**0.5
    
    avg_50k_dx = individual_50k_dxs.mean()
    std_50k_dx = individual_50k_dxs.std()/(1.0*len(individual_50k_dxs))**(0.5)
    std_50k_dx_measurement = (numpy.square(individual_50k_std_dxs).mean()/len(individual_50k_std_dxs))**0.5
    
    
    return individual_40k_dxs, individual_40k_std_dxs, avg_40k_dx, std_40k_dx, std_40k_dx_measurement, individual_50k_dxs, individual_50k_std_dxs, avg_50k_dx, std_50k_dx, std_50k_dx_measurement

def get_all_40k_50k_fitnesses():
    trajectories = parse_40k_50k_fitnesses()

    individual_40k_dxs = []
    individual_40k_std_dxs = []
    individual_50k_dxs = []
    individual_50k_std_dxs = []
    
    
    for line,i in zip(complete_lines,range(1,1+len(complete_lines))):
        x40 = trajectories[line][1][0]
        std_x40 = trajectories[line][2][0]
        x50 = trajectories[line][1][1]
        std_x50 = trajectories[line][2][1]
        x60 = trajectories[line][1][2]
        std_x60 = trajectories[line][2][2]
        
        dx_40k = x50-x40
        dx_50k = x60-x40
        
        individual_40k_dxs.append(dx_40k)
        individual_50k_dxs.append(dx_50k)
        
        std_40k_dx = (std_x50**2+std_x40**2)**0.5
        std_50k_dx = (std_x60**2+std_x40**2)**0.5
        
        individual_40k_std_dxs.append(std_40k_dx)
        individual_50k_std_dxs.append(std_50k_dx)
        
    individual_40k_dxs = numpy.array(individual_40k_dxs)
    individual_40k_std_dxs = numpy.array(individual_40k_std_dxs)
    individual_50k_dxs = numpy.array(individual_50k_dxs)
    individual_50k_std_dxs = numpy.array(individual_50k_std_dxs)
    
    avg_40k_dx = individual_40k_dxs.mean()
    std_40k_dx = individual_40k_dxs.std()/(1.0*len(individual_40k_dxs))**(0.5)
    std_40k_dx_measurement = (numpy.square(individual_40k_std_dxs).mean()/len(individual_40k_std_dxs))**0.5
    
    avg_50k_dx = individual_50k_dxs.mean()
    std_50k_dx = individual_50k_dxs.std()/(1.0*len(individual_50k_dxs))**(0.5)
    std_50k_dx_measurement = (numpy.square(individual_50k_std_dxs).mean()/len(individual_50k_std_dxs))**0.5
    
    
    return complete_lines,individual_40k_dxs, individual_40k_std_dxs, avg_40k_dx, std_40k_dx, std_40k_dx_measurement, individual_50k_dxs, individual_50k_std_dxs, avg_50k_dx, std_50k_dx, std_50k_dx_measurement

        
def get_mutator_fitness_trajectory():
    trajectories, line_data = parse_ancestor_fitnesses()

    n = len(complete_early_mutator_lines)

    avg_ts = numpy.zeros(41)
    avg_xs = numpy.zeros(41)
    avg_x2s = numpy.zeros(41)

    for line in complete_early_mutator_lines:
        ts = numpy.array(trajectories[line][0])*1.0
        xs = numpy.array(trajectories[line][1])
        avg_ts += ts
        avg_xs += xs
        avg_x2s += xs*xs
    
    avg_ts /= n
    avg_xs /= n
    avg_x2s /= n
    std_xs = numpy.sqrt((avg_x2s-avg_xs*avg_xs))*n/(n-1)
    return avg_ts, avg_xs, std_xs

def get_wiser_fitness_trajectory():
    trajectories, line_data = parse_ancestor_fitnesses()

    n = len(complete_nonmutator_lines)

    avg_ts = numpy.zeros(41)
    avg_xs = numpy.zeros(41)
    avg_x2s = numpy.zeros(41)

    for line in complete_nonmutator_lines:
        ts = numpy.array(trajectories[line][0])*1.0
        xs = numpy.log(numpy.array(trajectories[line][2]))
        avg_ts += ts
        avg_xs += xs
        avg_x2s += xs*xs
    
    avg_ts /= n
    avg_xs /= n
    avg_x2s /= n
    std_xs = numpy.sqrt((avg_x2s-avg_xs*avg_xs))*n/(n-1)
    return avg_ts, avg_xs, std_xs

 
def parse_mutation_data(): 
    muts =  {}
    wielgoss_data = parse_wielgoss_data()
    barrick_data = parse_barrick_data()

    for population,t in wielgoss_data.keys():
        if not t in muts:
            muts[t] = []
        for i in xrange(0,len(wielgoss_data[(population,t)])):
            clone_mutations = wielgoss_data[(population,t)][i]
            muts[t].append(len(clone_mutations))

    for population,t in barrick_data.keys():
        if not t in muts:
            muts[t] = []
        for i in xrange(0,len(barrick_data[(population,t)])):
            clone_mutations = barrick_data[(population,t)][i]
            muts[t].append(len(clone_mutations))

    mut_ts = numpy.array(sorted(muts.keys()))
    avg_muts = numpy.array([numpy.array(muts[t]).mean() for t in mut_ts])
    stderr_muts = numpy.array([numpy.array(muts[t]).std()/(1.0*len(muts[t]))**0.5 for t in mut_ts])
    return mut_ts, avg_muts, stderr_muts        
            
def parse_khan_fitnesses(filename="Khan_competition_data.txt"):
    fitnesses = {}

    file=open(filename,"r")
    line = file.readline()
    line = file.readline()
    line = file.readline()
    line = file.readline()
    line = file.readline()
    line = file.readline()
    line = file.readline()
    line = file.readline()
    line = file.readline()
    line = file.readline()
    line = file.readline()
    line = file.readline()
    line = file.readline()
    line = file.readline()
    line = file.readline()
    line = file.readline()
    
    for line in file:
        items = line.split()
        genotype = items[1]
        if genotype=="anc":
            genotype=""
            
        if not genotype in fitnesses:
            fitnesses[genotype] = []
        Ne0 = float(items[2])
        Nr0 = float(items[3])
        Nef = float(items[4])
        Nrf = float(items[5])
        X = calculate_X(Ne0,Nr0,Nef,Nrf,cycles=2)
        fitnesses[genotype].append(X)
        
    file.close()
    
    genotypes = fitnesses.keys()
    avg_fitnesses = numpy.array([numpy.array(fitnesses[genotype]).mean() for genotype in genotypes])
    std_fitnesses = numpy.array([numpy.array(fitnesses[genotype]).std()/(1.0*len(fitnesses[genotype]))**0.5 for genotype in genotypes])
    
    return genotypes, avg_fitnesses, std_fitnesses
 
wielgoss_ts = numpy.array([643.716, 1122.255, 1626.009, 2205.128, 5112.622, 10149.313, 15144.071, 20140.533, 25108.938, 30096.615, 35020.539, 39963.730])
wielgoss_ws = numpy.array([ 1.116, 1.169, 1.259, 1.268, 1.348, 1.476, 1.542, 1.551, 1.658, 1.657, 1.745, 1.632])
wielgoss_xs = numpy.array([w2x(w) for w in wielgoss_ws]) 
 
 
def parse_parallel_genes(filename):
    
    parallel_genes_file = open(filename,"r")
    line = parallel_genes_file.readline() # header
    parallel_genes = []
    for line in parallel_genes_file:
        items = line.split(",")
        parallel_genes.append( items[0].strip() )
        
    return parallel_genes
    
def parse_convergence_matrix(filename):
    
    convergence_matrix = {}
    
    convergence_matrix_file = open(filename,"r")
    # Header line
    line = convergence_matrix_file.readline()
    populations = [item.strip() for item in line.split(",")[2:]]
    
    for line in convergence_matrix_file:
        
        items = line.split(",")
        gene_name = items[0].strip()
        length = float(items[1])
        
        convergence_matrix[gene_name] = {'length':length, 'mutations': {population: [] for population in populations}}
        
        for population, item in zip(populations,items[2:]):
            if item.strip()=="":
                continue
                
            subitems = item.split(";")
            for subitem in subitems:
                subsubitems = subitem.split(":")
                mutation = (float(subsubitems[0]), long(subsubitems[1]), long(subsubitems[2]), float(subsubitems[3]))
                convergence_matrix[gene_name]['mutations'][population].append(mutation)
                

    return convergence_matrix
 
def old_parse_convergence_matrix(filename):
    
    convergence_matrix = []
    populations = []
    
    convergence_file = open(filename,"r")
    # header
    line = convergence_file.readline()
    gene_names = [item.strip() for item in line.split(",")[1:]]
    
    for line in convergence_file:    
        items = line.split(",")
        population = items[0]
        m = [float(item) for item in items[1:]]
        
        populations.append(population)
        convergence_matrix.append(m)
        
    convergence_matrix = numpy.array(convergence_matrix)
    return convergence_matrix, populations, gene_names


def parse_tenaillon2012_mutations(filename="additional_data/TenaillonEtal12_mutations.csv"):
    gene_data = parse_gene_list()
    repeat_data = parse_repeat_list()
    
    tenaillon_mutations = {}
    
    file = open(filename,"r")
    file.readline() # header
    for line in file:
        items = line.split(",")
        population = items[0].strip()
        location = long(items[1])
        tenaillon_var_type = items[2].strip()
        allele = items[3].strip()
        
        if population=='line108': # this was the mutator
            # skip
            continue
            
        if population not in tenaillon_mutations:
            tenaillon_mutations[population] = []

    
        # annotate mutation
        # first convert allele to our format
        if tenaillon_var_type == "deletion":
            allele = "indel;-%d" % long(allele.split("_")[0])
        elif tenaillon_var_type == "duplication":
            allele = "indel;+%d" % long(allele.split("_")[0])   
        elif tenaillon_var_type == 'insertion':
            allele = "indel;+%d" % long(allele.split("_")[0])  
        elif tenaillon_var_type == 'is_insertion':
            allele = "MOB_%s" % allele
        elif tenaillon_var_type == 'large_deletion':
            allele = "indel;-%d" % long(allele.split("_")[0])
        elif tenaillon_var_type == 'point_mutation':
            allele = allele
        else:
            # should never happen
            allele = 'indel;?'
        
        # then annotate the same way we do for our own data
        gene_name, var_type = annotate_variant(location, allele, gene_data, repeat_data)
    
        if not gene_name.startswith('ins'):
            # record! 
            
            tenaillon_mutations[population].append((location, gene_name, allele, var_type))
    
    return tenaillon_mutations
       
if __name__=='__main__':
    reference_sequence = parse_reference_genome()
    gene_data = parse_gene_list()
    repeat_data = parse_repeat_list()
    mask_data = parse_mask_list()
    gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = gene_data
    
    position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction = create_annotation_map(gene_data, repeat_data, mask_data)
    
    print "Total:", len(reference_sequence)
    print "Masked:", effective_gene_lengths['masked']
    print "Synonymous sites:", effective_gene_lengths['synonymous']
    print "Nonsynonymous sites:", effective_gene_lengths['nonsynonymous']
    print "Noncoding sites:", effective_gene_lengths['noncoding']
    
    print "Nonsynonymous:synonymous ratio:", effective_gene_lengths['nonsynonymous']/effective_gene_lengths['synonymous']
    print "Noncoding:synonymous ratio:", effective_gene_lengths['noncoding']/effective_gene_lengths['synonymous']
    print len(gene_names), "genes"
         
    operon_data = parse_operon_list(gene_data=gene_data,repeat_data=repeat_data,position_gene_map=position_gene_map)
    
    if 'arcA' in gene_names:
        idx = gene_names.index('arcA')
        #print start_positions[idx], end_positions[idx],strands[idx]
    else:
        #print 'arcA not in list'
        pass
    #mutations = parse_tenaillon_data()
    #print mutations
    
    for gene_name in gene_names:
        if gene_name.startswith('mut'):
            print gene_name, annotate_operon(gene_name,operon_data)
    