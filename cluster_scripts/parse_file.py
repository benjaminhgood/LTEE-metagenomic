import numpy
import sys
from math import fabs

base_table = {'A':'T','T':'A','G':'C','C':'G'}

codon_table = { 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'CGT': 'R', 'CGC': 'R', 'CGA':'R',
'CGG':'R', 'AGA':'R', 'AGG':'R', 'AAT':'N', 'AAC':'N', 'GAT':'D', 'GAC':'D', 'TGT':'C', 'TGC':'D', 'CAA':'Q', 'CAG':'Q', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'CAT':'H', 'CAC':'H', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'AAA':'K', 'AAG':'K', 'ATG':'M', 'TTT':'F', 'TTC':'F', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'TGG':'W', 'TAT':'Y', 'TAC':'Y', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'TAA':'!', 'TGA':'!', 'TAG':'!' }

def calculate_reverse_complement_sequence(dna_sequence):
    return "".join(base_table[base] for base in dna_sequence[::-1])
 
def calculate_codon_sequence(dna_sequence):
    return "".join(codon_table[dna_sequence[3*i:3*i+3]] for i in xrange(0,len(dna_sequence)/3))

def annotate_variant(location, allele, gene_data, repeat_data):
    
    gene_names, gene_start_positions, gene_end_positions, gene_sequences, strands = gene_data
    repeat_names, repeat_start_positions, repeat_end_positions = repeat_data
    
    #allele = allele.split('-->')[1].strip()
    #print allele
    # get gene
    gene_idxs = numpy.nonzero((location >= gene_start_positions) * (location <= gene_end_positions))[0]
    repeat_idxs = numpy.nonzero((location >= repeat_start_positions) * (location <= repeat_end_positions))[0]
    if len(gene_idxs) > 0:
        # if it is a gene
        gene_idx = gene_idxs[0]
        gene_name = gene_names[gene_idx]
        gene_start = gene_start_positions[gene_idx]
        gene_end = gene_end_positions[gene_idx]
        gene_sequence = gene_sequences[gene_idx]
        #print gene_sequence
        if allele.startswith('REL606'):
            var_type = 'jc'
        elif allele[0] == '+' or allele[0] == '-':
            indel_length = long(allele[1:].split(':')[0])
            if indel_length%3==0:
                var_type = 'non' 
            else:
                var_type = 'fs'
        else:
            gene_position = location-gene_start
            #print gene_idx,gene_start,gene_end
            
            new_gene_sequence = list(gene_sequence)
            new_gene_sequence[gene_position] = allele
            #print new_gene_sequence
            new_gene_sequence = "".join(new_gene_sequence)
             
            if strands[gene_idx] == 'reverse':
                gene_sequence = calculate_reverse_complement_sequence(gene_sequence)
                new_gene_sequence = calculate_reverse_complement_sequence(new_gene_sequence)
                
            codon_sequence = calculate_codon_sequence(gene_sequence)
            new_codon_sequence = calculate_codon_sequence(new_gene_sequence)
            if codon_sequence == new_codon_sequence:
                var_type = 'syn'
            else:
                var_type = 'non'
                
            if len(repeat_idxs) > 0:
                # a repeat sequence 
                repeat_idx = repeat_idxs[0]
                sys.stderr.write("Also a repeat sequence: %s %s\n" % (gene_name, repeat_names[repeat_idx]))
                
    else:
        var_type = ''
        if len(repeat_idxs) > 0:
            # a repeat sequence 
            repeat_idx = repeat_idxs[0]
            gene_name = repeat_names[repeat_idx]  
        else:
            gene_name = 'intergenic' 
    
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
        
    
def parse_reference_genome(filename="REL606.6.gbk"):
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
    
def print_reference_fasta(reference_sequence):
    print ">chr1"
    for i in xrange(0,len(reference_sequence),70):
        print reference_sequence[i:min([len(reference_sequence),i+70])]
    


def parse_gene_list(reference_sequence=None, filename="REL606.6.gbk"):

    if reference_sequence==None:
        reference_sequence=parse_reference_genome()

    gene_names = []
    start_positions = []
    end_positions = []
    gene_sequences = []
    strands = []

    file = open(filename,"r")
    line = file.readline()
    while line!="":
        items = line.split()
        feature = items[0]
        
        if feature == 'CDS':
            feature_location = items[1]
        
            while not line.strip().startswith('/gene'):
                line = file.readline()
            gene_name = line.split('=')[1].strip()[1:-1]
            
            while not line.strip().startswith('/translation'):
                line = file.readline()
            codon_sequences = [line.split('=')[1].strip()[1:-1]]
            
            line = file.readline()
            while not line.strip().startswith('gene'):
                codon_sequences.append(line.strip().rstrip('"'))
                line = file.readline()
            
            codon_sequence = "".join(codon_sequences)
            
            location_str = feature_location.lstrip("complement(").lstrip("join(").rstrip(")")
            location_strs = location_str.split(",")
            for location_str in location_strs:
            
                locations = [long(subitem) for subitem in location_str.split("..")]
            
                gene_start = locations[0]
                gene_end = locations[1]
                gene_sequence = reference_sequence[gene_start-1:gene_end]
                strand = 'forward'
                #print gene_start, gene_end, len(gene_sequence)
                
                if not len(gene_sequence)%3==0:
                    pass #print gene_name, gene_start, "Not a multiple of 3"
                else:            
                    if feature_location.startswith('complement'):
                        strand='reverse'
                        #gene_sequence = calculate_reverse_complement_sequence(gene_sequence)
            
                    #predicted_codon_sequence = calculate_codon_sequence(gene_sequence)
            
                    #if predicted_codon_sequence[:-1] != codon_sequence:
                    #    print "Not a good one!"
                        #print gene_sequence
                    #    print predicted_codon_sequence
                    #    print codon_sequence
                
                    start_positions.append(locations[0])
                    end_positions.append(locations[1])
                    gene_names.append(gene_name)
                    gene_sequences.append(gene_sequence)
                    strands.append(strand)
                
        line = file.readline()
    file.close()
    
    return gene_names, numpy.array(start_positions), numpy.array(end_positions), gene_sequences, strands
    
def parse_repeat_list(filename="REL606.6.gbk"):

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
        
            while not line.strip().startswith('/mobile_element'):
                line = file.readline()
            repeat_name = line.split('=')[1].strip()[1:-1]

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
        mutations.append((location, gene_name, ref_allele, alt_allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths))   
    file.close()
    
    # sort by position
    keys = [mutation[0] for mutation in mutations]
    keys, mutations = (list(t) for t in zip(*sorted(zip(keys, mutations))))
    return mutations    

def print_timecourse(mutations):

    times = mutations[0][5]
    print_strs = ['Position', 'Gene', 'Ref Allele', 'Alt Allele', 'Annotation']
    for t in zip(times):
            print_strs.append('AC:%d' % t)
            print_strs.append('DP:%d' % t)
            
    print ", ".join(print_strs)
    
    for location, gene_name, ref_allele, alt_allele, var_type, times, alts, depths in mutations:
        print_strs = [str(location),gene_name, ref_allele, alt_allele, var_type]
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

all_lines = []
for i in xrange(1,7):
    all_lines.append('m%d' % i)
    all_lines.append('p%d' % i)
 
complete_lines = ['m1', 'm4', 'm5', 'm6', 'p1', 'p2', 'p3', 'p4', 'p5']    
complete_nonmutator_lines = ['m5','m6','p1','p2','p4','p5']
complete_mutator_lines = ['m1','m4','p3']
complete_early_mutator_lines = ['m4','p3']
early_mutator_lines = ['m4','p3']

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
    return log(Nef*(100.0**(0))/Ne0)/log(Nrf*(100.0**(0))/Nr0)

def calculate_X(Ne0,Nr0,Nef,Nrf,cycles=1):
    return log( (Nef/Nrf) / (Ne0/Nr0) ) / (6.64*cycles)

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
    
def count_differences(mutation_list_1, mutation_list_2):
    unique_mutations = set()
    unique_mutations.update(mutation_list_1)
    unique_mutations.update(mutation_list_2)
    #print unique_mutations
    #print len(unique_mutations),len(mutation_list_1),len(mutation_list_2)
    return 2*len(unique_mutations)-len(mutation_list_1)-len(mutation_list_2) 

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
    
    
if __name__=='__main__':
    reference_sequence = parse_reference_genome()
    print_reference_fasta(reference_sequence)
