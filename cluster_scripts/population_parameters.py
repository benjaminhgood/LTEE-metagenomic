import numpy
import sys
from math import log10,exp

population_directory = 'data/breseq_output'

populations = ['m1','m2','m3','m4','m5','m6','p1','p2','p3','p4','p5','p6']

sample_names = {}
sample_times = {}

population_filenames = set()
clone_accessions = set()

#for filename in ['desai_population_list.txt']: #,'barrick_population_list.txt']:
#    population_file = open(filename,'r')
#    for line in population_file:
#        items = line.split()
#        sample_name = items[0]
#        populations = items[1].split(",")
#        timepoint = long(items[2])
#        for population in populations:
#            if population not in sample_names:
#                sample_names[population] = []
#                sample_times[population] = []

#            sample_names[population].append(sample_name)
#            sample_times[population].append(timepoint)

#    population_file.close()

population_file = open("population_samples.csv","r")
population_file.readline() # header
for line in population_file:
    items = line.split(",")
    sample_name = items[6].strip()
    flagged = ('1' in items[5])
    population_filenames.add(sample_name)
    if not flagged:
        populations = [subitem.strip() for subitem in items[1].split(";")]
        timepoint = float(items[2])
        for population in populations:
            if population not in sample_names:
                sample_names[population] = []
                sample_times[population] = []
            sample_names[population].append(sample_name)
            sample_times[population].append(timepoint)
population_file.close()

clone_file = open("clone_samples.csv","r")
clone_file.readline() # header
for line in clone_file:
    items = line.split(",")
    sample_name = items[4].strip()
    flagged = ('1' in items[3])
    clone_accessions.add(sample_name)
    if not flagged:
        populations = [subitem.strip() for subitem in items[1].split(";")]
        timepoint = 1000+float(items[2])
        for population in populations:
            if population not in sample_names:
                sample_names[population] = []
                sample_times[population] = []
            if len(sample_times[population])>0 and sample_times[population][-1]==timepoint:
                timepoint += 0.1
            sample_names[population].append(sample_name)
            sample_times[population].append(timepoint)
clone_file.close()

population_filenames = list(population_filenames)
population_filenames.sort()

clone_accessions = list(clone_accessions)
clone_accessions.sort()

# sort by timepoint
for population in sample_names.keys():
    sample_times[population], sample_names[population] = (list(x) for x in zip(*sorted(zip(sample_times[population],sample_names[population]))))

# straggler file
straggler_file = open('population_list_stragglers.txt','r')
unsorted_straggler_list = []
for line in straggler_file:
    if line[0]=='#':
        continue
    unsorted_straggler_list.append(line.strip())
straggler_file.close()

def print_parameters():
   for sample in sample_list:
      print sample

def parameters_to_string(parameters):
   return "_".join([str(param) for param in parameters])

def print_split_parameters(param_string):
    for item in param_string.split("_"):
        print item,


if __name__ == '__main__':
    if True:
        stragglers = False
        non_stragglers = False
        test = False
        just_clones = False
        clones = True
        if sys.argv[1] == "samples":
            populations = [item.strip() for item in sys.argv[2].split(",")]
            if len(sys.argv) > 3:
                if sys.argv[3] == 'stragglers':
                    stragglers = True
                if sys.argv[3] == 'non_stragglers':
                    non_stragglers = True
                if sys.argv[3] == 'test':
                    test = True
                if sys.argv[3] == 'clones':
                    just_clones = True
                if sys.argv[3] == 'non_clones':
                    clones=False
                if sys.argv[3] == 'non_straggler_non_clones':
                    non_stragglers = True
                    clones = False
 
            for population in populations:
                if population=='clones':
                    desired_samples = clone_accessions
                else:
                    desired_samples = sample_names[population]

                for sample_name in desired_samples:
                    if stragglers:
                        if sample_name in unsorted_straggler_list:
                            if clones:    
                                print sample_name
                            else:
                                if sample_name not in clone_accessions:
                                    print sample_name

                    elif non_stragglers == True:
                        if sample_name not in unsorted_straggler_list:
                            if clones:    
                                print sample_name
                            else:
                                if sample_name not in clone_accessions:
                                    print sample_name
                    elif test:
                        if sample_name == desired_samples[0]:
                            print sample_name
                    
                    elif just_clones:
                        if sample_name in clone_accessions:
                            print sample_name
                    else:
                        if clones:
                            print sample_name
                        else:
                            if sample_name not in clone_accessions:
                                print sample_name

        elif sys.argv[1] == 'times':
            population = sys.argv[2]
            for t in sample_times[population]:
                print t,
            print ""
        elif sys.argv[1] == 'directory':
            print population_directory
        elif sys.argv[1] == "runs":
            print num_runs
        elif sys.argv[1] == "split":
            print_split_parameters(sys.argv[2])  
        elif sys.argv[1] == "clone_accessions":
            for accession in clone_accessions:
                print accession
        elif sys.argv[1] == "clone_accession_index":
            idx = long(sys.argv[2])-1
            print clone_accessions[idx]
        elif sys.argv[1] == "population_filenames":
            for filename in population_filenames:
                print filename
        else:
	    print "Usage error!"
   
