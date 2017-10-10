import parse_file


# list of BioSample Attributes (columns) downloaded from the NCBI submission protal (Package Microbe; version 1.0)
# https://submit.ncbi.nlm.nih.gov/biosample/template/?package=Microbe.1.0&action=definition
column_header_str = "*sample_name	sample_title	bioproject_accession	*organism	strain	isolate	host	isolation_source	*collection_date	*geo_loc_name	*sample_type	altitude	biomaterial_provider	collected_by	culture_collection	depth	env_biome	genotype	host_tissue_sampled	identified_by	lab_host	lat_lon	mating_type	passage_history	samp_size	serotype	serovar	specimen_voucher	temp	description"

column_headers = column_header_str.split()

# Add extra column describing whether sample has been excluded from our analysis in the quality control step" 
column_headers.append("is_excluded")
column_headers.append("replicate")

sample_name_idx = 0
# Will be filled in below

bioproject_accession_idx = 2
bioproject_accession_str = "PRJNA380528"

organism_idx = 3
organsim_str = "Escherichia coli B str. REL606"

strain_idx = 4
# Will be filled in below

host_idx = 6
host_str = "not applicable"

isolation_source_idx = 7
isolation_source_str = "not applicable"

collection_date_idx = 8
collection_date_str = "not applicable"

geo_loc_idx = 9
geo_loc_str = "not applicable"

sample_type_idx = 10
sample_type_str = "mixed culture"

description_idx = len(column_headers)-3
# Will be filled in below
# Ara-6 30000 generation clone <- example description

flagged_idx = len(column_headers)-2
# Will be filled in below

replicate_idx = len(column_headers)-1

biosample_file_idx = 1
output_file = open("biosample_submission_%d.tsv" % biosample_file_idx, "w")
output_file.write("\t".join(column_headers))
  
population_file = open("population_samples.csv","r")
population_file.readline() # header
biosample_idx = 0
biosample_max_idx = 498

seen_sample_names = set([])

for line in population_file:

    biosample_idx += 1

    if biosample_idx >= biosample_max_idx:
        output_file.close()
        biosample_file_idx+=1
        output_file = open("biosample_submission_%d.tsv" % biosample_file_idx, "w")
        output_file.write("\t".join(column_headers))
        biosample_idx = 1

    items = line.split(",")
   
    strain_str = items[0].strip()
    
    populations = [subitem.strip() for subitem in items[1].split(";")]
    
    timepoint = float(items[2])
    timepoint_str = "generation %d" % (long(timepoint*1000))
    
    flagged = ('1' in items[5])
    
    
    sample_name = items[6].strip()
    flagged = ('1' in items[5])
    if flagged:
        flagged_str = "true"
    else:
        flagged_str = ""
    
    sample_name_str = strain_str
    replicate_str = "biological replicate 1"
    
    if sample_name_str in seen_sample_names:
    
        new_sample_name_str = sample_name_str
        idx = 1
        while new_sample_name_str in seen_sample_names:
            idx += 1
            new_sample_name_str = sample_name_str + (".%d" % idx)
        
        sample_name_str = new_sample_name_str
        replicate_str = "biological replicate %d" % idx
    
    
    
    population_str = ",".join([parse_file.get_pretty_name(population) for population in populations])
    if "BG" in strain_str:
        description_str = "%s clonal isolate" % (population_str)
    else:
        description_str = "%s %s mixed population sample" % (population_str, timepoint_str)
    
    
    column_strs = [" " for column in column_headers]
        
    column_strs[sample_name_idx] = sample_name_str
    
    column_strs[bioproject_accession_idx] = bioproject_accession_str 
    
    column_strs[organism_idx] = organsim_str 
    
    column_strs[strain_idx] = strain_str

    column_strs[host_idx] = host_str 
    column_strs[isolation_source_idx] = isolation_source_str 
    column_strs[collection_date_idx] = collection_date_str 
    column_strs[geo_loc_idx] = geo_loc_str 
    column_strs[sample_type_idx] = sample_type_str
    
    column_strs[description_idx] = description_str

    column_strs[flagged_idx] = flagged_str
    column_strs[replicate_idx] = replicate_str

    output_file.write("\n")
    output_file.write("\t".join(column_strs))
    
    seen_sample_names.add(sample_name_str)
