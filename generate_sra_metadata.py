import parse_file

# List of SRA metadata attributes (columns) downloaded from NCBI submission portal (ftp://ftp-trace.ncbi.nlm.nih.gov/sra/metadata_table/SRA_metadata_acc.tsv)
column_header_str = "biosample_accession	bioproject_accession	title	library_ID	design_description	library_strategy	library_source	library_selection	library_layout	platform	instrument_model	filetype	filename1	assembly"		
column_headers = column_header_str.split()
column_headers.append( "filename2" )
biosample_accession_idx = 0

bioproject_accession_idx = 1
bioproject_accession_str = "PRJNA380528"

title_idx = 2

library_id_idx = 3 # what to use here?

design_description_idx = 4 # will fill in with pop / generation

library_strategy_idx = 5
library_strategy_str = "WGS"

library_source_idx = 6
library_source_str = "GENOMIC" # could have also picked metagenomic, but not sure if that meant across species...

library_selection_idx = 7
library_selection_str = "RANDOM"

library_layout_idx = 8
library_layout_str = "Paired"

platform_idx = 9
platform_str = "ILLUMINA"

instrument_model_idx = 10
instrument_model_str = "Illumina HiSeq 2500"

filetype_idx = 11
filetype_str = "fastq"

filename_1_idx = 12 # filename of first read

assembly_idx = 13 # nothing!

filename_2_idx = 14 # filename of second read (if available)

# read in sample accessions / strain names from files
sample_accession_map = {}
for filename in ["biosample_submission_%d_results.tsv" % file_idx for file_idx in xrange(1,5)]:
    file = open(filename,"r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        sample_accession = items[0].strip()
        sample_title = items[2].strip()
        
        sample_accession_map[sample_title] = sample_accession
        
sra_filename_file = open("sra_filenames.txt", "w")

sra_file_idx = 1
output_file = open("sra_submission_%d.tsv" % sra_file_idx, "w")
output_file.write("\t".join(column_headers))
  
population_file = open("population_samples.csv","r")
population_file.readline() # header
sra_idx = 0
sra_max_idx = 495

seen_sample_names = set([])

for line in population_file:

    sra_idx += 1

    if sra_idx >= sra_max_idx:
        output_file.close()
        sra_file_idx+=1
        output_file = open("sra_submission_%d.tsv" % sra_file_idx, "w")
        output_file.write("\t".join(column_headers))
        sra_idx = 1

    items = line.split(",")
   
    strain_str = items[0].strip()
    
    populations = [subitem.strip() for subitem in items[1].split(";")]
    
    timepoint = float(items[2])
    timepoint_str = "generation %d" % (long(timepoint*1000))
    
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
    
    extraction_batch_str = items[3].strip()
    filename_stem_str = items[6].strip()
    
    library_id_str = filename_stem_str
    
    population_str = ",".join([parse_file.get_pretty_name(population) for population in populations])
    if "BG" in strain_str:
        description_str = "%s clonal isolate" % (population_str)
    else:
        description_str = "%s %s mixed population sample" % (population_str, timepoint_str)
    
    
    column_strs = [" " for column in column_headers]
    
    column_strs[biosample_accession_idx] = sample_accession_map[sample_name_str]
    
    column_strs[bioproject_accession_idx] = bioproject_accession_str
    
    column_strs[title_idx] = sample_name_str

    column_strs[library_id_idx] = library_id_str

    column_strs[design_description_idx] = description_str
    
    column_strs[library_strategy_idx] = library_strategy_str
    
    column_strs[library_source_idx] = library_source_str

    column_strs[library_selection_idx] = library_selection_str 

    column_strs[library_layout_idx] = library_layout_str
    
    column_strs[platform_idx] = platform_str

    column_strs[instrument_model_idx] = instrument_model_str 

    column_strs[filetype_idx] = filetype_str

    column_strs[filename_1_idx] = "%s.R1.fastq.gz" % filename_stem_str

    column_strs[filename_2_idx] = "%s.R2.fastq.gz" % filename_stem_str

    seen_sample_names.add(sample_name_str)
    
    # a hack because the DL_*** samples have an extra R0 sample
    if filename_stem_str.startswith("DL"):
    
        sra_filename_file.write(column_strs[filename_1_idx]+"\n")
        sra_filename_file.write(column_strs[filename_2_idx]+"\n")
        
    
        output_file.write("\n")
        output_file.write("\t".join(column_strs))
    
        sra_idx += 1
        
        column_strs[library_layout_idx] = 'Single'
        column_strs[library_id_idx] = library_id_str+".R0"
        
        column_strs[filename_1_idx] = "%s.R0.fastq.gz" % filename_stem_str
        column_strs[filename_2_idx] = " "

        sra_filename_file.write(column_strs[filename_1_idx]+"\n")

        output_file.write("\n")
        output_file.write("\t".join(column_strs))
    
    elif filename_stem_str.startswith("CL") :   
        
        # we have to append L1 and L2 to filename stem
        column_strs[filename_1_idx] = "%s.L1.R1.fastq.gz" % filename_stem_str
        column_strs[filename_2_idx] = "%s.L1.R2.fastq.gz" % filename_stem_str
        column_strs[library_id_idx] = filename_stem_str+".L1"
        
        output_file.write("\n")
        output_file.write("\t".join(column_strs))
        
        
        sra_filename_file.write(column_strs[filename_1_idx]+"\n")
        sra_filename_file.write(column_strs[filename_2_idx]+"\n")
        
        
        
        sra_idx += 1 
        
        # we have to append L1 and L2 to filename stem
        column_strs[filename_1_idx] = "%s.L2.R1.fastq.gz" % filename_stem_str
        column_strs[filename_2_idx] = "%s.L2.R2.fastq.gz" % filename_stem_str
        column_strs[library_id_idx] = filename_stem_str+".L2"
        
        
        sra_filename_file.write(column_strs[filename_1_idx]+"\n")
        sra_filename_file.write(column_strs[filename_2_idx]+"\n")
        
        
        output_file.write("\n")
        output_file.write("\t".join(column_strs))
        
    else:
        
        
        sra_filename_file.write(column_strs[filename_1_idx]+"\n")
        sra_filename_file.write(column_strs[filename_2_idx]+"\n")
        
        output_file.write("\n")
        output_file.write("\t".join(column_strs))
        
sra_filename_file.close()