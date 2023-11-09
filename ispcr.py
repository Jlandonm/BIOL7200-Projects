import subprocess
from operator import itemgetter
import tempfile
import os

def ispcr(primer_file: str, assembly_file: str, max_amplicon_size: int) -> str:

    blast_final = step_one(primer_file, assembly_file)
    paired_hits = step_two(blast_final, max_amplicon_size)
    sequences = step_three(paired_hits, assembly_file)
    return sequences

#Get and sort hits
def step_one(primer_file: str, assembly_file: str) -> list[list[str]]:
    outfmt = "'6 std qlen'"
    command = f'blastn -query {primer_file} -subject {assembly_file} -outfmt {outfmt} -task blastn-short'

    blast_out = run_blast(command)
    blast_final = filter_hits(blast_out)

    return blast_final

#Helper function to run blast
def run_blast(cmd):
    result = subprocess.run(args=cmd, capture_output=True, text=True, shell=True)
    blast_output = result.stdout.split("\n")
    blast_output.remove("")
    return blast_output

#Helper function to filter hits
def filter_hits(blast_in):
    blast_final = []
    #Filtering our hits here
    for line in blast_in:
        values = line.split("\t")
        values[8] = int(values[8])
        
        #If we have a hit
        if float(values[3]) >= (float(values[12]) * 0.8):
            blast_final.append(values)
    
    #In order to sort the output like needed, we need to convert the sequence position to int then convert back after we sort
    blast_final_sorted = sorted(blast_final, key=itemgetter(8))

    #Here we convert back to string
    for hit in blast_final_sorted:
        hit[8] = str(hit[8])
            
    return blast_final_sorted

#Pair the hits
def step_two(sorted_hits: list[list[str]], max_amplicon_size: int) -> list[tuple[list[str]]]:

    paired_hits = []
    #Loop through the hits
    for hit in sorted_hits:
        start_pos = int(hit[8])
        end_pos = int(hit[9])

        #Look for matches
        for hit2 in sorted_hits:
            start_pos2 = int(hit2[8])
            end_pos2 = int(hit2[9])

            # If the hits are facing each other and not overlapping
            if start_pos < end_pos and start_pos2 > end_pos2 and end_pos < end_pos2:
                # Making sure amplicon is the corect size
                if end_pos2 - end_pos < max_amplicon_size:
                    temp_tuple = (hit, hit2)
            # If the hits are facing each other and not overlapping (in other direction)
            elif start_pos > end_pos and start_pos2 < end_pos2 and end_pos > end_pos2:
                # Making sure amplicon is the corect size
                if end_pos - end_pos2 < max_amplicon_size:
                    temp_tuple = (hit2, hit)
            
            #Making sure the match isn't already in there
            if 'temp_tuple' in locals() and temp_tuple not in paired_hits:
                paired_hits.append(temp_tuple)

    return paired_hits

#Run Seqtk
def step_three(hit_pairs: list[tuple[list[str]]], assembly_file: str) -> str:
    bedfile = make_bedfile(hit_pairs)
    sequences = run_seqtk(bedfile, assembly_file)

    return sequences

#Making bedfile
def make_bedfile(hit_pairs: list[tuple[list[str]]]) -> str:
    bedfile = ""
    for hit in hit_pairs:
        id = hit[0][1]
        start_primer_pos = hit[0][9]
        end_primer_pos = str(int(hit[1][9]) - 1)
        bedfile = bedfile + f'{id}\t{start_primer_pos}\t{end_primer_pos}\n'
    bedfile = bedfile[:-1]

    return bedfile

#Running seqtk
def run_seqtk(bedfile: str, assembly_file: str) -> str:
    custom_filename = "temp_bedfile.bed"

    # Get the directory where temporary files are stored
    temp_dir = tempfile.gettempdir()

    # Create a temporary file with the custom name
    temp_file_path = os.path.join(temp_dir, custom_filename)

    # Open the temporary file in write mode
    with open(temp_file_path, 'w+') as temp_file:
        # Write your string to the temporary file
        temp_file.write(bedfile)
        #Go back to top of file
        temp_file.seek(0)

        #Initiate subprocess command using seqtk subseq
        command = f'seqtk subseq {assembly_file} {temp_file.name}'
        result = subprocess.run(args=command, capture_output=True, text=True, shell=True)

        #Get rid of trailing newline character
        seqtk_output = result.stdout
        seqtk_output = seqtk_output[:-1]

        return seqtk_output
    
