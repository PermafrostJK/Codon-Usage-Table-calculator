def FASTA_reader(filename):
    n=0
    sequence = ''
    f = open(filename, 'r')
    l = l=f.readlines()
    for i in l:
        if i[0] == '>':
            n+=1
        elif n<2:
            sequence+=i[:-1] #removes the '\n' from the end of each line
        else:
            break
    return sequence

def complement_finder(seq):
    reverse_seq = seq[::-1]
    comp_seq = ''
    d_codon_comp = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}
    for i in reverse_seq:
        comp_seq += (d_codon_comp[i])
    return comp_seq

def sequence_extractor(loc, genome): # this does not handle complementary sequences
    # loc is a list of location ranges in the form of "a..b"; 
    # genome is the complied genomic sequence obtained from FASTA
    seq = ''
    for i in loc: # here we expect the structure of i to be [a..b],  where a & b are the start and end nucleotide positions
        beginning = int(i.split('.')[0])-1 # the nucleotide at the 190 position in the genome should have index 189
        end = int(i.split('.')[-1]) # no -1 here since we want to include the final nucleotide specified in loc

        seq += genome[beginning:end] # somehow this step is generating non-integer indices
    return seq

def list_of_aa(codon_table='codon_table.csv'):
    aa_list = []
    with open(codon_table) as read_file:
        read_file.readline()  # skip the first line
        for line in read_file:
            split_line = line.strip().split(',')
            
            if not(split_line[2] in aa_list):
                aa_list.append(split_line[2])
    return(aa_list)

def codon_aa_table(codon_table='codon_table.csv'):
    codon_dict = {}
    with open(codon_table) as read_file:
        read_file.readline()  # skip the first line
        for line in read_file:
            split_line = line.strip().split(',')
            
            codon_dict[split_line[0]] = split_line[2]
    return(codon_dict)

def codon_counting_table(codon_table='codon_table.csv'):
    codon_dict = {}
    with open(codon_table) as read_file:
        read_file.readline()  # skip the first line
        for line in read_file:
            split_line = line.strip().split(',')
            
            codon_dict[split_line[0]] = 0
    return(codon_dict)

def codon_counter(seq, cc):
    a = 0
    b = 3
    while b <= len(seq): 
        codon = seq[a:b]
        cc[codon] += 1
        a += 3 
        b += 3    
    return cc

#==========================================================================================================================

f=open('/Users/jackfrostljx/Desktop/Genomics_and_Bioinformatics/Critical_Paper_Analysis/Data/firmicutes.gbff', 'r') 
l=f.readlines() # use // to indicate termination of file reading
f.close()
for i in l:
    if '//' in i:
        index = l.index(i)
        l = l[:index]
        break
l_loc, l_gene_loc, l_seq = [],[],[]
x, p = 0, 0
y = 1 
d_aa_codon= {}

cc = codon_counting_table(codon_table='codon_table.csv')
ca = codon_aa_table(codon_table='codon_table.csv')
l_aa = list_of_aa(codon_table='codon_table.csv')
genome=FASTA_reader('Ecoli_test.fna')

for i in l:
    if i[0:5]==' '*5 and i[5:9] == 'gene': # since we are dealing with prokaryotes, no need to extract CDS.
        l_loc.append([])
        l_loc[x].append(i)
        position = l.index(i)
        while l[position+y][0:21] ==' '*21 and l[position+y][21].isdigit(): 
        # does not seem to work to screen for location info stored in multiple lines
        # but seeing as it is prokaryotes which do not normally contain exons, can relatively safely ignore
            l_loc[x].append(l[position+y])
            y+=1
        x+=1
        y=1


for i in l_loc: # for each CDS 
    location = [] # list/list of lists 
    for j in i: # for each location line of that CDS
        if 'join' in j: # ['feature_name     join(a..b, c..d, e..f', 'm..n, x..y)']
            pass 
# there are no genes that have joined locations in the E. coli genome, so I omit this part
# of the code for simplicity

# likewise here we ommit the multi-line checking and simply treat each gene as having its location
# documented in one line in the gbff file alone
            # if 'complement' in j: # ['  feature_name   complement(a..b)']
            #     proto_loc = j.split()[-1] # extracts 'complement(a..b)'
            #     index_start = proto_loc.index('(')+1
            #     index_end = proto_loc.index(')')
            #     location.append(proto_loc[index_start:index_end])
            # else:
            #     proto_loc = j.split() 
            #     for k in proto_loc: # for each segment of that CDS
            #         if k[0:10] == '':
            #             location.append(k)
        else:
            if 'complement' in j: 
            # ['  feature_name   complement(a..b)'] # in future completed versions of this code 
            # the check for presence of 'complement' should take place before "for j in i:", since if there a multi-line
            # location is denoting the gene's presense on the complementary strand, then non-first row location information
            # needs to be treated the same way as the first row which contains the str 'complement'
                proto_loc = j.split()[-1] # extracts 'complement(a..b)'
                index_start = proto_loc.index('(')+1
                index_end = proto_loc.index(')')
                location.append(proto_loc[index_start:index_end])
                location.append('comp') 
            else:
                proto_loc = j.split() 
                for k in proto_loc: # for each segment of that CDS
                    if k[0].isdigit():
                        location.append(k)
    l_gene_loc.append(location)

for i in l_gene_loc:
    if len(i) == 1:
        seq = sequence_extractor(i, genome)
        codon_counter(seq, cc)
    elif len(i) == 2: 
# since only genes on the complementary strand have their location list appended with 'comp', 
# they are the only ones that have a list length of 2
        seq = sequence_extractor([i[0]], genome)
        comp_seq = complement_finder(seq)
        codon_counter(comp_seq, cc) 
    l_seq.append(seq)
for k,v in ca.items():
 	d_aa_codon[v] = d_aa_codon.get(v,[])+[k]

d_aa_count = {}
for k in cc:
    for i in l_aa:
        if k in d_aa_codon[i]:
            d_aa_count[i] = d_aa_count.get(i,0) + cc[k]
d_CUT ={} # codon usage table
for k in cc:
    try:
        d_CUT[k] = float(str(cc[k]/d_aa_count[ca[k]])[0:7]) # 4 decimal places
    except:
        d_CUT[k] = 0

total_coding_region = 0
repeat_tracker = [] #handles variant repeats in genes 
for i in l_seq:
    if not(l_gene_loc[l_seq.index(i)] in repeat_tracker):
        total_coding_region += len(i)
        repeat_tracker.append(l_gene_loc[l_seq.index(i)])
        for j in i:
            if j=="G" or j=='C':
                p+=1

genome_GC_content = str(100*p/len(genome)) + '%' 
percent_CDS = str(total_coding_region/len(genome)*100) + '%' 

print('Number of individual codon appearance: \n', cc,'\n') 
print('Amino acid - codon table:\n', d_aa_codon,'\n')
print('Number of each amino acid:\n', d_aa_count,'\n')
print('Codon Usage Table: \n', d_CUT,'\n')
print('GC content:', genome_GC_content, '\n'+'Pecent of genome that codes fo protein:', percent_CDS)

# add the feature which allows the addition of codons that are not used by the organism as 0 usage
# make sure the code works for random code snipets 