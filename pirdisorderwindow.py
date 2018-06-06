#For sequence creation and annotation
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

#For making and analyzing alignment object
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.SubsMat import FreqTable

#Necessary to construct frequency table
from Bio import Alphabet

#For constructing alignment diagram
from Bio.Graphics import GenomeDiagram

#For formatting genome diagram
from reportlab.lib import colors
from reportlab.lib.units import cm

#For naming output file
from os.path import basename

################################################################################
''' Define functions '''
################################################################################

#Calculate location in MSA based on location in unaligned sequence
def forwardlocationfinder(sequence,unalignedlocation):
    MSAcounter = 0
    residuecounter = 0
    for i in sequence:
        if i == "-":
            MSAcounter += 1
        elif i != "-":
            residuecounter += 1
            MSAcounter += 1
        if residuecounter >= unalignedlocation:
            break
    return MSAcounter

#Calculate location in unaligned sequence based on location in MSA
def reverselocationfinder(sequence,MSAlocation):
    location = 0
    for i in sequence[:MSAlocation]:
        if i !="-":
            location += 1
        else:
            continue
    return location

#Input a list of values, return all intervals of values below a threshold
def valueintervalfinder(scores, threshold):
    i = 0
    output = []
    while i < len(scores):
        searchforfirstchar = True
        searchforlastchar = True
        while searchforfirstchar:
            if scores[i] <= threshold:
                searchforfirstchar = False
                firstcharposition = i
                while searchforlastchar and i < len(scores):
                    if scores[i] > threshold:
                        searchforlastchar = False
                        lastcharposition = i
                        output.append((firstcharposition, lastcharposition))
                    elif i == len(scores) - 1:
                        lastcharposition = i
                        output.append((firstcharposition, lastcharposition))
                        i += 1
                    else:
                        i += 1
            elif i == len(scores) - 1:
                searchforfirstchar = False
                i += 1
            else:
                i += 1
    return output

#Input a string, return all intervals of a character
def stringintervalfinder(string, searchfor):
    i = 0
    output = []
    while i < len(string):
        searchforfirstD = True
        searchforlastD = True
        while searchforfirstD:
            if string[i] == searchfor:
                searchforfirstD = False
                firstDposition = i
                while searchforlastD and i < len(string):
                    if string[i] != searchfor:
                        searchforlastD = False
                        lastDposition = i
                        output.append((firstDposition, lastDposition))
                    elif i == len(string) - 1:
                        lastDposition = i
                        output.append((firstDposition, lastDposition))
                        i += 1
                    else:
                        i += 1
            elif i == len(string) - 1:
                searchforfirstD = False
                i += 1
            else:
                i += 1
    return output

#Input a string, return all intervals that are not a character
def reversestringintervalfinder(string, avoid):
    i = 0
    output = []
    while i < len(string):
        searchforfirstchar = True
        searchforlastchar = True
        while searchforfirstchar:
            if string[i] != "-":
                searchforfirstchar = False
                firstcharposition = i
                while searchforlastchar and i < len(string):
                    if string[i] == avoid:
                        searchforlastchar = False
                        lastcharposition = i
                        output.append((firstcharposition, lastcharposition))
                    elif i == len(string) - 1:
                        lastcharposition = i
                        output.append((firstcharposition, lastcharposition))
                        i += 1
                    else:
                        i += 1
            elif i == len(string) - 1:
                searchforfirstchar = False
                i += 1
            else:
                i += 1
    return output

################################################################################
''' Open disorder prediction file and pull out predictions '''
################################################################################

#Input base of filename (without "_protein_translation")
filename = str(input("Base file name (no extension): "))

#Read folder of results from SPOT-disorder
all_descriptions = []
protein_sequences = []
predictions = []
disorder_location = str(input("Folder of disorder results: "))
with open(disorder_location + "/list.desc") as file:
    for description in file:
        with open(disorder_location + "/" + str(description.split()[0]) +
            ".spotd", "rt") as disopred:
            all_descriptions.append(str(description[:-1]))
            protein_sequences.append("")
            predictions.append("")
            disopredtext = disopred.read().split()
            for i in disopredtext[5::4]:
                protein_sequences[-1] += i
            for i in disopredtext[7::4]:
                predictions[-1] += i

################################################################################
''' Align disorder predictions using MSA '''
################################################################################

#Read alignment file, align disorder predictions by sequence alignment
alignment = AlignIO.read(filename + "_protein_translation_MSA.fa", "fasta")
for value in range(len(alignment)):
    alignseq = alignment[value]
    predseq = predictions[value]
    modifiedpredseq = ""
    disordercount = 0
    counter = 0
    for i in alignseq:
        if i != "-":
            if disordercount < len(predseq):
                modifiedpredseq += predseq[disordercount]
                disordercount += 1
                counter +=1
            else:
                break
        elif i == "-":
            modifiedpredseq += "-"
            counter += 1
    predictions[value] = modifiedpredseq

################################################################################
''' Calculate conservation of sequences in alignment using sliding window '''
################################################################################

expect_freq = {'K':0.0635,'R':0.0536,'H':0.0224,'E':0.0666,'D':0.0539,
    'T':0.0515,'S':0.0885,'W':0.0125,'V':0.0676,'L':0.0936,'Y':0.0287,
    'M':0.0244,'I':0.0528,'G':0.0661,'A':0.0567,'N':0.0432,'C':0.0175,
    'Q':0.0345,'F':0.0426,'P':0.0487,'X':0
    }
summary_align = AlignInfo.SummaryInfo(alignment)
e_freq_table = FreqTable.FreqTable(expect_freq, FreqTable.FREQ,
    Alphabet.IUPAC.Alphabet)
pseudo_count = 0.1
chars_to_ignore = ['-']
log_base=2

window_size = int(input("Sliding window size for conservation: "))
if window_size == 1:
    half_window = 0
elif window_size % 2 == 1:
    half_window = int((window_size-1)/2)
else:
    raise ValueError("Window size should be an odd number")

conservation_scores = []
for j in range(half_window):
    window_score = summary_align.information_content(
        start=0,
        end=j+half_window+1,
        e_freq_table=e_freq_table,
        log_base=2,
        chars_to_ignore=chars_to_ignore,
        pseudo_count=pseudo_count)
    conservation_scores.append(round(float(window_score)/(j+half_window+1),2))
for j in range(half_window+1, len(alignment[0])-(half_window+1)):
    window_score = summary_align.information_content(
        start=j-half_window,
        end=j+half_window+1,
        e_freq_table=e_freq_table,
        log_base=2,
        chars_to_ignore=chars_to_ignore,
        pseudo_count=pseudo_count)
    conservation_scores.append(round(float(window_score)/(window_size),2))
for j in range(len(alignment[0])-half_window, len(alignment[0])):
    window_score = summary_align.information_content(
        start=j-half_window,
        end=len(alignment[0]),
        e_freq_table=e_freq_table,
        log_base=2,
        chars_to_ignore=chars_to_ignore,
        pseudo_count=pseudo_count)
    conservation_scores.append(round(
        float(window_score)/(len(alignment)-j+half_window),2))

first_quartile = conservation_scores[int(0.25*len(conservation_scores))]
median = conservation_scores[int(0.5*len(conservation_scores))]
third_quartile = conservation_scores[int(0.75*len(conservation_scores))]
print("Summary of conservation: \nFirst quartile: " + str(first_quartile) +
        "\nAverage: " + str(round(sum(conservation_scores)/
            len(conservation_scores),2)) +
        "\nMedian: " + str(median) +
        "\nThird quartile: " + str(third_quartile) +
        "\nMaximum: " + str(max(conservation_scores)))

score_threshold = float(input("Conservation setting (bits): "))
i = 0
conservationlocations = valueintervalfinder(conservation_scores,
    score_threshold)

################################################################################
''' Construct Genome Diagram '''
################################################################################

name = str(str(basename(__file__)).replace(".py","") + "_" +
    str(filename) + str(score_threshold))
gd_diagram = GenomeDiagram.Diagram(name)
max_len=0
for record in alignment:
    max_len = max(max_len, len(record))

#Set opacity for conservation mask (optional input setting)
alphasetting = 0.6
#float(input("Opacity of conservation mask: (1 is full, 0 is none) "))

athalianaindices = []
for value in range(len(alignment)):
    record = alignment[value]

    #Name each track by species, or full description otherwise
    if all_descriptions[value].find("[") != -1:
        name_start = all_descriptions[value].find("[")+1
    else:
        name_start = 0
    if all_descriptions[value].find("]") != -1:
        name_end = all_descriptions[value].find("]")
    else:
        name_end = len(all_descriptions[value])
    record.name = all_descriptions[value][name_start:name_end]

    #Save indices of sequences that come from A. thaliana
    if record.name == "Arabidopsis thaliana":
        athalianaindices.append(value)

    #Create track for this species (include scale labels if the first track)
    if value == 0:
        gd_track_for_features = gd_diagram.new_track(1,
                                name=record.name,
                                greytrack=1,
                                greytrack_labels=1,
                                greytrack_fontsize=16,
                                scale=True,
                                scale_largeticks=0.8,
                                scale_largetick_interval=100,
                                scale_smalltick_interval=20,
                                scale_fontsize=30,
                                scale_fontangle=0,
                                start=0,
                                end=max_len)
    elif value != 0:
        gd_track_for_features = gd_diagram.new_track(1,
                                name=record.name,
                                greytrack=1,
                                greytrack_labels=1,
                                greytrack_fontsize=16,
                                scale=True,
                                scale_largeticks=0.8,
                                scale_largetick_interval=100,
                                scale_smalltick_interval=20,
                                scale_largetick_labels=False,
                                start=0,
                                end=max_len)
    gd_feature_set = gd_track_for_features.new_set()

    #Determine intervals of disorder and alignment
    disorderlocations = stringintervalfinder(predictions[value], "D")
    alignmentlocations = reversestringintervalfinder(alignment[value], "-")

    #Add features to feature set
    for set,look in [(alignmentlocations,colors.blue),
                    (disorderlocations,colors.red),
                    (conservationlocations,colors.Color(0,0,0,
                        alpha=alphasetting))
        ]:
        for start, end in set:
            feature = SeqFeature(FeatureLocation(start, end))
            gd_feature_set.add_feature(feature, color=look)

################################################################################
''' Draw Genome Diagram '''
################################################################################

#Draw diagram and write to file (optional .svg file)
gd_diagram.draw(format="linear", pagesize=(45*cm, 80*cm), fragments=1,
                start=0, end=max_len, orientation="landscape")
gd_diagram.write(name + ".pdf", "PDF")
#gd_diagram.write(name + ".svg", "SVG")

################################################################################
''' Print conserved Arabidopsis sequence '''
################################################################################

#Convert from intervals to a list of all indices
notconservedresidues = []
for start, end in conservationlocations:
    for i in range(start,end):
        notconservedresidues.append(i)

#Determine which sequence is A. thaliana
if athalianaindices is not None:
    index = min(athalianaindices)
else:
    response = str(input("Is the first sequence in the alignment " +
        "the sequence of interest? (y/n) "))
    if response != "y":
        index = int(input("What is the index of the sequence of interest? " +
            "(indexing begins at 0) "))

#Show conserved residues in uppercase, not conserved residues in lowercase
#Include spacing between every 10 residues and line breaks between every 50
counter = 0
conservationstring = "Conservation in sequence: \n"
for i in range(len(alignment[index])):
    if alignment[index][i] != '-':
        if i in notconservedresidues:
            conservationstring += alignment[index][i].lower()
        elif i not in notconservedresidues:
            conservationstring += alignment[index][i].upper()
        counter += 1
        if (counter % 10) == 0:
            if (counter % 50) == 0:
                conservationstring += "\n"
            else:
                conservationstring += " "

print(conservationstring)

################################################################################
''' (Optional) Convert between sequence positions '''
################################################################################

#MSA to unaligned
response = str(input("Do you want to calculate residue locations from MSA " +
    "to unaligned sequence? (y/n) "))
while response == "y":
    selected = int(input("Location in MSA: "))
    print(reverselocationfinder(alignment[index],selected))
    response = str(input("Continue? (y/n) "))

#Aligned to MSA
response = str(input("Do you want to calculate residue locations from " +
    "unaligned to MSA sequence? (y/n) "))
while response == "y":
    selected = int(input("Location in unaligned sequence: "))
    print(forwardlocationfinder(alignment[index],selected))
    response = str(input("Continue? (y/n) "))
