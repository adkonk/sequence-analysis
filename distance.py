from Bio import SeqIO
import time
import matplotlib.pyplot as plt
import pickle
import operator
import math

all_letters = "ACDEFGHIKLMNPQRSTVWY"

start_time = time.time()

try:
    with open('distances.p', 'rb') as fp:
        distances = pickle.load(fp)
except:
    iterator =  SeqIO.parse('GCF_000001735.3_TAIR10_protein.fasta', 'fasta')

    cumulative_distances = {}
    for letter in all_letters:
        cumulative_distances[letter] = 0

    distances = {}
    for letter in all_letters:
        distances[letter] = []

    record_number = 1
    while True:
        try:
            record = next(iterator)
        except:
            break
        distance_by_letter = {}
        for i in all_letters:
            distance_by_letter[i] = []
            before = 0
            after = 0
            while after != -1:
                after = record.seq.find(i, before)
                distance_by_letter[i].append(after-before)
                before = after + 1
            if len(distance_by_letter[i]) != 1:
                distance_by_letter[i] = (sum(distance_by_letter[i][:-1])+len(distance_by_letter[i][:-1]))/len(distance_by_letter[i][:-1])
            else:
                distance_by_letter[i] = 0
            cumulative_distances[i] += distance_by_letter[i]
            if distance_by_letter[i] != 0:
                distances[i].append(distance_by_letter[i])
        print(record_number)
        record_number += 1

    with open('distances.p', 'wb') as fp:
        pickle.dump(distances, fp, protocol=pickle.HIGHEST_PROTOCOL)

for i in distances.keys():
    distances[i] = [math.log(j) for j in distances[i]]

data = []
for i in distances.keys():
    data.append(distances[i])

pir = "MEVAFPSSPPVKCIADERLTEISMDFTKLGFPEEDEEIKKLERSWSKLEESVEFNEDDEEEEEEFSFACVNGEGSPITADEAFEDGQIRPVFPLFNRDLLFEYENEDDKNDNVSVTDENRPRLRKLFVEDRNGNGDGEETEGSEKEPLGPYCSWTGGTVAEASPETCRKSNSTGFSKLWRFRDLVLRSNSDGRDAFVFLNNSNDKTRTRSSSSSSSTAAEENDKKVITEKKKGKEKTSTSSETKKKTTTTKSAHEKLYMRNRAMKEEVKHRSYLPYKQVGFFTNVNGLSRNIHPF"

distance_by_letter = {}
pir_distances = [0]
for i in distances.keys():
    distance_by_letter[i] = []
    before = 0
    after = 0
    while after != -1:
        after = pir.find(i, before)
        distance_by_letter[i].append(after-before)
        before = after + 1
    if len(distance_by_letter[i]) != 1:
        listlength = len(distance_by_letter[i][:-1])
        distance_by_letter[i] = (sum(distance_by_letter[i][:-1])+listlength)/listlength
    else:
        distance_by_letter[i] = 0
    pir_distances.append(math.log(distance_by_letter[i]))

print("--- %s seconds ---\n" % (time.time() - start_time))

plt.figure(1)
plt.boxplot(data,labels=distances.keys())

PIR, = plt.plot(pir_distances, marker='^',color='y',
                markeredgecolor='k',linestyle='',
                label='Average distance in PIR', markersize=15)

plt.title("Distances between residues over A. thaliana proteome")
plt.xlabel("Amino acid")
plt.ylabel("log(Average distance between residues)")
plt.legend(handles=[PIR])

plt.show()

percentiles = {}
for i in range(len(all_letters)):
    letter = all_letters[i]
    listlength = len(distances[letter])
    j = 0
    reflist = sorted(data[i])
    while pir_distances[i+1] >= reflist[j]:
        j += 1
    percentiles[all_letters[i]] = (j+1)/listlength

print("{:1} {:<12} {:<15}".format(" ","Amino acid", "Percentile"))
for k, v in percentiles.items():
    print("{:1} {:<12} {:<15}".format(" ",k, round(v,3)))

'''counter = 1
for letter in distances.keys():
    plt.subplot(4,5,counter)
    plt.hist(distances[letter], 50,facecolor='green', alpha=0.75)

    plt.title("Average distribution of " + letter)
    counter += 1

plt.show()'''
