# Import libraries up front
import json
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import csv
# From Table S13 in Plaisier et al., Cell Systems 2016
# These are Entrez IDs (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3013746/)
input = ['430', '1052', '1053', '1385', '84699', '9586', '1871', '1874', '144455', '79733', '1960', '1997', '2002', '2004', '80712', '2114', '2115', '2120', '51513', '2551', '2623', '2624', '2625', '9421', '3232', '10320', '3659', '3662', '3670', '91464', '3726', '10661', '11278', '128209', '10365', '9314', '1316', '51176', '9935', '23269', '4602', '4774', '4790', '7025', '9480', '5468', '5914', '5916', '3516', '5971', '864', '6257', '4093', '6659', '6660', '6662', '25803', '347853', '30009', '9496', '6929', '6925', '8463', '7022', '29842', '10155', '6935', '132625', '23051', '85416', '7707', '7764', '23528', '201516']

# Loading JSON file
# https://www.safaribooksonline.com/library/view/python-cookbook-3rd/9781449357337/ch06s02.html
# Example:
# import json
#
# # Reading data back
# with open('data.json', 'r') as f:
#      data = json.load(f)

# Reading TF regulator to TF target gene relationships into Python
# The json library we import takes care of most of the work
with open('tfbsDb_plus_and_minus_5000_entrez.json', 'r') as f:
     tfbsDb = json.load(f)

# Example set of keys in tfbsDb, they are Motif IDs (http://jaspar.genereg.net/search?q=Homo%20sapiens&collection=CORE&tax_group=vertebrates)
print(list(tfbsDb.keys())[0:5])

# Example set of values under a specific Motif ID, they are Entrez IDs
print(tfbsDb[list(tfbsDb.keys())[0]][0:5])

print(len(tfbsDb[list(tfbsDb.keys())[0]]))
# Read in humanTFs file
id2motif = dict()
motif2id = {}

with open('id_conversion/humanTFs_All.CSV','r') as inFile:
    # Use the readline() function to read in a single line
    # strip() gets rid of the newline character at the end of the line
    # split(',') splits up the line into columns based on commas
    # Tried using pandas,but dict conversion was more ocmplex and intensive than current code
    header = inFile.readline().strip().split(',')
    print(header)

    
    #id2motif = Entrez ID:Motif Name for Entrez ID,Motif name in HumanTfs 
    while 1:
        inLine = inFile.readline()
        if not inLine:
            break
        split = inLine.strip().split(',') 
        # TODO Fill out the id2moitf dictionary (key = Entrez ID, value = Motif Name)
        if not split[2] in id2motif:
            id2motif[split[2]] = []
        id2motif[split[2]].append(split[0])
        
        
    
        # TODO Fill out the motif2id dictionary (key = Motif Name, value = Entrez ID)
Family2Id = {}
Id2Family = {}
with open('id_conversion/tfFamilies.csv','r', encoding='iso-8859-1') as openFile:     # opening file 
    header = openFile.readline().strip().split(',')       # reading in first line of file as header
    #print(header)
    while 1:
        inLine = openFile.readline()
        if not inLine:
            break
        strip = inLine.strip().split(',')
        strip2 = strip[2].split(' ')         # strip2 to seperate Entrez id values 
        Family2Id[strip[0]] = strip2          # adding keys and values to Family2Id dictionary from file
        
        for IdList in strip2:
            Id2Family[IdList] = []
            Id2Family[IdList] = strip[0]

            
#print (Id2Family.keys())        
#print (Family2Id.values())

## To build a TF regulator to TF target gene network (constrained to TFs within the input list).
## This will require mapping from:
##     1. Input list of potential TF regulator Entrez Gene IDs (input)
##     2. List of Motif IDs for an Entrez Gene ID in the input list (either id2motif or motif2id)
##     3. TF target genes that are Entrez Gene IDs that are the values under a specific Motif ID in tfbsDb
##     4. Restrict TF target genes to only those in the input list
##     5. Add new entry to tfNetwor dictionary that has as the key the TF regulator and the values all the TF target genes
tfNetwork = {}
TFRegulator={}
for TranscriptionFactors in input:    # for loop iterates for each entrez id given in the input list
    if TranscriptionFactors in id2motif:    # if the transcription factors (entrez ids) are located in the id2motif dict
        for Motif in id2motif[TranscriptionFactors]:    # loop function that checks each motif that corresponds to entrez id in input list
            if Motif in tfbsDb: # checks the JSON file for the determined motifs
                targets = tfbsDb[Motif]          # assign targets based on the entrez ids that correspond to the motifs
                for eachTarget in targets:  # loop iterates for each target identified
                    if not TranscriptionFactors in tfNetwork: 
                        tfNetwork[TranscriptionFactors] = []
                    if eachTarget in input and not eachTarget in tfNetwork[TranscriptionFactors]:
                        tfNetwork[TranscriptionFactors].append(eachTarget)
    else:
        for TFFamilies in Family2Id: # for loop iterates for every family in tf families
            if TranscriptionFactors in Family2Id[TFFamilies]: #if there is a desired input entrez id in the family list
                for Id in Family2Id[TFFamilies]:    # for each id in the family list
                    if Id in id2motif: #check for id in id2motif
                        for Motifs in id2motif[Id]:    # loop function that checks motif2id in id2motif[eachTfReg]
                            if Motifs in tfbsDb:    # checks for motifs in JSON file
                                targets = tfbsDb[Motifs]          # assign targets from id2motif[eachTfReg]/eachTfreg
                                for eachTarget in targets:  #loop iterates for each target identified
                                    if not TranscriptionFactors in tfNetwork:
                                        tfNetwork[TranscriptionFactors] = []
                                    if eachTarget in input and not eachTarget in tfNetwork[TranscriptionFactors]:
                                        tfNetwork[TranscriptionFactors].append(eachTarget)

netConnections = []
for TFreg in tfNetwork:
    for TFtarg in tfNetwork[TFreg]:
        netConnections.append((TFreg,TFtarg))

        
#print(netConnections)

G = nx.DiGraph()
G.add_edges_from(netConnections)
#print(G)
pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), node_size = 500)
nx.draw_networkx_labels(G, pos)
nx.draw_networkx_edges(G, pos)
plt.show()

