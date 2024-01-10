import csv,sys,re, glob,os
from collections import defaultdict, OrderedDict
import pandas as pd
import itertools
from collections import Counter
from itertools import chain



class HOG:
    def __init__(self, ID):
        self.ID = ID
        self.subgenomes = []
        self.cultivars = set()
        self.members = defaultdict(dict)
        self.passport = []
        

    def addSubGenome(self, s):
        self.subgenomes.append(s)

    def addOGID(self, o):
        self.OGID = o

    def addClade(self, clade):
        self.clade = clade

    def addCultivar(self, c):
        self.cultivars.add(c)

    def addMembers(self,k,v,n):
        self.members[k][v] = n

    def addPassport(self, p):
        self.passport.append(p)



class OG:
    def __init__(self, ID):
        self.ID = ID
        self.subgenomes = []
        self.cultivars = set()
        self.members = defaultdict(dict)
        

    def addSubGenome(self, s):
        self.subgenomes.append(s)

    def addCultivar(self, c):
        self.cultivars.add(c)

    def addMembers(self,k,v,n):
        self.members[k][v] = n



def getMissingGenes(hogsfile,allgIDs):

    dSpecies = defaultdict()
    dSpeciesGenes = defaultdict(list)
    dAllGenes = defaultdict(list)
    dGenesMissing = defaultdict(list)
    x = 3
    
    with open(hogsfile, "r") as fdh:
        lines = csv.reader(fdh,delimiter="\t")
        header = next(lines)
#         print(header)
        for h in header:
            if h not in ['HOG', 'OG', 'Gene Tree Parent Clade']:
#                 print(h)                
                _sp = re.sub('GsRTD_','',h)
                _sp = re.sub("_pep",'',_sp)
                _sp = re.sub(".primary.protein","",_sp)
#                 _sp = h.split('_')[1]
#                 print(_sp)
                dSpecies[x] = _sp
                x += 1
        if header != None:
            for line in lines:
                for i in range(3,(len(line))):
                    if line[i] != '':
                        mem = line[i].split(', ')
                        for g in mem:
#                             print(dSpecies[i],g)
                            dSpeciesGenes[dSpecies[i]].append(g)

    with open(allgIDs,'r') as fdh:
        lines = csv.reader(fdh,delimiter=":")
        for line in lines:
#                 print(line)
                sp = re.sub('.primary.protein.fa','',line[0])
                gid = re.sub('>','',line[1])
#                 print(sp,gid)
                dAllGenes[sp].append(gid)

    for s in dSpecies:
        _sp = dSpecies[s]
        result = set(dAllGenes[_sp]) - set(dSpeciesGenes[_sp])
        dGenesMissing[_sp] = result
    
    return dGenesMissing



def parseHOGs(hogsfile):

    dSpecies = defaultdict()
    dGeneNumbers = defaultdict(list)
    dHOGs = defaultdict()
    ddHOGs = {}
    x = 3

    with open(hogsfile, "r") as fdh:
        lines = csv.reader(fdh,delimiter="\t")
        header = next(lines)
#         print(header)
        for h in header:
            if h not in ['HOG', 'OG', 'Gene Tree Parent Clade']:
#                 print(h)                
                _sp = re.sub('GsRTD_','',h)
                _sp = re.sub("_pep",'',_sp)
                _sp = re.sub(".primary.protein","",_sp)
#                 _sp = h.split('_')[1]
#                 print(_sp)
                dSpecies[x] = _sp
                x += 1
        if header != None:
            for line in lines:
                HOGID = line[0]
                OGID = line[1]
                clade = line[2]
                dGeneNumbers[HOGID] = [0 if "" == cell else len(cell.split(", ")) for cell in line[3:]]
                dHOGs[HOGID] = line[3:]
                hog = HOG(HOGID)

                for i in range(3,(len(line))):
                    hog.addMembers(dSpecies[i],'',line[i])
                    hog.addCultivar(dSpecies[i])
                ddHOGs[HOGID] = hog


    return dGeneNumbers,dHOGs,dSpecies,ddHOGs


def parsePassPort(passportdata,dSpecies):
    d = defaultdict(list)

    with open(passportdata, "r") as fdh:
        lines = csv.reader(fdh,delimiter=",")
        header = next(lines)
        if header != None:
            for line in lines:
#                 print(line[4],line[7],line[8])
#                 d[line[4]].append((line[7],line[8]))
                d[line[4]] = [line[7],line[8]]

#     print(len(d))

    return d

def parseGeneTags(genetags,dSpecies):
    d = defaultdict(list)

    with open(genetags, "r") as fdh:
        lines = csv.reader(fdh,delimiter=";")
        header = next(lines)
        if header != None:
            for line in lines:
                gID = re.sub('ID=','',line[0])
                m = re.sub('method=','',line[1])
#                 print(gID,m)


#     print(len(d))

    return d


def calcDistribution(dGeneNumbers,dGenesMissing):

    d = defaultdict(list)
    dd = {}
    c = 0
    for i in dGeneNumbers:
        for x in range(0,77):
            if dGeneNumbers[i].count(0) == x:
                d[x].append(sum(dGeneNumbers[i]))

    for k in d:
        dd[k] = sum(d[k])
        
    for i in dGenesMissing:
        c = c + len(dGenesMissing[i])

    dd[76]=c

    df = pd.DataFrame.from_dict(dd,orient='index',columns=['genes'])
    df.index.name = 'genomes'
    print(df)
    df.to_csv("barley.conserved_variable.GeneNumbers.distribution.csv")





def main():

    hogsfile = sys.argv[1]
    allgIDs = sys.argv[2]
    passPortData = '../of/230719_221008_BPGv1v2_passport_data.mod.csv'
    genetags = '../of/BPGv2.projected_gene_tags.txt'

    dMapping = {}
    dMethod = {}
    _sp = set()
    d = defaultdict(list)
    
    coreHOGs = {}
    scHOGs = {}
    shellHOGs = {}
    gtHOGs = {}
    cloud = []

    conservedHOGs = {}
    singleCopyConservedHOGs = {}
    variableHOGs = {}

    dStatus = defaultdict(set)

    dSingleCopySummary = defaultdict(list)
    dCoreSummary = defaultdict(list)
    dShellSummary = defaultdict(list)
    dGTSummary = defaultdict(list)
    dSummary = defaultdict(list)


    dGeneNumbers,dHOGs,dSpecies,ddHOGs = parseHOGs(hogsfile)
    dGenesMissing = getMissingGenes(hogsfile,allgIDs)

    dPassPortData = parsePassPort(passPortData,dSpecies)
    parseGeneTags(genetags,dSpecies)


    for h in ddHOGs:
        for c in ddHOGs[h].cultivars:
            ddHOGs[h].addPassport(dPassPortData[c][1])
#             print(c,dPassPortData[c][1])


#     for i in dSpecies:
#         print(i,dSpecies[i],dPassPortData[dSpecies[i]][1])


    for i in dGeneNumbers:
        if dGeneNumbers[i].count(0) == 0:
#             print(i,"core",dGeneNumbers[i],len(dGeneNumbers[i]),sum(dGeneNumbers[i]),dGeneNumbers[i].count(1),dGeneNumbers[i].count(0),len(dSpecies))
#             print(ddHOGs[i].ID,ddHOGs[i].cultivars,ddHOGs[i].passport)
#             print(i,"core",len(dGeneNumbers[i]),sum(dGeneNumbers[i]),dGeneNumbers[i].count(1),dGeneNumbers[i].count(0),len(dSpecies))
            conservedHOGs[i] = ddHOGs[i]
            coreHOGs[i] = ddHOGs[i]
        if dGeneNumbers[i].count(1) == len(dSpecies):
#             print(i,"single-copy core",dGeneNumbers[i],len(dGeneNumbers[i]),sum(dGeneNumbers[i]),dGeneNumbers[i].count(1),dGeneNumbers[i].count(0),len(dSpecies))
            scHOGs[i] = ddHOGs[i]
            singleCopyConservedHOGs[i] = ddHOGs[i]
#         if dGeneNumbers[i].count(0) in range(1,len(dSpecies)-1) and sum(dGeneNumbers[i]) != 1:
        if dGeneNumbers[i].count(0)  in range(1,len(dSpecies)-1) and sum(dGeneNumbers[i]) != 1:
#             print(i,"shell",dGeneNumbers[i],len(dGeneNumbers[i]),sum(dGeneNumbers[i]),dGeneNumbers[i].count(1),dGeneNumbers[i].count(0),len(dSpecies))
#             print(ddHOGs[i].ID,ddHOGs[i].cultivars,ddHOGs[i].passport)
#             print(i,)
            shellHOGs[i] = ddHOGs[i]
            variableHOGs[i] = ddHOGs[i]
        if dGeneNumbers[i].count(0) != 0 and sum(dGeneNumbers[i]) == 1:
# #             print(i,"genotype-specific",dGeneNumbers[i],len(dGeneNumbers[i]),sum(dGeneNumbers[i]),dGeneNumbers[i].count(1),dGeneNumbers[i].count(0),len(dSpecies))
#         if dGeneNumbers[i].count(0) in range(1,len(dSpecies)-1) and sum(dGeneNumbers[i]) == 1:
# #             print(i)
            gtHOGs[i] = ddHOGs[i]
            variableHOGs[i] = ddHOGs[i]

        if dGeneNumbers[i].count(0) == len(dSpecies)-1:
# #             print(i,"genotype-specific",dGeneNumbers[i],len(dGeneNumbers[i]),sum(dGeneNumbers[i]),dGeneNumbers[i].count(1),dGeneNumbers[i].count(0),len(dSpecies))
            gtHOGs[i] = ddHOGs[i]
            variableHOGs[i] = ddHOGs[i]

# #             print(i)

    for i in dGenesMissing:
        for g in dGenesMissing[i]:
            cloud.append(g)


    print("core HOGs:",len(coreHOGs))
    print("single-copy HOGs:", len(scHOGs))
    print("shell HOGs:", len(shellHOGs))
    print("genotype-specific HOGs:", len(gtHOGs))
    print("cloud:", len(cloud))    
    print("Total:",len(coreHOGs)+len(shellHOGs)+len(gtHOGs))
    print()
    print("conserved HOGs:", len(conservedHOGs))
    print("sinlge-copy conserved HOGs:", len(singleCopyConservedHOGs))
    print("variable HOGs:",len(variableHOGs))

    for i in shellHOGs:
        hog = shellHOGs[i]
        for v in range(0,len(dSpecies)):
            if dGeneNumbers[i][v] != 0:
                if hog.passport[v] == 'landrace':
#                     print(hog.ID,hog.passport,dGeneNumbers[i])
                    dStatus['landrace'].add(hog.ID)

    for i in shellHOGs:
        hog = shellHOGs[i]
        for v in range(0,len(dSpecies)):
            if dGeneNumbers[i][v] != 0:
                if hog.passport[v] == 'wild':
#                     print(hog.ID,hog.passport,dGeneNumbers[i])
                    dStatus['wild'].add(hog.ID)

    for i in shellHOGs:
        hog = shellHOGs[i]
        for v in range(0,len(dSpecies)):
            if dGeneNumbers[i][v] != 0:
                if hog.passport[v] == 'feral':
#                     print(hog.ID,hog.passport,dGeneNumbers[i])
                    dStatus['cultivar'].add(hog.ID)

    for i in shellHOGs:
        hog = shellHOGs[i]
        for v in range(0,len(dSpecies)):
            if dGeneNumbers[i][v] != 0:
                if hog.passport[v] == 'cultivar':
#                     print(hog.ID,hog.passport,dGeneNumbers[i])
                    dStatus['cultivar'].add(hog.ID)


    df = pd.DataFrame.from_dict(dStatus,orient='index')
#     print(df)
    df_transposed = df.T
    print(df_transposed)
    df_transposed.to_csv("shell_status_noferal.csv")


    for i in singleCopyConservedHOGs:
        for sp in ddHOGs[i].members:
            if ddHOGs[i].members[sp][''] != '': 
#                 print(i,sp,ddHOGs[i].members[sp],len(ddHOGs[i].members[sp][''].split(', ')))
                dSingleCopySummary[sp].append(len(ddHOGs[i].members[sp][''].split(', ')))
            else:
#                 print(i,sp,0)
                dSingleCopySummary[sp].append(0)

    for i in conservedHOGs:
        for sp in ddHOGs[i].members:
            if ddHOGs[i].members[sp][''] != '': 
#                 print(i,sp,ddHOGs[i].members[sp],len(ddHOGs[i].members[sp][''].split(', ')))
                dCoreSummary[sp].append(len(ddHOGs[i].members[sp][''].split(', ')))
            else:
#                 print(i,sp,0)
                dCoreSummary[sp].append(0)

    for i in variableHOGs:
        for sp in ddHOGs[i].members:
            if ddHOGs[i].members[sp][''] != '': 
#                 print(i,sp,ddHOGs[i].members[sp],len(ddHOGs[i].members[sp][''].split(', ')))
                dShellSummary[sp].append(len(ddHOGs[i].members[sp][''].split(', ')))
            else:
#                 print(i,sp,0)
                dShellSummary[sp].append(0)


    for i in dSingleCopySummary:
        dSummary[i] = (sum(dSingleCopySummary[i]),sum(dCoreSummary[i]),sum(dShellSummary[i]),len(dGenesMissing[i]),dPassPortData[i][1])

    df = pd.DataFrame.from_dict(dSummary,orient='index',columns=['single-copy','conserved','variable','cloud','status'])
    df.index.name = 'genotype'
    df['delta_conserved'] = df['conserved']-df['single-copy']
    df['variableTotal'] = df['cloud']+df['variable']
    df['Total'] = df['cloud']+df['conserved']+df['variable']
    df['conserved_pct'] = df['conserved'].div(df['Total']).mul(100)
    df['variable_pct'] = df['variableTotal'].div(df['Total']).mul(100)
#     df['cloud_pct'] = df['variableTotal'].div(df['Total']).mul(100)
#     df.loc['Mean',:] = df.mean(axis=0)
    print(df)
    df.to_csv('barley.conserved_variable.GeneNumbers_sorted.csv')


    calcDistribution(dGeneNumbers,dGenesMissing)

#     with open('core.HOGs.tsv','w') as fout:
#         hed = []
#         for i in dSpecies:
#             hed.append(dSpecies[i])
#         print('hog',*hed,sep='\t',file=fout)
#         for i in coreHOGs:
#             print(i,*dHOGs[i],sep='\t',file=fout)
# 
#     with open('single-copy.HOGs.tsv','w') as fout:
#         hed = []
#         for i in dSpecies:
#             hed.append(dSpecies[i])
#         print('hog',*hed,sep='\t',file=fout)
#         for i in scHOGs:
#             print(i,*dHOGs[i],sep='\t',file=fout)
# 
#     with open('shell.HOGs.tsv','w') as fout:
#         hed = []
#         for i in dSpecies:
#             hed.append(dSpecies[i])
#         print('hog',*hed,sep='\t',file=fout)
#         for i in shellHOGs:
#             print(i,*dHOGs[i],sep='\t',file=fout)
# 
#     with open('gt-specific.HOGs.tsv','w') as fout:
#         hed = []
#         for i in dSpecies:
#             hed.append(dSpecies[i])
#         print('hog',*hed,sep='\t',file=fout)
#         for i in gtHOGs:
#             print(i,*dHOGs[i],sep='\t',file=fout)
# 
#     with open('cloud.unassigned_genes.tsv','w') as fout:
#         for i in dGenesMissing:
#             print(i,*dGenesMissing[i],file=fout)
#             


      
if __name__ == "__main__":
    main()
