# modify the chromosome header line to include chrom id as chrN
# current: >CM021962.1 Canis lupus familiaris isolate Mischka breed German Shepherd chromosome 1, whole genome shotgun sequence
# new: >chr1 CM021962.1 Canis lupus familiaris isolate Mischka breed German Shepherd chromosome 1, whole genome shotgun sequence
# Zhanyang Zhu
# 6/20/2022

import re

OriFileName = "GCA_011100685.1_UU_Cfam_GSD_1.0_genomic.fna"
NewFileName = "Canis_familiaris.UU_Cfam_GSD_1.0.fa"
O = open(OriFileName, 'r')
N = open(NewFileName, 'w')

curLine = O.readline()
#
# search for header line that begins with > and include "chromosome N,"
# then insert ChrN after > sign:
#
foundChrN = True
printOn = True
numL = 0
while curLine != '': # end of file
  curLine = curLine.rstrip()
  if curLine.startswith(">"):
     # chromesome header line
     # >CM021992.1 Canis lupus familiaris isolate Mischka breed German Shepherd chromosome 31, whole genome shotgun sequence
     foundChrN = re.search(r'(?<= chromosome )[\dX]+(?=, )', curLine)
     foundMT = re.search(r'^>.* German Shepherd mitochondrion, complete sequence', curLine)
     printOn = True
     if foundChrN:
        # modify header:
        newChrID = ">chr" + foundChrN.group() + ' '
        curLine = curLine.replace(">", newChrID)
     elif foundMT: 
        # >CM022001.1 Canis lupus familiaris isolate Mischka breed German Shepherd mitochondrion, complete sequence, whole genome shotgun sequence
        # chrMT
        newChrID = ">chrMT "
        curLine = curLine.replace(">", newChrID)
     else:
        printOn = False
     if printOn:
        print(curLine)
  # copy the current line to output file:
  if printOn:
    print(curLine, file=N)
  # read a new line:
  curLine = O.readline()

# close original file and new file:
O.close()
N.close()