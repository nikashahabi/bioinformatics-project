import sys
sys.path.append("/opt/anaconda3/envs/viennarna/lib/python3.7/site-packages")
import RNA
import time
import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi
import forgi.graph.bulge_graph as fgb
import forgi.threedee.model.coarse_grain as ftmc


def rnaMutedStructure(rnaSeq, protein, rnaStart, proteinStart, frameshift, hashtable, hashtableAmino):
    mutationFile = open("mutationFile", "w")
    (ssInitial, mfeInitial) = RNA.fold(rnaSeq)
    initialRnaSeqData = [">" + rnaSeq + "\n", ssInitial + "\n", str(mfeInitial) + "\n"]
    mutationFile.writelines(initialRnaSeqData)
    ACGU = {"A", "C", "G", "U"}
    strOut = "replacing index {index} with {char}"
    strOutMfeDistance = "mfe of the initial sequence - mfe of the muted sequence = "
    strOutDistance = "distance = "
    rnaIndex, proteinStartIndex, rnaEndIndex = getIndex(rnaStart, proteinStart, frameshift)
    proteinIndex = proteinStartIndex
    while rnaIndex <= rnaEndIndex:
        for X in ACGU:
            if X != rnaSeq[rnaIndex] :
                newRNASeq = rnaSeq[:rnaIndex] + X + rnaSeq[rnaIndex + 1:]
                newProteinSeq = protein[:proteinIndex] + X + protein[proteinIndex + 1:]
                if checkIfProteinChanges(newProteinSeq, proteinIndex, proteinStartIndex, hashtable, hashtableAmino) is False:
                        (ss, mfe) = RNA.fold(newRNASeq)
                        mfeDistance = mfeInitial - mfe
                        dis = distance(ssInitial, ss)
                        # bpDistance = RNA.bp_distance(rnaSeq, newRNASeq)
                        seqRNAData = [">" + strOut.format(index=rnaIndex, char=X) + "\n",
                                      newRNASeq + "\n",
                                      str(ss) + "\n",
                                      str(mfe) + "\n",
                                      strOutMfeDistance + str(mfeDistance) + "\n",
                                      strOutDistance + str(dis) + "\n"]
                        mutationFile.writelines(seqRNAData)
                        # cg = ftmc.CoarseGrainRNA.from_dotbracket(dotbracket_str=ss, seq=newRNASeq)
                        # fvm.plot_rna(cg, text_kwargs={"fontweight": "bold", "size": 3}, lighten=0.7,
                        #              backbone_kwargs={"linewidth": 3})
                        # # plt.margins(x=0, y=-0.1)
                        # location = "figures/" + str(i) + "_" + X
                        # plt.savefig(location, dpi=600)
                        # plt.clf()
                else:
                    seqRNAData = [">" + strOut.format(index=rnaIndex, char=X) + "\n",
                                  newRNASeq + "\n",
                                  "This mutation changes the amino acid in the corresponding protein" + "\n"]
                    mutationFile.writelines(seqRNAData)

        proteinIndex += 1
        rnaIndex += 1




def dnaToRna(dnaSeq):
    rnaSeq = dnaSeq.replace("T", "U")
    rnaSeq = rnaSeq.replace("\n", "")
    return rnaSeq


def distance(ss1, ss2):
    def code(ss):
        stack = []
        codedStr = []
        codedStr[:0] = ss
        for i in range(len(ss)):
            if ss[i] == '(':
                stack.append(i)
            elif ss[i] == ')':
                popped = stack.pop(-1)
                codedStr[popped] = i + 1
                codedStr[i] = popped + 1
            elif ss[i] == '.':
                codedStr[i] = i + 1
            else:
                raise("structure not acceptable")
        return codedStr

    if len(ss1) != len(ss2):
        raise("two sequences do not have the same length")
    codedss1 = code(ss1)
    codedss2 = code(ss2)
    counter = 0
    for i in range(len(codedss1)):
        if codedss2[i] != codedss1[i]:
            counter += 1
    return counter


def checkIfProteinChanges(protein, proteinIndex, proteinStartIndex, hashtable, hashtableAmino):
    startCodon =  proteinStartIndex + 3*int((proteinIndex - proteinStartIndex)/3)
    currentAmino = hashtable[startCodon]
    mutedCodon = protein[startCodon: startCodon+3]
    mutedAmino = hashtableAmino[mutedCodon]
    if mutedAmino == currentAmino:
        return False
    return True


def createTranslateHash(protein, codonTable, frameshift=0):
    seq = protein.replace('\n', '').replace(' ', '')
    table = {}
    print(protein[frameshift: frameshift+3])
    for i in range(frameshift, len(seq), 3):
        codon = seq[i: i + 3]
        aminoAcid = codonTable.get(codon, '_')
        # if aminoAcid != '_':
        table[i] = aminoAcid
        print(aminoAcid, end="")
        # else:
        #     break
    return table


def getIndex(rnaStart, proteinStart, frameshift):
    proteinIndex = proteinStart + frameshift


start = time.time()
dnaSeq = """TTGCCTGTTTTTCCGCAACATAGTTACAGCTAAACATTTGCCCAAAACCATCTTCT
TATATATATATCTGCGATGGCGAGCCCAGCGGAAGGGATGTCCGCTTACTAATTCCGACA
CACCGGTTTAACCCCCCGGCGTGCTGGGGTGGGTGCCCCTAAGGGAGCGGGTTCTGTACT
TCCAGTAAGCGGCATTTGCGCTGTCATCGCCTTATCGAACCCGCTACTGAGATCATGTCC
TGAGTGGGTGAGTCGCACGCCCAATCGGCTGCA"""
rnaSeq = dnaToRna(dnaSeq)
codonTable = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
hashtable0 = createTranslateHash(dnaSeq, codonTable, frameshift=1)
print(hashtable0)
# hashtable1 = createTranslateHash(protein, codonTable, frameshift=1)
# hashtable2 = createTranslateHash(protein, codonTable, frameshift=2)
# rnaMutedStructure(rnaSeq, protein, rnaStart, proteinStart, frameshift, hashtable, codonTable)


# i = 102
# (ssInitial, mfeInitial) = RNA.fold(rnaSeq)
# ACGU = {"A", "C", "G"'', "U"}
# X = 'G'
# newRNASeq = rnaSeq[:i] + X + rnaSeq[i + 1:]
# print(newRNASeq)
# (ss, mfe) = RNA.fold(newRNASeq)
# ss = """...(((((((....((((((....((((....))))...............(((((...................)))))........))))))((..((((((..........))))))..))(((....((((((((....)))))).))..)))((....))(((((((((.(((.........((((.......))))..........))).)))))))))(((.((........)).))))))))))..(((.(((.....))))))."""
# cg = ftmc.CoarseGrainRNA.from_dotbracket(dotbracket_str=ss, seq=newRNASeq)
# fvm.plot_rna(cg)
# plt.show()
# location = "test"
# # plt.savefig(location, dpi=600)
# # plt.clf()
# end = time.time()



