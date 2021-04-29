import sys
sys.path.append("/opt/anaconda3/envs/viennarna/lib/python3.7/site-packages")
import RNA
import time


def rnaMutedStructure(rnaSeq):
    mutationFile = open("mutationFile", "w")
    (ssInitial, mfeInitial) = RNA.fold(rnaSeq)
    initialRnaSeqData = [">" + rnaSeq + "\n", ssInitial + "\n", str(mfeInitial) + "\n"]
    mutationFile.writelines(initialRnaSeqData)
    ACGU = {"A", "C", "G", "U"}
    strOut = "replacing index {index} with {char}"
    strOutMfeDistance = "mfe of the initial sequence - mfe of the muted sequence = "
    strOutDistance = "distance = "
    # strOutBPDistance = "base pair distance = "
    for i in range(len(rnaSeq)):
        for X in ACGU:
            if X != rnaSeq[i]:
                newRNASeq = rnaSeq[:i] + X + rnaSeq[i + 1:]
                (ss, mfe) = RNA.fold(newRNASeq)
                mfeDistance = mfeInitial - mfe
                dis = distance(ssInitial, ss)
                # bpDistance = RNA.bp_distance(rnaSeq, newRNASeq)
                seqRNAData = [">" + strOut.format(index=i, char=X) + "\n",
                              str(ss) + "\n",
                              str(mfe) + "\n",
                              strOutMfeDistance + str(mfeDistance) + "\n",
                              strOutDistance + str(dis) + "\n"]
                mutationFile.writelines(seqRNAData)


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



start = time.time()
dnaSeq = """TTTCTTGCCTGTTTTTCCGCAACATAGTTACAGCTAAACATTTGCCCAAAACCATCTTCT
TATATATATATCTGCGATGGCGAGCCCAGCGGAAGGGATGTCCGCTTACTAATTCCGACA
CACCGGTTTAACCCCCCGGCGTGCTGGGGTGGGTGCCCCTAAGGGAGCGGGTTCTGTACT
TCCAGTAAGCGGCATTTGCGCTGTCATCGCCTTATCGAACCCGCTACTGAGATCATGTCC
TGAGTGGGTGAGTCGCACGCCCAATCGGCTGCA"""
rnaSeq = dnaToRna(dnaSeq)
# rnaMutedStructure(rnaSeq)
end = time.time()
print(len(rnaSeq))
print(end - start)
print(RNA.fold(dnaToRna(dnaSeq)))
# print(type())


