import sys

def getSeq(fileName): #for fasta with one header
  seq = []
  file = open(fileName)
  next(file)
  for line in file:
    seq.extend(list(line.strip()))
  file.close()
  return seq

def reverseComplement(seq): #list
  seq.reverse()
  bases = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
  seq = [bases[base] for base in seq]
  return seq

def findStart(seq):
  for i in range(0, len(seq)-3, 3):
    if seq[i:i+3] == 'ATG':
      return i
  return -1

def findEnd(seq, startIndex):
  stopCodons = ['TAA', 'TAG', 'TGA']
  for i in range(startIndex+3, len(seq)-3, 3):
    if seq[i:i+3] in stopCodons:
      return i+2
  return -1

def orf(seq):
  begin = findStart(seq)
  if begin >= 0:
    end = findEnd(seq, begin)
    if end >= 0:
      return (begin, end)
  return False

def getRevOrfs(seq):
  orfs = []
  for i in range(2):
    if orf(seq[i:]):
      begin, end = orf(seq[i:])
      orfs.append((len(seq)-(end+i), len(seq)-(begin+i)))
  return orfs

#fileName = sys.argv[1]
fileName = 'Mycoplasma.fa'
#fileName = 'My.fa'
seq = getSeq(fileName)
revSeq = reverseComplement(seq)
seq = ''.join(seq)
orfs = getOrfs(seq)
orfs.append('Reverse Complement:')
orfs.extend(getRevOrfs(revSeq))
print(orfs)