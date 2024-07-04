from Bio import AlignIO
import os, sys
import csv
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from glob import glob


Ref="Reference.fasta"   #Reference sequence with which other sequences are to be aligned

with open("GapInfo.csv", 'w', newline='') as file:       #Saves info on which sequences have insertions, so that Ref sequence has alignment gaps

    writer = csv.writer(file)
    writer.writerow(["SN","Sequence","Insertion"])

    #Run for each IonCode fasta input....

    cnt=0      #Count for no. of sequences read

    for filename in glob('IonCode*.fasta'):

        Seq = str(filename)
        print()
        print(Seq)
        cnt += 1

        #Append Ref and Seq
        #Seq="IonCode_0101.fasta"
        c="Copy " + Ref + "+" + Seq + " ip.fasta > nul"     # '> nul' avoids system messages to display
        os.system(c)

        I="ip.fasta"        #Wuhan+Seq, Input to MAFFT
        O="op.fasta"        #MAFFT alignment Output

        #MAFFT Alignment
        mafft_cline=MafftCommandline(input=I)
        print(mafft_cline)
        stdout, stderr = mafft_cline()
        with open(O, 'w') as handle:
            handle.write(stdout)
            
        #Remove columns from op.fasta with gap in Wuhan sequence

        alignment = AlignIO.read(O, "fasta")

        OutputFile = "Aln_" + Seq       #Final output file for each sequence

        seqs=alignment[0]
        seqLen = len(seqs)

        edited=alignment[:,:0]

        c=0

        m="No"      #Insertions Yes/No in Seq i.e. Gaps Yes/No in Ref 

        for z in range(seqLen):
              if seqs[z] == '-':
                  c=c+1
                  if m == "No":
                      m="Yes"
              else:  
                  edited = edited + alignment[:,c:z+1]
                  c=z+1

        AlignIO.write(edited, OutputFile, "fasta")
        writer.writerow([cnt,Seq,m])

os.system("del ip.fasta")
os.system("del op.fasta")
os.system("Copy Aln*.fasta Final.fasta > nul")     #Final.fasta is the file to be considered
os.system("del Aln*.fasta")

print()       #Print blank line
print("No. of Sequences read = ",cnt)
print()

print()       #Print blank line
print("Final.fasta contains all the pairwise aligned sequences. Remove Reference sequence from this file!")
print()

print()       #Print blank line
print("Please refer GapInfo.csv to see information on insertions created in sequences during pairwise alignment with Reference sequence!")
print()

#Pause till any key pressed by the user
if sys.platform == 'win32':
    os.system('pause')
else:
    input('Press any key to continue...')

