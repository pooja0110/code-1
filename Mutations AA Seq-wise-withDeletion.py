#Data: Full length (29903 nt) alignment of COVID-19 sequences including first Wuhan seq as reference seq
#Program lists out the number of amino acid mutations and actual mutations in each sequence w.r.t. first Wuhan seq
#Ignores "-" as substitution

import sys
import csv
from Bio import AlignIO

InputFile = input("Enter sequence alignment filename: ")
alignment = AlignIO.read(InputFile, "fasta")

m=len(alignment)    #no. of seq in alignment, including first WUHAN as Ref_seq

# DNA codon table as Dictionary (keys and values)
codon   = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
           "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
           "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
           "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
           "TCT" : "S", "CCT" : "P", "ACT" : "T", "GCT" : "A",
           "TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "TCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
           "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "TAA" : "*", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "TAG" : "*", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
           "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "TGA" : "*", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
           }


nt="ACGT"                      #string of nucleotides



#START computation...

Tot_sub=0                      #Total no. of substitutions in entire genome



G=["nsp1","nsp2","nsp3","nsp4","nsp5","nsp6","nsp7","nsp8","nsp9","nsp10","nsp11","nsp13","nsp14","nsp15","nsp16"]       #names of genes
Gs=[266,806,2720,8555,10055,10973,11843,12092,12686,13025,13442,16237,18040,19621,20659]                                 #start of genes
Ge=[805,2719,8554,10054,10972,11842,12091,12685,13024,13441,13480,18039,19620,20658,21552]                               #end of genes, no stop codon
GL=len(G)                                                                                                                #length of G
F=["01","02","03","04","05","06","07","08","09","10","11","13","14","15","16"]                                           #File counter


for w in range(0,GL):

    gene=G[w]
    OutputFile = F[w]+"-"+gene+"-2.csv"
    s=Gs[w]                        #start of gene
    e=Ge[w]                        #end of gene, no stop codon
    sub_aln=alignment[:,s-1:e]     #slice of alignment from s to e, length=e-s+1
    n=len(sub_aln[0])              #nucleotide length of sub_aln

    Ref_seq=sub_aln[0]    #Wuhan reference seq in sub_aln

    print("")
    print("Gene = ",gene,s,e)    #print message on screen

    with open(OutputFile, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Sequence id", "#Mutations_"+gene,"Mutations_"+gene,"Deletions_"+gene])

        for j in range(1,m):     #Sequence j

             Mutations=""       #string of mutations
             Deletions=""       #string of deletions
             cnt=0              #count for seq-wise no. of mutations

             for i in range(0, n-(3+n%3)+1, 3):    #position i

                 Ref_aa=codon[ Ref_seq[i].upper()+Ref_seq[i+1].upper()+Ref_seq[i+2].upper() ]

                 seq=sub_aln[j]                        #actual seq no j

                 t1=seq[i].upper()                     #verify each nucleotide in codon if it is among A,C,G,T otherwise replace by "-"
                 if nt.find(t1) == -1:
                     t1="-"

                 t2=seq[i+1].upper()                   #verify each nucleotide in codon if it is among A,C,G,T otherwise replace by "-"
                 if nt.find(t2) == -1:
                     t2="-"

                 t3=seq[i+2].upper()                   #verify each nucleotide in codon if it is among A,C,G,T otherwise replace by "-"
                 if nt.find(t3) == -1:
                     t3="-"

                 aa=t1+t2+t3                           #nucleotides forming a codon in seq no j+1
 
                 if aa.find("-") >= 0:
                     aa="del"                          #amino acid="del" if any of three nucleotides is "-"
                 else:
                     aa=codon[aa]                      #amino acid when all nucleotides are present

                 if Ref_aa != aa and aa!="del":        #amino acid substitution found, except deletion
                     Tot_sub=Tot_sub+1
                     cnt=cnt+1
                     Mutations = Mutations + Ref_aa + str(int(i/3+1)) + aa + ","

                 if aa=="del":                         #Deletion found
                     Deletions = Deletions + Ref_aa + str(int(i/3+1)) + aa + ","

             writer.writerow([seq.id,cnt,Mutations,Deletions])




gene="nsp12_RdRp"
OutputFile = "12-"+gene+"-2.csv"
s1=13442                                                #start of gene
e1=13468                                                #end of gene, no stop codon
s2=13468                                                #start of gene
e2=16236                                                #end of gene, no stop codon
sub_aln=alignment[:,s1-1:e1] + alignment[:,s2-1:e2]     #slice of alignment from s1 to e1 and s2 to e2
n=len(sub_aln[0])                                       #nucleotide length of sub_aln

Ref_seq=sub_aln[0]    #Wuhan reference seq in sub_aln

print("")
print("Gene = ",gene,"  Join ",s1,"....",e1," & ",s2,"....",e2)      #print message on screen

with open(OutputFile, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Sequence id", "#Mutations_"+gene,"Mutations_"+gene,"Deletions_"+gene])
    
    for j in range(1,m):     #Sequence j

         Mutations=""       #string of mutations
         Deletions=""       #string of deletions
         cnt=0              #count for seq-wise no. of mutations

         for i in range(0, n-(3+n%3)+1, 3):    #position i

             Ref_aa=codon[ Ref_seq[i].upper()+Ref_seq[i+1].upper()+Ref_seq[i+2].upper() ]

             seq=sub_aln[j]                        #actual seq no j

             t1=seq[i].upper()                     #verify each nucleotide in codon if it is among A,C,G,T otherwise replace by "-"
             if nt.find(t1) == -1:
                 t1="-"

             t2=seq[i+1].upper()                   #verify each nucleotide in codon if it is among A,C,G,T otherwise replace by "-"
             if nt.find(t2) == -1:
                 t2="-"

             t3=seq[i+2].upper()                   #verify each nucleotide in codon if it is among A,C,G,T otherwise replace by "-"
             if nt.find(t3) == -1:
                 t3="-"

             aa=t1+t2+t3                           #nucleotides forming a codon in seq no j+1
 
             if aa.find("-") >= 0:
                 aa="del"                          #amino acid="del" if any of three nucleotides is "-"
             else:
                 aa=codon[aa]                      #amino acid when all nucleotides are present

             if Ref_aa != aa and aa!="del":        #amino acid substitution found, except deletion
                 Tot_sub=Tot_sub+1
                 cnt=cnt+1
                 Mutations = Mutations + Ref_aa + str(int(i/3+1)) + aa + ","

             if aa=="del":                         #Deletion found
                 Deletions = Deletions + Ref_aa + str(int(i/3+1)) + aa + ","
                     
         writer.writerow([seq.id,cnt,Mutations,Deletions])








G=["Spike","ORF3a","ORF4_E","ORF5_M","ORF6","ORF7a","ORF7b","ORF8","ORF9_N","ORF10"]
Gs=[21563,25393,26245,26523,27202,27394,27756,27894,28274,29558]
Ge=[25384,26220,26472,27191,27387,27759,27887,28259,29533,29674]            #includes stop codon
GL=len(G)
F=["17","18","19","20","21","22","23","24","25","26"]


for w in range(0,GL):

    gene=G[w]
    OutputFile = F[w]+"-"+gene+"-2.csv"
    s=Gs[w]                        #start of gene
    e=Ge[w]-3                      #end of gene, includes stop codon
    sub_aln=alignment[:,s-1:e]     #slice of alignment from s to e, length=e-s+1
    n=len(sub_aln[0])              #nucleotide length of sub_aln

    Ref_seq=sub_aln[0]    #Wuhan reference seq in sub_aln

    print("")
    print("Gene = ",gene,s,e)    #print message on screen

    with open(OutputFile, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Sequence id", "#Mutations_"+gene,"Mutations_"+gene,"Deletions_"+gene])
       
        for j in range(1,m):     #Sequence j

             Mutations=""       #string of mutations
             Deletions=""       #string of deletions
             cnt=0              #count for seq-wise no. of mutations

             for i in range(0, n-(3+n%3)+1, 3):    #position i

                 Ref_aa=codon[ Ref_seq[i].upper()+Ref_seq[i+1].upper()+Ref_seq[i+2].upper() ]

                 seq=sub_aln[j]                        #actual seq no j

                 t1=seq[i].upper()                     #verify each nucleotide in codon if it is among A,C,G,T otherwise replace by "-"
                 if nt.find(t1) == -1:
                     t1="-"

                 t2=seq[i+1].upper()                   #verify each nucleotide in codon if it is among A,C,G,T otherwise replace by "-"
                 if nt.find(t2) == -1:
                     t2="-"

                 t3=seq[i+2].upper()                   #verify each nucleotide in codon if it is among A,C,G,T otherwise replace by "-"
                 if nt.find(t3) == -1:
                     t3="-"

                 aa=t1+t2+t3                           #nucleotides forming a codon in seq no j+1
 
                 if aa.find("-") >= 0:
                     aa="del"                          #amino acid="del" if any of three nucleotides is "-"
                 else:
                     aa=codon[aa]                      #amino acid when all nucleotides are present

                 if Ref_aa != aa and aa!="del":        #amino acid substitution found, except deletion
                     Tot_sub=Tot_sub+1
                     cnt=cnt+1
                     Mutations = Mutations + Ref_aa + str(int(i/3+1)) + aa + ","

                 if aa=="del":                         #Deletion found
                     Deletions = Deletions + Ref_aa + str(int(i/3+1)) + aa + ","
                     
             writer.writerow([seq.id,cnt,Mutations,Deletions])





#Write Total no. of substitutions in separate file

with open("Tot_Sub-2.txt", 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Total no. of substitutions: "+str(Tot_sub)])

