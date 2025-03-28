#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import math


codon=["TTT","TTC","TTA","TTG","CTT","CTC","CTA","CTG","ATT","ATC",
       "ATA","ATG","GTT","GTC","GTA","GTG","TCT","TCC","TCA","TCG"
       ,"CCT","CCC","CCA","CCG","ACT","ACC","ACA","ACG","GCT","GCC","GCA","GCG"
       ,"TAT","TAC","TAA","TAG"
       ,"CAT","CAC","CAA","CAG","AAT","AAC","AAA","AAG","GAT","GAC","GAA","GAG"
       ,"TGT","TGC","TGA","TGG"
       ,"CGT","CGC","CGA","CGG","AGT","AGC","AGA","AGG","GGT","GGC","GGA","GGG"]

a_amine=["phe","phe","leu","leu","leu","leu","leu","leu","ile","ile",
       "ile","met","val","val","val","val","ser","ser","ser","ser"
       ,"pro","pro","pro","pro","thr","thr","thr","thr","ala","ala","ala","ala"
       ,"tyr","tyr","stop","stop"
       ,"his","his","gln","gln","asn","asn","lys","lys","asp","asp","glu","glu"
       ,"cys","cys","stop","trp"
       ,"arg","arg","arg","arg","ser","ser","arg","arg","gly","gly","gly","gly"]

bed="candidate.bed"
bed=open(bed, "r")
bed=bed.read()
bed=bed.split("\n")

dict_gene={}
gene_list=[]
for i in bed[1:]:
    if i!="":
        i=i.split("\t")
        dict_gene[i[3]]=[]
        gene_list.append(i[3])

transcrip=open("Alyrata_384_v2.1.transcript.fa","r")
transcrip=transcrip.read()
transcrip=transcrip.split('>')
for i in transcrip[1:]:
    if i!="":
        i=i.split("\n")
        name=i[0]
        name=name.split(".")
        name=name[0]
        if name in gene_list:
            seq=str(i[1])
            for j in i[2:]:
                if j!="":seq=seq+str(j)
            x=dict_gene[name]
            x.append(seq)
            dict_gene[i[3]]=x
            

fichier=open("biais_codon.csv", "w")
fichier.write("gene;codon;aa;count")

for i in gene_list:
    if i!="":
        x=dict_gene[i]
        codon_gene=[]
        for j in codon: codon_gene.append(0)
        for j in x:
            seq=str(j)
            while seq[0:3]!="ATG" and seq[0:3]!="":seq=str(seq[1:])
            while seq[0:3]!="TAA" and seq[0:3]!="TGA" and seq[0:3]!="TAG" and seq[0:3]!="":
                seq2=str(seq[0:3])
                for k in range(0,len(codon)):
                    if seq2==codon[k]:
                        x2=codon_gene[k]+1
                        codon_gene[k]=x2
                seq=str(seq[3:])
            seq2=str(seq[0:3])
            for k in range(0,len(codon)):
                if seq2==codon[k]:
                    x2=codon_gene[k]+1
                    codon_gene[k]=x2
        for j in range(0,len(codon)):
            fichier.write("\n"+str(i)+";"+str(codon[j])
                          +";"+str(a_amine[j])
                          +";"+str(codon_gene[j]))
fichier.close()
