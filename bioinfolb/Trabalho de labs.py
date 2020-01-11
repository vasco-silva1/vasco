# -*- coding: utf-8 -*-

from Bio import Entrez, SeqIO, SeqRecord, Seq
from Bio.Alphabet import IUPAC
from Bio.ExPASy import ScanProsite
from Bio import Medline
'''
Entrez.email="pg40958@alunos.uminho.pt"

handle=Entrez.esearch(db='nucleotide', term='Homo Sapiens[Orgn] AND FOXC2[Gene]',idtype='acc',retmax=283)
record=Entrez.read(handle)
handle.close()
​
print(record.keys())
print(record['Count'])
print(record['RetMax'])
​
idList= ','.join(record['IdList'][:5])
print(idList)
handle= Entrez.efetch(db='nucleotide',id=idList,retmode='xml')
records = Entrez.read(handle)
len(records)
print(records[0].keys())
print(records[0]['GBSeq_primary-accession'])
print(records[0]['GBSeq_other-seqids'])
print(records[0]['GBSeq_definition'])
print(records[0]['GBSeq_organism'])
'''
#querie com os nossos genes no organismo Homo Sapiens
handle=Entrez.egquery(term='Homo Sapiens[Orgn] AND FOXC2[Gene]')
record= Entrez.read(handle)
for row in record['eGQueryResult']:
    if row['DbName']=='nuccore':
        print(row['Count'])
        
#fazer search pelos abstract (esta a dar sempre o mesmo abstract)
handle = Entrez.esearch(db="pubmed", term="Homo Sapiens[Orgn] AND FOXC2[Gene]",retmode='text')
record= Entrez.read(handle)
id_lists=record['IdList']
handle.close()
handle = Entrez.efetch(db='pubmed',id=id_lists, rettype='Medline',retmode='text')
record= Medline.read(handle)
for rec in record:
    print(record['TI'])
    print()
    print(record['AB'])
    print()


handle=Entrez.esearch(db='nuccore', term='Homo Sapiens[Orgn] AND FOXC2[Gene]')
record= Entrez.read(handle)
gi_list=record['IdList']
print(gi_list)
gi_str=','.join(gi_list)
handle=Entrez.efetch(db='nuccore',id=gi_str,rettype='gb',retmode='text')
#text=handle.read()
#print(text)
##### Get -CDS-
records=SeqIO.parse(handle,'gb')
record=list(records)
print(record)
dic={}
record[0].annotations #Anotações
for rec in record:
    for feat in rec.features:
        if feat.type=='CDS':
            dic[rec.id]=feat.location.extract(rec).seq
print(dic['NM_005251.3'])

#tradução da sequencia
gene_teste=dic['NM_005251.3']            
print(gene_teste.translate())
a=gene_teste.translate()
    
# Get -ProteinACNumber-

handle=ScanProsite.scan(seq=a,db='sp')
result= ScanProsite.read(handle)
type(result)
result.n_match
len(result)
result.n_seq
res=result[2]
print(res['signature_ac'])
for i in result:
    print(i['signature_ac'])
    
file=open('teste.txt','w')
file.write(Entrez.efetch(db="gene", id=2303,retmode='xml',rettype='gb').read())
file.close()
### Nome do gene procurar no NCBI, bd  gene a sequencia
