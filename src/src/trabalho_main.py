# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 15:36:34 2019

@author: Bela
"""

from MySeq import MySeq
from MyBlast import MyBlast 
from upgma import UPGMA
from AlignSeq import AlignSeq
from MultipleAlign import MultipleAlign
from SubstMatrix import SubstMatrix
from Bio import Entrez, SeqIO, Phylo
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import ClustalwCommandline
import pickle
import re
import os



class GestorDNA:
   
    def __init__(self):
        self.dic={}
        self.seqs_alfabeto='ACTG'
        self.seqs_tipo='dna'
        self.name=None
        self.seq=None
        self.listaIds=[]
        self.blast_ids=[]
        self.cwd=os.getcwd().strip('src')
        self.file=None
        self.blastp_hits=[]
        self.blast_ids

    def getCWD(self):
        return self.cwd

    def getIds(self):
        return self.listaIds
    
    def setName(self,nome):
        self.name=nome
    
    def loaddic(self,filename):
        file=open(filename+'.txt','rb')
        self.name=filename
        self.dic=pickle.load(file)
        self.listIds()
        file.close()
    
    def addSeq(self,ident,seq):
        Id=ident
        sequen=seq
        if self.valSeq(seq):
            self.listaIds.append(Id)
            if Id not in self.dic.keys():
                self.dic[Id]={'Sequência':sequen}
                self.dic[Id]['Anotações']={}
            else:
                print('Identificador já utilizado')
        else:
            print('Sequência DNA inválida')
            
    def addAnnotation(self,tipo,anno,Id):
        if tipo not in self.dic[Id]['Anotações'].keys():
            self.dic[Id]['Anotações'][tipo]=[]
            self.dic[Id]['Anotações'][tipo].append(anno)
        else:
            self.dic[Id]['Anotações'][tipo].append(anno)
            
    def valSeq(self,seq):
        vs=MySeq(seq,'dna')
        if vs.validaER():
            return True
        else:
            return False
        
    def printbyid(self,Id):
        for key,val in self.dic[Id].items():
            print (key+':\n')
            if type(val)==dict:
                self.__printSubDic(val)
            elif type(val)==list:
                i=0
                for x in range(len(val)):
                    print(str(i)+'->'+val[x])
                    i+=1
                print()
            else:
                print(self.dic[Id][key])
            
            if type(val)==str:
                seq=MySeq(val,self.seqs_tipo)
                gc=self.conteudoGC(val)
                mostfreq=self.mostfreq(val)
                print('Comprimento: ', str(len(seq.seq)))
                print('Conteúdo G-C: ', gc)
                print('Frequência de cada Base Azotada: ' , mostfreq)
                
                
    def conteudoGC(self,mseq):
        G=mseq.count("G")
        C=mseq.count("C")
        conteudo=G+C
        conteudo=conteudo/len(mseq)*100
        #print("O conteúdo GC (percentagem) é :", conteudo)
        return conteudo

    def mostfreq(self,mseq):
        G=(mseq.count("G")/len(mseq))*100
        C=(mseq.count("C")/len(mseq))*100
        A=(mseq.count("A")/len(mseq))*100
        T=(mseq.count("T")/len(mseq))*100
        #print("A : ",A ,"T : ",T ,"G : ",G ,"C : ", C)
        return (("A : ",A ),("T : ",T ),("G : ",G ),("C : ", C))
        
    def printDic(self,dic=None):
        if dic==None:
            dic=self.dic
        else:
            dic=dic
        l=self.listaIds
        tid=None
        for x in range(len(l)):
            tid=l[x]
            print('>'+tid)
            for key,val in dic[tid].items():
                print (key+':')
                if type(val)==dict:
                    self.__printSubDic(val)
                elif type(val)==list:
                    i=0
                    for x in range(len(val)):
                        print(str(i)+'->'+val[x])
                        i+=1
                else:
                    print(dic[tid][key])
                    print()
                    
    def __printSubDic(self,subdic):
        for key, val in subdic.items():
            print(key+':')
            if type(val)==dict:
                self.printSubDic(val)
            elif type(val)==list:
                i=0
                for x in range(len(val)):
                    print(str(i)+'->'+val[x])
                    i+=1
            else:
                print(subdic[key][val])
                print()
        
    def listIds(self):
        for key in self.dic.keys():
            self.listaIds.append(key)
                
    def editAnnotations(self,Id,tipo,nr,new):
        self.dic[Id]['Anotações'][tipo][nr]=new
        
    def editSeq(self,Id,new):
        self.dic[Id]['Sequência']=new
        
    def search_seqs(self,query,match,mismatch,gap):
        db=[]
        sm=SubstMatrix()
        sm.createFromMatchPars(match,mismatch,"ATGC")
        aseq=AlignSeq(sm,gap)
        seq1=query
        maximo=-999
        imax=None
        for i in range(len(db)):
            seq2=db[i]
            sc=aseq.needlemanWunsch(seq1,seq2)
            if sc>maximo:
                maximo=sc
                imax=i
        else:
            pass
        
        idmax=db[imax]
        
        for Id in self.listaIds:
            if self.dic[Id]['Sequência']==idmax:
                print('Melhor match encontrado na sequência com Id %s' %Id)
    
    def search_pattern(self,pat,Id):
        seq=self.dic[Id]['Sequência']
        matches = re.finditer(pat,seq)
        res=[]
        for match in matches:
            res.append(match.span())
        if len(res)!=0:
            print('Foram encontradas %d ocorrências do padrão na sequência %s' % (len(res),Id))
            i=1
            for match in res:
                print('%d ª ocorrência, posição: '%i,match)
                i+=1
        else:
            print('Não foram encontradas ocorrências do padrâo na sequência %s' %Id)
            return None

    def search_patdb(self,pat):
        for Id in self.listaIds:
            self.search_pattern(pat,Id)

    
    def search_freq(self,sim=None,Id=None):
        if Id != None:
            seq=self.dic[Id]['Sequência']
            res=[]
            matches=re.finditer(sim,seq)
            for match in matches:
                res.append(match.span())
            freq=(len(res))
            print('O simbolo apresenta uma frequência de %f \%' %freq)
        else:
            for Id in self.listaIds:
                seq=self.dic[Id]['Sequência']
                res=[]
                matches=re.finditer(sim,seq)
                for match in matches:
                    res.append(match.span())
                freq=(len(res))
                print('O simbolo apresenta uma frequência de %f \%, na sequência %s' %(freq,Id))
        
        
    def arvore(self):
        bd=[]
        for Id in self.listaIds:
            bd.append(MySeq(self.dic[Id]['Sequência'],'dna'))
        sm = SubstMatrix()
        sm.createFromMatchPars(3, -1, self.seqs_alfabeto)
        alin=AlignSeq(sm,-1)
        up=UPGMA(bd,alin)
        arv=up.run()
        arv.printtree()
        
    def traducao(self,Id):
        seq=self.dic[Id]['Sequência']
        ms=MySeq(seq,self.seqs_tipo)
        res=ms.traduzSeq()
        print(res)
    
    def orfs(self,Ids):
        listaids=[]
        for id in Ids:
            listaids.append(self.dic[id]['Sequência'])
        seq=self.dic[id]['Sequência']
        ms=MySeq(seq,self.seqs_tipo)
        res=ms.orfs()
        for x in res:
            x.printseq()
    
    def saveOrfs(self,Id):
        try:
            seq=self.dic[Id]['Sequência']
            ms=MySeq(seq,self.seqs_tipo)
            res=ms.orfs()
            file=open(self.name +'Orfs' + str(Id) + '.txt','w')
            for orf in res:
                pickle.dump(orf.seq,file)
            file.close()
        except:
            print("erro ao gravar orfs")
            

            
            
    def searchonline(self,Id,email):
        try:
            Entrez.email=email
            if type(Id)==str:
                handle=Entrez.efetch(db='nucleotide',id=Id,rettype='fasta',retmode='text')
                record=SeqIO.read(handle,'fasta')
                handle.close()
                self.addSeq(Id,str(record.seq))
            elif type(Id)==list:
                for id_ in Id:
                    try:
                        handle=Entrez.efetch(db='nucleotide',id=id_,rettype='fasta',retmode='text')
                        record=SeqIO.read(handle,'fasta')
                        self.addSeq(id_,str(record.seq))
                        handle.close()
                    except:
                        print('%s não encontrado na base de dados'%id_)
        except:
            print('Erro ao executar procura online')
            
    def blast(self,ID,program='blastn',base='nt',evalue=0.05,align=50,matrix='BLOSUM62'):
        '''
        Função que corre o blast de sequencias de dna contra a base de dados pretendida, tendo como parametros default
        o programa-blastn
        a matrix- BLOSUM62
        o e-value=0.001
        o nº de alinhamentos maximos=100
        '''
        self.blast_ids=[]
        seq=seq=self.dic[ID]['Sequência']
        result_handle=NCBIWWW.qblast(program,'nt',seq) 
        hit_id=[]
        blast_records = NCBIXML.parse(result_handle)
        e_value_thresh = 0.001
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < e_value_thresh:
                        hit_id.append(alignment.hit_id)
                        print("****Alignment****")
                        print("sequence:", alignment.title)
                        print("length:", alignment.length)
                        print("e value:", hsp.expect)
                        print(hsp.query[0:75] + "...")
                        print(hsp.match[0:75] + "...")
                        print(hsp.sbjct[0:75] + "...")
        self.blast_ids=hit_id
        for i in range(len(self.blast_ids)):
            print(self.blast_ids[i])
                     
        
    def get_blast_hits(self):
        return self.blast_ids
    
    def addBlastSeqs(self,ids):
        hit_id = list(set(ids))                  
        for ind in hit_id:               
            handle = Entrez.efetch(db ="nucleotide", rettype ='fasta', retmode ="fasta", id=ind)
            seq_record = SeqIO.read(handle, 'fasta')
            self.addSeq(ind,str(seq_record.seq))
            handle.close()
            
    def blastp(self,ID,program='blastn',base='swissprot',evalue=0.05,align=50,matrix='BLOSUM62'):
        '''
        Função que corre o blast de sequencias de dna contra a base de dados pretendida, tendo como parametros default
        o programa-blastn
        a matrix- BLOSUM62
        o e-value=0.001
        o nº de alinhamentos maximos=100
        '''
        self.blastp_ids=[]
        seq=self.dic[ID]['Sequência']
        result_handle=NCBIWWW.qblast(program,'nt',seq) 
        hit_id=[]
        blast_records = NCBIXML.parse(result_handle)
        e_value_thresh = 0.001
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < e_value_thresh:
                        hit_id.append(alignment.hit_id)
                        print("****Alignment****")
                        print("sequence:", alignment.title)
                        print("length:", alignment.length)
                        print("e value:", hsp.expect)
                        print(hsp.query[0:75] + "...")
                        print(hsp.match[0:75] + "...")
                        print(hsp.sbjct[0:75] + "...")
        self.blastp_ids=hit_id
                     
        
    def get_blastp_hits(self):
        return self.blast_ids
    
    def saveBlastpSeqs(self,ids):
        hit_id = list(set(ids))      
        cwd=self.cwd+'Proteínas\\'            
        for ind in hit_id:               
            file=open(cwd+ind+'GBfile.txt','w')                    
            file.write(Entrez.efetch(db='swissprot',id=ind,rettype='gb',retmode='text').read())                    
            file.close()
        
    def make_seqs_ma_exe(self,filename,Ids=None):
        try:
            file=open(self.cwd+'AlinMuExe\\'+filename+'.fasta','w')
            if Ids==None:
                for id_ in self.listaIds:
                    file.write('>'+id_+'\n')
                    file.write(self.dic[id_]['Sequência']+'\n')
            else:
                for id_ in Ids:
                    file.write('>'+id_+'\n')
                    file.write(self.dic[id_]['Sequência']+'\n')
            file.close()
        except:
            os.mkdir(self.cwd+'AlinMuExe\\')
            self.make_seqs_ma_exe()
    
    def multiple_alignment(self,filename):
        '''
        Função que realiza o alinhamento multiplo com as melhores sequencias 
        do blast corrido anteriormente, e escolhidas pela função anterior 
        '''
        cmdline=ClustalwCommandline(r'C:\Program Files (x86)\ClustalW2\clustalw2',infile=self.cwd+'AlinMuExe\\'+filename+'.fasta')
        cmdline()
        
    def create_tree(self,filename):
        '''
        Função que cria a árvore filogenética com base no alinhamento múltiplo
        '''
        infile=(self.cwd+'AlinMuExe\\'+filename+'.dnd')
        tree = Phylo.read(infile, "newick")
        Phylo.draw_ascii(tree)
        return tree
    
    def RemindSeq(self):
        for ID in self.listaIds:
            print(ID+ '->'+ self.dic[ID]['Sequência'])
            
    def MultiAli(self,Seqs,match=1,mismatch=-1,gap=-1):
        listaseqs=[]
        for seq in Seqs:
            listaseqs.append(self.dic[seq]['Sequência'])
        sm=SubstMatrix()
        sm.createFromMatchPars(match,mismatch,"ATGC")
        aseq=AlignSeq(sm,gap)
        mult=MultipleAlign(listaseqs,aseq)
        resultado=mult.alignConsensus()
        return resultado
    
    def import_fasta(self,filename):
        fasta_sequences = SeqIO.parse(open(filename),'fasta')
        for seq in fasta_sequences:
            if self.addSeq(seq.id,str(seq.seq).upper()):
                self.addAnnotation('Nome',seq.name,seq.id)
        
    def save_dic(self):
        file=open(self.name+'.txt','wb')
        pickle.dump(self.dic,file)
        file.close()

"""
    
    def save_dic(self):
        cwd=self.cwd+'databases\\'
        file=open(cwd+self.name+'.txt','wb')
        pickle.dump(self.dic,file)
        file.close()
        
    """