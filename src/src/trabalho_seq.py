# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 15:36:54 2019

@author: Bela
"""

from abc import ABC

class sequencia(ABC):
    
    def __init__(self,sequencia):
        '''
        Cria uma nova sequencia

        '''
        self.sequencia=sequencia.upper()
        self.alfabeto='ACTG'
        self.propriedades={}
         
    def __eq__(self, obj):
        if isinstance(obj,sequencia):
            return self.sequencia == obj.sequencia
        else:
            return False
        
    def __str__(self):
        return self.sequencia
    
    def __getitem__(self,n):
        return self.sequencia[n]
    
    def __getslice__(self,i,j):
        return self.sequencia[i,j]
    
    def __len__(self):
        return len(self.sequencia)

    def add_prop(self,key,value):
        self.propriedades[key]=value
        
    def get_prop(self,key):
        if key not in self.propriedades.keys():
            raise IndexError()
        return self.propriedades[key]

    def getlenght(self):
        '''
        
        Imprime o comprimento da sequência
        '''
        return len(self.seq)
    
    def traduzCodao(self, cod):
        '''
        Dados os codões da sequencia de DNA faz a sua correspondencia
        com o aminoacido de acordo com o dicionario tc.

        '''
        tc = {"GCT": "A", "GCC": "A", "GCA": "A", "GCC": "A", "TGT": "C", "TGC": "C",
              "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TTT": "F", "TTC": "F",
              "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G", "CAT": "H", "CAC": "H",
              "ATA": "I", "ATT": "I", "ATC": "I",
              "AAA": "K", "AAG": "K",
              "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
              "ATG": "M", "AAT": "N", "AAC": "N",
              "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
              "CAA": "Q", "CAG": "Q",
              "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
              "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
              "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
              "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
              "TGG": "W",
              "TAT": "Y", "TAC": "Y",
              "TAA": "_", "TAG": "_", "TGA": "_"}
        if cod in tc:
            aa = tc[cod]
        else:
            aa = "X"  
        return aa
    
    def traducao(self, iniPos=0):
        '''
        Retorna
        -------
        A sequencia proteica que corresponde à sequencia de DNA

        '''
        seqM = self.sequencia
        seqAA = ""
        for pos in range(iniPos, len(seqM)-2, 3):
            cod = seqM[pos:pos+3]
            aa = self.traduzCodao(cod)
            seqAA = seqAA + aa
        return seqAA
    
    
    def valida(self):
        '''
        Valida a sequencia criada como uma sequencía de DNA
        '''
        for i in range(len(self.sequencia)):
            if self.sequencia[i]not in self.alfabeto:
                return ('Sequencia de DNA inválida')
        return ('Sequencia de DNA válida')
    
    def printsequencia(self):
        '''
        Imprime a sequência no terminal

        '''
        print(self.seq)
        
        
if __name__ == '__main__':
    seq = sequencia('ACGGGTCACACGT')
    print(seq.valida()) 
    print(seq.traducao())