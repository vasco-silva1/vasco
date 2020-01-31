# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 11:22:28 2020

@author: Bela
"""

from cmd import *
#from pyfiglet import Figlet
from trabalho_main import GestorDNA
import pickle
import sys

class GestorDNAShell(Cmd):
    intro = '''Interpretador de comandos para o gestor de sequências de DNA.\nEscreva 'start' para começar.
            '''
    #f=Figlet(font='slant')
    #print (f.renderText('Hello World'))
    prompt = 'GestorDNA> '
    
    def do_start(self,arg=None):
        while True:    
            print('''\n\n\n\n
                  --------Gestor de Sequências de DNA--------
                  
                  1.Criação ou importação de uma base de dados.
                  2.Manipulação de sequências.
                  3.Procura de padrões e similaridade.
                  4.Alinhamentos, Blast e Arvores Filogeneticas.
                  5.Alinhamento e Arvore Filogenetica com recurso a programas/packages externos às classes desenvolvidas nas aulas.
                  6.Tradução de sequência e tratamento de sequência proteica.
                  7.Guardar a instância do gestor.
                  8.Sair do programa.
                  
                  -------------------------------------------''')
            
            x=input('\n Escolha uma das opções apresentadas:')
            if x == "1":
                self.criacao()
            elif x == "2":
                self.manipulate()
            elif x == "3":
                self.search()
            elif x == "4":
                self.Align()
            elif x == "5":
                self.externos()
                pass
            elif x == "6":
                self.traducao()
                pass
            elif x == "7":
                self.saveDB()
            elif x == "8":
                sys.exit()
            else:
                print('\n Opção não válida, escolha outra vez')
                self.do_start()
            
                
            
    def criacao(self):
        while True:
            print('''\n\n\n\n
                  --------Gestor de Sequências de DNA--------
                  
                  1.Criar uma nova base de dados.
                  2.Importar uma base de dados criada anteriormente.
                  3.Adicionar sequências à base de dados manualmente.
                  4.Procurar sequências Online através de identificador de NCBI.
                  5.Importar sequência de ficheiro fasta.
                  6.Voltar ao menu anterior.
                  
                  -------------------------------------------''')
            x=input('\n Escolha uma das opções apresentadas:')
            if x=='1':
                self.create()
            elif x=='2':
                self.loadDB()
            elif x=='3':
                self.addSeq()
            elif x=='4':
                self.searchSeqsOnline()
            elif x=='5':
                self.importFasta()
            elif x=='6':
                self.do_start()
            else:
                print('\n Opção não válida, escolha outra vez')
           
    def create(self):
        try:
            name=input('Introduza um nome para a instância do gestor: ')
            eng.setName(name)                          
        except:
            print('Erro ao criar a instância do gestor.')
            
    def loadDB(self):
        try:
            name=input('Introduza o nome do ficheiro a carregar: ')
            eng.loaddic(name)
        except:
            print('Erro ao carregar a instância pretendida')
    
    def addSeq(self):
        try:
            ID=input('Introduza o identificador desejado: ').upper()
            seq=input('Introduza a sequência: ').upper()
            eng.addSeq(ID,seq)
        except:
            print('Erro ao adicionar sequência')
            
    def searchSeqsOnline(self):
        try:
            email=input('Introduza um email válido para efectuar a procurar na base de dados online: ')
            Ids= self.retrievelistaidsonline()
            eng.searchonline(Ids,email)
        except:
            print('Erro na procura online')   
    
    def importFasta(self):
        try:
            x=input('Introduza o nome do ficheiro fasta a carregar: ')
            filename=x+'.fasta'
            eng.import_fasta(filename)
        except:
            print('Erro ao carregar ficheiro fasta')

    def manipulate(self):
        while True:
             print('''\n\n\n\n
                  --------Gestor de Sequências de DNA--------
                  
                  1.Apresentar todas as sequências carregadas no gestor.
                  2.Apresentar informação sobre uma sequência
                  3.Procurar sequência por identificador.
                  4.Adicionar anotações a uma sequência.
                  5.Editar uma sequência.
                  6.Editar anotações de uma sequência.
                  7.Voltar ao menu anterior.
                  
                  -------------------------------------------''')
             x=input('\n Escolha uma das opções apresentadas:')
             if x=='1':
                 self.printGes()
             elif x=='2':
                 self.printSeq()
             elif x=='3':
                 self.searchid()
             elif x=='4':
                 self.addAnnotations()
             elif x=='5':
                 self.editSeq()
             elif x=='6':
                 self.editAnno()
             elif x=='7':
                 self.do_start()
             else:
                 print('\n Opção não válida, escolha outra vez')

    def printGes(self):
        try:
            eng.printDic()
        except:
            print('Erro ao apresentar o conteúdo do gestor')
    
    def printSeq(self):
        try:
            eng.RemindSeq()
            arg=input('Introduza o identificador da sequência: ').upper()
            if arg in eng.getIds():
                eng.printbyid(arg)
            else:
                print('Identificador não encontrado no gestor')
        except:
            print('Erro ao apresentar a sequência')

    def do_searchid(self):
        try:
            eng.RemindSeq()
            arg=input('Introduza o identificador: ').upper()
            if arg in eng.getIds():
                self.printSeq(arg)
            else:
                print('Identificador não encontrado no gestor')
        except:
            print('Erro ao procurar por identificador')
            
    def addAnno(self):
        try:
            eng.RemindSeq()
            ID=input('Introduza o identificador: ').upper()
            Anno=input('Introduza o texto da anotação: ')
            tipo=input('Introduza o tipo de anotação: ')
            if ID in eng.getIds():
                eng.addAnnotation(tipo,Anno,ID)
            else:
                print('Identificador não encontrado no gestor')
        except:
            print('Erro ao adicionar anotação')
            
    def editSeq(self):
        try:
            eng.RemindSeq()
            ID=input('Introduza o identificador: ').upper()
            seq=input('Introduza a nova sequência: ').upper()
            if ID in eng.getIds():
                eng.editSeq(ID,seq)
            else:
                print('Identificador não encontrado no gestor')
        except:
            print('Erro ao editar a sequência')
        
    def editAnno(self):
        try:
            eng.RemindSeq()
            ID=input('Introduza o identificador: ').upper()
            self.printSeq(ID)
            Anno=input('Introduza a nova anotação: ')
            tipo=input('Introduza o tipo de anotação a alterar: ')
            nr=input('Introduza o nr da anotação a alterar: ')
            if ID in eng.getIds():
                eng.editAnno(ID,tipo,nr,Anno)
            else:
                print('Identificador não encontrado no gestor')
        except:
            print('Erro ao editar a anotação')
            
    def search(self):
        while True:
            print('''\n\n\n\n
                  --------Gestor de Sequências de DNA--------
                  
                  1.Procurar sequência mais semelhante no gestor com uma sequência query.
                  2.Procurar um padrão ou subsequência contra uma sequência.
                  3.Procurar um padrão ou subsequência contra todas as sequências.
                  4.Procurar frequência de um símbolo ou padrão contra uma sequência.
                  5.Procurar frequência de um símbolo ou padrão contra todas as sequência.
                  6.Voltar ao menu anterior.
                  
                  -------------------------------------------''')
            x=input('\n Escolha uma das opções apresentadas:')
            if x=='1':
                self.search_seqs()
            elif x=='2':
                self.search_patt()
            elif x=='3':
                self.searc_pattDB()
            elif x=='4':
                self.search_freq()
            elif x=='5':
                self.search_freqDB()
            elif x=='6':
                self.do_start()
            else:
                print('\n Opção não válida, escolha outra vez')
    
    def search_seqs(self):
        try:
            query=input('Introduza a sequência query: ').upper()
            if eng.valSeq(query):
                eng.search_seqs(query)
            else:
                print('Sequência introduzida inválida')
        except:
            print('Erro ao procurar a sequência pretendida no gestor')
    
    def search_patt(self):
        try:
            eng.RemindSeq()
            Pat=input('Introduza o padrão a procurar como expressão regular: ')
            ID=input('Introduza o identificador: ')
            eng.search_pattern(Pat,ID)
        except:
            print('Erro ao procurar pelo padrão')
    
    def search_pattDB(self):
        try:
            Pat=input('Introduza o padrão a procurar como expressão regular: ')
            eng.search_patdb(Pat)
        except:
            print('Erro ao procurar pelo padrão')
    
    def search_freq(self):
        try:
            eng.RemindSeq()
            query=input('Introduza o simbolo ou subsequência: ')
            ID=input('Introduza o identificador: ')
            eng.search_freq(query,ID)
        except:
            print('Erro ao calcular a frequência')
    
    def search_freqDB(self):
        try:
            query=input('Introduza o simbolo ou subsequência: ')
            eng.search_freq(query)
        except:
            print('Erro ao calcular a frequência')
    
    def Align(self):
        while True:
            print('''\n\n\n\n
                  --------Gestor de Sequências de DNA--------
                  
                  1.Realizar alinhamento múltiplo.
                  2.Realizar Blast Online através do portal do NCBI.
                  3.Adicionar as sequências resultantes do Blast ao gestor.
                  4.Criar árvore filogenética usando UPGMA
                  5.Voltar ao menu anterior.
                  
                  -------------------------------------------''')
            x=input('\n Escolha uma das opções apresentadas:')
            if x=='1':
                pass
            elif x=='2':
                self.blast()
            elif x=='3':
                self.save_blast()
            elif x=='4':
                self.create_tree()
            elif x=='5':
                self.do_start()
            else:
                print('\n Opção não válida, escolha outra vez')
                
    def MultiAlign(self):
        try:
            eng.RemindSeq()
            Ids= self.retrievelistaseqs()
            match=int(input())
            mismatch=int(input())
            gap=int(input())
            eng.MultiAli(Ids,match,mismatch,gap)
        except:
            print('Erro ao correr o alinhamento múltiplo')
                
    def blast(self):
        try:
            eng.RemindSeq()
            ID=input('Introduza o identificador: ')
            db=input('Introduza a base de dados em que pretende correr o blast')
            eng.blast(ID,base=db)
        except:
            print('Erro ao realizar o blast online')
    
    def save_blast(self):
        try:
            print(eng.get_blast_hits)
            self.retrievelistaseqs()
            eng.addBlastSeqs()
        except:
            print('Erro ao adicionar sequências do blast')
    
    def create_tree(self):
        try:
            eng.arvore()
        except:
            print('Erro ao contruir a árvore')
    
    
    def externos(self):
        while True:
            print('''\n\n\n\n
                  --------Gestor de Sequências de DNA--------
                  
                  1.Realizar alinhamento múltiplo com recurso ao ClustalW
                  2.Criar árvore com recurso ao package Phylo do BioPython
                  3.Voltar ao menu anterior.
                  
                  -------------------------------------------''')
            x=input('\n Escolha uma das opções apresentadas:')
            if x=='1':
                self.AlignClustal()
            elif x=='2':
                self.CreateTreePhylo()
            elif x=='3':
                self.do_start()
            else:
                print('\n Opção não válida, escolha outra vez')
                
    def AlignClustal(self):
        try:
            eng.RemindSeq()
            nome= input('Introduza um nome para o ficheiros criados pelo ClustalW: ')
            Ids= self.retrievelistaseqs()
            if Ids == -1:
                eng.make_seqs_ma_exe(nome)
                eng.multiple_alignment(nome)
            else:
                eng.make_seqs_ma_exe(nome,Ids)
                eng.multiple_alignment(nome)
        except:
            print('Erro ao executar o alinhamento multiplo externo')
            
    def CreateTreePhylo(self):
        try:
            eng.create_tree()
        except:
            print('Erro ao criar árvore filogenética utilizando o Phylo')
            

    def retrievelistaseqs(self):
        i=0
        listaseqs=[]
        while True:
            i+=1
            item=input('''
                       Introduza o id da %dª sequência: 
                       Se pretender utilizar todas as sequências gravadas escreva "all"
                       '''%i).upper()
            if item=='':
                break
            elif item=='ALL':
                listaseqs=-1
            listaseqs.append(item)
        print("As sequências inseridas foram",listaseqs)
        
        return listaseqs
    
    def retrievelistaidsonline(self):
        i=0
        listaids=[]
        while True:
            i+=1
            item=input('Introduza o %dº identificador: '%i).upper()
            if item=='':
                break
            listaids.append(item)
        print("Os identificadores inseridos foram",listaids)
        
        return listaids
                
    def traducao(self):
        while True:
            print('''\n\n\n\n
                  --------Gestor de Sequências de DNA--------
                  
                  1.Efectuar tradução de sequência do gestor.
                  2.Realizar Blast protein Online através do portal do NCBI.
                  3.Gravar ficheiro das proteinas resultantes do blast.
                  4.Apresentar as orfs
                  5.Gravar orfs
                  6.Voltar ao menu anterior.
                  
                  -------------------------------------------''')
            x=input('\n Escolha uma das opções apresentadas:')
            if x=='1':
                self.traducao_seq()
            elif x=='2':
                self.blastp()
            elif x=='3':
                self.save_blast()
            elif x=='4':
                self.ORFS()
            elif x=='5':
                self.saveORFS
            elif x=='6':
                self.do_start()
            else:
                print('\n Opção não válida, escolha outra vez')
                
    def traducao_seq(self):
        try:
            eng.RemindSeq()
            ID=input().upper()
            eng.traducao(ID)
        except:
            print('Erro ao traduzir a sequência')

    def ORFS(self):
        try:
            eng.RemindSeq()
            Ids= self.retrievelistaseqs()
            eng.orfs(Ids)
        except:
            print("Erro na determinação das ORFS")
            
    def saveORFS(self):
        try:
            eng.saveOrfs()
        except:
            print('Erro ao gravar ')
        
            
    def blast(self):
        try:
            ID=input('Introduza o identificador: ')
            eng.blast(ID)
        except:
            print('Erro ao realizar o blast online')
    
    def save_blastp(self):
        try:
            print(eng.get_blastp_hits)
            ID=self.retrievelistaseqs()
            eng.saveBlastpSeqs(ID)
        except:
            print('Erro ao adicionar sequências do blast')
        
    def saveDB(self):
        try:
            eng.save_dic()
        except:
            print('Erro ao gravar a instância do gestor')
    
    def do_sair(self, arg=None):
        "Sair do programa Gestor de sequências de DNA: sair"
        print('Obrigado por ter utilizado o Gestor de sequências de DNA')
        global janela  # pois pretendo atribuir um valor a um identificador global
        if janela is not None:
                    del janela  # invoca o metodo destruidor de instancia __del__()
        return True

if __name__ == '__main__':
    eng = GestorDNA()
    print(eng.getCWD())
    janela = None
    sh = GestorDNAShell()
    sh.cmdloop()
