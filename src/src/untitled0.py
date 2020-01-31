# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 00:16:37 2020

@author: Bela
"""

from cmd import *
#from pyfiglet import Figlet
from trabalho_main import GestorDNA

def do_start(arg=None):
    while True:    
        print('''\n\n\n\n
              --------Gestor de Sequências de DNA--------
              
              1.Criação ou importação de uma base de dados.
              2.Manipulação de sequências.
              3.Procura de padrões e similaridade.
              4.Alinhamentos, Blast e Arvores Filogeneticas.
              5.
              6.
              7.Guardar a instância do gestor.
              8.Sair do programa.
              
              -------------------------------------------''')
        x=input('\n Escolha uma das opções apresentadas:')
        if x == '1':
            criação()
        elif x == '2':
            manipulate()
        elif x == '3':
            search()
        elif x == '4':
            Align()
        elif x == '5':
            menu5()
        elif x == '6':
            menu()
        elif x == '8':
            print("Para sair escreva 'sair'")
        else:
            print('\n Opção não válida, escolha outra vez')
            do_start()
        
def criação():
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
        x=int(input('\n Escolha uma das opções apresentadas:'))
        if x==1:
            create()
        elif x==2:
            loadDB()
        elif x==3:
            addSeq()
        elif x==4:
            searchSeqsOnline()
        elif x==5:
            importFasta()
        elif x==6:
            do_start()
        else:
            print('\n Opção não válida, escolha outra vez')
       
def create():
    try:
        name=input('Introduza um nome para a instância do gestor: ')
        eng.setName(name)                          
    except:
        print('Erro ao criar a instância do gestor.')
        
def loadBD():
    try:
        name=input('Introduza o nome do ficheiro .txt a carregar: ')
        eng.loaddic(name)
    except:
        print('Erro ao carregar a instância pretendida')

def addSeq():
    try:
        ID=input('Introduza o identificador desejado: ').upper()
        seq=input('Introduza a sequência: ').upper()
        eng.addSeq(ID,seq)
    except:
        print('Erro ao adicionar sequência')
        
def searchSeqsOnline():
    try:
        email=input('Introduza um email válido para efectuar a procurar na base de dados online: ')
        ids=input()('Introduza o(s) identificador(es) da(s) sequência(s) a procurar').split(' ')
        eng.searchonline(ids,email)
    except:
        print('Erro na procura online')   

def importFasta():
    pass

def manipulate():
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
         x=int(input('\n Escolha uma das opções apresentadas:'))
         if x==1:
             printGes()
         elif x==2:
             printSeq()
         elif x==3:
             searchid()
         elif x==4:
             addAnnotations()
         elif x==5:
             editSeq()
         elif x==6:
             editAnno()
         elif x==7:
             do_start()
         else:
             print('\n Opção não válida, escolha outra vez')

def printGes():
    try:
        eng.printDic()
    except:
        pass

def printSeq():
    try:
        arg=input('Introduza o identificador da sequência: ').upper()
        if arg in eng.getIds():
            eng.printbyid(arg)
        else:
            print('Identificador não encontrado no gestor')
    except:
        print('Erro ao apresentar a sequência')

def do_searchid():
    try:
        arg=input('Introduza o identificador: ').upper()
        if arg in eng.getIds():
            printSeq(arg)
        else:
            print('Identificador não encontrado no gestor')
    except:
        print('Erro ao procurar por identificador')
        
def addAnno():
    try:
        ID=input('Introduza o identificador: ').upper()
        Anno=input('Introduza o texto da anotação: ')
        tipo=input('Introduza o tipo de anotação: ')
        if ID in eng.getIds():
            eng.addAnnotation(tipo,Anno,ID)
        else:
            print('Identificador não encontrado no gestor')
    except:
        print('Erro ao adicionar anotação')
        
def editSeq():
    try:
        ID=input('Introduza o identificador: ').upper()
        seq=input('Introduza a nova sequência: ').upper()
        if ID in eng.getIds():
            eng.editSeq(ID,seq)
        else:
            print('Identificador não encontrado no gestor')
    except:
        print('Erro ao editar a sequência')
    
def editAnno():
    try:
        ID=input('Introduza o identificador: ').upper()
        printSeq(ID)
        Anno=input('Introduza a nova anotação: ')
        tipo=input('Introduza o tipo de anotação a alterar: ')
        nr=input('Introduza o nr da anotação a alterar: ')
        if ID in eng.getIds():
            eng.editAnno(ID,tipo,nr,Anno)
        else:
            print('Identificador não encontrado no gestor')
    except:
        print('Erro ao editar a anotação')
        
def search():
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
        x=int(input('\n Escolha uma das opções apresentadas:'))
        if x==1:
            search_seqs()
        elif x==2:
            search_patt()
        elif x==3:
            searc_pattDB()
        elif x==4:
            search_freq()
        elif x==5:
            search_freqDB()
        elif x==6:
            do_start()
        else:
            print('\n Opção não válida, escolha outra vez')

def search_seqs():
    try:
        query=input('Introduza a sequência query: ').upper()
        if eng.valSeq(query):
            eng.search_seqs(query)
        else:
            print('Sequência introduzida inválida')
    except:
        print('Erro ao procurar a sequência pretendida no gestor')

def search_patt():
    try:
        Pat=input('Introduza o padrão a procurar como expressão regular: ')
        ID=input('Introduza o identificador: ')
        eng.search_pattern(Pat,ID)
    except:
        print('Erro ao procurar pelo padrão')

def search_pattDB():
    try:
        Pat=input('Introduza o padrão a procurar como expressão regular: ')
        eng.search_patdb(Pat)
    except:
        print('Erro ao procurar pelo padrão')

def search_freq():
    try:
        query=input('Introduza o simbolo ou subsequência: ')
        ID=input('Introduza o identificador: ')
        eng.search_freq(query,ID)
    except:
        print('Erro ao calcular a frequência')

def search_freqDB():
    try:
        query=input('Introduza o simbolo ou subsequência: ')
        eng.search_freq(query)
    except:
        print('Erro ao calcular a frequência')

def Align():
    hits=[]
    while True:
        print('''\n\n\n\n
              --------Gestor de Sequências de DNA--------
              
              1.
              2.Realizar Blast Online através do portal do NCBI.
              3.Adicionar as sequências resultantes do Blast ao gestor.
              4.Criar árvore filogenética usando UPGMA
              6.Voltar ao menu anterior.
              
              -------------------------------------------''')
        x=int(input('\n Escolha uma das opções apresentadas:'))
        if x==1:
            pass
        elif x==2:
            search_patt()
        elif x==3:
            searc_pattDB()
        elif x==4:
            search_freq()
        elif x==5:
            search_freqDB()
        elif x==6:
            do_start()
        else:
            print('\n Opção não válida, escolha outra vez')


def saveDB():
    try:
        eng.save_dic()
    except:
        print('Erro ao gravar a instância do gestor')

def do_sair( arg=None):
    "Sair do programa Gestor de sequências de DNA: sair"
    print('Obrigado por ter utilizado o Gestor de sequências de DNA')
    global janela  # pois pretendo atribuir um valor a um identificador global
    if janela is not None:
                del janela  # invoca o metodo destruidor de instancia __del__()
    return True

if __name__ == '__main__':
    eng=GestorDNA()
    do_start()
