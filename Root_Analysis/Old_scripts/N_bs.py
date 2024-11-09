#!/usr/bin/env python
############################################################
#               Code for generating multiple jets plots
#               Gabriel V. Vian
#               MsC Student from Institute of Theoretical Physics -- IFT UNESP
############################################################
#	ATTENTION: CHANGE ACCORDING TO YOUR INSTALLATION
#
# this is needed if running outside Delphes directory

DelphesDirectory = "/home/gabriel/Desktop/MG5_aMC_v2_9_17/Delphes"
#/Desktop/MG5_aMC_v2_9_17$

############################################################

############################################################
#	INITIALIZATION BLOCK
#
import sys
import ROOT
import os
import math
import subprocess #to allow me to open the pictures and do not interrupt the code 


# Check for user input
if len(sys.argv) < 1:
    print("Please input the correct number of .ROOT files!")
    sys.exit(1)

# Saving the absolute path to the ".root" input files
inputFiles = sys.argv[1:]
inputFiles = [os.path.abspath(file) for file in inputFiles]

############################################################

############################################################
#	DELPHES/ROOT LIBRARIES BLOCK
#
# Save current directory
current_dir = os.getcwd()

# Change to Delphes directory
os.chdir(DelphesDirectory)
new_dir = os.getcwd()

# Loading Delphes libraries
ROOT.gSystem.Load("libDelphes")
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')

############################################################

############################################################
#	Criando diretórios importantes 
############################################################

# Diretório do script
script_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(script_dir, "N_jets")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
##############################################################
#Definindo pesos, pro hisgrama distribuir de acordo com a seção de choque:
##############################################################
#PESO = (SEÇÃO DE CHOQUE * LUMINOSIDADE)/(NUMERO TOTAL DE EVENTOS GERADOS (10000))

##############################################################

P_1= 7559.214707 #tt
P_2= 13.19608992 #ttbb
#P_3= 2.690111995 #ttbbbb
P_4= 2.683539397 #tth
P_5= 1.155650301 #ttz
P_6= 0.048070682 #tttt
P_7= 0.002389038 #signal
P_8= 0.001228458 #ttzh
P_9= 0.000411384 #ttzz

lista_de_pesos = [P_1, P_2, P_4, P_5, P_6, P_7, P_8, P_9]
############################################
#Bloco de leitura de cada arquivo .root
############################################################
for i, inputFile in enumerate(inputFiles):
    print(f"Reading file {inputFile}")
   
    # Obter o peso correspondente ao arquivo atual
    peso_atual = lista_de_pesos[i]
    
    # Criando cadeia de árvores ROOT
    chain = ROOT.TChain("Delphes")
    chain.Add(inputFile)

    # Criando objeto ExRootTreeReader para ler a cadeia de árvores
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

      # Criando histogramas para PT e ETA dos jato b 
    hist_jet_number = ROOT.TH1F("hist_jet_number", f"Number of b- Jets", 40, 0, 12)


############################################
#Bloco de leitura comandos em cada arquivo
############################################################





    # Loop sobre todos os eventos na árvore
    for event in range(numberOfEntries):
        # Obter o evento atual
        treeReader.ReadEntry(event)

        # Obtenha o galho dos jatos
        branchJet = treeReader.UseBranch("Jet")
        
        # Contar o número de jatos de |sabor| == 5 no evento atual
        num_jets_flavor5 = 0
        for jetIndex in range(branchJet.GetEntries()):
            jetObject = branchJet.At(jetIndex)
            if abs(jetObject.Flavor) == 5:
                num_jets_flavor5 += 1
        
        # Preencher o histograma com o número de jatos de |sabor| == 5 no evento atual
        hist_jet_number.Fill( num_jets_flavor5)

        print(f"Número de jatos b no evento {event}:", num_jets_flavor5)




    # Crie um canvas para os histogramas
    canvas = ROOT.TCanvas("canvas", "Canvas", 1200, 1000)
    canvas.SetLogy()

    #hist_jet_number.Sumw2(0)
  #  hist_jet_number.SetLineColor(ROOT.kBlue)
  #  hist_jet_number.SetLineWidth(2)
    hist_jet_number.Sumw2(0)
    hist_jet_number.Draw()


    image_path = os.path.join(output_dir, f"bjet_number_perevent_{os.path.basename(inputFile).replace('.root', '.png')}")
    canvas.SaveAs(image_path)


    # Abra os histogramas quando prontos
    subprocess.Popen(["xdg-open", image_path])
