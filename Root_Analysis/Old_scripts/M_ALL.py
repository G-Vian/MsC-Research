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
import itertools


# Check for user input
if len(sys.argv) < 9:
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
output_dir = os.path.join(script_dir, "M")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

##############################################################
##############################################################

# Defina a função para calcular a massa invariante
def invariant_mass(jet1, jet2):
    # Extrai as componentes de momento e energia para os jatos j1 e j2 --> Todas essas quantidades estão acessiveis no branch jet
    p1x, p1y, p1z, E1 = jet1.PT * ROOT.TMath.Cos(jet1.Phi), jet1.PT * ROOT.TMath.Sin(jet1.Phi), jet1.PT * ROOT.TMath.SinH(jet1.Eta), jet1.PT * ROOT.TMath.CosH(jet1.Eta)
    p2x, p2y, p2z, E2 = jet2.PT * ROOT.TMath.Cos(jet2.Phi), jet2.PT * ROOT.TMath.Sin(jet2.Phi), jet2.PT * ROOT.TMath.SinH(jet2.Eta), jet2.PT * ROOT.TMath.CosH(jet2.Eta)

    # Crie os vetores de momento de quatro dimensões para os jatos j1 e j2--> isso é um vetor de lorentz : https://root.cern.ch/doc/master/classTLorentzVector.html
    jet1_vec = ROOT.Math.PxPyPzEVector(p1x, p1y, p1z, E1)
    jet2_vec = ROOT.Math.PxPyPzEVector(p2x, p2y, p2z, E2)

    # Calcule a soma dos vetores de momento dos dois jatos
    total_vec = jet1_vec + jet2_vec

    # Calcule a massa invariante (total_vec.M() calcula a massa invariante usando a energia e o momento contidos no vetor de Lorentz total_vec, isso é uma definição do ROOT)
    mass = total_vec.M()

    return mass
##############################################################
##############################################################


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

##############################################################
#	LOADING TREES BLOCK
#
# Lista para armazenar os histogramas de PT e ETA de cada arquivo
histograms_HT = []
mean_values_HT = []
std_dev_values_HT = []

# Desativar as estatísticas padrão
ROOT.gStyle.SetOptStat(0)

# Ativar a exibição do valor médio e desvio padrão personalizado
ROOT.gStyle.SetOptFit(1111)

# Cores para cada processo
colors = [ROOT.kBlack, ROOT.kGray,  ROOT.kCyan, ROOT.kRed, ROOT.kMagenta, ROOT.kBlue, ROOT.kGreen+1, ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]



#########################every possible pair of all jets#########################


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
    histogram_HT = ROOT.TH1F("hist_M", "M; M (GeV/c^{²}) (all pair combinationw of jets); Number of events", 50, 0, 800)
     # Loop sobre os eventos
    for event in range(numberOfEntries):
        # Obter o evento atual
        treeReader.ReadEntry(event)
        
        jets = []
        
        # Obter o galho "Jet"
        branchJet = treeReader.UseBranch("Jet")
        
        # Certificar-se de que existem pelo menos dois jatos no evento
        if branchJet.GetEntries() < 2:
            continue
        
        for jet in branchJet:
            jets.append(jet)
        
        if jets:

            for jet in jets:
                combinacoes = itertools.combinations(jets, 2)
                for combination in combinacoes:
                    jet1, jet2 = combination
                    if jet1 != jet2: 
        
                        mass_invariant = invariant_mass(jet1, jet2)
                        print("Massa invariante dos jatos no evento", event, ": ", mass_invariant) 
    # Obter o número total de eventos

    total_events = histogram_HT.GetEntries()


    # Adicionando histogramas à lista
    histograms_HT.append(histogram_HT)

    # Calculando média e desvio padrão para PT
    mean_ht = histogram_HT.GetMean()
    std_dev_ht = histogram_HT.GetStdDev()

    mean_values_HT.append(mean_ht)
    std_dev_values_HT.append(std_dev_ht)
    
# Criando canvas para os histogramas de PT e ETA
canvas_HT = ROOT.TCanvas("canvas_HT", "HT Histogram", 1600, 1000)


# Calcular o valor máximo de todos os histogramas
max_y_ht = max([histogram_HT.GetMaximum() for histogram_HT in histograms_HT])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de ETA com base no valor máximo calculado
for histogram_HT in histograms_HT:
    histogram_HT.SetMaximum(max_y_ht * margin_factor)




canvas_HT.cd()
for i, histogram_HT in enumerate(histograms_HT):
    histogram_HT.SetLineColor(colors[i])
    histogram_HT.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_HT.SetLineWidth(2)
    histogram_HT.Draw("SAME")
    histogram_HT.Sumw2(0)
    histogram_HT.SetMinimum(0.0001)  # Definindo escala logarítmica
    canvas_HT.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Calculando média e desvio padrão
    mean_ht = histogram_HT.GetMean()
    std_dev_ht = histogram_HT.GetStdDev()

    # Criando legenda personalizada
    legend_HT = ROOT.TLegend(0.92, 1.02 - 0.12 * i, 0.95, 0.83 - 0.12 * i) # Posição fora da figura
    legend_HT.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_HT.SetEntrySeparation(0.0005) 
    legend_HT.SetBorderSize(0)  # Removendo borda da legenda
    legend_HT.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_HT.GetListOfFunctions().Remove(legend_HT)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda


# Adicionando o nome do arquivo à legenda

    legend_HT.AddEntry(0, "", "")
    legend_HT.AddEntry(histogram_HT, f"  {os.path.basename(inputFiles[i]).replace('.root', '')}", "f100")
    legend_HT.AddEntry(0, f"Mean: {mean_ht:.2f}", "")
    #legend_PT.AddEntry(0, "", "")
    legend_HT.AddEntry(0, f"Std Dev: {std_dev_ht:.2f}", "")
    legend_HT.AddEntry(0, "", "")

    legend_HT.Draw()


canvas_HT.Update()

  # Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 100  
    # Definir o limite superior para PT e ETA com base no valor máximo de y encontrado nao funcionou
histogram_HT.SetMaximum(max_y_ht * margin_factor)


# Salvando os histogramas em arquivos de imagem
output_file_path_HT = os.path.join(output_dir, "M_all.png")
canvas_HT.SaveAs(output_file_path_HT)
# Abra os histogramas quando prontos
subprocess.Popen(["xdg-open", output_file_path_HT])





