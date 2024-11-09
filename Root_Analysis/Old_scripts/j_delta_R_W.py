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
output_dir = os.path.join(script_dir, "j_dR")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
##############################################################
#Definindo pesos, pro hisgrama distribuir de acordo com a seção de choque:
##############################################################
#PESO = (SEÇÃO DE CHOQUE * LUMINOSIDADE)/(NUMERO TOTAL DE EVENTOS GERADOS (10000))

##############################################################

#P_1= 7559.214707 #tt
P_2= 13.19608992 #ttbb
#P_3= 2.690111995 #ttbbbb
P_4= 2.683539397 #tth
P_5= 1.155650301 #ttz
P_6= 0.048070682 #tttt
P_7= 0.002389038 #signal
P_8= 0.001228458 #ttzh
P_9= 0.000411384 #ttzz

#lista_de_pesos = [P_1, P_2, P_4, P_5, P_6, P_7, P_8, P_9]
lista_de_pesos = [P_2, P_4, P_5, P_6, P_8, P_9, P_7]

#########################################################################
######################## Delta_R (j1 j2) ##############################
#########################################################################
#	LOADING TREES BLOCK
#
ROOT.gROOT.SetBatch(True)

def calculate_invariant_mass(jet1, jet2):
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


def calculate_delta_r(jet1, jet2):
    delta_eta = jet1.Eta - jet2.Eta
    delta_phi = math.fabs(jet1.Phi - jet2.Phi)
    if delta_phi > math.pi:
        delta_phi = 2 * math.pi - delta_phi
    return math.sqrt(delta_eta**2 + delta_phi**2)


# Lista para armazenar os histogramas de PT e ETA de cada arquivo
histograms_PT = []
mean_values_PT = []
mean_values_ETA = []
std_dev_values_PT = []
std_dev_values_ETA = []
modal_values = []
N_entries = []
min_values = []
max_values = []


histograms_PT2 = []
mean_values_PT2 = []
mean_values_ETA2 = []
std_dev_values_PT2 = []
std_dev_values_ETA2 = []
modal_values2 = []
N_entries2 = []
min_values2 = []
max_values2 = []


# Desativar as estatísticas padrão
ROOT.gStyle.SetOptStat(0)

# Ativar a exibição do valor médio e desvio padrão personalizado
ROOT.gStyle.SetOptFit(1111)

# Cores para cada processo
colors =  [ ROOT.kMagenta, ROOT.kBlue,  ROOT.kGreen+1, ROOT.kBlack, ROOT.kGray, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]


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
    histogram_PT = ROOT.TH1F(f"hist_pt_{i}", "#Delta R of q-jets from W[1]; #Delta R ; Number of events (L = 200 fb^{-1})", 40, 0, 10)
    histogram_PT2 = ROOT.TH1F(f"hist_pt_{i}", "#Delta R of q-jets from W[2]; #Delta R ; Number of events (L = 200 fb^{-1})", 40, 0, 10)

    # histogram_ETA = ROOT.TH1F(f"hist_eta_{i}", "#eta jb[1]; #eta ; Number of events (L = 200 fb^{-1})", 40, -10, 10)

    # Loop sobre os eventos
    for entry in range(numberOfEntries):
        print(f"Processing event {entry} of file {inputFile}")

        # Lendo entrada específica
        treeReader.ReadEntry(entry)
        jets = treeReader.UseBranch("Jet")
        jets_list = []
        jet_pairs = []

        for jet in jets:
            # Verificando se a partícula é um jato b
            if abs(jet.Flavor) in [0, 1, 2, 3, 4, 21] and jet.PT >= 30 and abs(jet.Eta) <= 3:
                    jets_list.append(jet)
            
        combinations = itertools.combinations(jets_list, 2)
        
        for jet1, jet2 in combinations:
            if jet1 != jet2:
                flavor_pairs = {
                    (1, 2), (1, 3), (4, 2), (4, 3), 
                    (2, 1), (3, 1), (2, 4), (3, 4)
                }
                jet_flavors = (abs(jet1.Flavor), abs(jet2.Flavor))
                
                if jet_flavors in flavor_pairs:
                    invariant_mass = calculate_invariant_mass(jet1, jet2)
                    delta_r = calculate_delta_r(jet1, jet2)
                    pt_sum = jet1.PT + jet2.PT
                    
                    # Armazenando os pares de jatos e suas propriedades
                    jet_pairs.append({
                        'jets': (jet1, jet2),
                        'PT': (jet1.PT, jet2.PT),
                        'Eta': (jet1.Eta, jet2.Eta),
                        'Phi': (jet1.Phi, jet2.Phi),
                        'mass': invariant_mass,
                        'delta_r': delta_r,
                        'pt_sum': pt_sum
                    })

        # Selecionar os dois pares de jatos com a massa mais próxima do bóson W
        mass_W = 80.4  # GeV/c^2
        jet_pairs_sorted_by_mass = sorted(jet_pairs, key=lambda pair: abs(pair['mass'] - mass_W))

        closest_pairs = jet_pairs_sorted_by_mass[:2]

        # Preencher os histogramas com os valores de ΔR dos pares selecionados
        if len(closest_pairs) > 1:

            if closest_pairs[0]['pt_sum'] > closest_pairs[1]['pt_sum'] :
                histogram_PT.Fill(closest_pairs[0]['delta_r'])
                histogram_PT2.Fill(closest_pairs[1]['delta_r'])
            
            if closest_pairs[1]['pt_sum'] > closest_pairs[0]['pt_sum'] :
                histogram_PT.Fill(closest_pairs[1]['delta_r'])
                histogram_PT2.Fill(closest_pairs[0]['delta_r'])

    # Adicionando histogramas à lista
    histograms_PT.append(histogram_PT)
    histograms_PT2.append(histogram_PT2)

    # Calculando média e desvio padrão para PT
    # Calculando média e desvio padrão para PT
    mean_pt = histogram_PT.GetMean()
    std_dev_pt = histogram_PT.GetStdDev()
    # Número de entradas
    n_entries = histogram_PT.GetEntries()

    # Valor modal (bin com o máximo número de entradas)
    max_bin_content = 0
    modal_value = 0

    # Iterar sobre todos os bins do histograma
    for bin in range(1, histogram_PT.GetNbinsX() + 1):
        bin_content = histogram_PT.GetBinContent(bin)
        if bin_content > max_bin_content:
            max_bin_content = bin_content
            modal_value = histogram_PT.GetBinCenter(bin)


        min_value = float('inf')
        max_value = -float('inf')

    # Iterar sobre os bins do histograma para encontrar o valor mínimo e máximo do eixo x
    for bin in range(1, histogram_PT.GetNbinsX() + 1):
        bin_content = histogram_PT.GetBinContent(bin)
        # Verificar se o conteúdo do bin é diferente de zero
        if bin_content != 0:
            bin_x_value = histogram_PT.GetXaxis().GetBinLowEdge(bin)
            min_value = min(min_value, bin_x_value)
            max_value = max(max_value, histogram_PT.GetXaxis().GetBinUpEdge(bin))

    mean_values_PT.append(mean_pt)
    std_dev_values_PT.append(std_dev_pt)
    modal_values.append(modal_value)
    N_entries.append(n_entries)
    min_values.append(min_value)
    max_values.append(max_value)

#########################################################3
   # Calculando média e desvio padrão para PT
    # Calculando média e desvio padrão para PT
    mean_pt2 = histogram_PT2.GetMean()
    std_dev_pt2 = histogram_PT2.GetStdDev()
    # Número de entradas
    n_entries2 = histogram_PT2.GetEntries()

    # Valor modal (bin com o máximo número de entradas)
    max_bin_content2 = 0
    modal_value2 = 0

    # Iterar sobre todos os bins do histograma
    for bin in range(1, histogram_PT2.GetNbinsX() + 1):
        bin_content2 = histogram_PT2.GetBinContent(bin)
        if bin_content2 > max_bin_content2:
            max_bin_content2 = bin_content2
            modal_value2 = histogram_PT2.GetBinCenter(bin)


        min_value2 = float('inf')
        max_value2 = -float('inf')

    # Iterar sobre os bins do histograma para encontrar o valor mínimo e máximo do eixo x
    for bin in range(1, histogram_PT2.GetNbinsX() + 1):
        bin_content = histogram_PT2.GetBinContent(bin)
        # Verificar se o conteúdo do bin é diferente de zero
        if bin_content != 0:
            bin_x_value2 = histogram_PT2.GetXaxis().GetBinLowEdge(bin)
            min_value2 = min(min_value2, bin_x_value2)
            max_value2 = max(max_value2, histogram_PT2.GetXaxis().GetBinUpEdge(bin))

    mean_values_PT2.append(mean_pt2)
    std_dev_values_PT2.append(std_dev_pt2)
    modal_values2.append(modal_value2)
    N_entries2.append(n_entries2)
    min_values2.append(min_value2)
    max_values2.append(max_value2)

# Calcular a soma de todos os elementos, exceto o último, em N_entries
sum_except_last1 = sum(N_entries[:-1])

# Calcular a razão R
if sum_except_last1 != 0:
    R1 = N_entries[-1] / sum_except_last1
else:
    R1 = float('inf')  # Evitar divisão por zero



# Encontrar o menor valor de bin não nulo para definir o limite inferior de y
min_nonzero_bin_content = float('inf')
for histogram_PT in histograms_PT:
    for i in range(1, histogram_PT.GetNbinsX() + 1):
        bin_content = histogram_PT.GetBinContent(i)
        if bin_content > 0 and bin_content < min_nonzero_bin_content:
            min_nonzero_bin_content = bin_content
# Definir o limite inferior de y (use um fator de segurança se necessário)
y_min = min_nonzero_bin_content * 0.9 if min_nonzero_bin_content != float('inf') else 0.01


# Calcular a soma de todos os elementos, exceto o último, em N_entries
sum_except_last2 = sum(N_entries2[:-1])

# Calcular a razão R
if sum_except_last2 != 0:
    R2 = N_entries2[-1] / sum_except_last2
else:
    R2 = float('inf')  # Evitar divisão por zero



# Encontrar o menor valor de bin não nulo para definir o limite inferior de y
min_nonzero_bin_content2 = float('inf')

for histogram_PT2 in histograms_PT2:
    for i in range(1, histogram_PT2.GetNbinsX() + 1):
        bin_content = histogram_PT2.GetBinContent(i)
        if bin_content > 0 and bin_content < min_nonzero_bin_content:
            min_nonzero_bin_content = bin_content

# Definir o limite inferior de y (use um fator de segurança se necessário)
y_min2 = min_nonzero_bin_content2 * 0.9 if min_nonzero_bin_content2 != float('inf') else 0.01
    
# Criando canvas para os histogramas de PT e ETA
canvas_width = 2200  # Ajuste a largura do canvas
canvas_height = 1000  # Altura do canvas
canvas_PT = ROOT.TCanvas("canvas_PT", "PT Histogram", canvas_width, canvas_height)

# Ajustando margens do canvas para criar mais espaço à direita
canvas_PT.SetLeftMargin(0.1)
canvas_PT.SetRightMargin(0.15)
canvas_PT.SetBottomMargin(0.1)
canvas_PT.SetTopMargin(0.1)

# canvas_ETA = ROOT.TCanvas("canvas_ETA", "ETA Histogram", 1600, 1000)


# Criando canvas para os histogramas de PT e ETA
canvas_width = 2200  # Ajuste a largura do canvas
canvas_height = 1000  # Altura do canvas
canvas_PT = ROOT.TCanvas("canvas_PT", "PT Histogram", canvas_width, canvas_height)

# Ajustando margens do canvas para criar mais espaço à direita
canvas_PT.SetLeftMargin(0.1)
canvas_PT.SetRightMargin(0.15)
canvas_PT.SetBottomMargin(0.1)
canvas_PT.SetTopMargin(0.1)
# canvas_ETA = ROOT.TCanvas("canvas_ETA", "ETA Histogram", 1600, 1000)

# Calcular o valor máximo de todos os histogramas
max_y_pt = max([histogram_PT.GetMaximum() for histogram_PT in histograms_PT])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de PT com base no valor máximo calculado
for histogram_PT in histograms_PT:
    histogram_PT.SetMaximum(max_y_pt * margin_factor)
    histogram_PT.SetMinimum(y_min)

canvas_PT.cd()
for i, histogram_PT in enumerate(histograms_PT):
    histogram_PT.SetLineColor(colors[i])
    if i != len(histograms_PT) - 1: 
        histogram_PT.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_PT.SetLineWidth(2)
    histogram_PT.Draw("SAME")
    histogram_PT.Sumw2(0)
    histogram_PT.SetMinimum(0.0001)  # Definindo escala logarítmica
    canvas_PT.SetLogy()  # Ativando escala logarítmica no eixo Y

# Criando legenda personalizada
# Criando legenda personalizada
    legend_PT = ROOT.TLegend(0.85, 1.02 - 0.12 * i, 0.92, 0.83 - 0.12 * i)
    legend_PT.SetTextSize(0.02)
    legend_PT.SetEntrySeparation(0.0005)
    legend_PT.SetBorderSize(0)
    legend_PT.SetFillStyle(0)
    histogram_PT.GetListOfFunctions().Remove(legend_PT)

    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(histogram_PT, f"  {os.path.basename(inputFiles[i]).replace('100.root', '')}", "f100")
    legend_PT.AddEntry(0, f"Mean: {mean_values_PT[i]:.1f}, Std Dev: {std_dev_values_PT[i]:.1f}", "")
    legend_PT.AddEntry(0, f"Min: {min_values[i]:.1f}, Max: {max_values[i]:.1f}", "")
    legend_PT.AddEntry(0, f"Entries: {N_entries[i]:.1f}, Mod. V: {modal_values[i]:.1f}", "")

    legend_PT.AddEntry(0, "", "")
    legend_PT.Draw()

text = ROOT.TLatex()
text.SetNDC()
text.SetTextSize(0.03)
text.DrawLatex(0.85, 0.02, f"R: {R1:.3f}")
canvas_PT.Update()

output_file_path_PT = os.path.join(output_dir, "deltarw1.png")
canvas_PT.SaveAs(output_file_path_PT)

# Abrir a imagem gerada
##subprocess.Popen(["xdg-open", output_file_path_PT])
def clear_canvas(c):
    c.Clear()
    c.Update()

# Limpa o canvas
clear_canvas(canvas_PT)

   
# Criando canvas para os histogramas de PT e ETA
canvas_width = 2200  # Ajuste a largura do canvas
canvas_height = 1000  # Altura do canvas
canvas_PT = ROOT.TCanvas("canvas_PT", "PT Histogram", canvas_width, canvas_height)

# Ajustando margens do canvas para criar mais espaço à direita
canvas_PT.SetLeftMargin(0.1)
canvas_PT.SetRightMargin(0.15)
canvas_PT.SetBottomMargin(0.1)
canvas_PT.SetTopMargin(0.1)

# canvas_ETA = ROOT.TCanvas("canvas_ETA", "ETA Histogram", 1600, 1000)


# Criando canvas para os histogramas de PT e ETA
canvas_width = 2200  # Ajuste a largura do canvas
canvas_height = 1000  # Altura do canvas
canvas_PT = ROOT.TCanvas("canvas_PT", "PT Histogram", canvas_width, canvas_height)

# Ajustando margens do canvas para criar mais espaço à direita
canvas_PT.SetLeftMargin(0.1)
canvas_PT.SetRightMargin(0.15)
canvas_PT.SetBottomMargin(0.1)
canvas_PT.SetTopMargin(0.1)


# Calcular o valor máximo de todos os histogramas
max_y_pt2 = max([histogram_PT2.GetMaximum() for histogram_PT2 in histograms_PT2])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de PT com base no valor máximo calculado
for histogram_PT2 in histograms_PT2:
    histogram_PT2.SetMaximum(max_y_pt2 * margin_factor)
    histogram_PT2.SetMinimum(y_min2)

canvas_PT.cd()
for i, histogram_PT2 in enumerate(histograms_PT2):
    histogram_PT2.SetLineColor(colors[i])
    if i != len(histograms_PT2) - 1: 
        histogram_PT2.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_PT2.SetLineWidth(2)
    histogram_PT2.Draw("SAME")
    histogram_PT2.Sumw2(0)
    canvas_PT.SetLogy()  # Ativando escala logarítmica no eixo Y

# Criando legenda personalizada
# Criando legenda personalizada
    legend_PT = ROOT.TLegend(0.85, 1.02 - 0.12 * i, 0.92, 0.83 - 0.12 * i)
    legend_PT.SetTextSize(0.02)
    legend_PT.SetEntrySeparation(0.0005)
    legend_PT.SetBorderSize(0)
    legend_PT.SetFillStyle(0)
    histogram_PT.GetListOfFunctions().Remove(legend_PT)

    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(histogram_PT2, f"  {os.path.basename(inputFiles[i]).replace('100.root', '')}", "f100")
    legend_PT.AddEntry(0, f"Mean: {mean_values_PT2[i]:.1f}, Std Dev: {std_dev_values_PT2[i]:.1f}", "")
    legend_PT.AddEntry(0, f"Min: {min_values2[i]:.1f}, Max: {max_values2[i]:.1f}", "")
    legend_PT.AddEntry(0, f"Entries: {N_entries2[i]:.1f}, Mod. V: {modal_values2[i]:.1f}", "")

    legend_PT.AddEntry(0, "", "")
    legend_PT.Draw()

text = ROOT.TLatex()
text.SetNDC()
text.SetTextSize(0.03)
text.DrawLatex(0.85, 0.02, f"R: {R2:.3f}")
canvas_PT.Update()

output_file_path_PT = os.path.join(output_dir, "deltarw2.png")
canvas_PT.SaveAs(output_file_path_PT)





