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


def deltaR(jet1, jet2):
    # Extrai as coordenadas eta e phi dos objetos jet1 e jet2
    eta_jet1 = jet1.Eta.p
    eta_jet2 = jet2.Eta
    phi_jet1 = jet1.Phi
    phi_jet2 = jet2.Phi
    
    # Calcula as diferenças em eta e phi
    d_eta = eta_jet1 - eta_jet2
    d_phi = phi_jet1 - phi_jet2
    

    
    # Calcula o deltaR usando a fórmula padrão
    delta_R = math.sqrt(d_eta**2 + d_phi**2)
    
    return delta_R


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
    histogram_PT = ROOT.TH1F(f"hist_pt_{i}", "#Delta R j[1] and j[2]; #Delta R ; Number of events (L = 200 fb^{-1})", 40, 0, 10)
    # histogram_ETA = ROOT.TH1F(f"hist_eta_{i}", "#eta jb[1]; #eta ; Number of events (L = 200 fb^{-1})", 40, -10, 10)

    # Loop sobre os eventos
    for entry in range(numberOfEntries):
        print(f"Processing event {entry} of file {inputFile}")

        # Lendo entrada específica
        treeReader.ReadEntry(entry)
        jets = treeReader.UseBranch("Jet")

        # Inicializando listas para armazenar os PTs e ETAs dos jato b
        pt_b = []
        eta_b = []
        phi_b = []

        # Loop sobre as partículas
        for jet in jets:
            # Verificando se a partícula é um jato b
            if abs(jet.Flavor) in [0, 1, 2, 3, 4, 21] and jet.PT >= 30 and abs(jet.Eta) <= 3 :
                pt_b.append(jet.PT)
                eta_b.append(jet.Eta)
                phi_b.append(jet.Phi)

        # Verifica se há bottoms neste evento
        if len(pt_b) >= 2:
            # Ordena os PTs dos jato b em ordem decrescente
            pt_b.sort(reverse=True)
            idx1 = pt_b.index(pt_b[0])
            idx2 = pt_b.index(pt_b[1])

            # Calcula as diferenças em eta e phi
            d_eta = eta_b[idx1] - eta_b[idx2]
            d_phi = phi_b[idx1] - phi_b[idx2]
            if d_phi > math.pi:
                d_phi -= 2 * math.pi
            elif d_phi < -math.pi:
                d_phi += 2 * math.pi
    
            # Calcula o deltaR usando a fórmula padrão
            delta_R = math.sqrt(d_eta**2 + d_phi**2)
            
            # Preenche o histograma com o PT do jato b mais energético
            histogram_PT.Fill(delta_R, peso_atual)
   # Adicionando histogramas à lista
    histograms_PT.append(histogram_PT)

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


# Calcular a soma de todos os elementos, exceto o último, em N_entries
sum_except_last1 = sum(N_entries[:-1])

# Calcular a razão R
if sum_except_last1 != 0:
    R1 = N_entries[-1] / sum_except_last1
else:
    R1 = float('inf')  # Evitar divisão por zero



# Encontrar o menor valor de bin não nulo para definir o limite inferior de y
min_nonzero_bin_content = float('inf')

for histogram_list in [histograms_PT]:
    for hist in histogram_list:
        for i in range(1, hist.GetNbinsX() + 1):
            bin_content = hist.GetBinContent(i)
            if bin_content > 0 and bin_content < min_nonzero_bin_content:
                min_nonzero_bin_content = bin_content

# Definir o limite inferior de y (use um fator de segurança se necessário)
y_min = min_nonzero_bin_content * 0.9 if min_nonzero_bin_content != float('inf') else 0.01
    
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
        # Definindo escala logarítmica
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

output_file_path_PT = os.path.join(output_dir, "deltarj1j2.png")
canvas_PT.SaveAs(output_file_path_PT)

# Abrir a imagem gerada
##subprocess.Popen(["xdg-open", output_file_path_PT])


#########################################################################
######################## Delta_R (j3 j4) ##############################
#########################################################################


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
    histogram_PT = ROOT.TH1F(f"hist_pt_{i}", "#Delta R j[3] and j[4]; #Delta R ; Number of events (L = 200 fb^{-1})", 40, 0, 12)
    # histogram_ETA = ROOT.TH1F(f"hist_eta_{i}", "#eta jb[1]; #eta ; Number of events (L = 200 fb^{-1})", 40, -10, 10)

    # Loop sobre os eventos
    for entry in range(numberOfEntries):
        print(f"Processing event {entry} of file {inputFile}")

        # Lendo entrada específica
        treeReader.ReadEntry(entry)
        jets = treeReader.UseBranch("Jet")

        # Inicializando listas para armazenar os PTs e ETAs dos jato b
        pt_b = []
        eta_b = []
        phi_b = []

        # Loop sobre as partículas
        for jet in jets:
            # Verificando se a partícula é um jato b
            if abs(jet.Flavor) in [0, 1, 2, 3, 4, 21] and jet.PT >= 30 and abs(jet.Eta) <= 3:
                pt_b.append(jet.PT)
                eta_b.append(jet.Eta)
                phi_b.append(jet.Phi)

        # Verifica se há bottoms neste evento
        if len(pt_b) >= 4:
            # Ordena os PTs dos jato b em ordem decrescente
            pt_b.sort(reverse=True)
            idx1 = pt_b.index(pt_b[2])
            idx2 = pt_b.index(pt_b[3])

            # Calcula as diferenças em eta e phi
            d_eta = eta_b[idx1] - eta_b[idx2]
            d_phi = phi_b[idx1] - phi_b[idx2]
            if d_phi > math.pi:
                d_phi -= 2 * math.pi
            elif d_phi < -math.pi:
                d_phi += 2 * math.pi
    
            # Calcula o deltaR usando a fórmula padrão
            delta_R = math.sqrt(d_eta**2 + d_phi**2)
            
            # Preenche o histograma com o PT do jato b mais energético
            histogram_PT.Fill(delta_R, peso_atual)

    # Adicionando histogramas à lista
    histograms_PT.append(histogram_PT)

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


# Calcular a soma de todos os elementos, exceto o último, em N_entries
sum_except_last1 = sum(N_entries[:-1])

# Calcular a razão R
if sum_except_last1 != 0:
    R1 = N_entries[-1] / sum_except_last1
else:
    R1 = float('inf')  # Evitar divisão por zero



# Encontrar o menor valor de bin não nulo para definir o limite inferior de y
min_nonzero_bin_content = float('inf')

for histogram_list in [histograms_PT]:
    for hist in histogram_list:
        for i in range(1, hist.GetNbinsX() + 1):
            bin_content = hist.GetBinContent(i)
            if bin_content > 0 and bin_content < min_nonzero_bin_content:
                min_nonzero_bin_content = bin_content

# Definir o limite inferior de y (use um fator de segurança se necessário)
y_min = min_nonzero_bin_content * 0.9 if min_nonzero_bin_content != float('inf') else 0.01
    
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
        # Definindo escala logarítmica
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

# Salvando os histogramas em arquivos de imagem

output_file_path_PT = os.path.join(output_dir, "deltarj3j4.png")
canvas_PT.SaveAs(output_file_path_PT)

# Abrir a imagem gerada
##subprocess.Popen(["xdg-open", output_file_path_PT])



#########################################################################
######################## Delta_R (j1 j3) ##############################
#########################################################################


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
    histogram_PT = ROOT.TH1F(f"hist_pt_{i}", "#Delta R j[1] and j[3]; #Delta R ; Number of events (L = 200 fb^{-1})", 40, 0, 12)
    # histogram_ETA = ROOT.TH1F(f"hist_eta_{i}", "#eta jb[1]; #eta ; Number of events (L = 200 fb^{-1})", 40, -10, 10)

    # Loop sobre os eventos
    for entry in range(numberOfEntries):
        print(f"Processing event {entry} of file {inputFile}")

        # Lendo entrada específica
        treeReader.ReadEntry(entry)
        jets = treeReader.UseBranch("Jet")

        # Inicializando listas para armazenar os PTs e ETAs dos jato b
        pt_b = []
        eta_b = []
        phi_b = []

        # Loop sobre as partículas
        for jet in jets:
            # Verificando se a partícula é um jato b
            if abs(jet.Flavor) in [0, 1, 2, 3, 4, 21] and jet.PT >= 30 and abs(jet.Eta) <= 3:
                pt_b.append(jet.PT)
                eta_b.append(jet.Eta)
                phi_b.append(jet.Phi)

        # Verifica se há bottoms neste evento
        if len(pt_b) >= 3:
            # Ordena os PTs dos jato b em ordem decrescente
            pt_b.sort(reverse=True)
            idx1 = pt_b.index(pt_b[0])
            idx2 = pt_b.index(pt_b[2])

            # Calcula as diferenças em eta e phi
            d_eta = eta_b[idx1] - eta_b[idx2]
            d_phi = phi_b[idx1] - phi_b[idx2]
            if d_phi > math.pi:
                d_phi -= 2 * math.pi
            elif d_phi < -math.pi:
                d_phi += 2 * math.pi
    
            # Calcula o deltaR usando a fórmula padrão
            delta_R = math.sqrt(d_eta**2 + d_phi**2)
            
            # Preenche o histograma com o PT do jato b mais energético
            histogram_PT.Fill(delta_R, peso_atual)

    # Adicionando histogramas à lista
    histograms_PT.append(histogram_PT)

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


# Calcular a soma de todos os elementos, exceto o último, em N_entries
sum_except_last1 = sum(N_entries[:-1])

# Calcular a razão R
if sum_except_last1 != 0:
    R1 = N_entries[-1] / sum_except_last1
else:
    R1 = float('inf')  # Evitar divisão por zero



# Encontrar o menor valor de bin não nulo para definir o limite inferior de y
min_nonzero_bin_content = float('inf')

for histogram_list in [histograms_PT]:
    for hist in histogram_list:
        for i in range(1, hist.GetNbinsX() + 1):
            bin_content = hist.GetBinContent(i)
            if bin_content > 0 and bin_content < min_nonzero_bin_content:
                min_nonzero_bin_content = bin_content

# Definir o limite inferior de y (use um fator de segurança se necessário)
y_min = min_nonzero_bin_content * 0.9 if min_nonzero_bin_content != float('inf') else 0.01
    
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
        # Definindo escala logarítmica
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

# Salvando os histogramas em arquivos de imagem

output_file_path_PT = os.path.join(output_dir, "deltarj1j3.png")
canvas_PT.SaveAs(output_file_path_PT)

# Abrir a imagem gerada
##subprocess.Popen(["xdg-open", output_file_path_PT])





#########################################################################
######################## Delta_R (j1 j4) ##############################
#########################################################################


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
    histogram_PT = ROOT.TH1F(f"hist_pt_{i}", "#Delta R j[1] and j[4]; #Delta R ; Number of events (L = 200 fb^{-1})", 40, 0, 12)
    # histogram_ETA = ROOT.TH1F(f"hist_eta_{i}", "#eta jb[1]; #eta ; Number of events (L = 200 fb^{-1})", 40, -10, 10)

    # Loop sobre os eventos
    for entry in range(numberOfEntries):
        print(f"Processing event {entry} of file {inputFile}")

        # Lendo entrada específica
        treeReader.ReadEntry(entry)
        jets = treeReader.UseBranch("Jet")

        # Inicializando listas para armazenar os PTs e ETAs dos jato b
        pt_b = []
        eta_b = []
        phi_b = []

        # Loop sobre as partículas
        for jet in jets:
            # Verificando se a partícula é um jato b
            if abs(jet.Flavor) in [0, 1, 2, 3, 4, 21] and jet.PT >= 30 and abs(jet.Eta) <= 3:
                pt_b.append(jet.PT)
                eta_b.append(jet.Eta)
                phi_b.append(jet.Phi)

        # Verifica se há bottoms neste evento
        if len(pt_b) >= 4:
            # Ordena os PTs dos jato b em ordem decrescente
            pt_b.sort(reverse=True)
            idx1 = pt_b.index(pt_b[0])
            idx2 = pt_b.index(pt_b[3])

            # Calcula as diferenças em eta e phi
            d_eta = eta_b[idx1] - eta_b[idx2]
            d_phi = phi_b[idx1] - phi_b[idx2]
            if d_phi > math.pi:
                d_phi -= 2 * math.pi
            elif d_phi < -math.pi:
                d_phi += 2 * math.pi
    
            # Calcula o deltaR usando a fórmula padrão
            delta_R = math.sqrt(d_eta**2 + d_phi**2)
            
            # Preenche o histograma com o PT do jato b mais energético
            histogram_PT.Fill(delta_R, peso_atual)

    # Adicionando histogramas à lista
    histograms_PT.append(histogram_PT)

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


# Calcular a soma de todos os elementos, exceto o último, em N_entries
sum_except_last1 = sum(N_entries[:-1])

# Calcular a razão R
if sum_except_last1 != 0:
    R1 = N_entries[-1] / sum_except_last1
else:
    R1 = float('inf')  # Evitar divisão por zero



# Encontrar o menor valor de bin não nulo para definir o limite inferior de y
min_nonzero_bin_content = float('inf')

for histogram_list in [histograms_PT]:
    for hist in histogram_list:
        for i in range(1, hist.GetNbinsX() + 1):
            bin_content = hist.GetBinContent(i)
            if bin_content > 0 and bin_content < min_nonzero_bin_content:
                min_nonzero_bin_content = bin_content

# Definir o limite inferior de y (use um fator de segurança se necessário)
y_min = min_nonzero_bin_content * 0.9 if min_nonzero_bin_content != float('inf') else 0.01
    
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
        # Definindo escala logarítmica
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

# Salvando os histogramas em arquivos de imagem

output_file_path_PT = os.path.join(output_dir, "deltarj1j4.png")
canvas_PT.SaveAs(output_file_path_PT)

# Abrir a imagem gerada
##subprocess.Popen(["xdg-open", output_file_path_PT])




#########################################################################
######################## Delta_R (j2 j3) ##############################
#########################################################################


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
    histogram_PT = ROOT.TH1F(f"hist_pt_{i}", "#Delta R j[2] and j[3]; #Delta R ; Number of events (L = 200 fb^{-1})", 40, 0, 12)
    # histogram_ETA = ROOT.TH1F(f"hist_eta_{i}", "#eta jb[1]; #eta ; Number of events (L = 200 fb^{-1})", 40, -10, 10)

    # Loop sobre os eventos
    for entry in range(numberOfEntries):
        print(f"Processing event {entry} of file {inputFile}")

        # Lendo entrada específica
        treeReader.ReadEntry(entry)
        jets = treeReader.UseBranch("Jet")

        # Inicializando listas para armazenar os PTs e ETAs dos jato b
        pt_b = []
        eta_b = []
        phi_b = []

        # Loop sobre as partículas
        for jet in jets:
            # Verificando se a partícula é um jato b
            if abs(jet.Flavor) in [0, 1, 2, 3, 4, 21] and jet.PT >= 30 and abs(jet.Eta) <= 3:
                pt_b.append(jet.PT)
                eta_b.append(jet.Eta)
                phi_b.append(jet.Phi)

        # Verifica se há bottoms neste evento
        if len(pt_b) >= 3:
            # Ordena os PTs dos jato b em ordem decrescente
            pt_b.sort(reverse=True)
            idx1 = pt_b.index(pt_b[1])
            idx2 = pt_b.index(pt_b[2])

            # Calcula as diferenças em eta e phi
            d_eta = eta_b[idx1] - eta_b[idx2]
            d_phi = phi_b[idx1] - phi_b[idx2]
            if d_phi > math.pi:
                d_phi -= 2 * math.pi
            elif d_phi < -math.pi:
                d_phi += 2 * math.pi
    
            # Calcula o deltaR usando a fórmula padrão
            delta_R = math.sqrt(d_eta**2 + d_phi**2)
            
            # Preenche o histograma com o PT do jato b mais energético
            histogram_PT.Fill(delta_R, peso_atual)

    # Adicionando histogramas à lista
    histograms_PT.append(histogram_PT)

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


# Calcular a soma de todos os elementos, exceto o último, em N_entries
sum_except_last1 = sum(N_entries[:-1])

# Calcular a razão R
if sum_except_last1 != 0:
    R1 = N_entries[-1] / sum_except_last1
else:
    R1 = float('inf')  # Evitar divisão por zero



# Encontrar o menor valor de bin não nulo para definir o limite inferior de y
min_nonzero_bin_content = float('inf')

for histogram_list in [histograms_PT]:
    for hist in histogram_list:
        for i in range(1, hist.GetNbinsX() + 1):
            bin_content = hist.GetBinContent(i)
            if bin_content > 0 and bin_content < min_nonzero_bin_content:
                min_nonzero_bin_content = bin_content

# Definir o limite inferior de y (use um fator de segurança se necessário)
y_min = min_nonzero_bin_content * 0.9 if min_nonzero_bin_content != float('inf') else 0.01
    
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
        # Definindo escala logarítmica
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

# Salvando os histogramas em arquivos de imagem

output_file_path_PT = os.path.join(output_dir, "deltarj2j3.png")
canvas_PT.SaveAs(output_file_path_PT)

# Abrir a imagem gerada
##subprocess.Popen(["xdg-open", output_file_path_PT])



#########################################################################
######################## Delta_R (j2 j4) ##############################
#########################################################################


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
    histogram_PT = ROOT.TH1F(f"hist_pt_{i}", "#Delta R j[2] and j[4]; #Delta R ; Number of events (L = 200 fb^{-1})", 40, 0, 12)
    # histogram_ETA = ROOT.TH1F(f"hist_eta_{i}", "#eta jb[1]; #eta ; Number of events (L = 200 fb^{-1})", 40, -10, 10)

    # Loop sobre os eventos
    for entry in range(numberOfEntries):
        print(f"Processing event {entry} of file {inputFile}")

        # Lendo entrada específica
        treeReader.ReadEntry(entry)
        jets = treeReader.UseBranch("Jet")

        # Inicializando listas para armazenar os PTs e ETAs dos jato b
        pt_b = []
        eta_b = []
        phi_b = []

        # Loop sobre as partículas
        for jet in jets:
            # Verificando se a partícula é um jato b
            if abs(jet.Flavor) in [0, 1, 2, 3, 4, 21] and jet.PT >= 30 and abs(jet.Eta) <= 3:
                pt_b.append(jet.PT)
                eta_b.append(jet.Eta)
                phi_b.append(jet.Phi)

        # Verifica se há bottoms neste evento
        if len(pt_b) >= 4:
            # Ordena os PTs dos jato b em ordem decrescente
            pt_b.sort(reverse=True)
            idx1 = pt_b.index(pt_b[1])
            idx2 = pt_b.index(pt_b[3])

            # Calcula as diferenças em eta e phi
            d_eta = eta_b[idx1] - eta_b[idx2]
            d_phi = phi_b[idx1] - phi_b[idx2]
            if d_phi > math.pi:
                d_phi -= 2 * math.pi
            elif d_phi < -math.pi:
                d_phi += 2 * math.pi
    
            # Calcula o deltaR usando a fórmula padrão
            delta_R = math.sqrt(d_eta**2 + d_phi**2)
            
            # Preenche o histograma com o PT do jato b mais energético
            histogram_PT.Fill(delta_R, peso_atual)

    # Adicionando histogramas à lista
    histograms_PT.append(histogram_PT)

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


# Calcular a soma de todos os elementos, exceto o último, em N_entries
sum_except_last1 = sum(N_entries[:-1])

# Calcular a razão R
if sum_except_last1 != 0:
    R1 = N_entries[-1] / sum_except_last1
else:
    R1 = float('inf')  # Evitar divisão por zero



# Encontrar o menor valor de bin não nulo para definir o limite inferior de y
min_nonzero_bin_content = float('inf')

for histogram_list in [histograms_PT]:
    for hist in histogram_list:
        for i in range(1, hist.GetNbinsX() + 1):
            bin_content = hist.GetBinContent(i)
            if bin_content > 0 and bin_content < min_nonzero_bin_content:
                min_nonzero_bin_content = bin_content

# Definir o limite inferior de y (use um fator de segurança se necessário)
y_min = min_nonzero_bin_content * 0.9 if min_nonzero_bin_content != float('inf') else 0.01
    
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
        # Definindo escala logarítmica
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

# Salvando os histogramas em arquivos de imagem
output_file_path_PT = os.path.join(output_dir, "deltarj2j4.png")
canvas_PT.SaveAs(output_file_path_PT)

# Abrir a imagem gerada
#subprocess.Popen(["xdg-open", output_file_path_PT])