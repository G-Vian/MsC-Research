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


def Invariant_mass(jet1, jet2, jet3, jet4):
    # Extrai as componentes de momento e energia para os jatos j1 e j2 --> Todas essas quantidades estão acessiveis no branch jet
    p1x, p1y, p1z, E1 = jet1.PT * ROOT.TMath.Cos(jet1.Phi), jet1.PT * ROOT.TMath.Sin(jet1.Phi), jet1.PT * ROOT.TMath.SinH(jet1.Eta), jet1.PT * ROOT.TMath.CosH(jet1.Eta)
    p2x, p2y, p2z, E2 = jet2.PT * ROOT.TMath.Cos(jet2.Phi), jet2.PT * ROOT.TMath.Sin(jet2.Phi), jet2.PT * ROOT.TMath.SinH(jet2.Eta), jet2.PT * ROOT.TMath.CosH(jet2.Eta)
    p3x, p3y, p3z, E3 = jet3.PT * ROOT.TMath.Cos(jet3.Phi), jet3.PT * ROOT.TMath.Sin(jet3.Phi), jet3.PT * ROOT.TMath.SinH(jet3.Eta), jet3.PT * ROOT.TMath.CosH(jet3.Eta)
    p4x, p4y, p4z, E4 = jet4.PT * ROOT.TMath.Cos(jet4.Phi), jet4.PT * ROOT.TMath.Sin(jet4.Phi), jet4.PT * ROOT.TMath.SinH(jet4.Eta), jet4.PT * ROOT.TMath.CosH(jet4.Eta)

    # Crie os vetores de momento de quatro dimensões para os jatos j1 e j2--> isso é um vetor de lorentz : https://root.cern.ch/doc/master/classTLorentzVector.html
    jet1_vec = ROOT.Math.PxPyPzEVector(p1x, p1y, p1z, E1)
    jet2_vec = ROOT.Math.PxPyPzEVector(p2x, p2y, p2z, E2)
    jet3_vec = ROOT.Math.PxPyPzEVector(p3x, p3y, p3z, E3)
    jet4_vec = ROOT.Math.PxPyPzEVector(p4x, p4y, p4z, E4)

    # Calcule a soma dos vetores de momento dos dois jatos
    Total_vec = jet1_vec + jet2_vec + jet3_vec + jet4_vec

    # Calcule a massa invariante (total_vec.M() calcula a massa invariante usando a energia e o momento contidos no vetor de Lorentz total_vec, isso é uma definição do ROOT)
    Mass = Total_vec.M()

    return Mass
##############################################################
##############################################################

def deltaR(jet1, jet2):
    # Extrai as coordenadas eta e phi dos objetos jet1 e jet2
    eta_jet1 = jet1.Eta
    eta_jet2 = jet2.Eta
    phi_jet1 = jet1.Phi
    phi_jet2 = jet2.Phi
    
    # Calcula as diferenças em eta e phi
    d_eta = eta_jet1 - eta_jet2
    d_phi = phi_jet1 - phi_jet2
    
    # Normaliza a diferença de phi para o intervalo [-pi, pi]
    #if d_phi > math.pi:
    #    d_phi -= 2*math.pi
    #elif d_phi < -math.pi:
    #    d_phi += 2*math.pi
    
    # Calcula o deltaR usando a fórmula padrão
    delta_R = math.sqrt(d_eta**2 + d_phi**2)
    
    return delta_R


##############################################################
#Definindo pesos, pro hisgrama distribuir de acordo com a seção de choque:
##############################################################
#PESO = (SEÇÃO DE CHOQUE * LUMINOSIDADE)/(NUMERO TOTAL DE EVENTOS GERADOS (10000))

##############################################################
# Pesos originais
P_2 = 13.19608992  # ttbb
P_4 = 2.683539397  # tth
P_5 = 1.155650301  # ttz
P_6 = 0.048070682  # tttt
P_7 = 0.002389038  # signal
P_8 = 0.00122845   # ttzh
P_9 = 0.000411384  # ttzz

# Calcular a soma total dos pesos
soma_total = P_2 + P_4 + P_5 + P_6 + P_7 + P_8 + P_9

# Normalizar os pesos
P_2_normalizado = P_2 / soma_total
P_4_normalizado = P_4 / soma_total
P_5_normalizado = P_5 / soma_total
P_6_normalizado = P_6 / soma_total
P_7_normalizado = P_7 / soma_total
P_8_normalizado = P_8 / soma_total
P_9_normalizado = P_9 / soma_total

# Multiplicar por um fator comum para ajustar a escala de visualização
fator_escala = 1000

P_2_ajustado = P_2_normalizado * fator_escala
P_4_ajustado = P_4_normalizado * fator_escala
P_5_ajustado = P_5_normalizado * fator_escala
P_6_ajustado = P_6_normalizado * fator_escala
P_7_ajustado = P_7_normalizado * fator_escala
P_8_ajustado = P_8_normalizado * fator_escala
P_9_ajustado = P_9_normalizado * fator_escala

# Lista de pesos ajustados
lista_de_pesos = [P_2_ajustado, P_4_ajustado, P_5_ajustado, P_6_ajustado, P_7_ajustado, P_8_ajustado, P_9_ajustado]


##############################################################
#	LOADING TREES BLOCK
#
# Lista para armazenar os histogramas
histograms_1 = []
histograms_2 = []
mean_values_1 = []
std_dev_values_1 = []
mean_values_2 = []
std_dev_values_2 = []

modal_values_1 = []
modal_values_2 = []
N_entries1 = []
min_values1 = []
max_values1= []

N_entries2 = []
min_values2 = []
max_values2 = []


# Desativar as estatísticas padrão
ROOT.gStyle.SetOptStat(0)

# Ativar a exibição do valor médio e desvio padrão personalizado
ROOT.gStyle.SetOptFit(1111)

# Cores para cada processo
colors =  [ ROOT.kMagenta, ROOT.kBlue,  ROOT.kGreen+1, ROOT.kBlack, ROOT.kGray, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]



#########################M dos bottoms com cut#########################
delta_h = 5
mh = 91


for i, inputFile in enumerate(inputFiles):
    print(f"Reading file {inputFile}")

    # Histogramas globais para todos os eventos
    hist1 = ROOT.TH1F("hist1", "Invariant mass of z[1]", 50, 80, 100)
    hist2 = ROOT.TH1F("hist2", "Invariant mass of z[2]", 50, 80, 100)

    # Obter o peso correspondente ao arquivo atual
    peso_atual = lista_de_pesos[i]
    peso_atual = peso_atual
    
    # Criando cadeia de árvores ROOT
    chain = ROOT.TChain("Delphes")
    chain.Add(inputFile)

    # Criando objeto ExRootTreeReader para ler a cadeia de árvores
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

    for event in range(numberOfEntries):
        # Obter o evento atual
        treeReader.ReadEntry(event)

        # Obter o galho "Jet"
        branchJet = treeReader.UseBranch("Jet")

        # Certificar-se de que existem pelo menos dois jatos no evento
        if branchJet.GetEntries() < 2:
            continue

        jets = [jet for jet in branchJet if abs(jet.Flavor) == 5 and jet.PT >= 30 and abs(jet.Eta) <= 2.5]

        #Certificar-se de que existem pelo menos quatro jatos b
        if len(jets) < 4:
            # Selecionar jatos b que têm PT menor que 30 ou maior que 20
            additional_jets = [jet for jet in branchJet if abs(jet.Flavor) == 5 and 20 <= jet.PT < 30] # and abs(jet.Eta) <= 2.5]

            # Ordenar essa lista com base na proximidade do PT ao valor 30
            additional_jets.sort(key=lambda jet: abs(jet.PT - 30))

            # Adicionar esses jatos à lista principal até que a condição len(jets) < 4 seja satisfeita
            while len(jets) < 4 and additional_jets:
                jets.append(additional_jets.pop(0))

        # Caso ainda não haja quatro jatos b, continuar para o próximo evento
        if len(jets) < 4:
            continue
        
        # Calcular combinações de jatos b e suas massas invariantes
        combinacoes = itertools.combinations(jets, 2)
        massas_invariantes = [(invariant_mass(jet1, jet2), jet1, jet2) for jet1, jet2 in combinacoes]

        # Verificar se há pelo menos duas massas invariantes
        if len(massas_invariantes) < 2:
            continue

        # Inicializar variáveis para armazenar as melhores combinações
        melhor_massa1 = None
        melhor_massa2 = None
        menor_diferenca_total = float('inf')

        # Encontrar a melhor combinação de pares de jatos
        for (massa1, jet1a, jet1b), (massa2, jet2a, jet2b) in itertools.combinations(massas_invariantes, 2):
            if len({jet1a, jet1b, jet2a, jet2b}) == 4:  # Certificar-se de que os jatos são diferentes
                diferenca_total = abs(massa1 - mh) + abs(massa2 - mh)
                if diferenca_total < menor_diferenca_total:
                    menor_diferenca_total = diferenca_total
                    melhor_massa1 = massa1
                    melhor_massa2 = massa2

        # Verificar se encontramos uma combinação válida
        if melhor_massa1 is not None and melhor_massa2 is not None:
            # Preencher histogramas com as massas mais próximas de mh
            hist1.Fill(melhor_massa1, peso_atual)
            hist2.Fill(melhor_massa2, peso_atual)
                # Imprimir os resultados para verificação
            print(f"   Par de jatos 1: Massa Invariante = {massa1}, Jet 1: PT = {jet1a.PT}, Eta = {jet1a.Eta}, Phi = {jet1a.Phi}, Jet 2: PT = {jet1b.PT}, Eta = {jet1b.Eta}, Phi = {jet1b.Phi}")
            print(f"   Par de jatos 2: Massa Invariante = {massa2}, Jet 1: PT = {jet2a.PT}, Eta = {jet2a.Eta}, Phi = {jet2a.Phi}, Jet 2: PT = {jet2b.PT}, Eta = {jet2b.Eta}, Phi = {jet2b.Phi}")
        massas_invariantes.clear()
    # Normalizar os histogramas após preencher todos os eventos
    if hist1.Integral() > 0:
        hist1.Scale(1.0 / hist1.Integral())
    if hist2.Integral() > 0:
        hist2.Scale(1.0 / hist2.Integral())
    # Adicionando histogramas à lista
    histograms_1.append(hist1)
    histograms_2.append(hist2)

    # Calculando média e desvio padrão 
    mean_1 = hist1.GetMean()
    std_dev_1 = hist1.GetStdDev()
    mean_2 = hist2.GetMean()
    std_dev_2 = hist2.GetStdDev()

    mean_values_1.append(mean_1)
    std_dev_values_1.append(std_dev_1)
    mean_values_2.append(mean_2)
    std_dev_values_2.append(std_dev_2)

    # Número de entradas
    n_entries1 = hist1.GetEntries()

    # Valor modal (bin com o máximo número de entradas)
    max_bin_content = 0
    modal_value1 = 0

    # Iterar sobre todos os bins do histograma
    for bin in range(1, hist1.GetNbinsX() + 1):
        bin_content = hist1.GetBinContent(bin)
        if bin_content > max_bin_content:
            max_bin_content = bin_content
            modal_value1 = hist1.GetBinCenter(bin)

    min_value1 = float('inf')
    max_value1 = -float('inf')

    # Iterar sobre os bins do histograma para encontrar o valor mínimo e máximo do eixo x
    for bin in range(1, hist1.GetNbinsX() + 1):
        bin_content = hist1.GetBinContent(bin)
        # Verificar se o conteúdo do bin é diferente de zero
        if bin_content != 0:
            bin_x_value = hist1.GetXaxis().GetBinLowEdge(bin)
            min_value1 = min(min_value1, bin_x_value)
            max_value1 = max(max_value1, hist1.GetXaxis().GetBinUpEdge(bin))

    modal_values_1.append(modal_value1)
    N_entries1.append(n_entries1)
    min_values1.append(min_value1)
    max_values1.append(max_value1)

    # Número de entradas
    n_entries2 = hist2.GetEntries()

    # Valor modal (bin com o máximo número de entradas)
    max_bin_content = 0
    modal_value2 = 0

    # Iterar sobre todos os bins do histograma
    for bin in range(1, hist2.GetNbinsX() + 1):
        bin_content = hist2.GetBinContent(bin)
        if bin_content > max_bin_content:
            max_bin_content = bin_content
            modal_value2 = hist2.GetBinCenter(bin)

    min_value2 = float('inf')
    max_value2 = -float('inf')

    # Iterar sobre os bins do histograma para encontrar o valor mínimo e máximo do eixo x
    for bin in range(1, hist2.GetNbinsX() + 1):
        bin_content = hist2.GetBinContent(bin)
        # Verificar se o conteúdo do bin é diferente de zero
        if bin_content != 0:
            bin_x_value = hist2.GetXaxis().GetBinLowEdge(bin)
            min_value2 = min(min_value2, bin_x_value)
            max_value2 = max(max_value2, hist2.GetXaxis().GetBinUpEdge(bin))

    modal_values_2.append(modal_value2)
    N_entries2.append(n_entries2)
    min_values2.append(min_value2)
    max_values2.append(max_value2)




# Criando canvas para os histogramas de PT e ETA
## Criando canvas para os histogramas de PT e ETA
canvas_width = 2200  # Ajuste a largura do canvas
canvas_height = 1000  # Altura do canvas
canvas_1 = ROOT.TCanvas("canvas_1", "h1", canvas_width, canvas_height)
canvas_2 = ROOT.TCanvas("canvas_2", "h2", canvas_width, canvas_height)

# Ajustando margens do canvas para criar mais espaço à direita
canvas_1.SetLeftMargin(0.1)
canvas_1.SetRightMargin(0.15)
canvas_1.SetBottomMargin(0.1)
canvas_1.SetTopMargin(0.1)



# Ajustando margens do canvas para criar mais espaço à direita
canvas_2.SetLeftMargin(0.1)
canvas_2.SetRightMargin(0.15)
canvas_2.SetBottomMargin(0.1)
canvas_2.SetTopMargin(0.1)



# Calcular o valor máximo de todos os histogramas
max_y_1 = max([hist1.GetMaximum() for hist1 in histograms_1])
max_y_2 = max([hist2.GetMaximum() for hist2 in histograms_2])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas com base no valor máximo calculado
for hist1 in histograms_1:
    hist1.SetMaximum(max_y_1 * margin_factor)
for hist2 in histograms_2:
    hist2.SetMaximum(max_y_2 * margin_factor)
# Desenhando os histogramas no canvas_1
canvas_1.cd()
for i, hist1 in enumerate(histograms_1):
    hist1.SetLineColor(colors[i])
    
    # Verifica se é o último histograma
    if i != len(histograms_1) - 1:
        hist1.SetFillColor(colors[i])
    
    hist1.SetLineWidth(2)
    hist1.Sumw2(0)

    hist1.Draw("SAME")
    
    # Criando legenda personalizada
    legend_PT = ROOT.TLegend(0.85, 1.02 - 0.12 * i, 0.92, 0.83 - 0.12 * i)
    legend_PT.SetTextSize(0.02)
    legend_PT.SetEntrySeparation(0.0005)
    legend_PT.SetBorderSize(0)
    legend_PT.SetFillStyle(0)
    hist1.GetListOfFunctions().Remove(legend_PT)

    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(hist1, f"  {os.path.basename(inputFiles[i]).replace('100.root', '')}", "f100")
    legend_PT.AddEntry(0, f"Mean: {mean_values_1[i]:.1f}, Std Dev: {std_dev_values_1[i]:.1f}", "")
    legend_PT.AddEntry(0, f"Min: {min_values1[i]:.1f}, Max: {max_values1[i]:.1f}", "")
    legend_PT.AddEntry(0, f"Entries: {N_entries1[i]:.1f}, Mod. V: {modal_values_1[i]:.1f}", "")

    legend_PT.AddEntry(0, "", "")
    legend_PT.Draw()

canvas_1.Update()

# Salvando o histograma em um arquivo de imagem
output_file_path_1 = os.path.join(output_dir, "mz1_v1.png")
canvas_1.SaveAs(output_file_path_1)
subprocess.Popen(["xdg-open", output_file_path_1])


# Desenhando os histogramas no canvas_2
canvas_2.cd()
for i, hist2 in enumerate(histograms_2):
    hist2.SetLineColor(colors[i])
    
    # Verifica se é o último histograma
    if i != len(histograms_2) - 1:
        hist2.SetFillColor(colors[i])
    
    hist2.SetLineWidth(2)
    hist2.Sumw2(0)

    hist2.Draw("SAME")
    
    # Criando legenda personalizada
    legend_PT = ROOT.TLegend(0.85, 1.02 - 0.12 * i, 0.92, 0.83 - 0.12 * i)
    legend_PT.SetTextSize(0.02)
    legend_PT.SetEntrySeparation(0.0005)
    legend_PT.SetBorderSize(0)
    legend_PT.SetFillStyle(0)
    hist1.GetListOfFunctions().Remove(legend_PT)

    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(hist2, f"  {os.path.basename(inputFiles[i]).replace('100.root', '')}", "f100")
    legend_PT.AddEntry(0, f"Mean: {mean_values_2[i]:.1f}, Std Dev: {std_dev_values_2[i]:.1f}", "")
    legend_PT.AddEntry(0, f"Min: {min_values2[i]:.1f}, Max: {max_values2[i]:.1f}", "")
    legend_PT.AddEntry(0, f"Entries: {N_entries2[i]:.1f}, Mod. V: {modal_values_2[i]:.1f}", "")

    legend_PT.AddEntry(0, "", "")
    legend_PT.Draw()

canvas_2.Update()
# Salvando o histograma em um arquivo de imagem
output_file_path_2 = os.path.join(output_dir, "mz2_v1.png")
canvas_2.SaveAs(output_file_path_2)
subprocess.Popen(["xdg-open", output_file_path_2])
#######