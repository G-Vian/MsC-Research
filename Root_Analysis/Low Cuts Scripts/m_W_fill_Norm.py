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
output_file_path = os.path.join(output_dir, "processing_log.txt")

# Excluir o arquivo de log antigo, se existir
if os.path.exists(output_file_path):
    os.remove(output_file_path)
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


#Essa função olha para cada jato de lista "additional_jets" e calcula sua massa invariante com cada jato da lista "jets". Então, ela escolhe aquele jato da lista "additional_jets"
# que, combinado com os jatos da lista "jets", possui a massa invariante mais próxima da massa do Higgs

def select_best_jet(additional_jets, jets, mh):
    best_jet = None
    best_mass_diff = float('inf')
    
    for additional_jet in additional_jets:
        for jet in jets:
            mass = invariant_mass(additional_jet, jet)
            mass_diff = abs(mass - mh)
            if mass_diff < best_mass_diff:
                best_mass_diff = mass_diff
                best_jet = additional_jet
    
    return best_jet




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



##############################################################
P_1  = 75.59214707  # tt
P_2  = 1.319608992  # ttbb
P_3  = 0.4078219283 # ttwh
P_4  = 0.2690111995 # ttbbbb
P_5  = 0.2683539397 # tth
P_6  = 0.2267598608 # ttwz
P_7  = 0.1155650301 # ttz
P_8  = 0.003341408  # tttt
P_9  = 0.0028875448 # ttww
P_10 = 0.0003030022 # tttw
P_11 = 0.0002389038 # tthh
P_12 = 0.0001228458 # ttzh
P_13 = 0.0000411384 # ttzz

lista_de_pesos = [
    P_3,   # ttwh
    P_6,   # ttwz
    P_4,   # ttbbbb
    P_7,   # ttz
    P_8,   # tttt
    P_10,  # tttw
    P_12,  # ttzh
    P_13,  # ttzz
    P_11   # tthh (signal)
]

colors = [
    ROOT.kOrange+4, # ttwh
    ROOT.kAzure+3,  # ttwz
    ROOT.kBlue,       # ttbbbb
    ROOT.kOrange+7, # ttz
    ROOT.kTeal-3,   # tttt
    ROOT.kYellow,     # tttw
    ROOT.kBlack,      # ttzh
    ROOT.kCyan,       # ttzz
    ROOT.kRed         # tthh (signal)
]
#########################M dos bottoms com cut#########################
delta_h = 25
mh = 80
   


with open(output_file_path, "a") as output_file:

    for i, inputFile in enumerate(inputFiles):
        print(f"Reading file {inputFile}")
    # Histogramas globais para todos os eventos
        hist1 = ROOT.TH1F("hist1", "Invariant mass of the most energetic W;  M_{j_{q} j_{q}} (GeV); Number of events (Normalized to one)", 20, 60, 100)
        hist2 = ROOT.TH1F("hist2", "Invariant mass of the 2nd most energetic W;  M_{j_{q} j_{q}} (GeV); Number of events (Normalized to one)", 20, 60, 100)
    
        # Obter o peso correspondente ao arquivo atual
        peso_atual = lista_de_pesos[i]
        
        # Criando cadeia de árvores ROOT
        chain = ROOT.TChain("Delphes")
        chain.Add(inputFile)

        # Criando objeto ExRootTreeReader para ler a cadeia de árvores
        treeReader = ROOT.ExRootTreeReader(chain)
        numberOfEntries = treeReader.GetEntries()


    # A definição dos galhos deve ser feito uma única vez antes de iniciar o loop. Não é necessário (nem recomendado) chamar UseBranch a cada iteração do loop, pois isso pode causar sobrecarga.

            # OBTENDO OS GALHOS
        branchJet = treeReader.UseBranch("Jet")
        branchMissingET = treeReader.UseBranch("MissingET")
        branchHT = treeReader.UseBranch("ScalarHT")
        branchelec = treeReader.UseBranch("Electron")
        branchmu = treeReader.UseBranch("Muon")


        # Contadores para acompanhar quantos eventos foram selecionados
        total_events = 0
        events_passed_MET_HT_nJets = 0
        events_passed_final_selection = 0


        # Loop sobre os eventos
        for entry in range(numberOfEntries):
            print(f"Processing event {entry} of file {inputFile}")
            total_events += 1  # Incrementa o contador de eventos totais

            # Lendo entrada específica
            treeReader.ReadEntry(entry)


            #OBTENDO AS GRANDEZAS DE CADA GALHO
                
            missingET = branchMissingET.At(0)  # Assume-se que há apenas um objeto de MET por evento
            MET_value = missingET.MET        
            nJets = branchJet.GetEntries()
            ht_object = branchHT.At(0)  # Assume-se que há apenas um objeto de HT por evento
            HT_value = ht_object.HT        #jet_PT = branchJet.PT #Não precisa pra esses objetos
            #jet_ETA = branchJet.ETA
            #elec_PT = branchelec.PT
            #mu_PT = branchmu.PT
            #elec_ETA = branchelec.ETA
            #mu_ETA = branchmu.ETA
            #LISTA DE JATOS E LEPTONS SELECIONADOS QUE DEVE SER RESETADA PRA CADA EVENTO NOVO
            selected_jets = []
            q_jets = []
            b_jets = []
            electrons = []
            muons = []
            
            ###CUT FLOW###

            #APLICANDO CORTES (VERIFICANDO SE O EVENTO SATISFAZ OS CRITÉRIOS)
            #PRIMEIRA CONDICIONAL
            if   nJets >= 0 and  HT_value > 500 : 
                events_passed_MET_HT_nJets += 1 #pra eu contar quantos eventos passaram
                selected_jets = [jet for jet in branchJet if jet.PT >= 40 and abs(jet.Eta) <= 2.5 ]        
                q_jets = [jet for jet in selected_jets if abs(jet.Flavor)  in [0, 1, 2, 3, 4]  ] 
                b_jets = [jet for jet in selected_jets if abs(jet.Flavor) == 5]
                electrons = [electron for electron in branchelec if electron.PT > 20 and abs(electron.Eta) <= 2.5  ] 
                muons = [muon for muon in branchmu if muon.PT > 20 and abs(muon.Eta) <= 2.5 ] 

            #SEGUNDA CONDICIONAL, DENTRO DA PRIMEIRA, QUE OBSERVA O NUMERO DE OBJETOS
                if len(q_jets) >= 0 and len(b_jets) >= 0 and len(muons) == 0 and len(electrons) == 0 : 
                    events_passed_final_selection += 1  # Incrementa o contador para eventos que passam essa seleção final

        
                    # Calcular combinações de jatos b e suas massas invariantes
                    combinations = itertools.combinations(q_jets, 2)
                    massas_invariantes = []
                    
                    for jet1, jet2 in combinations:
                        if jet1 != jet2:
                            flavor_pairs = {
                                (1, 2), (1, 3), (4, 2), (4, 3), 
                                (2, 1), (3, 1), (2, 4), (3, 4)
                            }
                            jet_flavors = (abs(jet1.Flavor), abs(jet2.Flavor))
                            
                            if jet_flavors in flavor_pairs:
                                mass = invariant_mass(jet1, jet2)
                                massas_invariantes.append((mass, jet1, jet2))
                    # Lista para armazenar valores de X e suas combinações correspondentes
                    valores_X = []

                    # Gerar todas as combinações possíveis de dois pares de jatos diferentes
                    pares_de_pares = list(itertools.combinations(massas_invariantes, 2))    
                    for (massa1, jet1a, jet1b), (massa2, jet2a, jet2b) in pares_de_pares:
                        # Certificar-se de que os jatos são diferentes
                        if len({jet1a, jet1b, jet2a, jet2b}) == 4:
                            # Calcular X
                            X = ((massa1 - mh)**2 / delta_h**2) + ((massa2 - mh)**2 / delta_h**2)
                            # Armazenar o valor de X e os pares correspondentes
                            valores_X.append((X, (jet1a, jet1b), (jet2a, jet2b)))

                    # Encontrar a combinação que minimiza X
                    if valores_X:
                        min_X_combination = min(valores_X, key=lambda x: x[0])

                        # Verificar se as combinações satisfazem os critérios adicionais
                        X_min, (jet1a, jet1b), (jet2a, jet2b) = min_X_combination
                        massa1 = invariant_mass(jet1a, jet1b)
                        massa2 = invariant_mass(jet2a, jet2b)
                        PT_1 = jet1a.PT + jet1b.PT
                        PT_2 = jet2a.PT + jet2b.PT
                        if abs(massa1 - mh) <= delta_h and abs(massa2 - mh) <= delta_h and PT_1 > PT_2:
                            # Preencher histogramas com os pesos
                            hist1.Fill(massa1, peso_atual)
                            hist2.Fill(massa2, peso_atual)
                        if abs(massa1 - mh) <= delta_h and abs(massa2 - mh) <= delta_h and PT_1 < PT_2:
                            hist1.Fill(massa2, peso_atual)
                            hist2.Fill(massa1, peso_atual)
                            # Imprimir os resultados para verificação
                       
                    massas_invariantes.clear()



            else:
                    continue  
        # Salvar os resultados em um arquivo de texto para o arquivo de entrada atual
        output_file.write(f"Results for file: {inputFile}\n")
        output_file.write(f"Total de eventos processados: {total_events}\n")
        output_file.write(f"Eventos que passaram MET >= 40, HT > 500 e nJets >= 6: {events_passed_MET_HT_nJets}\n")
        output_file.write(f"Eventos que passaram a seleção final: {events_passed_final_selection}\n")
        output_file.write("-" * 40 + "\n")
        

        n_entries1 = peso_atual*events_passed_final_selection

        n_entries2 = peso_atual*events_passed_final_selection
       
      

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
                max_value2 = max(max_value2, hist1.GetXaxis().GetBinUpEdge(bin))


        modal_values_2.append(modal_value2)
        N_entries2.append(n_entries2)
        min_values2.append(min_value2)
        max_values2.append(max_value2)



# Calcular a soma de todos os elementos, exceto o último, em N_entries1
sum_except_last1 = sum(N_entries1[:-1])

# Calcular a razão R
if sum_except_last1 != 0:
    R1 = N_entries1[-1] / sum_except_last1
else:
    R1 = float('inf')  # Evitar divisão por zero

# Formatar R1 em notação científica com tolerância a valores pequenos

if R1 != "Infinity":
    R1_str = f"{R1:.2e}"  # Formatar R1 em notação científica padrão
    base, exponent = R1_str.split('e')
    R1 = f"R: {float(base):.2f} \\times 10^{{{int(exponent)}}}"  # Formatar para LaTeX
else:
    R1 = "R: \\infty"  # Caso infinito


# Calcular a soma de todos os elementos, exceto o último, em N_entries1
sum_except_last2 = sum(N_entries2[:-1])

# Calcular a razão R
if sum_except_last2 != 0:
    R2 = N_entries2[-1] / sum_except_last2
else:
    R2 = float('inf')  # Evitar divisão por zero

# Formatar R1 em notação científica com tolerância a valores pequenos

if R2 != "Infinity":
    R2_str = f"{R2:.2e}"  # Formatar R1 em notação científica padrão
    base, exponent = R2_str.split('e')
    R2 = f"R: {float(base):.2f} \\times 10^{{{int(exponent)}}}"  # Formatar para LaTeX
else:
    R2 = "R: \\infty"  # Caso infinito


# Criando canvas para os histogramas de PT e ETA
## Criando canvas para os histogramas de PT e ETA
canvas_width = 2200  # Ajuste a largura do canvas
canvas_height = 1300  # Altura do canvas
canvas_1 = ROOT.TCanvas("canvas_1", "h1", canvas_width, canvas_height)
canvas_2 = ROOT.TCanvas("canvas_2", "h2", canvas_width, canvas_height)

# Ajustando margens do canvas para criar mais espaço à direita
canvas_1.SetLeftMargin(0.1)
canvas_1.SetRightMargin(0.12)
  
canvas_1.SetBottomMargin(0.15)
canvas_1.SetTopMargin(0.1)



# Ajustando margens do canvas para criar mais espaço à direita
canvas_2.SetLeftMargin(0.1)
canvas_2.SetRightMargin(0.12)
  
canvas_2.SetBottomMargin(0.15)  
canvas_2.SetTopMargin(0.1)



# Calcular o valor máximo de todos os histogramas
max_y_1 = max([hist1.GetMaximum() for hist1 in histograms_1])
max_y_2 = max([hist2.GetMaximum() for hist2 in histograms_2])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.05  

# Definir o limite superior para todos os histogramas com base no valor máximo calculado
for hist1 in histograms_1:
    hist1.SetMaximum(max_y_1 * margin_factor)
for hist2 in histograms_2:
    hist2.SetMaximum(max_y_2 * margin_factor)

# Desenhando os histogramas no canvas_1

# Desenhando os histogramas no canvas_1
canvas_1.cd()
for i, hist1 in enumerate(histograms_1):
    hist1.SetLineColor(colors[i])
    hist1.SetLineWidth(5)
    
    # Verifica se é o último histograma
    if i != len(histograms_1) - 1:
        hist1.SetFillColor(colors[i])
        hist1.SetLineWidth(2)
    hist1.Sumw2(0)

    hist1.Draw("SAME")
    
    # Criando legenda personalizada
    n_entries = len(histograms_1)  # Número de entradas na legenda

    # Posição inicial (topo) da legenda
    y1 = 0.956

    # Posição final (base) da legenda
    y2 = 0.01
  


    legend_PT = ROOT.TLegend(0.90, y1 - i * (y1 - y2) / n_entries, 1.0 , y1 - (i + 1) * (y1 - y2) / n_entries)
    legend_PT.SetTextSize(0.023)  # NOVO PADRÃO
    legend_PT.SetEntrySeparation(0.5)  # NOVO PADRÃO
    legend_PT.SetBorderSize(0)
    legend_PT.SetFillStyle(0)
    hist1.GetListOfFunctions().Remove(legend_PT)

    # Adicionando apenas o nome do processo com a cor
    legend_PT.AddEntry(hist1, f"  {os.path.basename(inputFiles[i]).replace('100.root', '')}", "f")

    # Desenhando a legenda
    legend_PT.Draw()



text = ROOT.TLatex()
text.SetNDC()
text.SetTextSize(0.04)  # Ajuste o tamanho do texto para menor
text.DrawLatex(0.615, 0.83, R1)


canvas_1.Update()

# Salvando o histograma em um arquivo de imagem
output_file_path_1 = os.path.join(output_dir, "mW1_fill_norm.png")
canvas_1.SaveAs(output_file_path_1)
subprocess.Popen(["xdg-open", output_file_path_1])


# Desenhando os histogramas no canvas_2
canvas_2.cd()
for i, hist2 in enumerate(histograms_2):
    hist2.SetLineColor(colors[i])
    hist2.SetLineWidth(5)
    
    # Verifica se é o último histograma
    if i != len(histograms_2) - 1:
        hist2.SetFillColor(colors[i])
    
        hist2.SetLineWidth(2)

    hist2.Sumw2(0)

    hist2.Draw("SAME")
    
    # Criando legenda personalizada
    n_entries = len(histograms_2)  # Número de entradas na legenda

    # Posição inicial (topo) da legenda
    y1 = 0.956

    # Posição final (base) da legenda
    y2 = 0.01
  


    legend_PT = ROOT.TLegend(0.90, y1 - i * (y1 - y2) / n_entries, 1.0 , y1 - (i + 1) * (y1 - y2) / n_entries)
    legend_PT.SetTextSize(0.023)  # NOVO PADRÃO
    legend_PT.SetEntrySeparation(0.5)  # NOVO PADRÃO
    legend_PT.SetBorderSize(0)
    legend_PT.SetFillStyle(0)
    hist2.GetListOfFunctions().Remove(legend_PT)

    # Adicionando apenas o nome do processo com a cor
    legend_PT.AddEntry(hist2, f"  {os.path.basename(inputFiles[i]).replace('100.root', '')}", "f")

    # Desenhando a legenda
    legend_PT.Draw()


text = ROOT.TLatex()
text.SetNDC()
text.SetTextSize(0.04)  # Ajuste o tamanho do texto para menor
text.DrawLatex(0.615, 0.83, R2)

canvas_2.Update()
# Salvando o histograma em um arquivo de imagem
output_file_path_2 = os.path.join(output_dir, "mW2_fill_norm.png")
canvas_2.SaveAs(output_file_path_2)
subprocess.Popen(["xdg-open", output_file_path_2])
#######