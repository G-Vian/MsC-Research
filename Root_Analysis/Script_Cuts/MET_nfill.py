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
output_dir = os.path.join(script_dir, "MET_Nfill")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
output_file_path = os.path.join(output_dir, "processing_log.txt")

# Excluir o arquivo de log antigo, se existir
if os.path.exists(output_file_path):
    os.remove(output_file_path)

##############################################################
#Definindo pesos, pro hisgrama distribuir de acordo com a seção de choque:
##############################################################
#PESO = (SEÇÃO DE CHOQUE * LUMINOSIDADE)/(NUMERO TOTAL DE EVENTOS GERADOS (10000))

##############################################################



P_1  = 7559.214707 #tt
P_2  = 13.19608992 # ttbb
P_3  = 4.078219283 # ttwh
P_4  = 2.690111995 #ttbbbb
P_5 = 2.683539397 #tth
P_6 = 2.267598608 # ttwz
P_7  = 1.155650301 # ttz
P_8  = 0.03341408 # tttt
P_9  = 0.028875448 # ttww
P_10  = 0.003030022 #tttw
P_11  = 0.002389038 #tthh
P_12  = 0.001228458 #ttzh
P_13  = 0.000411384 # ttzz


lista_de_pesos = [P_1, P_2, P_3, P_4, P_5, P_6, P_7, P_8, P_9, P_10 , P_12, P_13, P_11]

############################################################
######################## H_T ##############################
############################################################
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
colors =  [ ROOT.kMagenta+3, ROOT.kMagenta, ROOT.kOrange+4, ROOT.kBlue,  ROOT.kGreen+3, ROOT.kAzure+3,  ROOT.kOrange+7, ROOT.kTeal-3, ROOT.kViolet, ROOT.kYellow, ROOT.kBlack, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]
with open(output_file_path, "a") as output_file:

    for i, inputFile in enumerate(inputFiles):
        print(f"Reading file {inputFile}")
        histogram_HT = ROOT.TH1F("hist_ht", "Missing transverse energy ; MET (GeV); Number of events (L = 200 fb^{-1})", 80, 0, 500)

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
            if   nJets >= 6 and HT_value > 500 : 
                events_passed_MET_HT_nJets += 1 #pra eu contar quantos eventos passaram
                selected_jets = [jet for jet in branchJet if jet.PT >= 40 and abs(jet.Eta) <= 2.5 ]        
                q_jets = [jet for jet in selected_jets if abs(jet.Flavor)  in [0, 1, 2, 3, 4]  ] 
                b_jets = [jet for jet in selected_jets if abs(jet.Flavor) == 5]
                electrons = [electron for electron in branchelec if electron.PT > 20 and abs(electron.Eta) <= 2.5  ] 
                muons = [muon for muon in branchmu if muon.PT > 20 and abs(muon.Eta) <= 2.5 ] 

            #SEGUNDA CONDICIONAL, DENTRO DA PRIMEIRA, QUE OBSERVA O NUMERO DE OBJETOS
                if len(q_jets) >= 2 and len(b_jets) >= 4 and len(muons) == 0 and len(electrons) == 0 : 
                    events_passed_final_selection += 1  # Incrementa o contador para eventos que passam essa seleção final

                    histogram_HT.Fill(MET_value, peso_atual)



            else:
                    continue  # Se a primeira condicional nao for satisfeita, vá para o próximo event
        print(f"Total de eventos processados: {total_events}")
        print(f"Eventos que passaram MET >= 40, HT > 500 e nJets >= 6: {events_passed_MET_HT_nJets}")
        print(f"Eventos que passaram a seleção final: {events_passed_final_selection}")
        # Salvar os resultados em um arquivo de texto para o arquivo de entrada atual
        output_file.write(f"Results for file: {inputFile}\n")
        output_file.write(f"Total de eventos processados: {total_events}\n")
        output_file.write(f"Eventos que passaram MET >= 40, HT > 500 e nJets >= 6: {events_passed_MET_HT_nJets}\n")
        output_file.write(f"Eventos que passaram a seleção final: {events_passed_final_selection}\n")
        output_file.write("-" * 40 + "\n")

        # Adicionando histogramas à lista
        histograms_1.append(histogram_HT)

        # Calculando média e desvio padrão 
        mean_1 = histogram_HT.GetMean()
        std_dev_1 = histogram_HT.GetStdDev()
        
        
        mean_values_1.append(mean_1)
        std_dev_values_1.append(std_dev_1)
        
            # Número de entradas
        n_entries1 = histogram_HT.GetEntries()

            # Valor modal (bin com o máximo número de entradas)
        max_bin_content = 0
        modal_value1 = 0

        # Iterar sobre todos os bins do histograma
        for bin in range(1, histogram_HT.GetNbinsX() + 1):
            bin_content = histogram_HT.GetBinContent(bin)
            if bin_content > max_bin_content:
                max_bin_content = bin_content
                modal_value1 = histogram_HT.GetBinCenter(bin)


            min_value1 = float('inf')
            max_value1 = -float('inf')

            # Iterar sobre os bins do histograma para encontrar o valor mínimo e máximo do eixo x
        for bin in range(1, histogram_HT.GetNbinsX() + 1):
            bin_content = histogram_HT.GetBinContent(bin)
            # Verificar se o conteúdo do bin é diferente de zero
            if bin_content != 0:
                bin_x_value = histogram_HT.GetXaxis().GetBinLowEdge(bin)
                min_value1 = min(min_value1, bin_x_value)
                max_value1 = max(max_value1, histogram_HT.GetXaxis().GetBinUpEdge(bin))


        modal_values_1.append(modal_value1)
        N_entries1.append(n_entries1)
        min_values1.append(min_value1)
        max_values1.append(max_value1)

    



# Calcular a soma de todos os elementos, exceto o último, em N_entries1
sum_except_last1 = sum(N_entries1[:-1])

# Calcular a razão R
if sum_except_last1 != 0:
    R1 = N_entries1[-1] / sum_except_last1
else:
    R1 = float('inf')  # Evitar divisão por zero



# Encontrar o menor valor de bin não nulo para definir o limite inferior de y
min_nonzero_bin_content = float('inf')

for hist in [histogram_HT]:
    for i in range(1, hist.GetNbinsX() + 1):
        bin_content = hist.GetBinContent(i)
        if bin_content > 0 and bin_content < min_nonzero_bin_content:
            min_nonzero_bin_content = bin_content

# Definir o limite inferior de y (use um fator de segurança se necessário)
y_min = min_nonzero_bin_content * 0.9 if min_nonzero_bin_content != float('inf') else 0.01




# Criando canvas para os histogramas de PT e ETA
## Criando canvas para os histogramas de PT e ETA
canvas_width = 2200  # Ajuste a largura do canvas
canvas_height = 1300  # Altura do canvas
canvas_HT = ROOT.TCanvas("canvas_HT", "h1", canvas_width, canvas_height)

# Ajustando margens do canvas para criar mais espaço à direita
canvas_HT.SetLeftMargin(0.1)
canvas_HT.SetRightMargin(0.19)
canvas_HT.SetBottomMargin(0.15)
canvas_HT.SetTopMargin(0.1)




# Calcular o valor máximo de todos os histogramas
max_y_1 = max([histogram_HT.GetMaximum() for histogram_HT in histograms_1])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas com base no valor máximo calculado
for histogram_HT in histograms_1:
    histogram_HT.SetMaximum(max_y_1 * margin_factor)





canvas_HT.cd()
for i, histogram_HT in enumerate(histograms_1):
    histogram_HT.SetLineColor(colors[i])

    # Verifica se é o último histograma
    histogram_HT.SetLineWidth(5)

    if i != len(histograms_1) - 1:
        histogram_HT.SetLineWidth(2)
    histogram_HT.Sumw2(0)
    histogram_HT.SetMinimum(1e-1)

    histogram_HT.Draw("SAME")
    canvas_HT.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Criando legenda personalizada

#################################################################
    n_entries = len(histograms_1)  # Número de entradas na legenda

    # Posição inicial (topo) da legenda
    y1 = 0.956

    # Posição final (base) da legenda
    y2 = 0.01


    # Distribuindo uniformemente
    legend_PT = ROOT.TLegend(0.85, y1 - i * (y1 - y2) / n_entries, 0.92, y1 - (i + 1) * (y1 - y2) / n_entries)

    legend_PT.SetTextSize(0.023) #NOVO PADRÃO
    legend_PT.SetEntrySeparation(0.8) #NOVO PADRÃO
    legend_PT.SetBorderSize(0)
    legend_PT.SetFillStyle(0)
    histogram_HT.GetListOfFunctions().Remove(legend_PT)
################################################################## #NOVO PADRÃO
    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(histogram_HT, f"  {os.path.basename(inputFiles[i]).replace('100.root', '')}", "f100")
    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(0, f"Mean: {mean_values_1[i]:.1f}, Std Dev: {std_dev_values_1[i]:.1f}", "")
    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(0, f"Min: {min_values1[i]:.1f}, Max: {max_values1[i]:.1f}", "")
    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(0, f"Entr.: {N_entries1[i]:.1f}, Mod. V: {modal_values_1[i]:.1f}", "")
    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(0, "", "")
    legend_PT.Draw()



text = ROOT.TLatex()
text.SetNDC()
text.SetTextSize(0.03)
text.DrawLatex(0.80, 0.01, f"R: {R1:.3f}")

canvas_HT.Update()



# Salvando os histogramas em arquivos de imagem
output_file_path_HT = os.path.join(output_dir, "MET.png")
canvas_HT.SaveAs(output_file_path_HT)
# Abra os histogramas quando prontos
#subprocess.Popen(["xdg-open", output_file_path_HT])

