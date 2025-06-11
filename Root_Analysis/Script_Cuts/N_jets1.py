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
output_dir = os.path.join(script_dir, "N_jets1")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
output_file_path = os.path.join(output_dir, "processing_log2.txt")

# Excluir o arquivo de log antigo, se existir
if os.path.exists(output_file_path):
    os.remove(output_file_path)

##############################################################
#Definindo pesos, pro hisgrama distribuir de acordo com a seção de choque:
##############################################################
#PESO = (SEÇÃO DE CHOQUE * LUMINOSIDADE)/(NUMERO TOTAL DE EVENTOS GERADOS (10000))

##############################################################
P_1  = 75.59214707  # tt
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



lista_de_pesos = [P_1, P_2, P_3, P_6, P_5, P_4 , P_7, P_8, P_9, P_10 , P_12, P_13, P_11] #sem cut
colors =  [ ROOT.kMagenta+3, ROOT.kMagenta, ROOT.kOrange+4, ROOT.kAzure+3, ROOT.kGreen+3, ROOT.kBlue,  ROOT.kOrange+7, ROOT.kTeal-3, ROOT.kViolet, ROOT.kYellow, ROOT.kBlack, ROOT.kCyan,  ROOT.kRed ]  # sem corte

#lista_de_pesos =  [P_1, P_3, P_5, P_6, P_2, P_4, P_7, P_8, P_12, P_10, P_13, P_9, P_11]
#colors =  [  ROOT.kMagenta+3, ROOT.kOrange+4, ROOT.kGreen+3, ROOT.kAzure+3,  ROOT.kMagenta,  ROOT.kBlue, ROOT.kOrange+7, ROOT.kTeal-3,  ROOT.kBlack, ROOT.kYellow,  ROOT.kCyan, ROOT.kViolet, ROOT.kRed ] 

############################################
# Lista para armazenar os histogramas
histograms_jets = []
histograms_b = []
histograms_q = []
n_entries_jets = []
n_entries_bjets = []
n_entries_qjets = []

############################################
#Bloco de leitura de cada arquivo .root
############################################################


# Desativar as estatísticas padrão
ROOT.gStyle.SetOptStat(0)

# Ativar a exibição do valor médio e desvio padrão personalizado
ROOT.gStyle.SetOptFit(1111)


with open(output_file_path, "a") as output_file:

    for i, inputFile in enumerate(inputFiles):
        print(f"Reading file {inputFile}")
        hist_jet_number = ROOT.TH1F("hist_jet_number", f"Number of jets per event", 12, 5, 16)
        hist_b_number = ROOT.TH1F("hist_bjet_number", f"Number of b-tagged jets per event", 8, 3, 10)
        hist_q_number = ROOT.TH1F("hist_qjet_number", f"Number of light-quark jets per event", 11, 0, 12)

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
            if nJets >= 6 and HT_value > 500 : 
                events_passed_MET_HT_nJets += 1 #pra eu contar quantos eventos passaram
                selected_jets = [jet for jet in branchJet if jet.PT >= 40 and abs(jet.Eta) <= 2.5 ]        
                q_jets = [jet for jet in selected_jets if abs(jet.Flavor)  in [0, 1, 2, 3, 4]  ] 
                b_jets = [jet for jet in selected_jets if abs(jet.Flavor) == 5]
                electrons = [electron for electron in branchelec if electron.PT > 20 and abs(electron.Eta) <= 2.5  ] 
                muons = [muon for muon in branchmu if muon.PT > 20 and abs(muon.Eta) <= 2.5 ] 

            #SEGUNDA CONDICIONAL, DENTRO DA PRIMEIRA, QUE OBSERVA O NUMERO DE OBJETOS
                if len(q_jets) >= 2 and len(b_jets) >= 4 and len(muons) == 0 and len(electrons) == 0 : 
                    events_passed_final_selection += 1  # Incrementa o contador para eventos que passam essa seleção final
                    num_jets = 0
                    num_bjets = 0
                    num_qjets= 0
                    for jetObject in q_jets:
                        num_jets += 1
                        num_qjets += 1
                    for jetObject in b_jets:
                        num_jets += 1
                        num_bjets += 1
                    # Preencher o histograma com o número de jatos de |sabor| == 5 no evento atual
                    hist_jet_number.Fill(num_jets, peso_atual)
                    hist_b_number.Fill(num_bjets, peso_atual)
                    hist_q_number.Fill(num_qjets, peso_atual)

         



            else:
                    continue  # Se a primeira condicional nao for satisfeita, vá para o próximo event

        # Salvar os resultados em um arquivo de texto para o arquivo de entrada atual
        output_file.write(f"Results for file: {inputFile}\n")
        output_file.write(f"Total de eventos processados: {total_events}\n")
        output_file.write(f"Eventos que passaram MET >= 40, HT > 500 e nJets >= 6: {events_passed_MET_HT_nJets}\n")
        output_file.write(f"Eventos que passaram a seleção final: {events_passed_final_selection}\n")
        output_file.write("-" * 40 + "\n")

    

        histograms_jets.append(hist_jet_number)
        histograms_b.append(hist_b_number)
        histograms_q.append(hist_q_number)

        # Número de entradas que sera usado pra calular o R e será o numero de eventos contados N apresentados na legenda do histograma.
        # Só considera-se numero de eventos acima ou igual a 1 e aqueles que são fracionarios são arredondados pra baixo
        # Não é o numero de eventos selecionados, é o numero de eventos selecionados e pesados com valor acima ou iguao a um e que seja inteiro.
        n_entriesj = 0
        for bin in range(1, hist_jet_number.GetNbinsX() + 1):
            bin_content = hist_jet_number.GetBinContent(bin)
            if bin_content >= 1:
                n_entriesj += int(bin_content)  # Arredonda para baixo truncando o valor
        
        n_entriesb = 0
        for bin in range(1, hist_b_number.GetNbinsX() + 1):
            bin_content = hist_b_number.GetBinContent(bin)
            if bin_content >= 1:
                n_entriesb += int(bin_content)  # Arredonda para baixo truncando o valor
        
        n_entriesq = 0
        for bin in range(1, hist_q_number.GetNbinsX() + 1):
            bin_content = hist_q_number.GetBinContent(bin)
            if bin_content >= 1:
                n_entriesq += int(bin_content)  # Arredonda para baixo truncando o valor
        n_entries_jets.append(n_entriesj)
        n_entries_bjets.append(n_entriesb)
        n_entries_qjets.append(n_entriesq)


 # Criando canvas para os histogramas de PT e ETA
## Criando canvas para os histogramas de PT e ETA
canvas_width = 2200  # Ajuste a largura do canvas
canvas_height = 1300  # Altura do canvas
canvas_jet = ROOT.TCanvas("cnv_jet", "h1", canvas_width, canvas_height)
canvas_bjet = ROOT.TCanvas("cnv_bjet", "h2", canvas_width, canvas_height)
canvas_qjet = ROOT.TCanvas("cnv_qjet", "h3", canvas_width, canvas_height)

# Ajustando margens do canvas para criar mais espaço à direita
canvas_jet.SetLeftMargin(0.1)
canvas_jet.SetRightMargin(0.12)
canvas_jet.SetBottomMargin(0.15)
canvas_jet.SetTopMargin(0.1)

canvas_bjet.SetLeftMargin(0.1)
canvas_bjet.SetRightMargin(0.12)
canvas_bjet.SetBottomMargin(0.15)
canvas_bjet.SetTopMargin(0.1)

canvas_qjet.SetLeftMargin(0.1)
canvas_qjet.SetRightMargin(0.12)
canvas_qjet.SetBottomMargin(0.15)
canvas_qjet.SetTopMargin(0.1)



canvas_jet.cd()
for i, hist_jet_number in enumerate(histograms_jets):
    hist_jet_number.SetLineColor(colors[i])
    hist_jet_number.SetLineWidth(5)

    # Verifica se é o último histograma
    if i != len(histograms_jets) - 1:
        hist_jet_number.SetFillColor(colors[i])
        hist_jet_number.SetLineWidth(2)

    hist_jet_number.Sumw2(0)
    hist_jet_number.SetMinimum(1e-1)

    hist_jet_number.Draw("SAME")
    canvas_jet.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Criando legenda personalizada
#################################################################
    n_entries = len(histograms_jets)  # Número de entradas na legenda

    # Posição inicial (topo) da legenda
    y1 = 0.956

    # Posição final (base) da legenda
    y2 = 0.01


    # Distribuindo uniformemente
    legend_PT = ROOT.TLegend(0.90, y1 - i * (y1 - y2) / n_entries, 1.0 , y1 - (i + 1) * (y1 - y2) / n_entries)
    legend_PT.SetTextSize(0.023)  # NOVO PADRÃO
    legend_PT.SetEntrySeparation(0.5)  # NOVO PADRÃO
    legend_PT.SetBorderSize(0)
    legend_PT.SetFillStyle(0)
    hist_jet_number.GetListOfFunctions().Remove(legend_PT)

    # Adicionando apenas o nome do processo com a cor
    legend_PT.AddEntry(hist_jet_number, f"  {os.path.basename(inputFiles[i]).replace('100.root', '')}", "f")

    # Desenhando a legenda
    legend_PT.Draw()


canvas_jet.Update()



# Salvando os histogramas em arquivos de imagem
output_file_path_HT = os.path.join(output_dir, "Njets.png")
canvas_jet.SaveAs(output_file_path_HT)
# Abra os histogramas quando prontos
subprocess.Popen(["xdg-open", output_file_path_HT])




canvas_bjet.cd()
for i, hist_bjet_number in enumerate(histograms_b):
    hist_bjet_number.SetLineColor(colors[i])
    hist_bjet_number.SetLineWidth(5)

    # Verifica se é o último histograma
    if i != len(histograms_b) - 1:
        hist_bjet_number.SetFillColor(colors[i])
        hist_bjet_number.SetLineWidth(2)

    hist_bjet_number.Sumw2(0)
    hist_bjet_number.SetMinimum(1e-1)

    hist_bjet_number.Draw("SAME")
    canvas_bjet.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Criando legenda personalizada
#################################################################
    n_entries = len(histograms_b)  # Número de entradas na legenda

    # Posição inicial (topo) da legenda
    y1 = 0.956

    # Posição final (base) da legenda
    y2 = 0.01


    # Distribuindo uniformemente
    legend_PT = ROOT.TLegend(0.90, y1 - i * (y1 - y2) / n_entries, 1.0 , y1 - (i + 1) * (y1 - y2) / n_entries)
    legend_PT.SetTextSize(0.023)  # NOVO PADRÃO
    legend_PT.SetEntrySeparation(0.5)  # NOVO PADRÃO
    legend_PT.SetBorderSize(0)
    legend_PT.SetFillStyle(0)
    hist_bjet_number.GetListOfFunctions().Remove(legend_PT)

    # Adicionando apenas o nome do processo com a cor
    legend_PT.AddEntry(hist_bjet_number, f"  {os.path.basename(inputFiles[i]).replace('100.root', '')}", "f")

    # Desenhando a legenda
    legend_PT.Draw()

canvas_bjet.Update()



# Salvando os histogramas em arquivos de imagem
output_file_path_HT2 = os.path.join(output_dir, "Nbjets.png")
canvas_bjet.SaveAs(output_file_path_HT2)
# Abra os histogramas quando prontos
subprocess.Popen(["xdg-open", output_file_path_HT2])






canvas_qjet.cd()
for i, hist_qjet_number in enumerate(histograms_q):
    hist_qjet_number.SetLineColor(colors[i])
    hist_qjet_number.SetLineWidth(5)

    # Verifica se é o último histograma
    if i != len(histograms_q) - 1:
        hist_qjet_number.SetFillColor(colors[i])
        hist_qjet_number.SetLineWidth(2)

    hist_qjet_number.Sumw2(0)
    hist_qjet_number.SetMinimum(1e-1)

    hist_qjet_number.Draw("SAME")
    canvas_qjet.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Criando legenda personalizada
#################################################################
    n_entries = len(histograms_q)  # Número de entradas na legenda

    # Posição inicial (topo) da legenda
    y1 = 0.956

    # Posição final (base) da legenda
    y2 = 0.01


    # Distribuindo uniformemente
    legend_PT = ROOT.TLegend(0.90, y1 - i * (y1 - y2) / n_entries, 1.0 , y1 - (i + 1) * (y1 - y2) / n_entries)
    legend_PT.SetTextSize(0.023)  # NOVO PADRÃO
    legend_PT.SetEntrySeparation(0.5)  # NOVO PADRÃO
    legend_PT.SetBorderSize(0)
    legend_PT.SetFillStyle(0)
    hist_qjet_number.GetListOfFunctions().Remove(legend_PT)

    # Adicionando apenas o nome do processo com a cor
    legend_PT.AddEntry(hist_qjet_number, f"  {os.path.basename(inputFiles[i]).replace('100.root', '')}", "f")

    # Desenhando a legenda
    legend_PT.Draw()


canvas_qjet.Update()



# Salvando os histogramas em arquivos de imagem
output_file_path_HT3 = os.path.join(output_dir, "Nqjets.png")
canvas_qjet.SaveAs(output_file_path_HT3)
# Abra os histogramas quando prontos
subprocess.Popen(["xdg-open", output_file_path_HT3])