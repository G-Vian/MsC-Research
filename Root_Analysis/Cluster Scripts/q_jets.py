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
output_dir = os.path.join(script_dir, "q_jets")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)



output_file_path1 = os.path.join(output_dir, "processing_log1.txt")

# Excluir o arquivo de log antigo, se existir
if os.path.exists(output_file_path1):
    os.remove(output_file_path1)


output_file_path2 = os.path.join(output_dir, "processing_log2.txt")

# Excluir o arquivo de log antigo, se existir
if os.path.exists(output_file_path2):
    os.remove(output_file_path2)


output_file_path3 = os.path.join(output_dir, "processing_log3.txt")

# Excluir o arquivo de log antigo, se existir
if os.path.exists(output_file_path3):
    os.remove(output_file_path3)


output_file_path4 = os.path.join(output_dir, "processing_log4.txt")

# Excluir o arquivo de log antigo, se existir
if os.path.exists(output_file_path4):
    os.remove(output_file_path4)
##############################################################
#Definindo pesos, pro hisgrama distribuir de acordo com a seção de choque:
##############################################################
#PESO = (SEÇÃO DE CHOQUE * LUMINOSIDADE)/(NUMERO TOTAL DE EVENTOS GERADOS (10000))

##############################################################
##############################################################
ROOT.gROOT.SetBatch(True) #NAo abre mais as imagens em formato root quando terminar


P_1 = 75.59214707  # tt
P_2 = 1.319608992  # ttbb
P_3 = 0.4078219283  # ttwh
P_4 = 0.2690111995  # ttbbbb
P_5 = 0.2683539397  # tth
P_6 = 0.2267598608  # ttwz
P_7 = 0.1155650301  # ttz
P_8 = 0.003341408  # tttt
P_9 = 0.0028875448  # ttww
P_10 = 0.0003030022  # tttw
P_11 = 0.0002389038  # tthh
P_12 = 0.0001228458  # ttzh
P_13 = 0.0000411384  # ttzz

lista_de_pesos = [
    P_1,  # tt
    P_2,  # ttbb
    P_4,  # ttbbbb
    P_5,  # tth
    P_3,  # ttwh
    P_12,  # ttzh
    P_7,  # ttz
    P_13,  # ttzz
    P_6,  # ttwz
    P_9,  # ttww
    P_8,  # tttt
    P_10,  # tttw
    P_11  # tthh
]
colors =  [  ROOT.kMagenta+3, ROOT.kOrange+4, ROOT.kGreen+3, ROOT.kAzure+3,  ROOT.kMagenta,  ROOT.kBlue, ROOT.kOrange+7, ROOT.kTeal-3,  ROOT.kBlack, ROOT.kYellow,  ROOT.kCyan, ROOT.kViolet, ROOT.kRed ] 
############################################################
######################## jq[1] ##############################
############################################################
#	LOADING TREES BLOCK
#
# Lista para armazenar os histogramas de PT e ETA de cada arquivo
modal_values_1 = []
modal_values_2 = []
N_entries1 = []
histograms_PT = []
histograms_ETA = []
mean_values_PT = []
mean_values_ETA = []
std_dev_values_PT = []
std_dev_values_ETA = []
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

with open(output_file_path1, "a") as output_file1:

    for i, inputFile in enumerate(inputFiles):
        print(f"Reading file {inputFile}")
    
        # Obter o peso correspondente ao arquivo atual
        peso_atual = lista_de_pesos[i]
        
        # Criando cadeia de árvores ROOT
        chain = ROOT.TChain("Delphes")
        chain.Add(inputFile)

        # Criando histogramas para PT e ETA dos jato b 
        histogram_PT = ROOT.TH1F("hist_pt", "P_{T} of the most energetic light-quark jet; P_{T} (GeV/c); Number of events (L = 200 fb^{-1})", 20,0, 1400)
        histogram_ETA = ROOT.TH1F("hist_pt", "#eta jq[1]; #eta ; Number of events (L = 200 fb^{-1})", 70, -5,5)
        
        # Criando objeto ExRootTreeReader para ler a cadeia de árvores
        treeReader = ROOT.ExRootTreeReader(chain)
        numberOfEntries = treeReader.GetEntries()
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

                   # Inicializando listas para armazenar os PTs e ETAs dos jato b
                    pt_b = []
                    eta_b = []
                    # Loop sobre as partículas
                    for jet in q_jets:
                        # Verificando se a partícula é um jato b
                        if jet in q_jets:
                            pt_b.append(jet.PT)
                            eta_b.append(jet.Eta)



                    # Verifica se há bottoms neste evento
                    if pt_b:
                        # Ordena os PTs dos jato b em ordem decrescente
                        pt_b.sort(reverse=True)

                        # Preenche o histograma com o PT do jato b mais energético
                        histogram_PT.Fill(pt_b[0], peso_atual)

                        # Preenche o histograma com o ETA do jato b mais energético
                        histogram_ETA.Fill(eta_b[pt_b.index(pt_b[0])], peso_atual)



            else:
                    continue  # Se a primeira condicional nao for satisfeita, vá para o próximo event
                
        # Obter o número total de eventos

        total_events_hist = histogram_PT.GetEntries()
        print(f"Total de eventos processados: {total_events}")
        print(f"Eventos que passaram MET >= 40, HT > 500 e nJets >= 6: {events_passed_MET_HT_nJets}")
        print(f"Eventos que passaram a seleção final: {events_passed_final_selection}")
        # Salvar os resultados em um arquivo de texto para o arquivo de entrada atual
        output_file1.write(f"Results for file: {inputFile}\n")
        output_file1.write(f"Total de eventos processados: {total_events}\n")
        output_file1.write(f"Eventos que passaram MET >= 40, HT > 500 e nJets >= 6: {events_passed_MET_HT_nJets}\n")
        output_file1.write(f"Eventos que passaram a seleção final: {events_passed_final_selection}\n")
        output_file1.write(f"Eventos que passaram a seleção final: {events_passed_final_selection}\n")
        output_file1.write(f"Eventos que entraram no histograma: {total_events_hist}\n")

        output_file1.write("-" * 40 + "\n")
                


                
        # Obter o número total de eventos

        total_events = histogram_ETA.GetEntries()


        # Adicionando histogramas à lista
        histograms_PT.append(histogram_PT)
        histograms_ETA.append(histogram_ETA)

    
    # Calculando média e desvio padrão 
        mean_1 = histogram_PT.GetMean()
        std_dev_1 = histogram_PT.GetStdDev()
        mean_2 = histogram_ETA.GetMean()
        std_dev_2 = histogram_ETA.GetStdDev()
    
    
        mean_values_PT.append(mean_1)
        std_dev_values_PT.append(std_dev_1)
        mean_values_ETA.append(mean_2)
        std_dev_values_ETA.append(std_dev_2)
        n_entries1 = 0  # Inicializa a variável
        n_entries1 = peso_atual*events_passed_final_selection

        # Valor modal (bin com o máximo número de entradas)
        max_bin_content = 0
        modal_value1 = 0

        # Iterar sobre todos os bins do histograma
        for bin in range(1, histogram_PT.GetNbinsX() + 1):
            bin_content = histogram_PT.GetBinContent(bin)
            if bin_content > max_bin_content:
                max_bin_content = bin_content
                modal_value1 = histogram_PT.GetBinCenter(bin)


            min_value1 = float('inf')
            max_value1 = -float('inf')

        # Iterar sobre os bins do histograma para encontrar o valor mínimo e máximo do eixo x
        for bin in range(1, histogram_PT.GetNbinsX() + 1):
            bin_content = histogram_PT.GetBinContent(bin)
            # Verificar se o conteúdo do bin é diferente de zero
            if bin_content != 0:
                bin_x_value = histogram_PT.GetXaxis().GetBinLowEdge(bin)
                min_value1 = min(min_value1, bin_x_value)
                max_value1 = max(max_value1, histogram_PT.GetXaxis().GetBinUpEdge(bin))


        modal_values_1.append(modal_value1)
        N_entries1.append(n_entries1)
        min_values1.append(min_value1)
        max_values1.append(max_value1)

        n_entries2 = peso_atual*events_passed_final_selection

        # Valor modal (bin com o máximo número de entradas)
        max_bin_content = 0
        modal_value2 = 0

        # Iterar sobre todos os bins do histograma
        for bin in range(1, histogram_ETA.GetNbinsX() + 1):
            bin_content = histogram_ETA.GetBinContent(bin)
            if bin_content > max_bin_content:
                max_bin_content = bin_content
                modal_value2 = histogram_ETA.GetBinCenter(bin)


            min_value2 = float('inf')
            max_value2 = -float('inf')

        # Iterar sobre os bins do histograma para encontrar o valor mínimo e máximo do eixo x
        for bin in range(1, histogram_ETA.GetNbinsX() + 1):
            bin_content = histogram_ETA.GetBinContent(bin)
            # Verificar se o conteúdo do bin é diferente de zero
            if bin_content != 0:
                bin_x_value = histogram_ETA.GetXaxis().GetBinLowEdge(bin)
                min_value2 = min(min_value2, bin_x_value)
                max_value2 = max(max_value2, histogram_PT.GetXaxis().GetBinUpEdge(bin))


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

    # Criação dos histogramas combinados (Código Novo)
    hist_tt_nb = histograms_PT[0].Clone("hist_tt_nb")
    hist_tth_x = histograms_PT[3].Clone("hist_tth_x")
    hist_tt_nV = histograms_PT[6].Clone("hist_tt_nV")
    hist_ttt_x = histograms_PT[10].Clone("hist_ttt_x")

    # Somando histogramas conforme as categorias (Código Novo)
    for i in range(1, 3):
        hist_tt_nb.Add(histograms_PT[i])

    for i in range(4, 6):
        hist_tth_x.Add(histograms_PT[i])

    for i in range(7, 10):
        hist_tt_nV.Add(histograms_PT[i])

    for i in range(11, 12):
        hist_ttt_x.Add(histograms_PT[i])

    # Criando canvas para os histogramas (Código Novo)
    canvas_PT = ROOT.TCanvas("canvas_HT", "h1", 2200, 1300)
    canvas_PT.SetLeftMargin(0.1)
    canvas_PT.SetRightMargin(0.22)
    canvas_PT.SetBottomMargin(0.15)
    canvas_PT.SetTopMargin(0.1)
    canvas_PT.SetLogy()

    # Definir cores para cada categoria (Código Novo)
    colors_dict = {
        "tt+nb": ROOT.kBlack,
        "tth+X": ROOT.kBlue,
        "tt+nV": ROOT.kGreen,
        "ttt+X": ROOT.kMagenta,
        "Signal": ROOT.kRed
    }

    # Ajustando histogramas combinados (Código Novo)
    for hist, label, color in [
        (hist_tt_nb, "tt+nb", colors_dict["tt+nb"]),
        (hist_tth_x, "tth+X", colors_dict["tth+X"]),
        (hist_tt_nV, "tt+nV", colors_dict["tt+nV"]),
        (hist_ttt_x, "ttt+X", colors_dict["ttt+X"])
    ]:
        hist.SetLineColor(color)
        hist.SetFillColor(color)
        hist.SetLineWidth(2)
        hist.Sumw2(0)  # Desativar incertezas
        hist.Draw("SAME")

    # Último histograma sem preenchimento (Código Novo)
    histograms_PT[-1].SetLineColor(colors_dict["Signal"])
    histograms_PT[-1].SetLineWidth(5)
    histograms_PT[-1].Sumw2(0)  # Desativar incertezas
    histograms_PT[-1].Draw("SAME")

    # Criando legenda (Código Novo)
    legend = ROOT.TLegend(0.79, 0.75, 0.95, 0.94)
    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    # Adicionando entradas na legenda (Código Novo)
    # Adicionando entradas na legenda (Código Novo)
    for label, hist in [
        ("tt + tt2b + tt4b", hist_tt_nb),
        ("tth + tthZ + tthW", hist_tth_x),
        ("ttWW + ttZZ + ttWZ+ ttZ", hist_tt_nV),
        ("tttW + tttt", hist_ttt_x),
        ("Signal", histograms_PT[-1])
    ]:
        legend.AddEntry(hist, label, "f")



    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.04)  # Ajuste o tamanho do texto para menor
    text.DrawLatex(0.615, 0.83, R1)


    legend.Draw()

canvas_PT.Update()
output_file_path_PT = os.path.join(output_dir, "jq[1]_pt.png")
canvas_PT.SaveAs(output_file_path_PT)


############################################################
######################## jq[2] ##############################
############################################################
#	LOADING TREES BLOCK
#
# Lista para armazenar os histogramas de PT e ETA de cada arquivo
modal_values_1 = []
modal_values_2 = []
N_entries1 = []
histograms_PT = []
histograms_ETA = []
mean_values_PT = []
mean_values_ETA = []
std_dev_values_PT = []
std_dev_values_ETA = []
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

# Criando histogramas para PT e ETA dos jato b




with open(output_file_path2, "a") as output_file2:

    for i, inputFile in enumerate(inputFiles):
        print(f"Reading file {inputFile}")
    
        # Obter o peso correspondente ao arquivo atual
        peso_atual = lista_de_pesos[i]
        
        # Criando cadeia de árvores ROOT
        chain = ROOT.TChain("Delphes")
        chain.Add(inputFile)

        # Criando histogramas para PT e ETA dos jato b 
        histogram_PT = ROOT.TH1F("hist_pt", "P_{T} of the 2nd most energetic light-quark jet; P_{T} (GeV/c); Number of events (L = 200 fb^{-1})", 20,0, 800)
        histogram_ETA = ROOT.TH1F("hist_pt", "#eta jq[2]; #eta ; Number of events (L = 200 fb^{-1})", 70, -5,5)
        
        # Criando objeto ExRootTreeReader para ler a cadeia de árvores
        treeReader = ROOT.ExRootTreeReader(chain)
        numberOfEntries = treeReader.GetEntries()
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

                   # Inicializando listas para armazenar os PTs e ETAs dos jato b
                    pt_b = []
                    eta_b = []
                    # Loop sobre as partículas
                    for jet in q_jets:
                        # Verificando se a partícula é um jato b
                        if jet in q_jets:
                            pt_b.append(jet.PT)
                            eta_b.append(jet.Eta)

                    if pt_b and len(pt_b) >= 2 :
                    # Ordena os PTs dos jato b em ordem decrescente
                        pt_b.sort(reverse=True)

                        
            
                        # Preenche o histograma com o PT do segundo jato b mais energético
                        histogram_PT.Fill(pt_b[1], peso_atual)

                        # Preenche o histograma com o ETA do segundo jato b mais energético
                        histogram_ETA.Fill(eta_b[pt_b.index(pt_b[1])], peso_atual)


            else:
                    continue  # Se a primeira condicional nao for satisfeita, vá para o próximo event
                
        # Obter o número total de eventos

        total_events_hist = histogram_PT.GetEntries()
        print(f"Total de eventos processados: {total_events}")
        print(f"Eventos que passaram MET >= 40, HT > 500 e nJets >= 6: {events_passed_MET_HT_nJets}")
        print(f"Eventos que passaram a seleção final: {events_passed_final_selection}")
        # Salvar os resultados em um arquivo de texto para o arquivo de entrada atual
        output_file2.write(f"Results for file: {inputFile}\n")
        output_file2.write(f"Total de eventos processados: {total_events}\n")
        output_file2.write(f"Eventos que passaram MET >= 40, HT > 500 e nJets >= 6: {events_passed_MET_HT_nJets}\n")
        output_file2.write(f"Eventos que passaram a seleção final: {events_passed_final_selection}\n")
        output_file2.write(f"Eventos que passaram a seleção final: {events_passed_final_selection}\n")
        output_file2.write(f"Eventos que entraram no histograma: {total_events_hist}\n")

        output_file2.write("-" * 40 + "\n")
        # Verifica se há bottoms neste evento


            
        total_events = histogram_ETA.GetEntries()


        # Adicionando histogramas à lista
        histograms_PT.append(histogram_PT)
        histograms_ETA.append(histogram_ETA)

    # Calculando média e desvio padrão 
        mean_1 = histogram_PT.GetMean()
        std_dev_1 = histogram_PT.GetStdDev()
        mean_2 = histogram_ETA.GetMean()
        std_dev_2 = histogram_ETA.GetStdDev()
    
    
        mean_values_PT.append(mean_1)
        std_dev_values_PT.append(std_dev_1)
        mean_values_ETA.append(mean_2)
        std_dev_values_ETA.append(std_dev_2)
        
        n_entries1 = 0  # Inicializa a variável
        n_entries1 = peso_atual*events_passed_final_selection


        # Valor modal (bin com o máximo número de entradas)
        max_bin_content = 0
        modal_value1 = 0

        # Iterar sobre todos os bins do histograma
        for bin in range(1, histogram_PT.GetNbinsX() + 1):
            bin_content = histogram_PT.GetBinContent(bin)
            if bin_content > max_bin_content:
                max_bin_content = bin_content
                modal_value1 = histogram_PT.GetBinCenter(bin)


            min_value1 = float('inf')
            max_value1 = -float('inf')

        # Iterar sobre os bins do histograma para encontrar o valor mínimo e máximo do eixo x
        for bin in range(1, histogram_PT.GetNbinsX() + 1):
            bin_content = histogram_PT.GetBinContent(bin)
            # Verificar se o conteúdo do bin é diferente de zero
            if bin_content != 0:
                bin_x_value = histogram_PT.GetXaxis().GetBinLowEdge(bin)
                min_value1 = min(min_value1, bin_x_value)
                max_value1 = max(max_value1, histogram_PT.GetXaxis().GetBinUpEdge(bin))


        modal_values_1.append(modal_value1)
        N_entries1.append(n_entries1)
        min_values1.append(min_value1)
        max_values1.append(max_value1)

        n_entries2 = peso_atual*events_passed_final_selection

        # Valor modal (bin com o máximo número de entradas)
        max_bin_content = 0
        modal_value2 = 0

        # Iterar sobre todos os bins do histograma
        for bin in range(1, histogram_ETA.GetNbinsX() + 1):
            bin_content = histogram_ETA.GetBinContent(bin)
            if bin_content > max_bin_content:
                max_bin_content = bin_content
                modal_value2 = histogram_ETA.GetBinCenter(bin)


            min_value2 = float('inf')
            max_value2 = -float('inf')

        # Iterar sobre os bins do histograma para encontrar o valor mínimo e máximo do eixo x
        for bin in range(1, histogram_ETA.GetNbinsX() + 1):
            bin_content = histogram_ETA.GetBinContent(bin)
            # Verificar se o conteúdo do bin é diferente de zero
            if bin_content != 0:
                bin_x_value = histogram_ETA.GetXaxis().GetBinLowEdge(bin)
                min_value2 = min(min_value2, bin_x_value)
                max_value2 = max(max_value2, histogram_PT.GetXaxis().GetBinUpEdge(bin))


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

    # Criação dos histogramas combinados (Código Novo)
    hist_tt_nb = histograms_PT[0].Clone("hist_tt_nb")
    hist_tth_x = histograms_PT[3].Clone("hist_tth_x")
    hist_tt_nV = histograms_PT[6].Clone("hist_tt_nV")
    hist_ttt_x = histograms_PT[10].Clone("hist_ttt_x")

    # Somando histogramas conforme as categorias (Código Novo)
    for i in range(1, 3):
        hist_tt_nb.Add(histograms_PT[i])

    for i in range(4, 6):
        hist_tth_x.Add(histograms_PT[i])

    for i in range(7, 10):
        hist_tt_nV.Add(histograms_PT[i])

    for i in range(11, 12):
        hist_ttt_x.Add(histograms_PT[i])

    # Criando canvas para os histogramas (Código Novo)
    canvas_PT = ROOT.TCanvas("canvas_HT", "h1", 2200, 1300)
    canvas_PT.SetLeftMargin(0.1)
    canvas_PT.SetRightMargin(0.22)
    canvas_PT.SetBottomMargin(0.15)
    canvas_PT.SetTopMargin(0.1)
    canvas_PT.SetLogy()

    # Definir cores para cada categoria (Código Novo)
    colors_dict = {
        "tt+nb": ROOT.kBlack,
        "tth+X": ROOT.kBlue,
        "tt+nV": ROOT.kGreen,
        "ttt+X": ROOT.kMagenta,
        "Signal": ROOT.kRed
    }

    # Ajustando histogramas combinados (Código Novo)
    for hist, label, color in [
        (hist_tt_nb, "tt+nb", colors_dict["tt+nb"]),
        (hist_tth_x, "tth+X", colors_dict["tth+X"]),
        (hist_tt_nV, "tt+nV", colors_dict["tt+nV"]),
        (hist_ttt_x, "ttt+X", colors_dict["ttt+X"])
    ]:
        hist.SetLineColor(color)
        hist.SetFillColor(color)
        hist.SetLineWidth(2)
        hist.Sumw2(0)  # Desativar incertezas
        hist.Draw("SAME")

    # Último histograma sem preenchimento (Código Novo)
    histograms_PT[-1].SetLineColor(colors_dict["Signal"])
    histograms_PT[-1].SetLineWidth(5)
    histograms_PT[-1].Sumw2(0)  # Desativar incertezas
    histograms_PT[-1].Draw("SAME")

    # Criando legenda (Código Novo)
    legend = ROOT.TLegend(0.79, 0.75, 0.95, 0.94)
    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    # Adicionando entradas na legenda (Código Novo)
    # Adicionando entradas na legenda (Código Novo)
    for label, hist in [
        ("tt + tt2b + tt4b", hist_tt_nb),
        ("tth + tthZ + tthW", hist_tth_x),
        ("ttWW + ttZZ + ttWZ+ ttZ", hist_tt_nV),
        ("tttW + tttt", hist_ttt_x),
        ("Signal", histograms_PT[-1])
    ]:
        legend.AddEntry(hist, label, "f")



    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.04)  # Ajuste o tamanho do texto para menor
    text.DrawLatex(0.615, 0.83, R1)


    legend.Draw()

canvas_PT.Update()
output_file_path_PT = os.path.join(output_dir, "jq[2]_pt.png")
canvas_PT.SaveAs(output_file_path_PT)
############################################################
######################## jq[3} ##############################
############################################################

#	LOADING TREES BLOCK
#
# Lista para armazenar os histogramas de PT e ETA de cada arquivo
modal_values_1 = []
modal_values_2 = []
N_entries1 = []
histograms_PT = []
histograms_ETA = []
mean_values_PT = []
mean_values_ETA = []
std_dev_values_PT = []
std_dev_values_ETA = []
min_values1 = []
max_values1= []

N_entries2 = []
min_values2 = []
max_values2 = []

# Desativar as estatísticas padrão
ROOT.gStyle.SetOptStat(0)

# Ativar a exibição do valor médio e desvio padrão personalizado
ROOT.gStyle.SetOptFit(1111)


# Criando histogramas para PT e ETA dos jato b




with open(output_file_path3, "a") as output_file3:

    for i, inputFile in enumerate(inputFiles):
        print(f"Reading file {inputFile}")
    
        # Obter o peso correspondente ao arquivo atual
        peso_atual = lista_de_pesos[i]
        
        # Criando cadeia de árvores ROOT
        chain = ROOT.TChain("Delphes")
        chain.Add(inputFile)

        # Criando histogramas para PT e ETA dos jato b 
        histogram_PT = ROOT.TH1F("hist_pt", "P_{T} of the 3rd most energetic light-quark jet; P_{T} (GeV/c); Number of events (L = 200 fb^{-1})", 20,0, 500)
        histogram_ETA = ROOT.TH1F("hist_pt", "#eta jq[3]; #eta ; Number of events (L = 200 fb^{-1})", 70, -5,5)
        
        # Criando objeto ExRootTreeReader para ler a cadeia de árvores
        treeReader = ROOT.ExRootTreeReader(chain)
        numberOfEntries = treeReader.GetEntries()
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

                   # Inicializando listas para armazenar os PTs e ETAs dos jato b
                    pt_b = []
                    eta_b = []
                    # Loop sobre as partículas
                    for jet in q_jets:
                        # Verificando se a partícula é um jato b
                        if jet in q_jets:
                            pt_b.append(jet.PT)
                            eta_b.append(jet.Eta)

                    if pt_b and len(pt_b) >= 3 :
                    # Ordena os PTs dos jato b em ordem decrescente
                        pt_b.sort(reverse=True)

                        
            
                        # Preenche o histograma com o PT do teceiro jato b mais energético
                        histogram_PT.Fill(pt_b[2], peso_atual)

                        # Preenche o histograma com o ETA do terceiro jato b mais energético
                        histogram_ETA.Fill(eta_b[pt_b.index(pt_b[2])], peso_atual)

            else:
                    continue  # Se a primeira condicional nao for satisfeita, vá para o próximo event
                
        # Obter o número total de eventos

        total_events_hist = histogram_PT.GetEntries()
        print(f"Total de eventos processados: {total_events}")
        print(f"Eventos que passaram MET >= 40, HT > 500 e nJets >= 6: {events_passed_MET_HT_nJets}")
        print(f"Eventos que passaram a seleção final: {events_passed_final_selection}")
        # Salvar os resultados em um arquivo de texto para o arquivo de entrada atual
        output_file3.write(f"Results for file: {inputFile}\n")
        output_file3.write(f"Total de eventos processados: {total_events}\n")
        output_file3.write(f"Eventos que passaram MET >= 40, HT > 500 e nJets >= 6: {events_passed_MET_HT_nJets}\n")
        output_file3.write(f"Eventos que passaram a seleção final: {events_passed_final_selection}\n")
        output_file3.write(f"Eventos que passaram a seleção final: {events_passed_final_selection}\n")
        output_file3.write(f"Eventos que entraram no histograma: {total_events_hist}\n")

        output_file3.write("-" * 40 + "\n")



        # Obter o número total de eventos

        total_events = histogram_ETA.GetEntries()


        # Adicionando histogramas à lista
        histograms_PT.append(histogram_PT)
        histograms_ETA.append(histogram_ETA)
    # Calculando média e desvio padrão 
        mean_1 = histogram_PT.GetMean()
        std_dev_1 = histogram_PT.GetStdDev()
        mean_2 = histogram_ETA.GetMean()
        std_dev_2 = histogram_ETA.GetStdDev()
    
    
        mean_values_PT.append(mean_1)
        std_dev_values_PT.append(std_dev_1)
        mean_values_ETA.append(mean_2)
        std_dev_values_ETA.append(std_dev_2)

        n_entries1 = peso_atual*events_passed_final_selection


        # Valor modal (bin com o máximo número de entradas)
        max_bin_content = 0
        modal_value1 = 0

        # Iterar sobre todos os bins do histograma
        for bin in range(1, histogram_PT.GetNbinsX() + 1):
            bin_content = histogram_PT.GetBinContent(bin)
            if bin_content > max_bin_content:
                max_bin_content = bin_content
                modal_value1 = histogram_PT.GetBinCenter(bin)


            min_value1 = float('inf')
            max_value1 = -float('inf')

        # Iterar sobre os bins do histograma para encontrar o valor mínimo e máximo do eixo x
        for bin in range(1, histogram_PT.GetNbinsX() + 1):
            bin_content = histogram_PT.GetBinContent(bin)
            # Verificar se o conteúdo do bin é diferente de zero
            if bin_content != 0:
                bin_x_value = histogram_PT.GetXaxis().GetBinLowEdge(bin)
                min_value1 = min(min_value1, bin_x_value)
                max_value1 = max(max_value1, histogram_PT.GetXaxis().GetBinUpEdge(bin))


        modal_values_1.append(modal_value1)
        N_entries1.append(n_entries1)
        min_values1.append(min_value1)
        max_values1.append(max_value1)

        n_entries2 = peso_atual*events_passed_final_selection

        # Valor modal (bin com o máximo número de entradas)
        max_bin_content = 0
        modal_value2 = 0

        # Iterar sobre todos os bins do histograma
        for bin in range(1, histogram_ETA.GetNbinsX() + 1):
            bin_content = histogram_ETA.GetBinContent(bin)
            if bin_content > max_bin_content:
                max_bin_content = bin_content
                modal_value2 = histogram_ETA.GetBinCenter(bin)


            min_value2 = float('inf')
            max_value2 = -float('inf')

        # Iterar sobre os bins do histograma para encontrar o valor mínimo e máximo do eixo x
        for bin in range(1, histogram_ETA.GetNbinsX() + 1):
            bin_content = histogram_ETA.GetBinContent(bin)
            # Verificar se o conteúdo do bin é diferente de zero
            if bin_content != 0:
                bin_x_value = histogram_ETA.GetXaxis().GetBinLowEdge(bin)
                min_value2 = min(min_value2, bin_x_value)
                max_value2 = max(max_value2, histogram_PT.GetXaxis().GetBinUpEdge(bin))


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

    # Criação dos histogramas combinados (Código Novo)
    hist_tt_nb = histograms_PT[0].Clone("hist_tt_nb")
    hist_tth_x = histograms_PT[3].Clone("hist_tth_x")
    hist_tt_nV = histograms_PT[6].Clone("hist_tt_nV")
    hist_ttt_x = histograms_PT[10].Clone("hist_ttt_x")

    # Somando histogramas conforme as categorias (Código Novo)
    for i in range(1, 3):
        hist_tt_nb.Add(histograms_PT[i])

    for i in range(4, 6):
        hist_tth_x.Add(histograms_PT[i])

    for i in range(7, 10):
        hist_tt_nV.Add(histograms_PT[i])

    for i in range(11, 12):
        hist_ttt_x.Add(histograms_PT[i])

    # Criando canvas para os histogramas (Código Novo)
    canvas_PT = ROOT.TCanvas("canvas_HT", "h1", 2200, 1300)
    canvas_PT.SetLeftMargin(0.1)
    canvas_PT.SetRightMargin(0.22)
    canvas_PT.SetBottomMargin(0.15)
    canvas_PT.SetTopMargin(0.1)
    canvas_PT.SetLogy()

    # Definir cores para cada categoria (Código Novo)
    colors_dict = {
        "tt+nb": ROOT.kBlack,
        "tth+X": ROOT.kBlue,
        "tt+nV": ROOT.kGreen,
        "ttt+X": ROOT.kMagenta,
        "Signal": ROOT.kRed
    }

    # Ajustando histogramas combinados (Código Novo)
    for hist, label, color in [
        (hist_tt_nb, "tt+nb", colors_dict["tt+nb"]),
        (hist_tth_x, "tth+X", colors_dict["tth+X"]),
        (hist_tt_nV, "tt+nV", colors_dict["tt+nV"]),
        (hist_ttt_x, "ttt+X", colors_dict["ttt+X"])
    ]:
        hist.SetLineColor(color)
        hist.SetFillColor(color)
        hist.SetLineWidth(2)
        hist.Sumw2(0)  # Desativar incertezas
        hist.Draw("SAME")

    # Último histograma sem preenchimento (Código Novo)
    histograms_PT[-1].SetLineColor(colors_dict["Signal"])
    histograms_PT[-1].SetLineWidth(5)
    histograms_PT[-1].Sumw2(0)  # Desativar incertezas
    histograms_PT[-1].Draw("SAME")

    # Criando legenda (Código Novo)
    legend = ROOT.TLegend(0.79, 0.75, 0.95, 0.94)
    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    # Adicionando entradas na legenda (Código Novo)
    # Adicionando entradas na legenda (Código Novo)
    for label, hist in [
        ("tt + tt2b + tt4b", hist_tt_nb),
        ("tth + tthZ + tthW", hist_tth_x),
        ("ttWW + ttZZ + ttWZ+ ttZ", hist_tt_nV),
        ("tttW + tttt", hist_ttt_x),
        ("Signal", histograms_PT[-1])
    ]:
        legend.AddEntry(hist, label, "f")



    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.04)  # Ajuste o tamanho do texto para menor
    text.DrawLatex(0.615, 0.83, R1)


    legend.Draw()

canvas_PT.Update()
output_file_path_PT = os.path.join(output_dir, "jq[3]_pt.png")
canvas_PT.SaveAs(output_file_path_PT)


############################################################
######################## jq[4] ##############################
############################################################
#	LOADING TREES BLOCK
#
# Lista para armazenar os histogramas de PT e ETA de cada arquivo
modal_values_1 = []
modal_values_2 = []
N_entries1 = []
histograms_PT = []
histograms_ETA = []
mean_values_PT = []
mean_values_ETA = []
std_dev_values_PT = []
std_dev_values_ETA = []
min_values1 = []
max_values1= []

N_entries2 = []
min_values2 = []
max_values2 = []

# Desativar as estatísticas padrão
ROOT.gStyle.SetOptStat(0)

# Ativar a exibição do valor médio e desvio padrão personalizado
ROOT.gStyle.SetOptFit(1111)




with open(output_file_path3, "a") as output_file3:

    for i, inputFile in enumerate(inputFiles):
        print(f"Reading file {inputFile}")
    
        # Obter o peso correspondente ao arquivo atual
        peso_atual = lista_de_pesos[i]
        
        # Criando cadeia de árvores ROOT
        chain = ROOT.TChain("Delphes")
        chain.Add(inputFile)

        # Criando histogramas para PT e ETA dos jato b 
        histogram_PT = ROOT.TH1F("hist_pt", "P_{T} of the 4th most energetic light-quark jet; P_{T} (GeV/c); Number of events (L = 200 fb^{-1})", 20,0, 250)
        histogram_ETA = ROOT.TH1F("hist_pt", "#eta jq[4]; #eta ; Number of events (L = 200 fb^{-1})", 70, -5,5)
        
        # Criando objeto ExRootTreeReader para ler a cadeia de árvores
        treeReader = ROOT.ExRootTreeReader(chain)
        numberOfEntries = treeReader.GetEntries()
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

                   # Inicializando listas para armazenar os PTs e ETAs dos jato b
                    pt_b = []
                    eta_b = []
                    # Loop sobre as partículas
                    for jet in q_jets:
                        # Verificando se a partícula é um jato b
                        if jet in q_jets:
                            pt_b.append(jet.PT)
                            eta_b.append(jet.Eta)

                                # Verifica se há bottoms neste evento
                    if pt_b and len(pt_b) >= 4 :
                    # Ordena os PTs dos jato b em ordem decrescente
                        pt_b.sort(reverse=True)

                        
            
                        # Preenche o histograma com o PT do 4º jato b mais energético
                        histogram_PT.Fill(pt_b[3], peso_atual)

                        # Preenche o histograma com o ETA do 4º jato b mais energético
                        histogram_ETA.Fill(eta_b[pt_b.index(pt_b[3])], peso_atual)
            else:
                    continue  # Se a primeira condicional nao for satisfeita, vá para o próximo event
                
        # Obter o número total de eventos

        total_events_hist = histogram_PT.GetEntries()
        print(f"Total de eventos processados: {total_events}")
        print(f"Eventos que passaram MET >= 40, HT > 500 e nJets >= 6: {events_passed_MET_HT_nJets}")
        print(f"Eventos que passaram a seleção final: {events_passed_final_selection}")
        # Salvar os resultados em um arquivo de texto para o arquivo de entrada atual
        output_file3.write(f"Results for file: {inputFile}\n")
        output_file3.write(f"Total de eventos processados: {total_events}\n")
        output_file3.write(f"Eventos que passaram MET >= 40, HT > 500 e nJets >= 6: {events_passed_MET_HT_nJets}\n")
        output_file3.write(f"Eventos que passaram a seleção final: {events_passed_final_selection}\n")
        output_file3.write(f"Eventos que passaram a seleção final: {events_passed_final_selection}\n")
        output_file3.write(f"Eventos que entraram no histograma: {total_events_hist}\n")

        output_file3.write("-" * 40 + "\n")
        
        
        # Obter o número total de eventos

        total_events = histogram_ETA.GetEntries()


        # Adicionando histogramas à lista
        histograms_PT.append(histogram_PT)
        histograms_ETA.append(histogram_ETA)
    # Calculando média e desvio padrão 
        mean_1 = histogram_PT.GetMean()
        std_dev_1 = histogram_PT.GetStdDev()
        mean_2 = histogram_ETA.GetMean()
        std_dev_2 = histogram_ETA.GetStdDev()
    
    
        mean_values_PT.append(mean_1)
        std_dev_values_PT.append(std_dev_1)
        mean_values_ETA.append(mean_2)
        std_dev_values_ETA.append(std_dev_2)

        n_entries1 = peso_atual*events_passed_final_selection


        # Valor modal (bin com o máximo número de entradas)
        max_bin_content = 0
        modal_value1 = 0

        # Iterar sobre todos os bins do histograma
        for bin in range(1, histogram_PT.GetNbinsX() + 1):
            bin_content = histogram_PT.GetBinContent(bin)
            if bin_content > max_bin_content:
                max_bin_content = bin_content
                modal_value1 = histogram_PT.GetBinCenter(bin)


            min_value1 = float('inf')
            max_value1 = -float('inf')

        # Iterar sobre os bins do histograma para encontrar o valor mínimo e máximo do eixo x
        for bin in range(1, histogram_PT.GetNbinsX() + 1):
            bin_content = histogram_PT.GetBinContent(bin)
            # Verificar se o conteúdo do bin é diferente de zero
            if bin_content != 0:
                bin_x_value = histogram_PT.GetXaxis().GetBinLowEdge(bin)
                min_value1 = min(min_value1, bin_x_value)
                max_value1 = max(max_value1, histogram_PT.GetXaxis().GetBinUpEdge(bin))


        modal_values_1.append(modal_value1)
        N_entries1.append(n_entries1)
        min_values1.append(min_value1)
        max_values1.append(max_value1)

        n_entries2 = peso_atual*events_passed_final_selection

        # Valor modal (bin com o máximo número de entradas)
        max_bin_content = 0
        modal_value2 = 0

        # Iterar sobre todos os bins do histograma
        for bin in range(1, histogram_ETA.GetNbinsX() + 1):
            bin_content = histogram_ETA.GetBinContent(bin)
            if bin_content > max_bin_content:
                max_bin_content = bin_content
                modal_value2 = histogram_ETA.GetBinCenter(bin)


            min_value2 = float('inf')
            max_value2 = -float('inf')

        # Iterar sobre os bins do histograma para encontrar o valor mínimo e máximo do eixo x
        for bin in range(1, histogram_ETA.GetNbinsX() + 1):
            bin_content = histogram_ETA.GetBinContent(bin)
            # Verificar se o conteúdo do bin é diferente de zero
            if bin_content != 0:
                bin_x_value = histogram_ETA.GetXaxis().GetBinLowEdge(bin)
                min_value2 = min(min_value2, bin_x_value)
                max_value2 = max(max_value2, histogram_PT.GetXaxis().GetBinUpEdge(bin))


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

    # Criação dos histogramas combinados (Código Novo)
    hist_tt_nb = histograms_PT[0].Clone("hist_tt_nb")
    hist_tth_x = histograms_PT[3].Clone("hist_tth_x")
    hist_tt_nV = histograms_PT[6].Clone("hist_tt_nV")
    hist_ttt_x = histograms_PT[10].Clone("hist_ttt_x")

    # Somando histogramas conforme as categorias (Código Novo)
    for i in range(1, 3):
        hist_tt_nb.Add(histograms_PT[i])

    for i in range(4, 6):
        hist_tth_x.Add(histograms_PT[i])

    for i in range(7, 10):
        hist_tt_nV.Add(histograms_PT[i])

    for i in range(11, 12):
        hist_ttt_x.Add(histograms_PT[i])

    # Criando canvas para os histogramas (Código Novo)
    canvas_PT = ROOT.TCanvas("canvas_HT", "h1", 2200, 1300)
    canvas_PT.SetLeftMargin(0.1)
    canvas_PT.SetRightMargin(0.22)
    canvas_PT.SetBottomMargin(0.15)
    canvas_PT.SetTopMargin(0.1)
    canvas_PT.SetLogy()

    # Definir cores para cada categoria (Código Novo)
    colors_dict = {
        "tt+nb": ROOT.kBlack,
        "tth+X": ROOT.kBlue,
        "tt+nV": ROOT.kGreen,
        "ttt+X": ROOT.kMagenta,
        "Signal": ROOT.kRed
    }

    # Ajustando histogramas combinados (Código Novo)
    for hist, label, color in [
        (hist_tt_nb, "tt+nb", colors_dict["tt+nb"]),
        (hist_tth_x, "tth+X", colors_dict["tth+X"]),
        (hist_tt_nV, "tt+nV", colors_dict["tt+nV"]),
        (hist_ttt_x, "ttt+X", colors_dict["ttt+X"])
    ]:
        hist.SetLineColor(color)
        hist.SetFillColor(color)
        hist.SetLineWidth(2)
        hist.Sumw2(0)  # Desativar incertezas
        hist.Draw("SAME")

    # Último histograma sem preenchimento (Código Novo)
    histograms_PT[-1].SetLineColor(colors_dict["Signal"])
    histograms_PT[-1].SetLineWidth(5)
    histograms_PT[-1].Sumw2(0)  # Desativar incertezas
    histograms_PT[-1].Draw("SAME")

    # Criando legenda (Código Novo)
    legend = ROOT.TLegend(0.79, 0.75, 0.95, 0.94)
    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    # Adicionando entradas na legenda (Código Novo)
    # Adicionando entradas na legenda (Código Novo)
    for label, hist in [
        ("tt + tt2b + tt4b", hist_tt_nb),
        ("tth + tthZ + tthW", hist_tth_x),
        ("ttWW + ttZZ + ttWZ+ ttZ", hist_tt_nV),
        ("tttW + tttt", hist_ttt_x),
        ("Signal", histograms_PT[-1])
    ]:
        legend.AddEntry(hist, label, "f")



    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.04)  # Ajuste o tamanho do texto para menor
    text.DrawLatex(0.615, 0.83, R1)


    legend.Draw()

canvas_PT.Update()
output_file_path_PT = os.path.join(output_dir, "jq[4]_pt.png")
canvas_PT.SaveAs(output_file_path_PT)