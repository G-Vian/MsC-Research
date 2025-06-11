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
        hist_jet_number = ROOT.TH1F("hist_jet_number", "Number of jets per event; Number of Jets; Number of events (normalized to one)", 16, -1, 16)
        hist_b_number = ROOT.TH1F("hist_bjet_number", "Number of b-tagged jets per event; Number of b-quark Jets; Number of events (normalized to one)", 10, -1, 10)
        hist_q_number = ROOT.TH1F("hist_qjet_number", "Number of light-quark jets per event; Number of light-quark Jets; Number of events (normalized to one)", 11, -1, 11)

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
            if nJets >= 0 and HT_value > 500 : 
                events_passed_MET_HT_nJets += 1 #pra eu contar quantos eventos passaram
                selected_jets = [jet for jet in branchJet if jet.PT >= 40 and abs(jet.Eta) <= 2.5 ]        
                q_jets = [jet for jet in selected_jets if abs(jet.Flavor)  in [0, 1, 2, 3, 4]  ] 
                b_jets = [jet for jet in selected_jets if abs(jet.Flavor) == 5]
                electrons = [electron for electron in branchelec if electron.PT > 20 and abs(electron.Eta) <= 2.5  ] 
                muons = [muon for muon in branchmu if muon.PT > 20 and abs(muon.Eta) <= 2.5 ] 

            #SEGUNDA CONDICIONAL, DENTRO DA PRIMEIRA, QUE OBSERVA O NUMERO DE OBJETOS
                if len(q_jets) >= 0 and len(b_jets) >= 0 and len(muons) == 0 and len(electrons) == 0 : 
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

    
        n_entriesj = peso_atual*hist_jet_number.GetEntries()
        
        
        n_entriesb = peso_atual*hist_b_number.GetEntries()
        
        
        n_entriesq = peso_atual*hist_q_number.GetEntries()




        n_entries_jets.append(n_entriesj)
        n_entries_bjets.append(n_entriesb)
        n_entries_qjets.append(n_entriesq)


        # Adicionando histogramas à lista
        histograms_jets.append(hist_jet_number)
        histograms_b.append(hist_b_number)
        histograms_q.append(hist_q_number)

    # Criação dos histogramas combinados (Código Novo)
    hist_tt_nb1 = histograms_jets[0].Clone("hist_tt_nb")
    hist_tth_x1 = histograms_jets[3].Clone("hist_tth_x")
    hist_tt_nV1 = histograms_jets[6].Clone("hist_tt_nV")
    hist_ttt_x1 = histograms_jets[10].Clone("hist_ttt_x")

    # Somando histogramas conforme as categorias (Código Novo)
    for i in range(1, 3):
        hist_tt_nb1.Add(histograms_jets[i])

    for i in range(4, 6):
        hist_tth_x1.Add(histograms_jets[i])

    for i in range(7, 10):
        hist_tt_nV1.Add(histograms_jets[i])

    for i in range(11, 12):
        hist_ttt_x1.Add(histograms_jets[i])


    # Criando canvas para os histogramas (Código Novo)
    canvas_HT = ROOT.TCanvas("canvas_HT", "h1", 2200, 1300)
    canvas_HT.SetLeftMargin(0.1)
    canvas_HT.SetRightMargin(0.19)
    canvas_HT.SetBottomMargin(0.15)
    canvas_HT.SetTopMargin(0.1)
    canvas_HT.SetLogy()

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
        (hist_tt_nb1, "tt+nb", colors_dict["tt+nb"]),
        (hist_tth_x1, "tth+X", colors_dict["tth+X"]),
        (hist_tt_nV1, "tt+nV", colors_dict["tt+nV"]),
        (hist_ttt_x1, "ttt+X", colors_dict["ttt+X"])
    ]:
        hist.SetLineColor(color)
        hist.SetFillColor(color)
        hist.SetLineWidth(2)
        hist.Sumw2(0)  # Desativar incertezas
        hist.SetMinimum(0.1)  # Garante que o eixo Y comece do zero        
        hist.Draw("SAME")

    # Último histograma sem preenchimento (Código Novo)
    histograms_jets[-1].SetLineColor(colors_dict["Signal"])
    histograms_jets[-1].SetLineWidth(5)
    histograms_jets[-1].Sumw2(0)  # Desativar incertezas
    histograms_jets[-1].Draw("SAME")

    # Criando legenda (Código Novo)
    legend = ROOT.TLegend(0.81, 0.75, 0.95, 0.94)
    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    # Adicionando entradas na legenda (Código Novo)
    # Adicionando entradas na legenda (Código Novo)
    for label, hist in [
        ("tt + tt2b + tt4b", hist_tt_nb1),
        ("tth + tthZ + tthW", hist_tth_x1),
        ("ttWW + ttZZ + ttWZ+ ttZ", hist_tt_nV1),
        ("tttW + tttt", hist_ttt_x1),
        ("Signal", histograms_jets[-1])
    ]:
        legend.AddEntry(hist, label, "f")


 


    legend.Draw()


    # Atualizando e salvando o histograma (Código Novo)
    canvas_HT.Update()
    output_file_path = os.path.join(output_dir, "NJC.png")
    canvas_HT.SaveAs(output_file_path)
    subprocess.Popen(["xdg-open", output_file_path])


    # Criação dos histogramas combinados (Código Novo)
    hist_tt_nb2 = histograms_b[0].Clone("hist_tt_nb")
    hist_tth_x2 = histograms_b[3].Clone("hist_tth_x")
    hist_tt_nV2 = histograms_b[6].Clone("hist_tt_nV")
    hist_ttt_x2 = histograms_b[10].Clone("hist_ttt_x")

    # Somando histogramas conforme as categorias (Código Novo)
    for i in range(1, 3):
        hist_tt_nb2.Add(histograms_b[i])

    for i in range(4, 6):
        hist_tth_x2.Add(histograms_b[i])

    for i in range(7, 10):
        hist_tt_nV2.Add(histograms_b[i])

    for i in range(11, 12):
        hist_ttt_x2.Add(histograms_b[i])



    # Criando canvas para os histogramas (Código Novo)
    canvas_HT = ROOT.TCanvas("canvas_HT", "h1", 2200, 1300)
    canvas_HT.SetLeftMargin(0.1)
    canvas_HT.SetRightMargin(0.19)
    canvas_HT.SetBottomMargin(0.15)
    canvas_HT.SetTopMargin(0.1)
    canvas_HT.SetLogy()

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
        (hist_tt_nb2, "tt+nb", colors_dict["tt+nb"]),
        (hist_tth_x2, "tth+X", colors_dict["tth+X"]),
        (hist_tt_nV2, "tt+nV", colors_dict["tt+nV"]),
        (hist_ttt_x2, "ttt+X", colors_dict["ttt+X"])
    ]:
        hist.SetLineColor(color)
        hist.SetFillColor(color)
        hist.SetLineWidth(2)
        hist.Sumw2(0)  # Desativar incertezas
        hist.SetMinimum(0.1)  # Garante que o eixo Y comece do zero        
        hist.Draw("SAME")

    # Último histograma sem preenchimento (Código Novo)
    histograms_b[-1].SetLineColor(colors_dict["Signal"])
    histograms_b[-1].SetLineWidth(5)
    histograms_b[-1].Sumw2(0)  # Desativar incertezas
    histograms_b[-1].Draw("SAME")

    # Criando legenda (Código Novo)
    legend = ROOT.TLegend(0.81, 0.75, 0.95, 0.94)
    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    # Adicionando entradas na legenda (Código Novo)
    # Adicionando entradas na legenda (Código Novo)
    for label, hist in [
        ("tt + tt2b + tt4b", hist_tt_nb2),
        ("tth + tthZ + tthW", hist_tth_x2),
        ("ttWW + ttZZ + ttWZ+ ttZ", hist_tt_nV2),
        ("tttW + tttt", hist_ttt_x2),
        ("Signal", histograms_b[-1])
    ]:
        legend.AddEntry(hist, label, "f")


 


    legend.Draw()


    # Atualizando e salvando o histograma (Código Novo)
    canvas_HT.Update()
    output_file_path = os.path.join(output_dir, "NJBC.png")
    canvas_HT.SaveAs(output_file_path)
    subprocess.Popen(["xdg-open", output_file_path])





    # Criação dos histogramas combinados (Código Novo)
    hist_tt_nb3 = histograms_q[0].Clone("hist_tt_nb")
    hist_tth_x3 = histograms_q[3].Clone("hist_tth_x")
    hist_tt_nV3 = histograms_q[6].Clone("hist_tt_nV")
    hist_ttt_x3 = histograms_q[10].Clone("hist_ttt_x")

    # Somando histogramas conforme as categorias (Código Novo)
    for i in range(1, 3):
        hist_tt_nb3.Add(histograms_q[i])

    for i in range(4, 6):
        hist_tth_x3.Add(histograms_q[i])

    for i in range(7, 10):
        hist_tt_nV3.Add(histograms_q[i])

    for i in range(11, 12):
        hist_ttt_x3.Add(histograms_q[i])




    # Criando canvas para os histogramas (Código Novo)
    canvas_HT = ROOT.TCanvas("canvas_HT", "h1", 2200, 1300)
    canvas_HT.SetLeftMargin(0.1)
    canvas_HT.SetRightMargin(0.19)
    canvas_HT.SetBottomMargin(0.15)
    canvas_HT.SetTopMargin(0.1)
    canvas_HT.SetLogy()

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
        (hist_tt_nb3, "tt+nb", colors_dict["tt+nb"]),
        (hist_tth_x3, "tth+X", colors_dict["tth+X"]),
        (hist_tt_nV3, "tt+nV", colors_dict["tt+nV"]),
        (hist_ttt_x3, "ttt+X", colors_dict["ttt+X"])
    ]:
        hist.SetLineColor(color)
        hist.SetFillColor(color)
        hist.SetLineWidth(2)
        hist.Sumw2(0)  # Desativar incertezas
        hist.SetMinimum(0.1)  # Garante que o eixo Y comece do zero
        hist.Draw("SAME")

    # Último histograma sem preenchimento (Código Novo)
    histograms_q[-1].SetLineColor(colors_dict["Signal"])
    histograms_q[-1].SetLineWidth(5)
    histograms_q[-1].Sumw2(0)  # Desativar incertezas
    histograms_q[-1].Draw("SAME")

    # Criando legenda (Código Novo)
    legend = ROOT.TLegend(0.81, 0.75, 0.95, 0.94)
    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    # Adicionando entradas na legenda (Código Novo)
    # Adicionando entradas na legenda (Código Novo)
    for label, hist in [
        ("tt + tt2b + tt4b", hist_tt_nb3),
        ("tth + tthZ + tthW", hist_tth_x3),
        ("ttWW + ttZZ + ttWZ+ ttZ", hist_tt_nV3),
        ("tttW + tttt", hist_ttt_x3),
        ("Signal", histograms_q[-1])
    ]:
        legend.AddEntry(hist, label, "f")


 


    legend.Draw()


    # Atualizando e salvando o histograma (Código Novo)
    canvas_HT.Update()
    output_file_path = os.path.join(output_dir, "NJQC.png")
    canvas_HT.SaveAs(output_file_path)
    subprocess.Popen(["xdg-open", output_file_path])


