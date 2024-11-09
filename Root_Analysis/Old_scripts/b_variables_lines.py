#!/usr/bin/env python



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


# Check for user input
if len(sys.argv) < 8:
    print("Usage: python3 templateAnalysis.py input_file1.root input_file2.root input_file3.root ")
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
output_dir = os.path.join(script_dir, "b_variables_lines")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
##############################################################
#Definindo pesos, pro hisgrama distribuir de acordo com a seção de choque:
##############################################################
#PESO = (SEÇÃO DE CHOQUE * LUMINOSIDADE)/(NUMERO TOTAL DE EVENTOS GERADOS (10000))

##############################################################

P_1= 7559.214707 #tt
P_2= 13.19608992#*F_6  #ttbb
P_3= 2.690111995 #ttbbbb
P_4= 2.683539397 #tth
P_5= 1.155650301#*F_4 #ttz
P_6= 0.048070682#*F_5 #tttt
P_8= 0.00122845#*F_1 #ttzh
P_9= 0.000411384#*F_3 #ttzz
P_7 = 0.002389038#*F_2 #signal
                         
lista_de_pesos = [ P_1, P_2, P_4, P_5, P_6, P_7, P_8, P_9]

############################################################
######################## b[1] ##############################
############################################################
#	LOADING TREES BLOCK
#
# Lista para armazenar os histogramas de PT e ETA de cada arquivo
histograms_PT = []
histograms_ETA = []
mean_values_PT = []
mean_values_ETA = []
std_dev_values_PT = []
std_dev_values_ETA = []

# Desativar as estatísticas padrão
ROOT.gStyle.SetOptStat(0)

# Ativar a exibição do valor médio e desvio padrão personalizado
ROOT.gStyle.SetOptFit(1111)

# Cores para cada processo
#colors =  [ ROOT.kOrange+4, ROOT.kMagenta, ROOT.kBlue,  ROOT.kGreen+1, ROOT.kBlack, ROOT.kGray, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]

colors =  [ ROOT.kOrange+4, ROOT.kMagenta, ROOT.kBlue,  ROOT.kGreen+1, ROOT.kBlack, ROOT.kGray, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]
# Criando histogramas para PT e ETA dos quarks b

#histogram_PT = ROOT.TH1F("hist_pt", "P_{T} b[1]; P_{T} (GeV/c); Number of events", 40, -40, 1000)

#histogram_ETA = ROOT.TH1F("hist_pt", "#eta b[1]; #eta ; Number of events", 40, -8, 8)


# Loop sobre os arquivos .rootcolors =  [ ROOT.kOrange+4, ROOT.kMagenta, ROOT.kBlue,  ROOT.kGreen+1, ROOT.kBlack, ROOT.kGray, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]

luminosidade_desejada = 200  # Em fb^-1


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

      # Criando histogramas para PT e ETA dos quarks b ?????::::::::::::::::::::::::::::::
    histogram_PT = ROOT.TH1F("hist_pt", "P_{T} b[1]; P_{T} (GeV/c); Number of events", 40, -50, 1400)
    histogram_ETA = ROOT.TH1F("hist_pt", "#eta b[1]; #eta ; Number of events", 40, -8, 8)

     # Loop sobre os eventos
    for entry in range(numberOfEntries):
        print(f"Processing event {entry} of file {inputFile}")

        # Lendo entrada específica
        treeReader.ReadEntry(entry)
        particles = treeReader.UseBranch("Particle")

        # Inicializando listas para armazenar os PTs e ETAs dos quarks b
        pt_b_quarks = []
        eta_b_quarks = []

        # Loop sobre as partículas
        for particle in particles:
            # Verificando se a partícula é um quark b
            if particle.PID == 5 and particle.Status == 23 :  # Quark b tem PID = 5
                pt_b_quarks.append(particle.PT)
                eta_b_quarks.append(particle.Eta)



        # Verifica se há bottoms neste evento
        if pt_b_quarks:
            # Ordena os PTs dos quarks b em ordem decrescente
            pt_b_quarks.sort(reverse=True)

            # Preenche o histograma com o PT do quark b mais energético
            histogram_PT.Fill(pt_b_quarks[0], peso_atual)

            # Preenche o histograma com o ETA do quark b mais energético
            histogram_ETA.Fill(eta_b_quarks[pt_b_quarks.index(pt_b_quarks[0])], peso_atual)
	    # Imprime o PT e o ETA de cada partícula com o maior PT
            for pt, eta in zip(pt_b_quarks, eta_b_quarks):
    	        if pt == pt_b_quarks[0]:
    	        	print(f"Highest PT particle - PT: {pt}, ETA: {eta}")
            
    # Obter o número total de eventos

    total_events = histogram_ETA.GetEntries()

# Normalizar o histograma ETA
    if total_events > 0:
       num_bins = histogram_ETA.GetNbinsX()
       for i in range(1, num_bins + 1):
           bin_content = histogram_ETA.GetBinContent(i)
           histogram_ETA.SetBinContent(i, bin_content / total_events)

    # Adicionando histogramas à lista ?????::::::::::::::::::::::::::::::
    histograms_PT.append(histogram_PT)
    histograms_ETA.append(histogram_ETA)

    # Calculando média e desvio padrão para PT
    mean_pt = histogram_PT.GetMean()
    std_dev_pt = histogram_PT.GetStdDev()

    # Calculando média e desvio padrão para ETA
    mean_eta = histogram_ETA.GetMean()
    std_dev_eta = histogram_ETA.GetStdDev()

    mean_values_PT.append(mean_pt)
    mean_values_ETA.append(mean_eta)
    std_dev_values_PT.append(std_dev_pt)
    std_dev_values_ETA.append(std_dev_eta)
    
# Criando canvas para os histogramas de PT e ETA
canvas_PT = ROOT.TCanvas("canvas_PT", "PT Histogram", 1400, 800)
canvas_ETA = ROOT.TCanvas("canvas_ETA", "ETA Histogram", 1400, 800)



# Calcular o valor máximo de todos os histogramas
max_y_pt = max([histogram_PT.GetMaximum() for histogram_PT in histograms_PT])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de ETA com base no valor máximo calculado
for histogram_PT in histograms_PT:
    histogram_PT.SetMaximum(max_y_pt * margin_factor)

# Calcular o valor máximo de todos os histogramas
max_y_eta = max([histogram_ETA.GetMaximum() for histogram_ETA in histograms_ETA])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de ETA com base no valor máximo calculado
for histogram_ETA in histograms_ETA:
    histogram_ETA.SetMaximum(max_y_eta * margin_factor)



canvas_PT.cd()
for i, histogram_PT in enumerate(histograms_PT):
    histogram_PT.SetLineColor(colors[i])
    histogram_PT.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_PT.SetLineWidth(2)
    histogram_PT.Draw("SAME")
    histogram_PT.Sumw2(0)
    histogram_PT.SetMinimum(0.0001)  # Definindo escala logarítmica
    canvas_PT.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Calculando média e desvio padrão
    mean_eta = histogram_PT.GetMean()
    std_dev_eta = histogram_PT.GetStdDev()

    # Criando legenda personalizada
    legend_PT = ROOT.TLegend(0.92, 0.95 - 0.13 * i, 0.95, 0.80 - 0.13 * i)  # Posição fora da figura
    legend_PT.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_PT.SetEntrySeparation(0.07) 
    legend_PT.SetBorderSize(0)  # Removendo borda da legenda
    legend_PT.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_PT.GetListOfFunctions().Remove(legend_PT)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda


# Adicionando o nome do arquivo à legenda


    legend_PT.AddEntry(histogram_PT, f"  {os.path.basename(inputFiles[i]).replace('.root', '')}", "f100")
    legend_PT.AddEntry(0, f"Mean: {mean_eta:.2f}", "")
    legend_PT.AddEntry(0, f"Std Dev: {std_dev_eta:.2f}", "")
    legend_PT.Draw()


canvas_PT.Update()

  # Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 100  
    # Definir o limite superior para PT e ETA com base no valor máximo de y encontrado nao funcionou
histogram_PT.SetMaximum(max_y_pt * margin_factor)



# Configurando e desenhando histogramas de ETA
canvas_ETA.cd()
for i, histogram_ETA in enumerate(histograms_ETA):
    histogram_ETA.SetLineColor(colors[i])
    histogram_ETA.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_ETA.SetLineWidth(2)
    histogram_ETA.Draw("SAME")
    histogram_ETA.Sumw2(0)
    histogram_ETA.SetMinimum(0.0001)  # Definindo escala logarítmica
    canvas_ETA.SetLogy()  # Ativando escala logarítmica no eixo Y



    # Calculando média e desvio padrão
    mean_eta = histogram_ETA.GetMean()
    std_dev_eta = histogram_ETA.GetStdDev()

    # Criando legenda personalizada
    legend_ETA = ROOT.TLegend(0.92, 0.95 - 0.13 * i, 0.95, 0.80 - 0.13 * i)  # Posição fora da figura
    legend_ETA.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_ETA.SetEntrySeparation(0.07) 
    legend_ETA.SetBorderSize(0)  # Removendo borda da legenda
    legend_ETA.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_ETA.GetListOfFunctions().Remove(legend_ETA)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda
    legend_ETA.AddEntry(histogram_ETA, f"  {os.path.basename(inputFiles[i]).replace('.root', '')}", "f100")
    legend_ETA.AddEntry(0, f"Mean: {mean_eta:.2f}", "")
    legend_ETA.AddEntry(0, f"Std Dev: {std_dev_eta:.2f}", "")
    legend_ETA.Draw()


canvas_ETA.Update()




# Salvando os histogramas em arquivos de imagem
output_file_path_PT = os.path.join(output_dir, "b[1]_pt.png")
canvas_PT.SaveAs(output_file_path_PT)

output_file_path_ETA = os.path.join(output_dir, "b[1]_eta.png")
canvas_ETA.SaveAs(output_file_path_ETA)

# Aguardando até que os canvas sejam fechados pelo usuário
#canvas_PT.WaitPrimitive()
#canvas_ETA.WaitPrimitive()


############################################################
######################## b~[1] ##############################
############################################################
#	LOADING TREES BLOCK
#
# Lista para armazenar os histogramas de PT e ETA de cada arquivo
histograms_PT = []
histograms_ETA = []
mean_values_PT = []
mean_values_ETA = []
std_dev_values_PT = []
std_dev_values_ETA = []

# Desativar as estatísticas padrão
ROOT.gStyle.SetOptStat(0)

# Ativar a exibição do valor médio e desvio padrão personalizado
ROOT.gStyle.SetOptFit(1111)

# Cores para cada processo
colors =  [ ROOT.kOrange+4, ROOT.kMagenta, ROOT.kBlue,  ROOT.kGreen+1, ROOT.kBlack, ROOT.kGray, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]



# Loop sobre os arquivos .rootcolors =  [ ROOT.kOrange+4, ROOT.kMagenta, ROOT.kBlue,  ROOT.kGreen+1, ROOT.kBlack, ROOT.kGray, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]

for i, inputFile in enumerate(inputFiles):
    print(f"Reading file {inputFile}")

    # Criando cadeia de árvores ROOT
    chain = ROOT.TChain("Delphes")
    chain.Add(inputFile)

    # Criando objeto ExRootTreeReader para ler a cadeia de árvores
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

    # Criando histogramas para PT e ETA dos quarks b
    histogram_PT = ROOT.TH1F("hist_pt", "P_{T} b~[1]; P_{T} (GeV/c); Number of events", 40, -40, 800)
    histogram_ETA = ROOT.TH1F("hist_pt", "#eta b~[1]; #eta ; Number of events", 40, -8, 8)

    # Loop sobre os eventos
    for entry in range(numberOfEntries):
        print(f"Processing event {entry} of file {inputFile}")

        # Lendo entrada específica
        treeReader.ReadEntry(entry)
        particles = treeReader.UseBranch("Particle")

        # Inicializando listas para armazenar os PTs e ETAs dos quarks b
        pt_b_quarks = []
        eta_b_quarks = []

        # Loop sobre as partículas
        for particle in particles:
            # Verificando se a partícula é um quark b
            if particle.PID == -5 and particle.Status == 23:  # Quark b tem PID = 5
                pt_b_quarks.append(particle.PT)
                eta_b_quarks.append(particle.Eta)

        # Verifica se há bottoms neste evento
        if pt_b_quarks:
            # Ordena os PTs dos quarks b em ordem decrescente
            pt_b_quarks.sort(reverse=True)

            # Preenche o histograma com o PT do quark b mais energético
            histogram_PT.Fill(pt_b_quarks[0])

            # Preenche o histograma com o ETA do quark b mais energético
            histogram_ETA.Fill(eta_b_quarks[pt_b_quarks.index(pt_b_quarks[0])])
	    
            for pt, eta in zip(pt_b_quarks, eta_b_quarks):
    	        if pt == pt_b_quarks[0] and eta == eta_b_quarks[pt_b_quarks.index(pt_b_quarks[0])]:
    	            print(f"Highest PT antiparticle - PT: {pt}, ETA: {eta}")

    total_events = histogram_PT.GetEntries()

# Normalizar o histograma PT
    if total_events > 0:
       num_bins = histogram_PT.GetNbinsX()
       for i in range(1, num_bins + 1):
           bin_content = histogram_PT.GetBinContent(i)
           histogram_PT.SetBinContent(i, bin_content / total_events)        
    # Obter o número total de eventos
    total_events = histogram_PT.GetEntries()

# Normalizar o histograma ETA
    if total_events > 0:
       num_bins = histogram_ETA.GetNbinsX()
       for i in range(1, num_bins + 1):
           bin_content = histogram_ETA.GetBinContent(i)
           histogram_ETA.SetBinContent(i, bin_content / total_events)


    # Adicionando histogramas à lista
    histograms_PT.append(histogram_PT)
    histograms_ETA.append(histogram_ETA)

    # Calculando média e desvio padrão para PT
    mean_pt = histogram_PT.GetMean()
    std_dev_pt = histogram_PT.GetStdDev()

    # Calculando média e desvio padrão para ETA
    mean_eta = histogram_ETA.GetMean()
    std_dev_eta = histogram_ETA.GetStdDev()

    mean_values_PT.append(mean_pt)
    mean_values_ETA.append(mean_eta)
    std_dev_values_PT.append(std_dev_pt)
    std_dev_values_ETA.append(std_dev_eta)
    
    
# Criando canvas para os histogramas de PT e ETA


# Calcular o valor máximo de todos os histogramas
max_y_pt = max([histogram_PT.GetMaximum() for histogram_PT in histograms_PT])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de ETA com base no valor máximo calculado
for histogram_PT in histograms_PT:
    histogram_PT.SetMaximum(max_y_pt * margin_factor)

# Calcular o valor máximo de todos os histogramas
max_y_eta = max([histogram_ETA.GetMaximum() for histogram_ETA in histograms_ETA])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de ETA com base no valor máximo calculado
for histogram_ETA in histograms_ETA:
    histogram_ETA.SetMaximum(max_y_eta * margin_factor)

canvas_PT = ROOT.TCanvas("canvas_PT", "PT Histogram", 1400, 800)
canvas_ETA = ROOT.TCanvas("canvas_ETA", "ETA Histogram", 1400, 800)
# Configurando e desenhando histogramas de ETA
canvas_PT.cd()
for i, histogram_PT in enumerate(histograms_PT):
    histogram_PT.SetLineColor(colors[i])
#    #histogram_PT.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_PT.SetLineWidth(2)
    histogram_PT.Draw("SAME")
    histogram_PT.SetMinimum(0.0001)  # Definindo escala logarítmica
    canvas_PT.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Calculando média e desvio padrão
    mean_eta = histogram_PT.GetMean()
    std_dev_eta = histogram_PT.GetStdDev()

    # Criando legenda personalizada
    legend_PT = ROOT.TLegend(0.92, 0.95 - 0.13 * i, 0.95, 0.80 - 0.13 * i)  # Posição fora da figura
    legend_PT.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_PT.SetEntrySeparation(0.07) 
    legend_PT.SetBorderSize(0)  # Removendo borda da legenda
    legend_PT.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_PT.GetListOfFunctions().Remove(legend_PT)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda


# Adicionando o nome do arquivo à legenda


    legend_PT.AddEntry(histogram_PT, f"  {os.path.basename(inputFiles[i]).replace('.root', '')}", "f100")
    legend_PT.AddEntry(0, f"Mean: {mean_eta:.2f}", "")
    legend_PT.AddEntry(0, f"Std Dev: {std_dev_eta:.2f}", "")
    legend_PT.Draw()


canvas_PT.Update()
     # Encontrar o valor máximo de y para PT e ETA
max_y_pt = histogram_PT.GetMaximum()
  # Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2 
    # Definir o limite superior para PT e ETA com base no valor máximo de y encontrado
histogram_ETA.SetMaximum(max_y_eta * margin_factor)


# Configurando e desenhando histogramas de ETA
canvas_ETA.cd()
for i, histogram_ETA in enumerate(histograms_ETA):
    histogram_ETA.SetLineColor(colors[i])
#    #histogram_ETA.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_ETA.SetLineWidth(2)
    histogram_ETA.Draw("SAME")
    histogram_ETA.SetMinimum(0.0001)  # Definindo escala logarítmica
    canvas_ETA.SetLogy()  # Ativando escala logarítmica no eixo Y



    # Calculando média e desvio padrão
    mean_eta = histogram_ETA.GetMean()
    std_dev_eta = histogram_ETA.GetStdDev()

    # Criando legenda personalizada
    legend_ETA = ROOT.TLegend(0.92, 0.95 - 0.13 * i, 0.95, 0.80 - 0.13 * i)  # Posição fora da figura
    legend_ETA.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_ETA.SetEntrySeparation(0.07) 
    legend_ETA.SetBorderSize(0)  # Removendo borda da legenda
    legend_ETA.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_ETA.GetListOfFunctions().Remove(legend_ETA)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda
    legend_ETA.AddEntry(histogram_ETA, f"  {os.path.basename(inputFiles[i]).replace('.root', '')}", "f100")
    legend_ETA.AddEntry(0, f"Mean: {mean_eta:.2f}", "")
    legend_ETA.AddEntry(0, f"Std Dev: {std_dev_eta:.2f}", "")
    legend_ETA.Draw()


canvas_ETA.Update()



# Salvando os histogramas em arquivos de imagem
output_file_path_PT = os.path.join(output_dir, "b~[1]_pt.png")
canvas_PT.SaveAs(output_file_path_PT)

output_file_path_ETA = os.path.join(output_dir, "b~[1]_eta.png")
canvas_ETA.SaveAs(output_file_path_ETA)

# Aguardando até que os canvas sejam fechados pelo usuário
#canvas_PT.WaitPrimitive()
#canvas_ETA.WaitPrimitive()


############################################################
######################## b[2] ##############################
############################################################
#	LOADING TREES BLOCK
#
# Lista para armazenar os histogramas de PT e ETA de cada arquivo
histograms_PT = []
histograms_ETA = []
mean_values_PT = []
mean_values_ETA = []
std_dev_values_PT = []
std_dev_values_ETA = []
# Lista para armazenar os bottoms com segundo maior PT de cada evento
second_highest_bottom_pts = []
second_highest_bottom_etas = []
second_highest_bottom_phis = []

# Desativar as estatísticas padrão
ROOT.gStyle.SetOptStat(0)

# Ativar a exibição do valor médio e desvio padrão personalizado
ROOT.gStyle.SetOptFit(1111)

# Cores para cada processo
colors =  [ ROOT.kOrange+4, ROOT.kMagenta, ROOT.kBlue,  ROOT.kGreen+1, ROOT.kBlack, ROOT.kGray, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]



# Loop sobre os arquivos .rootcolors =  [ ROOT.kOrange+4, ROOT.kMagenta, ROOT.kBlue,  ROOT.kGreen+1, ROOT.kBlack, ROOT.kGray, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]


for i, inputFile in enumerate(inputFiles):
    print(f"Reading file {inputFile}")

    # Criando cadeia de árvores ROOT
    chain = ROOT.TChain("Delphes")
    chain.Add(inputFile)

    # Criando objeto ExRootTreeReader para ler a cadeia de árvores
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

  # Criando histogramas para PT e ETA dos quarks b
    histogram_PT = ROOT.TH1F("hist_pt", "P_{T} b[2]; P_{T} (GeV/c); Number of events", 40, -40, 1000)
    histogram_ETA = ROOT.TH1F("hist_pt", "#eta b[2]; #eta ; Number of events", 40, -10, 10)

     # Loop sobre os eventos
    for entry in range(numberOfEntries):
        print(f"Processing event {entry} of file {inputFile}")

        # Lendo entrada específica
        treeReader.ReadEntry(entry)
        particles = treeReader.UseBranch("Particle")

        # Inicializando listas para armazenar os PTs e ETAs dos quarks b
        pt_b_quarks = []
        eta_b_quarks = []

        # Loop sobre as partículas
        for particle in particles:
            # Verificando se a partícula é um quark b
            if particle.PID == 5 and particle.Status == 23 :  # Quark b tem PID = 5
                pt_b_quarks.append(particle.PT)
                eta_b_quarks.append(particle.Eta)



        # Verifica se há bottoms neste evento
        if pt_b_quarks and len(pt_b_quarks) >= 2:
            # Ordena os PTs dos quarks b em ordem decrescente
            pt_b_quarks.sort(reverse=True)

            # Preenche o histograma com o PT do segundo quark b mais energético
            histogram_PT.Fill(pt_b_quarks[1])

            # Preenche o histograma com o ETA do segundo quark b mais energético
            histogram_ETA.Fill(eta_b_quarks[pt_b_quarks.index(pt_b_quarks[1])])
	    # Imprime o PT e o ETA de cada partícula com o segundo maior PT
            for pt, eta in zip(pt_b_quarks, eta_b_quarks):
    	        if pt == pt_b_quarks[1] and eta == eta_b_quarks[pt_b_quarks.index(pt_b_quarks[1])]:
      	            print(f"Second highest PT particle - PT: {pt}, ETA: {eta}")
        

    # Adicionando histogramas à lista
    histograms_PT.append(histogram_PT)
    histograms_ETA.append(histogram_ETA)

    # Calculando média e desvio padrão para PT
    mean_pt = histogram_PT.GetMean()
    std_dev_pt = histogram_PT.GetStdDev()

    # Calculando média e desvio padrão para ETA
    mean_eta = histogram_ETA.GetMean()
    std_dev_eta = histogram_ETA.GetStdDev()

    mean_values_PT.append(mean_pt)
    mean_values_ETA.append(mean_eta)
    std_dev_values_PT.append(std_dev_pt)
    std_dev_values_ETA.append(std_dev_eta)
    # Criando canvas para os histogramas de PT e ETA
canvas_PT = ROOT.TCanvas("canvas_PT", "PT Histogram", 1400, 800)
canvas_ETA = ROOT.TCanvas("canvas_ETA", "ETA Histogram", 1400, 800)
# Configurando e desenhando histogramas de ETA



# Calcular o valor máximo de todos os histogramas
max_y_pt = max([histogram_PT.GetMaximum() for histogram_PT in histograms_PT])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de ETA com base no valor máximo calculado
for histogram_PT in histograms_PT:
    histogram_PT.SetMaximum(max_y_pt * margin_factor)

# Calcular o valor máximo de todos os histogramas
max_y_eta = max([histogram_ETA.GetMaximum() for histogram_ETA in histograms_ETA])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de ETA com base no valor máximo calculado
for histogram_ETA in histograms_ETA:
    histogram_ETA.SetMaximum(max_y_eta * margin_factor)
    
    
canvas_PT.cd()
for i, histogram_PT in enumerate(histograms_PT):
    histogram_PT.SetLineColor(colors[i])
    #histogram_PT.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_PT.SetLineWidth(2)
    histogram_PT.Draw("SAME")
    histogram_PT.SetMinimum(0.0001)  # Definindo escala logarítmica
    canvas_PT.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Calculando média e desvio padrão
    mean_eta = histogram_PT.GetMean()
    std_dev_eta = histogram_PT.GetStdDev()

    # Criando legenda personalizada
    legend_PT = ROOT.TLegend(0.92, 0.95 - 0.13 * i, 0.95, 0.80 - 0.13 * i)  # Posição fora da figura
    legend_PT.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_PT.SetEntrySeparation(0.07) 
    legend_PT.SetBorderSize(0)  # Removendo borda da legenda
    legend_PT.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_PT.GetListOfFunctions().Remove(legend_PT)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda


# Adicionando o nome do arquivo à legenda


    legend_PT.AddEntry(histogram_PT, f"  {os.path.basename(inputFiles[i]).replace('.root', '')}", "f100")
    legend_PT.AddEntry(0, f"Mean: {mean_eta:.2f}", "")
    legend_PT.AddEntry(0, f"Std Dev: {std_dev_eta:.2f}", "")
    legend_PT.Draw()

canvas_PT.Update()


# Configurando e desenhando histogramas de ETA
canvas_ETA.cd()
for i, histogram_ETA in enumerate(histograms_ETA):
    histogram_ETA.SetLineColor(colors[i])
    #histogram_ETA.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_ETA.SetLineWidth(2)
    histogram_ETA.Draw("SAME")
    histogram_ETA.SetMinimum(0.1)  # Definindo escala logarítmica
    canvas_ETA.SetLogy()  # Ativando escala logarítmica no eixo Y



    # Calculando média e desvio padrão
    mean_eta = histogram_ETA.GetMean()
    std_dev_eta = histogram_ETA.GetStdDev()

    # Criando legenda personalizada
    legend_ETA = ROOT.TLegend(0.92, 0.95 - 0.13 * i, 0.95, 0.80 - 0.13 * i)  # Posição fora da figura
    legend_ETA.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_ETA.SetEntrySeparation(0.07) 
    legend_ETA.SetBorderSize(0)  # Removendo borda da legenda
    legend_ETA.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_ETA.GetListOfFunctions().Remove(legend_ETA)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda
    legend_ETA.AddEntry(histogram_ETA, f"  {os.path.basename(inputFiles[i]).replace('.root', '')}", "f100")
    legend_ETA.AddEntry(0, f"Mean: {mean_eta:.2f}", "")
    legend_ETA.AddEntry(0, f"Std Dev: {std_dev_eta:.2f}", "")
    legend_ETA.Draw()

canvas_ETA.Update()



# Salvando os histogramas em arquivos de imagem
output_file_path_PT = os.path.join(output_dir, "b[2]_pt.png")
canvas_PT.SaveAs(output_file_path_PT)

output_file_path_ETA = os.path.join(output_dir, "b[2]_eta.png")
canvas_ETA.SaveAs(output_file_path_ETA)

# Aguardando até que os canvas sejam fechados pelo usuário
#canvas_PT.WaitPrimitive()
#canvas_ETA.WaitPrimitive()


############################################################
######################## b~[2] ##############################
############################################################
#	LOADING TREES BLOCK
#
# Lista para armazenar os histogramas de PT e ETA de cada arquivo
histograms_PT = []
histograms_ETA = []
mean_values_PT = []
mean_values_ETA = []
std_dev_values_PT = []
std_dev_values_ETA = []
# Lista para armazenar os bottoms com segundo maior PT de cada evento
second_highest_bottom_pts = []
second_highest_bottom_etas = []
second_highest_bottom_phis = []

# Desativar as estatísticas padrão
ROOT.gStyle.SetOptStat(0)

# Ativar a exibição do valor médio e desvio padrão personalizado
ROOT.gStyle.SetOptFit(1111)

# Cores para cada processo
colors =  [ ROOT.kOrange+4, ROOT.kMagenta, ROOT.kBlue,  ROOT.kGreen+1, ROOT.kBlack, ROOT.kGray, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]



# Loop sobre os arquivos .rootcolors =  [ ROOT.kOrange+4, ROOT.kMagenta, ROOT.kBlue,  ROOT.kGreen+1, ROOT.kBlack, ROOT.kGray, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]


for i, inputFile in enumerate(inputFiles):
    print(f"Reading file {inputFile}")

    # Criando cadeia de árvores ROOT
    chain = ROOT.TChain("Delphes")
    chain.Add(inputFile)

    # Criando objeto ExRootTreeReader para ler a cadeia de árvores
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

  # Criando histogramas para PT e ETA dos quarks b
    histogram_PT = ROOT.TH1F("hist_pt", "P_{T} b~[2]; P_{T} (GeV/c); Number of events", 40, -40, 1000)
    histogram_ETA = ROOT.TH1F("hist_pt", "#eta b~[2]; #eta ; Number of events", 40, -10, 10)

     # Loop sobre os eventos
    for entry in range(numberOfEntries):
        print(f"Processing event {entry} of file {inputFile}")

        # Lendo entrada específica
        treeReader.ReadEntry(entry)
        particles = treeReader.UseBranch("Particle")

        # Inicializando listas para armazenar os PTs e ETAs dos quarks b
        pt_b_quarks = []
        eta_b_quarks = []

        # Loop sobre as partículas
        for particle in particles:
            # Verificando se a partícula é um quark b
            if particle.PID == -5 and particle.Status == 23 :  # Quark b tem PID = 5
                pt_b_quarks.append(particle.PT)
                eta_b_quarks.append(particle.Eta)



        # Verifica se há bottoms neste evento
        if pt_b_quarks and len(pt_b_quarks) >= 2:
            # Ordena os PTs dos quarks b em ordem decrescente
            pt_b_quarks.sort(reverse=True)

            # Preenche o histograma com o PT do segundo quark b mais energético
            histogram_PT.Fill(pt_b_quarks[1])

            # Preenche o histograma com o ETA do segundo quark b mais energético
            histogram_ETA.Fill(eta_b_quarks[pt_b_quarks.index(pt_b_quarks[1])])
	    # Imprime o PT e o ETA de cada partícula com o segundo maior PT
            for pt, eta in zip(pt_b_quarks, eta_b_quarks):
    	        if pt == pt_b_quarks[1] and eta == eta_b_quarks[pt_b_quarks.index(pt_b_quarks[1])]:
      	            print(f"Second highest PT antiparticle - PT: {pt}, ETA: {eta}")
        

    # Adicionando histogramas à lista
    histograms_PT.append(histogram_PT)
    histograms_ETA.append(histogram_ETA)

    # Calculando média e desvio padrão para PT
    mean_pt = histogram_PT.GetMean()
    std_dev_pt = histogram_PT.GetStdDev()

    # Calculando média e desvio padrão para ETA
    mean_eta = histogram_ETA.GetMean()
    std_dev_eta = histogram_ETA.GetStdDev()

    mean_values_PT.append(mean_pt)
    mean_values_ETA.append(mean_eta)
    std_dev_values_PT.append(std_dev_pt)
    std_dev_values_ETA.append(std_dev_eta)
    # Criando canvas para os histogramas de PT e ETA
canvas_PT = ROOT.TCanvas("canvas_PT", "PT Histogram", 1400, 800)
canvas_ETA = ROOT.TCanvas("canvas_ETA", "ETA Histogram", 1400, 800)


# Calcular o valor máximo de todos os histogramas
max_y_pt = max([histogram_PT.GetMaximum() for histogram_PT in histograms_PT])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de ETA com base no valor máximo calculado
for histogram_PT in histograms_PT:
    histogram_PT.SetMaximum(max_y_pt * margin_factor)

# Calcular o valor máximo de todos os histogramas
max_y_eta = max([histogram_ETA.GetMaximum() for histogram_ETA in histograms_ETA])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de ETA com base no valor máximo calculado
for histogram_ETA in histograms_ETA:
    histogram_ETA.SetMaximum(max_y_eta * margin_factor)




canvas_PT.cd()
for i, histogram_PT in enumerate(histograms_PT):
    histogram_PT.SetLineColor(colors[i])
    #histogram_PT.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_PT.SetLineWidth(2)
    histogram_PT.Draw("SAME")
    histogram_PT.SetMinimum(0.001)  # Definindo escala logarítmica
    canvas_PT.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Calculando média e desvio padrão
    mean_eta = histogram_PT.GetMean()
    std_dev_eta = histogram_PT.GetStdDev()

    # Criando legenda personalizada
    legend_PT = ROOT.TLegend(0.92, 0.95 - 0.13 * i, 0.95, 0.80 - 0.13 * i)  # Posição fora da figura
    legend_PT.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_PT.SetEntrySeparation(0.07) 
    legend_PT.SetBorderSize(0)  # Removendo borda da legenda
    legend_PT.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_PT.GetListOfFunctions().Remove(legend_PT)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda


# Adicionando o nome do arquivo à legenda


    legend_PT.AddEntry(histogram_PT, f"  {os.path.basename(inputFiles[i]).replace('.root', '')}", "f100")
    legend_PT.AddEntry(0, f"Mean: {mean_eta:.2f}", "")
    legend_PT.AddEntry(0, f"Std Dev: {std_dev_eta:.2f}", "")
    legend_PT.Draw()


canvas_PT.Update()



# Configurando e desenhando histogramas de ETA
canvas_ETA.cd()
for i, histogram_ETA in enumerate(histograms_ETA):
    histogram_ETA.SetLineColor(colors[i])
    #histogram_ETA.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_ETA.SetLineWidth(2)
    histogram_ETA.Draw("SAME")
    histogram_ETA.SetMinimum(0.1)  # Definindo escala logarítmica
    canvas_ETA.SetLogy()  # Ativando escala logarítmica no eixo Y



    # Calculando média e desvio padrão
    mean_eta = histogram_ETA.GetMean()
    std_dev_eta = histogram_ETA.GetStdDev()

    # Criando legenda personalizada
    legend_ETA = ROOT.TLegend(0.92, 0.95 - 0.13 * i, 0.95, 0.80 - 0.13 * i)  # Posição fora da figura
    legend_ETA.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_ETA.SetEntrySeparation(0.07) 
    legend_ETA.SetBorderSize(0)  # Removendo borda da legenda
    legend_ETA.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_ETA.GetListOfFunctions().Remove(legend_ETA)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda
    legend_ETA.AddEntry(histogram_ETA, f"  {os.path.basename(inputFiles[i]).replace('.root', '')}", "f100")
    legend_ETA.AddEntry(0, f"Mean: {mean_eta:.2f}", "")
    legend_ETA.AddEntry(0, f"Std Dev: {std_dev_eta:.2f}", "")
    legend_ETA.Draw()

canvas_ETA.Update()

# Salvando os histogramas em arquivos de imagem
output_file_path_PT = os.path.join(output_dir, "b~[2]_pt.png")
canvas_PT.SaveAs(output_file_path_PT)

output_file_path_ETA = os.path.join(output_dir, "b~[2]_eta.png")
canvas_ETA.SaveAs(output_file_path_ETA)

# Aguardando até que os canvas sejam fechados pelo usuário
#canvas_PT.WaitPrimitive()
#canvas_ETA.WaitPrimitive()


############################################################
######################## b[3] ##############################
############################################################
#	LOADING TREES BLOCK
#
# Lista para armazenar os histogramas de PT e ETA de cada arquivo
histograms_PT = []
histograms_ETA = []
mean_values_PT = []
mean_values_ETA = []
std_dev_values_PT = []
std_dev_values_ETA = []

# Desativar as estatísticas padrão
ROOT.gStyle.SetOptStat(0)

# Ativar a exibição do valor médio e desvio padrão personalizado
ROOT.gStyle.SetOptFit(1111)

# Cores para cada processo
colors =  [ ROOT.kOrange+4, ROOT.kMagenta, ROOT.kBlue,  ROOT.kGreen+1, ROOT.kBlack, ROOT.kGray, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]



# Loop sobre os arquivos .rootcolors =  [ ROOT.kOrange+4, ROOT.kMagenta, ROOT.kBlue,  ROOT.kGreen+1, ROOT.kBlack, ROOT.kGray, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]


for i, inputFile in enumerate(inputFiles):
    print(f"Reading file {inputFile}")

    # Criando cadeia de árvores ROOT
    chain = ROOT.TChain("Delphes")
    chain.Add(inputFile)

    # Criando objeto ExRootTreeReader para ler a cadeia de árvores
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

  # Criando histogramas para PT e ETA dos quarks b
    histogram_PT = ROOT.TH1F("hist_pt", "P_{T} b[3]; P_{T} (GeV/c); Number of events", 40, -40, 1000)

    histogram_ETA = ROOT.TH1F("hist_pt", "#eta b[3]; #eta ; Number of events", 40, -10, 10)

     # Loop sobre os eventos
    for entry in range(numberOfEntries):
        print(f"Processing event {entry} of file {inputFile}")

        # Lendo entrada específica
        treeReader.ReadEntry(entry)
        particles = treeReader.UseBranch("Particle")

        # Inicializando listas para armazenar os PTs e ETAs dos quarks b
        pt_b_quarks = []
        eta_b_quarks = []

        # Loop sobre as partículas
        for particle in particles:
            # Verificando se a partícula é um quark b
            if particle.PID == 5 and particle.Status == 23 :  # Quark b tem PID = 5
                pt_b_quarks.append(particle.PT)
                eta_b_quarks.append(particle.Eta)



        # Verifica se há bottoms neste evento
        if pt_b_quarks and len(pt_b_quarks) >= 3:
            # Ordena os PTs dos quarks b em ordem decrescente
            pt_b_quarks.sort(reverse=True)

            # Preenche o histograma com o PT do quark b menos energético
            histogram_PT.Fill(pt_b_quarks[2])

            # Preenche o histograma com o ETA do quark b menos energético
            histogram_ETA.Fill(eta_b_quarks[pt_b_quarks.index(pt_b_quarks[2])])
	    # Imprime o PT e o ETA de cada partícula com o maior PT
            for pt, eta in zip(pt_b_quarks, eta_b_quarks):
    	        if pt == pt_b_quarks[2] and eta == eta_b_quarks[pt_b_quarks.index(pt_b_quarks[2])]:
      	            print(f"Third highest PT particle - PT: {pt}, ETA: {eta}")
            
    
    # Adicionando histogramas à lista
    histograms_PT.append(histogram_PT)
    histograms_ETA.append(histogram_ETA)

    # Calculando média e desvio padrão para PT
    mean_pt = histogram_PT.GetMean()
    std_dev_pt = histogram_PT.GetStdDev()

    # Calculando média e desvio padrão para ETA
    mean_eta = histogram_ETA.GetMean()
    std_dev_eta = histogram_ETA.GetStdDev()

    mean_values_PT.append(mean_pt)
    mean_values_ETA.append(mean_eta)
    std_dev_values_PT.append(std_dev_pt)
    std_dev_values_ETA.append(std_dev_eta)
    
        # Criando canvas para os histogramas de PT e ETA
canvas_PT = ROOT.TCanvas("canvas_PT", "PT Histogram", 1400, 800)
canvas_ETA = ROOT.TCanvas("canvas_ETA", "ETA Histogram", 1400, 800)



# Calcular o valor máximo de todos os histogramas
max_y_pt = max([histogram_PT.GetMaximum() for histogram_PT in histograms_PT])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de ETA com base no valor máximo calculado
for histogram_PT in histograms_PT:
    histogram_PT.SetMaximum(max_y_pt * margin_factor)

# Calcular o valor máximo de todos os histogramas
max_y_eta = max([histogram_ETA.GetMaximum() for histogram_ETA in histograms_ETA])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de ETA com base no valor máximo calculado
for histogram_ETA in histograms_ETA:
    histogram_ETA.SetMaximum(max_y_eta * margin_factor)


canvas_PT.cd()
for i, histogram_PT in enumerate(histograms_PT):
    histogram_PT.SetLineColor(colors[i])
    #histogram_PT.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_PT.SetLineWidth(2)
    histogram_PT.Draw("SAME")
    histogram_PT.SetMinimum(0.001)  # Definindo escala logarítmica

    canvas_PT.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Calculando média e desvio padrão
    mean_eta = histogram_PT.GetMean()
    std_dev_eta = histogram_PT.GetStdDev()

    # Criando legenda personalizada
    legend_PT = ROOT.TLegend(0.92, 0.95 - 0.13 * i, 0.95, 0.80 - 0.13 * i)  # Posição fora da figura
    legend_PT.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_PT.SetEntrySeparation(0.07) 
    legend_PT.SetBorderSize(0)  # Removendo borda da legenda
    legend_PT.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_PT.GetListOfFunctions().Remove(legend_PT)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda


# Adicionando o nome do arquivo à legenda


    legend_PT.AddEntry(histogram_PT, f"  {os.path.basename(inputFiles[i]).replace('.root', '')}", "f100")
    legend_PT.AddEntry(0, f"Mean: {mean_eta:.2f}", "")
    legend_PT.AddEntry(0, f"Std Dev: {std_dev_eta:.2f}", "")
    legend_PT.Draw()


canvas_PT.Update()



# Configurando e desenhando histogramas de ETA
canvas_ETA.cd()
for i, histogram_ETA in enumerate(histograms_ETA):
    histogram_ETA.SetLineColor(colors[i])
    #histogram_ETA.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_ETA.SetLineWidth(2)
    histogram_ETA.Draw("SAME")
    histogram_ETA.SetMinimum(0.0001)  # Definindo escala logarítmica
    canvas_ETA.SetLogy()  # Ativando escala logarítmica no eixo Y



    # Calculando média e desvio padrão
    mean_eta = histogram_ETA.GetMean()
    std_dev_eta = histogram_ETA.GetStdDev()

    # Criando legenda personalizada
    legend_ETA = ROOT.TLegend(0.92, 0.95 - 0.13 * i, 0.95, 0.80 - 0.13 * i)  # Posição fora da figura
    legend_ETA.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_ETA.SetEntrySeparation(0.07) 
    legend_ETA.SetBorderSize(0)  # Removendo borda da legenda
    legend_ETA.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_ETA.GetListOfFunctions().Remove(legend_ETA)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda
    legend_ETA.AddEntry(histogram_ETA, f"  {os.path.basename(inputFiles[i]).replace('.root', '')}", "f100")
    legend_ETA.AddEntry(0, f"Mean: {mean_eta:.2f}", "")
    legend_ETA.AddEntry(0, f"Std Dev: {std_dev_eta:.2f}", "")
    legend_ETA.Draw()

canvas_ETA.Update()# Salvando os histogramas em arquivos de imagem
output_file_path_PT = os.path.join(output_dir, "b[3]_pt.png")
canvas_PT.SaveAs(output_file_path_PT)

output_file_path_ETA = os.path.join(output_dir, "b[3]_eta.png")
canvas_ETA.SaveAs(output_file_path_ETA)

# Aguardando até que os canvas sejam fechados pelo usuário
#canvas_PT.WaitPrimitive()
#canvas_ETA.WaitPrimitive()



############################################################
######################## b~[3] ##############################
############################################################
#	LOADING TREES BLOCK
#
# Lista para armazenar os histogramas de PT e ETA de cada arquivo
histograms_PT = []
histograms_ETA = []
mean_values_PT = []
mean_values_ETA = []
std_dev_values_PT = []
std_dev_values_ETA = []

# Desativar as estatísticas padrão
ROOT.gStyle.SetOptStat(0)

# Ativar a exibição do valor médio e desvio padrão personalizado
ROOT.gStyle.SetOptFit(1111)

# Cores para cada processo
colors =  [ ROOT.kOrange+4, ROOT.kMagenta, ROOT.kBlue,  ROOT.kGreen+1, ROOT.kBlack, ROOT.kGray, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]



# Loop sobre os arquivos .rootcolors =  [ ROOT.kOrange+4, ROOT.kMagenta, ROOT.kBlue,  ROOT.kGreen+1, ROOT.kBlack, ROOT.kGray, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]


for i, inputFile in enumerate(inputFiles):
    print(f"Reading file {inputFile}")

    # Criando cadeia de árvores ROOT
    chain = ROOT.TChain("Delphes")
    chain.Add(inputFile)

    # Criando objeto ExRootTreeReader para ler a cadeia de árvores
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

  # Criando histogramas para PT e ETA dos quarks b
    histogram_PT = ROOT.TH1F("hist_pt", "P_{T} b~[3]; P_{T} (GeV/c); Number of events", 40, -40, 1000)
    histogram_ETA = ROOT.TH1F("hist_pt", "#eta b~[3]; #eta ; Number of events", 40, -10, 10)


     # Loop sobre os eventos
    for entry in range(numberOfEntries):
        print(f"Processing event {entry} of file {inputFile}")

        # Lendo entrada específica
        treeReader.ReadEntry(entry)
        particles = treeReader.UseBranch("Particle")

        # Inicializando listas para armazenar os PTs e ETAs dos quarks b
        pt_b_quarks = []
        eta_b_quarks = []

        # Loop sobre as partículas
        for particle in particles:
            # Verificando se a partícula é um quark b
            if particle.PID == -5 and particle.Status == 23 :  # Quark b tem PID = 5
                pt_b_quarks.append(particle.PT)
                eta_b_quarks.append(particle.Eta)



        # Verifica se há bottoms neste evento
        if pt_b_quarks and len(pt_b_quarks) >= 3:
            # Ordena os PTs dos quarks b em ordem decrescente
            pt_b_quarks.sort(reverse=True)

            # Preenche o histograma com o PT do quark b menos energético
            histogram_PT.Fill(pt_b_quarks[2])

            # Preenche o histograma com o ETA do quark b menos energético
            histogram_ETA.Fill(eta_b_quarks[pt_b_quarks.index(pt_b_quarks[2])])
	    # Imprime o PT e o ETA de cada partícula com o menor PT
            for pt, eta in zip(pt_b_quarks, eta_b_quarks):
    	        if pt == pt_b_quarks[2] and eta == eta_b_quarks[pt_b_quarks.index(pt_b_quarks[2])] :
    	        	print(f"Third highest PT particle - PT: {pt}, ETA: {eta}")
        



    # Adicionando histogramas à lista
    histograms_PT.append(histogram_PT)
    histograms_ETA.append(histogram_ETA)

    # Calculando média e desvio padrão para PT
    mean_pt = histogram_PT.GetMean()
    std_dev_pt = histogram_PT.GetStdDev()

    # Calculando média e desvio padrão para ETA
    mean_eta = histogram_ETA.GetMean()
    std_dev_eta = histogram_ETA.GetStdDev()

    mean_values_PT.append(mean_pt)
    mean_values_ETA.append(mean_eta)
    std_dev_values_PT.append(std_dev_pt)
    std_dev_values_ETA.append(std_dev_eta)
    
    # Criando canvas para os histogramas de PT e ETA
canvas_PT = ROOT.TCanvas("canvas_PT", "PT Histogram", 1400, 800)
canvas_ETA = ROOT.TCanvas("canvas_ETA", "ETA Histogram", 1400, 800)


# Calcular o valor máximo de todos os histogramas
max_y_pt = max([histogram_PT.GetMaximum() for histogram_PT in histograms_PT])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de ETA com base no valor máximo calculado
for histogram_PT in histograms_PT:
    histogram_PT.SetMaximum(max_y_pt * margin_factor)

# Calcular o valor máximo de todos os histogramas
max_y_eta = max([histogram_ETA.GetMaximum() for histogram_ETA in histograms_ETA])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de ETA com base no valor máximo calculado
for histogram_ETA in histograms_ETA:
    histogram_ETA.SetMaximum(max_y_eta * margin_factor)



canvas_PT.cd()
for i, histogram_PT in enumerate(histograms_PT):
    histogram_PT.SetLineColor(colors[i])
    #histogram_PT.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_PT.SetLineWidth(2)
    histogram_PT.Draw("SAME")
    histogram_PT.SetMinimum(0.001)  # Definindo escala logarítmica
    canvas_PT.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Calculando média e desvio padrão
    mean_eta = histogram_PT.GetMean()
    std_dev_eta = histogram_PT.GetStdDev()

    # Criando legenda personalizada
    legend_PT = ROOT.TLegend(0.92, 0.95 - 0.13 * i, 0.95, 0.80 - 0.13 * i)  # Posição fora da figura
    legend_PT.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_PT.SetEntrySeparation(0.07) 
    legend_PT.SetBorderSize(0)  # Removendo borda da legenda
    legend_PT.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_PT.GetListOfFunctions().Remove(legend_PT)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda


# Adicionando o nome do arquivo à legenda


    legend_PT.AddEntry(histogram_PT, f"  {os.path.basename(inputFiles[i]).replace('.root', '')}", "f100")
    legend_PT.AddEntry(0, f"Mean: {mean_eta:.2f}", "")
    legend_PT.AddEntry(0, f"Std Dev: {std_dev_eta:.2f}", "")
    legend_PT.Draw()


canvas_PT.Update()



# Configurando e desenhando histogramas de ETA
canvas_ETA.cd()
for i, histogram_ETA in enumerate(histograms_ETA):
    histogram_ETA.SetLineColor(colors[i])
    #histogram_ETA.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_ETA.SetLineWidth(2)
    histogram_ETA.Draw("SAME")
    histogram_ETA.SetMinimum(0.0001)  # Definindo escala logarítmica
    canvas_ETA.SetLogy()  # Ativando escala logarítmica no eixo Y



    # Calculando média e desvio padrão
    mean_eta = histogram_ETA.GetMean()
    std_dev_eta = histogram_ETA.GetStdDev()

    # Criando legenda personalizada
    legend_ETA = ROOT.TLegend(0.92, 0.95 - 0.13 * i, 0.95, 0.80 - 0.13 * i)  # Posição fora da figura
    legend_ETA.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_ETA.SetEntrySeparation(0.07) 
    legend_ETA.SetBorderSize(0)  # Removendo borda da legenda
    legend_ETA.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_ETA.GetListOfFunctions().Remove(legend_ETA)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda
    legend_ETA.AddEntry(histogram_ETA, f"  {os.path.basename(inputFiles[i]).replace('.root', '')}", "f100")
    legend_ETA.AddEntry(0, f"Mean: {mean_eta:.2f}", "")
    legend_ETA.AddEntry(0, f"Std Dev: {std_dev_eta:.2f}", "")
    legend_ETA.Draw()

canvas_ETA.Update()# Salvando os histogramas em arquivos de imagem
output_file_path_PT = os.path.join(output_dir, "b~[3]_pt.png")
canvas_PT.SaveAs(output_file_path_PT)

output_file_path_ETA = os.path.join(output_dir, "b~[3]_eta.png")
canvas_ETA.SaveAs(output_file_path_ETA)

# Aguardando até que os canvas sejam fechados pelo usuário
canvas_PT.WaitPrimitive()
canvas_ETA.WaitPrimitive()


