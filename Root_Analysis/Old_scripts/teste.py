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
if len(sys.argv) < 3:
    print("Usage: python3 templateAnalysis.py input_file1.root input_file2.root")
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
output_dir = os.path.join(script_dir, "b_variables")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

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
colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen, ROOT.kCoral, ROOT.kMagenta, ROOT.kCyan, ROOT.kBlack]


# Loop sobre os arquivos .rootcolors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen, ROOT.kCoral, ROOT.kMagenta, ROOT.kCyan, ROOT.kBlack]

for i, inputFile in enumerate(inputFiles):
    print(f"Reading file {inputFile}")

    # Criando cadeia de árvores ROOT
    chain = ROOT.TChain("Delphes")
    chain.Add(inputFile)

    # Criando objeto ExRootTreeReader para ler a cadeia de árvores
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

  # Criando histogramas para PT e ETA dos quarks b
    histogram_PT = ROOT.TH1F("hist_pt", "P_{T} b[1]; P_{T} (GeV/c); Number of events", 40, 0, 600)

    histogram_ETA = ROOT.TH1F("hist_pt", "#eta b[1]; #eta ; Number of events", 40, -10, 10)

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
            if particle.PID == 5 and particle.Status == 23 :    # Quark b tem PID = 5
                pt_b_quarks.append(particle.PT)
                eta_b_quarks.append(particle.Eta)

        # Selecionando o quark b com maior PT e preenchendo os histogramas
        if pt_b_quarks:
            max_pt_b = max(pt_b_quarks)
            histogram_PT.Fill(max_pt_b)

        if eta_b_quarks:
            max_eta_b = max(eta_b_quarks)
            print(f"PT do quark mais energético : {max_pt_b}, ETA do quark mais energético {max_eta_b}.")
            histogram_ETA.Fill(max_eta_b)
	histogram_ETA.Scale(1.0 / histogram_ETA.Integral())

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
canvas_PT = ROOT.TCanvas("canvas_PT", "PT Histogram", 800, 600)
canvas_ETA = ROOT.TCanvas("canvas_ETA", "ETA Histogram", 800, 600)

# Configurando e desenhando histogramas de PT
canvas_PT.cd()
for i, histogram_PT in enumerate(histograms_PT):
    histogram_PT.SetLineColor(colors[i])
    histogram_PT.SetLineWidth(2)
    histogram_PT.Draw("SAME")
    histogram_PT.SetMinimum(0.1)  # Definindo escala logarítmica
    canvas_PT.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Calculando média e desvio padrão
    mean_pt = histogram_PT.GetMean()
    std_dev_pt = histogram_PT.GetStdDev()

    # Criando legenda personalizada
    legend_PT = ROOT.TLegend(0.9, 0.9 - 0.1 * i, 0.95, 0.65 - 0.1 * i)  # Posição fora da figura
    legend_PT.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_PT.SetBorderSize(0)  # Removendo borda da legenda
    legend_PT.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_PT.GetListOfFunctions().Remove(legend_PT)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda
    legend_PT.AddEntry(histogram_PT, f"{os.path.basename(inputFiles[i]).replace('.root', '')}", "l")
    legend_PT.AddEntry(0, f"Mean: {mean_pt:.2f}", "")
    legend_PT.AddEntry(0, f"Std Dev: {std_dev_pt:.2f}", "")
    legend_PT.Draw()

canvas_PT.Update()


# Configurando e desenhando histogramas de ETA
canvas_ETA.cd()
for i, histogram_ETA in enumerate(histograms_ETA):
    histogram_ETA.SetLineColor(colors[i])
    histogram_ETA.SetLineWidth(2)
    histogram_ETA.Draw("SAME")
    histogram_ETA.SetMinimum(0.1)  # Definindo escala logarítmica
    canvas_ETA.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Calculando média e desvio padrão
    mean_eta = histogram_ETA.GetMean()
    std_dev_eta = histogram_ETA.GetStdDev()

    # Criando legenda personalizada
    legend_ETA = ROOT.TLegend(0.9, 0.9 - 0.1 * i, 0.95, 0.65 - 0.1 * i)  # Posição fora da figura
    legend_ETA.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_ETA.SetBorderSize(0)  # Removendo borda da legenda
    legend_ETA.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_ETA.GetListOfFunctions().Remove(legend_ETA)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda
    legend_ETA.AddEntry(histogram_ETA, f"{os.path.basename(inputFiles[i]).replace('.root', '')}", "l")
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
canvas_PT.WaitPrimitive()
canvas_ETA.WaitPrimitive()


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
colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen, ROOT.kCoral, ROOT.kMagenta, ROOT.kCyan, ROOT.kBlack]


# Loop sobre os arquivos .rootcolors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen, ROOT.kCoral, ROOT.kMagenta, ROOT.kCyan, ROOT.kBlack]

for i, inputFile in enumerate(inputFiles):
    print(f"Reading file {inputFile}")

    # Criando cadeia de árvores ROOT
    chain = ROOT.TChain("Delphes")
    chain.Add(inputFile)

    # Criando objeto ExRootTreeReader para ler a cadeia de árvores
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

  # Criando histogramas para PT e ETA dos quarks b
    histogram_PT = ROOT.TH1F("hist_pt", "P_{T} b~[1]; P_{T} (GeV/c); Number of events", 40, 0, 1000)
    histogram_ETA = ROOT.TH1F("hist_pt", "#eta b~[1]; #eta ; Number of events", 40, -10, 10)

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
            if particle.PID == -5 and particle.Status == 23 :  # Quark antib tem PID = -5
                pt_b_quarks.append(particle.PT)
                eta_b_quarks.append(particle.Eta)

        # Selecionando o quark b com maior PT e preenchendo os histogramas
        if pt_b_quarks:
            max_pt_b = max(pt_b_quarks)
            histogram_PT.Fill(max_pt_b)
    histogram_PT.Scale(1.0 / histogram_PT.Integral())    	
        if eta_b_quarks:
            max_eta_b = max(eta_b_quarks)
            print(f"PT do antiquark b mais energético : {max_pt_b}, ETA do antiquark b mais energético {max_eta_b}.")
            histogram_ETA.Fill(max_eta_b)
    histogram_ETA.Scale(1.0 / histogram_ETA.Integral())
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
canvas_PT = ROOT.TCanvas("canvas_PT", "PT Histogram", 800, 600)
canvas_ETA = ROOT.TCanvas("canvas_ETA", "ETA Histogram", 800, 600)

# Configurando e desenhando histogramas de PT
canvas_PT.cd()
for i, histogram_PT in enumerate(histograms_PT):
    histogram_PT.SetLineColor(colors[i])
    histogram_PT.SetLineWidth(2)
    histogram_PT.Draw("SAME")
    histogram_PT.SetMinimum(0.1)  # Definindo escala logarítmica
    canvas_PT.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Calculando média e desvio padrão
    mean_pt = histogram_PT.GetMean()
    std_dev_pt = histogram_PT.GetStdDev()

    # Criando legenda personalizada
    legend_PT = ROOT.TLegend(0.9, 0.9 - 0.1 * i, 0.95, 0.65 - 0.1 * i)  # Posição fora da figura
    legend_PT.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_PT.SetBorderSize(0)  # Removendo borda da legenda
    legend_PT.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_PT.GetListOfFunctions().Remove(legend_PT)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda
    legend_PT.AddEntry(histogram_PT, f"{os.path.basename(inputFiles[i]).replace('.root', '')}", "l")
    legend_PT.AddEntry(0, f"Mean: {mean_pt:.2f}", "")
    legend_PT.AddEntry(0, f"Std Dev: {std_dev_pt:.2f}", "")
    legend_PT.Draw()

canvas_PT.Update()


# Configurando e desenhando histogramas de ETA
canvas_ETA.cd()
for i, histogram_ETA in enumerate(histograms_ETA):
    histogram_ETA.SetLineColor(colors[i])
    histogram_ETA.SetLineWidth(2)
    histogram_ETA.Draw("SAME")
    histogram_ETA.SetMinimum(0.1)  # Definindo escala logarítmica
    canvas_ETA.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Calculando média e desvio padrão
    mean_eta = histogram_ETA.GetMean()
    std_dev_eta = histogram_ETA.GetStdDev()

    # Criando legenda personalizada
    legend_ETA = ROOT.TLegend(0.9, 0.9 - 0.1 * i, 0.95, 0.65 - 0.1 * i)  # Posição fora da figura
    legend_ETA.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_ETA.SetBorderSize(0)  # Removendo borda da legenda
    legend_ETA.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_ETA.GetListOfFunctions().Remove(legend_ETA)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda
    legend_ETA.AddEntry(histogram_ETA, f"{os.path.basename(inputFiles[i]).replace('.root', '')}", "l")
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
canvas_PT.WaitPrimitive()
canvas_ETA.WaitPrimitive()


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
colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen, ROOT.kCoral, ROOT.kMagenta, ROOT.kCyan, ROOT.kBlack]


# Loop sobre os arquivos .rootcolors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen, ROOT.kCoral, ROOT.kMagenta, ROOT.kCyan, ROOT.kBlack]

for i, inputFile in enumerate(inputFiles):
    print(f"Reading file {inputFile}")

    # Criando cadeia de árvores ROOT
    chain = ROOT.TChain("Delphes")
    chain.Add(inputFile)

    # Criando objeto ExRootTreeReader para ler a cadeia de árvores
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

  # Criando histogramas para PT e ETA dos quarks b
    histogram_PT = ROOT.TH1F("hist_pt", "P_{T} b[2]; P_{T} (GeV/c); Number of events", 40, 0, 1000)
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
        if pt_b_quarks:
            # Ordena os PTs dos quarks b em ordem decrescente
            pt_b_quarks.sort(reverse=True)

            # Preenche o histograma com o PT do segundo quark b mais energético
            histogram_PT.Fill(pt_b_quarks[1])

            # Preenche o histograma com o ETA do segundo quark b mais energético
            histogram_ETA.Fill(eta_b_quarks[pt_b_quarks.index(pt_b_quarks[1])])
	    # Imprime o PT e o ETA de cada partícula com o segundo maior PT
	    for pt, eta in zip(pt_b_quarks, eta_b_quarks):
		if pt == pt_b_quarks[1]:
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
canvas_PT = ROOT.TCanvas("canvas_PT", "PT Histogram", 800, 600)
canvas_ETA = ROOT.TCanvas("canvas_ETA", "ETA Histogram", 800, 600)

# Configurando e desenhando histogramas de PT
canvas_PT.cd()
for i, histogram_PT in enumerate(histograms_PT):
    histogram_PT.SetLineColor(colors[i])
    histogram_PT.SetLineWidth(2)
    histogram_PT.Draw("SAME")
    histogram_PT.SetMinimum(0.1)  # Definindo escala logarítmica
    canvas_PT.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Calculando média e desvio padrão
    mean_pt = histogram_PT.GetMean()
    std_dev_pt = histogram_PT.GetStdDev()

    # Criando legenda personalizada
    legend_PT = ROOT.TLegend(0.9, 0.9 - 0.1 * i, 0.95, 0.65 - 0.1 * i)  # Posição fora da figura
    legend_PT.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_PT.SetBorderSize(0)  # Removendo borda da legenda
    legend_PT.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_PT.GetListOfFunctions().Remove(legend_PT)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda
    legend_PT.AddEntry(histogram_PT, f"{os.path.basename(inputFiles[i]).replace('.root', '')}", "l")
    legend_PT.AddEntry(0, f"Mean: {mean_pt:.2f}", "")
    legend_PT.AddEntry(0, f"Std Dev: {std_dev_pt:.2f}", "")
    legend_PT.Draw()

canvas_PT.Update()


# Configurando e desenhando histogramas de ETA
canvas_ETA.cd()
for i, histogram_ETA in enumerate(histograms_ETA):
    histogram_ETA.SetLineColor(colors[i])
    histogram_ETA.SetLineWidth(2)
    histogram_ETA.Draw("SAME")
    histogram_ETA.SetMinimum(0.1)  # Definindo escala logarítmica
    canvas_ETA.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Calculando média e desvio padrão
    mean_eta = histogram_ETA.GetMean()
    std_dev_eta = histogram_ETA.GetStdDev()

    # Criando legenda personalizada
    legend_ETA = ROOT.TLegend(0.9, 0.9 - 0.1 * i, 0.95, 0.65 - 0.1 * i)  # Posição fora da figura
    legend_ETA.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_ETA.SetBorderSize(0)  # Removendo borda da legenda
    legend_ETA.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_ETA.GetListOfFunctions().Remove(legend_ETA)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda
    legend_ETA.AddEntry(histogram_ETA, f"{os.path.basename(inputFiles[i]).replace('.root', '')}", "l")
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
canvas_PT.WaitPrimitive()
canvas_ETA.WaitPrimitive()

