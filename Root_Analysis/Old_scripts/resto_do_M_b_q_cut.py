 
           #     histogram_HT.Fill(massa_selecionada_1, peso_atual)
           #     histogram_HT.Fill(massa_selecionada_2, peso_atual)
           #     # Limpar a lista de massas invariantes para o próximo evento
           #     massas_invariantes.clear()

    
        

            
    # Obter o número total de eventos

    total_events = histogram_HT.GetEntries()


    # Adicionando histogramas à lista
    histograms_HT.append(histogram_HT)

    # Calculando média e desvio padrão para PT
    mean_ht = histogram_HT.GetMean()
    std_dev_ht = histogram_HT.GetStdDev()

    mean_values_HT.append(mean_ht)
    std_dev_values_HT.append(std_dev_ht)
    
# Criando canvas para os histogramas de PT e ETA
canvas_HT = ROOT.TCanvas("canvas_HT", "HT Histogram", 1600, 1000)


# Calcular o valor máximo de todos os histogramas
max_y_ht = max([histogram_HT.GetMaximum() for histogram_HT in histograms_HT])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de ETA com base no valor máximo calculado
for histogram_HT in histograms_HT:
    histogram_HT.SetMaximum(max_y_ht * margin_factor)




canvas_HT.cd()
for i, histogram_HT in enumerate(histograms_HT):
    histogram_HT.SetLineColor(colors[i])
    histogram_HT.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_HT.SetLineWidth(2)
    histogram_HT.Draw("SAME")
    histogram_HT.Sumw2(0)
    histogram_HT.SetMinimum(0.0001)  # Definindo escala logarítmica
    canvas_HT.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Calculando média e desvio padrão
    mean_ht = histogram_HT.GetMean()
    std_dev_ht = histogram_HT.GetStdDev()

    # Criando legenda personalizada
    legend_HT = ROOT.TLegend(0.92, 1.02 - 0.12 * i, 0.95, 0.83 - 0.12 * i) # Posição fora da figura
    legend_HT.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_HT.SetEntrySeparation(0.0005) 
    legend_HT.SetBorderSize(0)  # Removendo borda da legenda
    legend_HT.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_HT.GetListOfFunctions().Remove(legend_HT)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda


# Adicionando o nome do arquivo à legenda

    legend_HT.AddEntry(0, "", "")
    legend_HT.AddEntry(histogram_HT, f"  {os.path.basename(inputFiles[i]).replace('.root', '')}", "f100")
    legend_HT.AddEntry(0, f"Mean: {mean_ht:.1f}", "")
    #legend_PT.AddEntry(0, "", "")
    legend_HT.AddEntry(0, f"Std Dev: {std_dev_ht:.1f}", "")
    legend_HT.AddEntry(0, "", "")

    legend_HT.Draw()


canvas_HT.Update()

  # Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 100  
    # Definir o limite superior para PT e ETA com base no valor máximo de y encontrado nao funcionou
histogram_HT.SetMaximum(max_y_ht * margin_factor)


# Salvando os histogramas em arquivos de imagem
output_file_path_HT = os.path.join(output_dir, "M_b_nocut.png")
canvas_HT.SaveAs(output_file_path_HT)
# Abra os histogramas quando prontos
subprocess.Popen(["xdg-open", output_file_path_HT])










##################M dos light jets sem cut################################

mh = 125 
massas_invariantes=[]
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
    histogram_HT = ROOT.TH1F("histogram", "M of pairs of light quark jets; M (GeV/c^{²}); Number of events", 50, 0, 1500)

    # Loop sobre todos os eventos na árvore
    for event in range(numberOfEntries):
        # Obter o evento atual
        treeReader.ReadEntry(event)
        
        jets = []
        
        # Obter o galho "Jet"
        branchJet = treeReader.UseBranch("Jet")
        
        # Certificar-se de que existem pelo menos dois jatos no evento
        if branchJet.GetEntries() < 2:
            continue
        
        for jet in branchJet:
            if abs(jet.Flavor) in [0, 1, 2, 3, 4]:
                jets.append(jet)
        
        if jets:

            for jet in jets:
                combinacoes = itertools.combinations(jets, 2)
                for combination in combinacoes:
                    jet1, jet2 = combination
                    if jet1 != jet2  and ((abs(jet1.Flavor) == 1 and abs(jet2.Flavor) == 2) or (abs(jet1.Flavor) == 1 and abs(jet2.Flavor) == 3) or (abs(jet1.Flavor) == 4 and abs(jet2.Flavor) == 2) or (abs(jet1.Flavor) == 4 and abs(jet2.Flavor) == 3) or (abs(jet1.Flavor) == 1 and abs(jet2.Flavor) == 2) or (abs(jet1.Flavor)== 1 and abs(jet2.Flavor) == 3) or (abs(jet1.Flavor) == 4 and abs(jet2.Flavor) == 2) or ( abs(jet1.Flavor) == 4 and abs(jet2.Flavor) == 3) or (abs(jet1.Flavor) == 2 and abs(jet2.Flavor) == 1) or (abs(jet1.Flavor) == 3 and abs(jet2.Flavor) == 1) or (abs(jet1.Flavor) == 2 and abs(jet2.Flavor) == 4) or (abs(jet1.Flavor) == 3 and abs(jet2.Flavor) == 4) or (abs(jet1.Flavor) == 2 and abs(jet2.Flavor) == 1) or (abs(jet1.Flavor) == 3 and abs(jet2.Flavor) == 1) or (abs(jet1.Flavor) == 2 and abs(jet2.Flavor) == 4) or (abs(jet1.Flavor) == 3 and abs(jet2.Flavor) == 4)) :
        
                        mass_invariant = invariant_mass(jet1, jet2)
                        print("Massa invariante dos jatos no evento", event, ": ", mass_invariant) 
                        histogram_HT.Fill(mass_invariant, peso_atual)
        
    # Obter o número total de eventos

    total_events = histogram_HT.GetEntries()


    # Adicionando histogramas à lista
    histograms_HT.append(histogram_HT)

    # Calculando média e desvio padrão para PT
    mean_ht = histogram_HT.GetMean()
    std_dev_ht = histogram_HT.GetStdDev()

    mean_values_HT.append(mean_ht)
    std_dev_values_HT.append(std_dev_ht)
    
# Criando canvas para os histogramas de PT e ETA
canvas_HT = ROOT.TCanvas("canvas_HT", "HT Histogram", 1600, 1000)


# Calcular o valor máximo de todos os histogramas
max_y_ht = max([histogram_HT.GetMaximum() for histogram_HT in histograms_HT])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de ETA com base no valor máximo calculado
for histogram_HT in histograms_HT:
    histogram_HT.SetMaximum(max_y_ht * margin_factor)




canvas_HT.cd()
for i, histogram_HT in enumerate(histograms_HT):
    histogram_HT.SetLineColor(colors[i])
    histogram_HT.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_HT.SetLineWidth(2)
    histogram_HT.Draw("SAME")
    histogram_HT.Sumw2(0)
    histogram_HT.SetMinimum(0.0001)  # Definindo escala logarítmica
    canvas_HT.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Calculando média e desvio padrão
    mean_ht = histogram_HT.GetMean()
    std_dev_ht = histogram_HT.GetStdDev()

    # Criando legenda personalizada
    legend_HT = ROOT.TLegend(0.92, 1.02 - 0.12 * i, 0.95, 0.83 - 0.12 * i) # Posição fora da figura
    legend_HT.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_HT.SetEntrySeparation(0.0005) 
    legend_HT.SetBorderSize(0)  # Removendo borda da legenda
    legend_HT.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_HT.GetListOfFunctions().Remove(legend_HT)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda


# Adicionando o nome do arquivo à legenda

    legend_HT.AddEntry(0, "", "")
    legend_HT.AddEntry(histogram_HT, f"  {os.path.basename(inputFiles[i]).replace('.root', '')}", "f100")
    legend_HT.AddEntry(0, f"Mean: {mean_ht:.2f}", "")
    #legend_PT.AddEntry(0, "", "")
    legend_HT.AddEntry(0, f"Std Dev: {std_dev_ht:.2f}", "")
    legend_HT.AddEntry(0, "", "")

    legend_HT.Draw()


canvas_HT.Update()

  # Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 100  
    # Definir o limite superior para PT e ETA com base no valor máximo de y encontrado nao funcionou
histogram_HT.SetMaximum(max_y_ht * margin_factor)


# Salvando os histogramas em arquivos de imagem
output_file_path_HT = os.path.join(output_dir, "M_js_nocut.png")
canvas_HT.SaveAs(output_file_path_HT)
# Abra os histogramas quando prontos
subprocess.Popen(["xdg-open", output_file_path_HT])





##################M dos light jets com cut################################

mh = 80.37 
massas_invariantes=[]
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
    histogram_HT = ROOT.TH1F("histogram", "M light quark jets pairs closer to m_{W}; M (GeV/c^{²}); Number of events", 50, 0, 1500)

    for event in range(numberOfEntries):
        # Obter o evento atual
        treeReader.ReadEntry(event)
        
        jets = []
        
        # Obter o galho "Jet"
        branchJet = treeReader.UseBranch("Jet")
        
        # Certificar-se de que existem pelo menos dois jatos no evento
        if branchJet.GetEntries() < 2:
            continue
        
        for jet in branchJet:
            if abs(jet.Flavor) in [0, 1, 2, 3, 4]:
                jets.append(jet)
        
        if jets:

            for jet in jets:
                combinacoes = itertools.combinations(jets, 2)
                for combination in combinacoes:
                    jet1, jet2 = combination
                    if jet1 != jet2 and (
        (abs(jet1.Flavor) == 1 and abs(jet2.Flavor) == 2) or 
        (abs(jet1.Flavor) == 1 and abs(jet2.Flavor) == 3) or 
        (abs(jet1.Flavor) == 4 and abs(jet2.Flavor) == 2) or 
        (abs(jet1.Flavor) == 4 and abs(jet2.Flavor) == 3) or 
        (abs(jet1.Flavor) == 1 and abs(jet2.Flavor) == 2) or 
        (abs(jet1.Flavor) == 1 and abs(jet2.Flavor) == 3) or 
        (abs(jet1.Flavor) == 4 and abs(jet2.Flavor) == 2) or 
        (abs(jet1.Flavor) == 4 and abs(jet2.Flavor) == 3) or 
        (abs(jet1.Flavor) == 2 and abs(jet2.Flavor) == 1) or 
        (abs(jet1.Flavor) == 3 and abs(jet2.Flavor) == 1) or 
        (abs(jet1.Flavor) == 2 and abs(jet2.Flavor) == 4) or 
        (abs(jet1.Flavor) == 3 and abs(jet2.Flavor) == 4) or 
        (abs(jet1.Flavor) == 2 and abs(jet2.Flavor) == 1) or 
        (abs(jet1.Flavor) == 3 and abs(jet2.Flavor) == 1) or 
        (abs(jet1.Flavor) == 2 and abs(jet2.Flavor) == 4) or 
        (abs(jet1.Flavor) == 3 and abs(jet2.Flavor) == 4)
    ):
        
                        mass_invariant = invariant_mass(jet1, jet2)
                        print("Massa invariante dos jatos no evento", event, ": ", mass_invariant) 
                        massas_invariantes.append(mass_invariant)
        
        
        
        # Verificar se há massas invariantes suficientes para preencher o histograma
            if len(massas_invariantes) >= 2:
                # Ordenar as massas invariantes em ordem crescente de diferença absoluta com a massa do bóson W
                massas_invariantes.sort(key=lambda x: abs(x - mw))
                # Selecionar as duas massas invariantes mais próximas da massa do bóson W
                massa_selecionada_1 = massas_invariantes[0]
                massa_selecionada_2 = massas_invariantes[1]
                # Preencher o histograma com as duas massas selecionadas
                histogram_HT.Fill(massa_selecionada_1, 0.0041)
                histogram_HT.Fill(massa_selecionada_2, 0.0041)
                # Limpar a lista de massas invariantes para o próximo evento
                massas_invariantes.clear()
        
    # Obter o número total de eventos

    total_events = histogram_HT.GetEntries()


    # Adicionando histogramas à lista
    histograms_HT.append(histogram_HT)

    # Calculando média e desvio padrão para PT
    mean_ht = histogram_HT.GetMean()
    std_dev_ht = histogram_HT.GetStdDev()

    mean_values_HT.append(mean_ht)
    std_dev_values_HT.append(std_dev_ht)
    
# Criando canvas para os histogramas de PT e ETA
canvas_HT = ROOT.TCanvas("canvas_HT", "HT Histogram", 1600, 1000)


# Calcular o valor máximo de todos os histogramas
max_y_ht = max([histogram_HT.GetMaximum() for histogram_HT in histograms_HT])

# Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 1.2  

# Definir o limite superior para todos os histogramas de ETA com base no valor máximo calculado
for histogram_HT in histograms_HT:
    histogram_HT.SetMaximum(max_y_ht * margin_factor)




canvas_HT.cd()
for i, histogram_HT in enumerate(histograms_HT):
    histogram_HT.SetLineColor(colors[i])
    histogram_HT.SetFillColor(colors[i])  # Definindo a mesma cor de linha como cor de preenchimento
    histogram_HT.SetLineWidth(2)
    histogram_HT.Draw("SAME")
    histogram_HT.Sumw2(0)
    histogram_HT.SetMinimum(0.0001)  # Definindo escala logarítmica
    canvas_HT.SetLogy()  # Ativando escala logarítmica no eixo Y

    # Calculando média e desvio padrão
    mean_ht = histogram_HT.GetMean()
    std_dev_ht = histogram_HT.GetStdDev()

    # Criando legenda personalizada
    legend_HT = ROOT.TLegend(0.92, 1.02 - 0.12 * i, 0.95, 0.83 - 0.12 * i) # Posição fora da figura
    legend_HT.SetTextSize(0.02)  # Definindo tamanho do texto da legenda
    legend_HT.SetEntrySeparation(0.0005) 
    legend_HT.SetBorderSize(0)  # Removendo borda da legenda
    legend_HT.SetFillStyle(0)  # Removendo fundo da legenda

    # Desativando a posição padrão da legenda
    histogram_HT.GetListOfFunctions().Remove(legend_HT)

    # Adicionando nome do arquivo, cor e valores médios/desvios padrão à legenda


# Adicionando o nome do arquivo à legenda

    legend_HT.AddEntry(0, "", "")
    legend_HT.AddEntry(histogram_HT, f"  {os.path.basename(inputFiles[i]).replace('.root', '')}", "f100")
    legend_HT.AddEntry(0, f"Mean: {mean_ht:.2f}", "")
    #legend_PT.AddEntry(0, "", "")
    legend_HT.AddEntry(0, f"Std Dev: {std_dev_ht:.2f}", "")
    legend_HT.AddEntry(0, "", "")

    legend_HT.Draw()


canvas_HT.Update()

  # Definir o fator de margem (por exemplo, 1.2 para 20% de margem)
margin_factor = 100  
    # Definir o limite superior para PT e ETA com base no valor máximo de y encontrado nao funcionou
histogram_HT.SetMaximum(max_y_ht * margin_factor)


# Salvando os histogramas em arquivos de imagem
output_file_path_HT = os.path.join(output_dir, "M_js_cut.png")
canvas_HT.SaveAs(output_file_path_HT)
# Abra os histogramas quando prontos
subprocess.Popen(["xdg-open", output_file_path_HT])











