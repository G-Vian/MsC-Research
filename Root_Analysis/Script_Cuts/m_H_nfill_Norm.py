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
import subprocess # to allow me to open the pictures and not interrupt the code 
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
#	Creating important directories 
############################################################

# Script directory
script_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(script_dir, "M")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
output_file_path = os.path.join(output_dir, "processing_log.txt")

# Delete the old log file, if it exists
if os.path.exists(output_file_path):
    os.remove(output_file_path)
##############################################################
##############################################################

# Define the function to compute the invariant mass
def invariant_mass(jet1, jet2):
    # Extract momentum and energy components for jets j1 and j2 --> All these quantities are accessible in the jet branch
    p1x, p1y, p1z, E1 = jet1.PT * ROOT.TMath.Cos(jet1.Phi), jet1.PT * ROOT.TMath.Sin(jet1.Phi), jet1.PT * ROOT.TMath.SinH(jet1.Eta), jet1.PT * ROOT.TMath.CosH(jet1.Eta)
    p2x, p2y, p2z, E2 = jet2.PT * ROOT.TMath.Cos(jet2.Phi), jet2.PT * ROOT.TMath.Sin(jet2.Phi), jet2.PT * ROOT.TMath.SinH(jet2.Eta), jet2.PT * ROOT.TMath.CosH(jet2.Eta)

    # Create four-momentum vectors for jets j1 and j2 --> this is a Lorentz vector: https://root.cern.ch/doc/master/classTLorentzVector.html
    jet1_vec = ROOT.Math.PxPyPzEVector(p1x, p1y, p1z, E1)
    jet2_vec = ROOT.Math.PxPyPzEVector(p2x, p2y, p2z, E2)

    # Compute the sum of the two jets' momentum vectors
    total_vec = jet1_vec + jet2_vec

    # Compute the invariant mass (total_vec.M() computes the invariant mass using the energy and momentum contained in the Lorentz vector total_vec, this is defined in ROOT)
    mass = total_vec.M()

    return mass


# This function looks at each jet in the "additional_jets" list and computes its invariant mass with each jet in the "jets" list. Then, it chooses the jet from "additional_jets"
# that, when combined with jets from the "jets" list, results in an invariant mass closest to the Higgs mass.

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
# List to store histograms
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


# Disable default statistics
ROOT.gStyle.SetOptStat(0)

# Enable display of mean value and custom standard deviation
ROOT.gStyle.SetOptFit(1111)

# Colors for each process
# with tt , ttbb & ttWW
#colors =  [ ROOT.kMagenta+3, ROOT.kMagenta, ROOT.kOrange+4, ROOT.kBlue,  ROOT.kGreen+3, ROOT.kAzure+3,  ROOT.kOrange+7, ROOT.kTeal-3, ROOT.kViolet, ROOT.kYellow, ROOT.kBlack, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]
# without tt , ttbb & ttWW
colors =  [   ROOT.kOrange+4, ROOT.kBlue,  ROOT.kGreen+3, ROOT.kAzure+3,  ROOT.kOrange+7, ROOT.kTeal-3, ROOT.kYellow, ROOT.kBlack, ROOT.kCyan,  ROOT.kRed ]  #ROOT.kBlue,ROOT.kOrange+4, ROOT.kViolet+5, ROOT.kAzure-7]


#########################M of bottoms with cut#########################
delta_h = 40
mh = 125 

   

with open(output_file_path, "a") as output_file:

    for i, inputFile in enumerate(inputFiles):
        print(f"Reading file {inputFile}")
     #   histogram_HT = ROOT.TH1F("hist_ht", "Missing transverse energy ; MET (GeV); Number of events (L = 200 fb^{-1})", 80, 0, 500)
    # Global histograms for all events
        hist1 = ROOT.TH1F("hist1", "Invariant mass of h[1];  M_{j_b j_b} (GeV); Number of events (Normalized to one)", 70, 70, 170)
        hist2 = ROOT.TH1F("hist2", "Invariant mass of h[2];  M_{j_b j_b} (GeV); Number of events (Normalized to one)", 70, 70, 170)
    
        # Get the weight corresponding to the current file
        peso_atual = 1
        
        # Create ROOT tree chain
        chain = ROOT.TChain("Delphes")
        chain.Add(inputFile)

        # Create ExRootTreeReader object to read the tree chain
        treeReader = ROOT.ExRootTreeReader(chain)
        numberOfEntries = treeReader.GetEntries()


        # Branch definitions must be done only once before starting the loop. It is not necessary (and not recommended) to call UseBranch in each iteration of the loop, as it can cause overhead.

        # GETTING BRANCHES
        branchJet = treeReader.UseBranch("Jet")
        branchMissingET = treeReader.UseBranch("MissingET")
        branchHT = treeReader.UseBranch("ScalarHT")
        branchelec = treeReader.UseBranch("Electron")
        branchmu = treeReader.UseBranch("Muon")


        # Counters to track how many events were selected
        total_events = 0
        events_passed_MET_HT_nJets = 0
        events_passed_final_selection = 0


        # Loop over events
        for entry in range(numberOfEntries):
            print(f"Processing event {entry} of file {inputFile}")
            total_events += 1  # Increment total event counter

            # Read specific entry
            treeReader.ReadEntry(entry)


            # GETTING QUANTITIES FROM EACH BRANCH
                
            missingET = branchMissingET.At(0)  # It is assumed that there is only one MET object per event
            MET_value = missingET.MET        
            nJets = branchJet.GetEntries()
            ht_object = branchHT.At(0)  # It is assumed that there is only one HT object per event
            HT_value = ht_object.HT        #jet_PT = branchJet.PT #Not needed for these objects
            #jet_ETA = branchJet.ETA
            #elec_PT = branchelec.PT
            #mu_PT = branchmu.PT
            #elec_ETA = branchelec.ETA
            #mu_ETA = branchmu.ETA
            #LIST OF SELECTED JETS AND LEPTONS THAT SHOULD BE RESET FOR EACH NEW EVENT
            selected_jets = []
            q_jets = []
            b_jets = []
            electrons = []
            muons = []
            
            ###CUT FLOW###

            # APPLYING CUTS (CHECKING IF THE EVENT SATISFIES THE CRITERIA)
            # FIRST CONDITIONAL
            if   nJets >= 6 and HT_value > 500 : 
                events_passed_MET_HT_nJets += 1 # to count how many events passed
                selected_jets = [jet for jet in branchJet if jet.PT >= 40 and abs(jet.Eta) <= 2.5 ]        
                q_jets = [jet for jet in selected_jets if abs(jet.Flavor)  in [0, 1, 2, 3, 4]  ] 
                b_jets = [jet for jet in selected_jets if abs(jet.Flavor) == 5]
                electrons = [electron for electron in branchelec if electron.PT > 20 and abs(electron.Eta) <= 2.5  ] 
                muons = [muon for muon in branchmu if muon.PT > 20 and abs(muon.Eta) <= 2.5 ] 

                # SECOND CONDITIONAL, INSIDE THE FIRST, THAT CHECKS THE NUMBER OF OBJECTS
                if len(q_jets) >= 2 and len(b_jets) >= 4 and len(muons) == 0 and len(electrons) == 0 : 
                    events_passed_final_selection += 1  # Increment counter for events that pass this final selection

                    # Compute combinations of b-jets and their invariant masses
                    combinacoes = itertools.combinations(b_jets, 2)
                    massas_invariantes = [(invariant_mass(jet1, jet2), jet1, jet2) for jet1, jet2 in combinacoes]

                    # List to store X values and their corresponding combinations
                    valores_X = []

                    # Generate all possible combinations of two different jet pairs
                    pares_de_pares = list(itertools.combinations(massas_invariantes, 2))    
                    for (massa1, jet1a, jet1b), (massa2, jet2a, jet2b) in pares_de_pares:
                        # Make sure the jets are different
                        if len({jet1a, jet1b, jet2a, jet2b}) == 4:
                            # Compute X
                            X = ((massa1 - mh)**2 / delta_h**2) + ((massa2 - mh)**2 / delta_h**2)
                            # Store the X value and corresponding pairs
                            valores_X.append((X, (jet1a, jet1b), (jet2a, jet2b)))

                    # Find the combination that minimizes X
                    if valores_X:
                        min_X_combination = min(valores_X, key=lambda x: x[0])

                        # Check if the combinations satisfy the additional criteria
                        X_min, (jet1a, jet1b), (jet2a, jet2b) = min_X_combination
                        massa1 = invariant_mass(jet1a, jet1b)
                        massa2 = invariant_mass(jet2a, jet2b)
                        PT_1 = jet1a.PT + jet1b.PT
                        PT_2 = jet2a.PT + jet2b.PT
                        if abs(massa1 - mh) <= delta_h and abs(massa2 - mh) <= delta_h and PT_1 > PT_2:
                            # Fill histograms with weights
                            hist1.Fill(massa1, peso_atual)
                            hist2.Fill(massa2, peso_atual)
                        if abs(massa1 - mh) <= delta_h and abs(massa2 - mh) <= delta_h and PT_1 < PT_2:
                            hist1.Fill(massa2, peso_atual)
                            hist2.Fill(massa1, peso_atual)
                        # Print results for verification
                        print(f"File {inputFile}, event {entry}: X = {X_min}")
                        print(f"   Jet pair 1: Invariant Mass = {massa1}, Jet 1: PT = {jet1a.PT}, Eta = {jet1a.Eta}, Phi = {jet1a.Phi}, Jet 2: PT = {jet1b.PT}, Eta = {jet1b.Eta}, Phi = {jet1b.Phi}")
                        print(f"   Jet pair 2: Invariant Mass = {massa2}, Jet 1: PT = {jet2a.PT}, Eta = {jet2a.Eta}, Phi = {jet2a.Phi}, Jet 2: PT = {jet2b.PT}, Eta = {jet2b.Eta}, Phi = {jet2b.Phi}")
                    massas_invariantes.clear()



            else:
                    continue  # If the first condition is not satisfied, go to the next event
        print(f"Total events processed: {total_events}")
        print(f"Events passing MET >= 40, HT > 500 and nJets >= 6: {events_passed_MET_HT_nJets}")
        print(f"Events passing the final selection: {events_passed_final_selection}")
        # Save the results in a text file for the current input file
        output_file.write(f"Results for file: {inputFile}\n")
        output_file.write(f"Total de eventos processados: {total_events}\n")
        output_file.write(f"Eventos que passaram MET >= 40, HT > 500 e nJets >= 6: {events_passed_MET_HT_nJets}\n")
        output_file.write(f"Eventos que passaram a seleção final: {events_passed_final_selection}\n")
        output_file.write("-" * 40 + "\n")

        # Normalize histograms after filling all events
        if hist1.Integral() > 0:
            hist1.Scale(1.0 / hist1.Integral())
        if hist2.Integral() > 0:
            hist2.Scale(1.0 / hist2.Integral())

        # Add histograms to the list
        histograms_1.append(hist1)
        histograms_2.append(hist2)

        # Calculate mean and standard deviation
        mean_1 = hist1.GetMean()
        std_dev_1 = hist1.GetStdDev()
        mean_2 = hist2.GetMean()
        std_dev_2 = hist2.GetStdDev()

        mean_values_1.append(mean_1)
        std_dev_values_1.append(std_dev_1)
        mean_values_2.append(mean_2)
        std_dev_values_2.append(std_dev_2)

        # Number of entries
        n_entries1 = hist1.GetEntries()

        # Modal value (bin with the highest number of entries)
        max_bin_content = 0
        modal_value1 = 0

        # Loop over all bins of the histogram
        for bin in range(1, hist1.GetNbinsX() + 1):
            bin_content = hist1.GetBinContent(bin)
            if bin_content > max_bin_content:
                max_bin_content = bin_content
                modal_value1 = hist1.GetBinCenter(bin)

            min_value1 = float('inf')
            max_value1 = -float('inf')

        # Loop over histogram bins to find min and max x-axis values
        for bin in range(1, hist1.GetNbinsX() + 1):
            bin_content = hist1.GetBinContent(bin)
            # Check if bin content is non-zero
            if bin_content != 0:
                bin_x_value = hist1.GetXaxis().GetBinLowEdge(bin)
                min_value1 = min(min_value1, bin_x_value)
                max_value1 = max(max_value1, hist1.GetXaxis().GetBinUpEdge(bin))

        modal_values_1.append(modal_value1)
        N_entries1.append(n_entries1)
        min_values1.append(min_value1)
        max_values1.append(max_value1)

        # Number of entries
        n_entries2 = hist2.GetEntries()

        # Modal value (bin with the highest number of entries)
        max_bin_content = 0
        modal_value2 = 0

        # Loop over all bins of the histogram
        for bin in range(1, hist2.GetNbinsX() + 1):
            bin_content = hist2.GetBinContent(bin)
            if bin_content > max_bin_content:
                max_bin_content = bin_content
                modal_value2 = hist2.GetBinCenter(bin)

            min_value2 = float('inf')
            max_value2 = -float('inf')

        # Loop over histogram bins to find min and max x-axis values
        for bin in range(1, hist2.GetNbinsX() + 1):
            bin_content = hist2.GetBinContent(bin)
            # Check if bin content is non-zero
            if bin_content != 0:
                bin_x_value = hist2.GetXaxis().GetBinLowEdge(bin)
                min_value2 = min(min_value2, bin_x_value)
                max_value2 = max(max_value2, hist1.GetXaxis().GetBinUpEdge(bin))

        modal_values_2.append(modal_value2)
        N_entries2.append(n_entries2)
        min_values2.append(min_value2)
        max_values2.append(max_value2)

# Sum all elements except the last one in N_entries1
sum_except_last1 = sum(N_entries1[:-1])

# Calculate ratio R
if sum_except_last1 != 0:
    R1 = N_entries1[-1] / sum_except_last1
else:
    R1 = float('inf')  # Avoid division by zero

# Sum all elements except the last one in N_entries2
sum_except_last2 = sum(N_entries2[:-1])

# Calculate ratio R
if sum_except_last2 != 0:
    R2 = N_entries2[-1] / sum_except_last2
else:
    R2 = float('inf')  # Avoid division by zero

# Creating canvas for PT and ETA histograms
canvas_width = 2200
canvas_height = 1300
canvas_1 = ROOT.TCanvas("canvas_1", "h1", canvas_width, canvas_height)
canvas_2 = ROOT.TCanvas("canvas_2", "h2", canvas_width, canvas_height)

# Adjust canvas margins to create more space on the right
canvas_1.SetLeftMargin(0.1)
canvas_1.SetRightMargin(0.19)
canvas_1.SetBottomMargin(0.15)
canvas_1.SetTopMargin(0.1)

# Adjust canvas margins to create more space on the right
canvas_2.SetLeftMargin(0.1)
canvas_2.SetRightMargin(0.19)
canvas_2.SetBottomMargin(0.15)
canvas_2.SetTopMargin(0.1)

# Compute max value among all histograms
max_y_1 = max([hist1.GetMaximum() for hist1 in histograms_1])
max_y_2 = max([hist2.GetMaximum() for hist2 in histograms_2])

# Set margin factor (e.g., 1.2 for 20% margin)
margin_factor = 1.2

# Set upper limit for all histograms based on computed max
for hist1 in histograms_1:
    hist1.SetMaximum(max_y_1 * margin_factor)
for hist2 in histograms_2:
    hist2.SetMaximum(max_y_2 * margin_factor)

# Drawing histograms on canvas_1
canvas_1.cd()
for i, hist1 in enumerate(histograms_1):
    hist1.SetLineColor(colors[i])
    hist1.SetLineWidth(2)
    hist1.Sumw2(0)
    hist1.Draw("SAME")

    # Create custom legend
    n_entries = len(histograms_1)  # Number of legend entries
    y1 = 0.956
    y2 = 0.01

    # Distribute uniformly
    legend_PT = ROOT.TLegend(0.835, y1 - i * (y1 - y2) / n_entries, 0.92, y1 - (i + 1) * (y1 - y2) / n_entries)
    legend_PT.SetTextSize(0.023)
    legend_PT.SetEntrySeparation(0.8)
    legend_PT.SetBorderSize(0)
    legend_PT.SetFillStyle(0)
    hist1.GetListOfFunctions().Remove(legend_PT)

    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(hist1, f"  {os.path.basename(inputFiles[i]).replace('100.root', '')}", "f100")
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

canvas_1.Update()

# Save histogram to image file
output_file_path_1 = os.path.join(output_dir, "mh1_Nfill_norm.png")
canvas_1.SaveAs(output_file_path_1)
subprocess.Popen(["xdg-open", output_file_path_1])

# Drawing histograms on canvas_2
canvas_2.cd()
for i, hist2 in enumerate(histograms_2):
    hist2.SetLineColor(colors[i])
    hist2.SetLineWidth(2)
    hist2.Sumw2(0)
    hist2.Draw("SAME")

    # Create custom legend
    n_entries = len(histograms_2)
    y1 = 0.956
    y2 = 0.01

    legend_PT = ROOT.TLegend(0.835, y1 - i * (y1 - y2) / n_entries, 0.92, y1 - (i + 1) * (y1 - y2) / n_entries)
    legend_PT.SetTextSize(0.023)
    legend_PT.SetEntrySeparation(0.8)
    legend_PT.SetBorderSize(0)
    legend_PT.SetFillStyle(0)
    hist2.GetListOfFunctions().Remove(legend_PT)

    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(hist2, f"  {os.path.basename(inputFiles[i]).replace('100.root', '')}", "f100")
    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(0, f"Mean: {mean_values_2[i]:.1f}, Std Dev: {std_dev_values_2[i]:.1f}", "")
    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(0, f"Min: {min_values2[i]:.1f}, Max: {max_values2[i]:.1f}", "")
    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(0, f"Entr.: {N_entries2[i]:.1f}, Mod. V: {modal_values_2[i]:.1f}", "")
    legend_PT.AddEntry(0, "", "")
    legend_PT.AddEntry(0, "", "")
    legend_PT.Draw()

text = ROOT.TLatex()
text.SetNDC()
text.SetTextSize(0.03)
text.DrawLatex(0.80, 0.01, f"R: {R1:.3f}")

canvas_2.Update()

# Save histogram to image file
output_file_path_2 = os.path.join(output_dir, "mh2_Nfill_norm.png")
canvas_2.SaveAs(output_file_path_2)
subprocess.Popen(["xdg-open", output_file_path_2])
