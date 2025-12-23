#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np
import awkward as ak
import uproot
import ROOT

# Cross sections in pb for each QCD bin
XSEC_MAP = {
    "QCD_Bin-PT-15to20": 885700000.0,
    "QCD_Bin-PT-20to30": 415700000.0,
    "QCD_Bin-PT-30to50": 112300000.0,
    "QCD_Bin-PT-50to80": 16730000.0,
    "QCD_Bin-PT-80to120": 2506000.0,
    "QCD_Bin-PT-120to170": 439800.0,
    "QCD_Bin-PT-170to300": 113300.0,
    "QCD_Bin-PT-300to470": 7581.0,
    "QCD_Bin-PT-470to600": 623.3,
    "QCD_Bin-PT-600to800": 178.7,
    "QCD_Bin-PT-800to1000": 30.62,
    "QCD_Bin-PT-1000to1500": 9.306,
    "QCD_Bin-PT-1500to2000": 0.5015,
    "QCD_Bin-PT-2000to2500": 0.04264,
    "QCD_Bin-PT-2500to3000": 0.004454,
    "QCD_Bin-PT-3000": 0.0005539,
}

LUMINOSITY = 112700.0  # pb^-1

def get_xsec(dirname):
    """Extract cross section from directory name"""
    for key in XSEC_MAP:
        if key in dirname:
            return XSEC_MAP[key]
    print(f"[WARNING] No cross section found for {dirname}, using 1.0")
    return 1.0

def calculate_weight(n_entries, xsec, lumi):
    return 1
    """Calculate event weight"""
    if n_entries == 0:
        return 0.0
    return (xsec * lumi) / n_entries

def flat(x):
    """Flatten awkward array to numpy"""
    return ak.flatten(x, axis=None).to_numpy()

def calculate_deltaR(eta1, phi1, eta2, phi2):
    """Calculate deltaR between two objects"""
    deta = eta1 - eta2
    dphi = np.arctan2(np.sin(phi1 - phi2), np.cos(phi1 - phi2))
    return np.sqrt(deta**2 + dphi**2)

def match_sv_to_jet(sv_eta, sv_phi, jet_eta, jet_phi, dr_threshold=0.4):
    """
    Match SVs to jets using deltaR < 0.4
    Returns: (matched_mask, matched_jet_indices)
    """
    # Broadcast to calculate all deltaR values
    # sv_eta: shape (n_events, n_sv)
    # jet_eta: shape (n_events, n_jets)
    
    # Expand dimensions for broadcasting
    sv_eta_exp = sv_eta[:, :, np.newaxis]  # (n_events, n_sv, 1)
    sv_phi_exp = sv_phi[:, :, np.newaxis]
    jet_eta_exp = jet_eta[:, np.newaxis, :]  # (n_events, 1, n_jets)
    jet_phi_exp = jet_phi[:, np.newaxis, :]
    
    # Calculate deltaR for all SV-jet pairs
    dr = calculate_deltaR(sv_eta_exp, sv_phi_exp, jet_eta_exp, jet_phi_exp)
    
    # Find minimum deltaR for each SV
    min_dr = ak.min(dr, axis=2)
    matched_mask = min_dr < dr_threshold
    
    # Find which jet each SV matched to
    matched_jet_idx = ak.argmin(dr, axis=2)
    
    return matched_mask, matched_jet_idx

def get_jet_flavor_for_sv(matched_mask, matched_jet_idx, jet_partonFlavour):
    """
    Get the parton flavour of the matched jet for each SV
    Returns array with same shape as SVs, -999 for unmatched
    """
    # Initialize with -999 (unmatched)
    sv_jet_flavour = ak.where(
        matched_mask,
        jet_partonFlavour[ak.local_index(matched_jet_idx, axis=1), matched_jet_idx],
        -999
    )
    return sv_jet_flavour

def process_file(input_file, output_file, dirname):
    """Process a single ROOT file and create histograms"""
    print(f"[INFO] Processing {input_file}")
    
    xsec = get_xsec(dirname)
    
    # Open input file
    with uproot.open(input_file, handler=uproot.source.xrootd.XRootDSource, timeout=300) as f:
        tree = f["analyzer/tree"]
        n_entries = tree.num_entries
        weight = calculate_weight(n_entries, xsec, LUMINOSITY)
        print(f"[INFO] Events: {n_entries}, XSec: {xsec} pb, Weight: {weight}")
        
        # Read data
        try:
            data = tree.arrays(
                ["jet_pt", "jet_eta", "jet_phi", "jet_partonFlavour",
                 "sv_mass", "sv_pt", "sv_eta", "sv_phi", "sv_lxy",
                 "lep_pt", "lep_eta", "lep_phi", 
                 "lep_isGlobalMuon"],
                library="ak"
            )
        except Exception as e:
            print(f"[ERROR] Failed to read branches: {e}")
            raise
    
    # Determine lepton type for each SV
    # lep_isGlobalMuon is (n_events, n_sv, 2)
    # If both tracks have isGlobalMuon >= 0, it's dimuon
    # If both tracks have isGlobalMuon == -1, it's dielectron
    is_muon_track = data["lep_isGlobalMuon"] >= 0
    is_dimuon = ak.all(is_muon_track, axis=2)  # (n_events, n_sv)
    is_dielectron = ak.all(~is_muon_track, axis=2)  # (n_events, n_sv)
    
    # Quality cuts - at least one SV per event
    mask_nonempty = ak.num(data["sv_pt"], axis=1) > 0
    data = data[mask_nonempty]
    is_dimuon = is_dimuon[mask_nonempty]
    is_dielectron = is_dielectron[mask_nonempty]
    
    n_events_pass = len(data)
    print(f"[INFO] Events passing cuts: {n_events_pass}/{n_entries}")
    
    if n_events_pass == 0:
        print("[WARNING] No events passed cuts")
        return
    
    # Match SVs to jets
    matched_to_jet, matched_jet_idx = match_sv_to_jet(
        data["sv_eta"], data["sv_phi"],
        data["jet_eta"], data["jet_phi"]
    )
    
    # Get jet flavour for each SV
    sv_jet_flavour = get_jet_flavor_for_sv(
        matched_to_jet, matched_jet_idx, data["jet_partonFlavour"]
    )
    
    # Create masks for jet flavours
    in_b_jet = matched_to_jet & (abs(sv_jet_flavour) == 5)
    in_c_jet = matched_to_jet & (abs(sv_jet_flavour) == 4)
    in_light_jet = matched_to_jet & ((abs(sv_jet_flavour) < 4) | (sv_jet_flavour == 21))
    outside_jet = ~matched_to_jet
    
    # Create output ROOT file
    out_file = ROOT.TFile(output_file, "RECREATE")
    
    # Lxy bins
    lxy_bins = [(0, 0.2), (0.2, 1), (1, 2.4), (2.4, 3.1), (3.1, 7), (7, 10)]
    
    # Lepton types
    lepton_types = ['muon', 'electron']
    
    # Jet categories
    jet_categories = ['in_b', 'in_c', 'in_light', 'outside']
    
    histograms = {}
    
    # Create histograms
    # Structure: h_invm_{lepton}_{jet_cat}_lxy{i}
    for lep_type in lepton_types:
        for jet_cat in jet_categories:
            for i, (lxy_min, lxy_max) in enumerate(lxy_bins):
                hname = f"h_invm_{lep_type}_{jet_cat}_lxy{i}"
                title = f"Invariant Mass ({lep_type}, {jet_cat}, L_{{xy}}=[{lxy_min},{lxy_max}] cm)"
                histograms[hname] = ROOT.TH1F(
                    hname,
                    f"{title};M(ll) [GeV];Events",
                    100, 0, 10
                )
                histograms[hname].Sumw2()
    
    # Also create inclusive histograms (all lxy)
    for lep_type in lepton_types:
        for jet_cat in jet_categories:
            hname = f"h_invm_{lep_type}_{jet_cat}_inclusive"
            title = f"Invariant Mass ({lep_type}, {jet_cat})"
            histograms[hname] = ROOT.TH1F(
                hname,
                f"{title};M(ll) [GeV];Events",
                100, 0, 10
            )
            histograms[hname].Sumw2()
    
    # Fill histograms
    for lep_type in lepton_types:
        # Select lepton type
        if lep_type == 'muon':
            lep_mask = is_dimuon
        else:  # electron
            lep_mask = is_dielectron
        
        for jet_cat in jet_categories:
            # Select jet category
            if jet_cat == 'in_b':
                jet_mask = in_b_jet
            elif jet_cat == 'in_c':
                jet_mask = in_c_jet
            elif jet_cat == 'in_light':
                jet_mask = in_light_jet
            else:  # outside
                jet_mask = outside_jet
            
            # Combined mask
            combined_mask = lep_mask & jet_mask
            
            if ak.sum(combined_mask) == 0:
                continue
            
            # Get invariant masses for this category
            invm_values = data["sv_mass"][combined_mask]
            lxy_values = data["sv_lxy"][combined_mask]
            
            # Fill inclusive histogram
            hname = f"h_invm_{lep_type}_{jet_cat}_inclusive"
            for val in flat(invm_values):
                #histograms[hname].Fill(val, weight)
                histograms[hname].Fill(val)
            
            # Fill lxy-binned histograms
            for i, (lxy_min, lxy_max) in enumerate(lxy_bins):
                lxy_bin_mask = (lxy_values >= lxy_min) & (lxy_values < lxy_max)
                invm_in_bin = invm_values[lxy_bin_mask]
                
                hname = f"h_invm_{lep_type}_{jet_cat}_lxy{i}"
                for val in flat(invm_in_bin):
                    #histograms[hname].Fill(val, weight)
                    histograms[hname].Fill(val)
    
    # Print some statistics
    print(f"[INFO] Total SVs: {ak.sum(ak.num(data['sv_mass'], axis=1))}")
    print(f"[INFO] Dimuon SVs: {ak.sum(ak.sum(is_dimuon, axis=1))}")
    print(f"[INFO] Dielectron SVs: {ak.sum(ak.sum(is_dielectron, axis=1))}")
    print(f"[INFO] SVs in b-jets: {ak.sum(ak.sum(in_b_jet, axis=1))}")
    print(f"[INFO] SVs in c-jets: {ak.sum(ak.sum(in_c_jet, axis=1))}")
    print(f"[INFO] SVs in light-jets: {ak.sum(ak.sum(in_light_jet, axis=1))}")
    print(f"[INFO] SVs outside jets: {ak.sum(ak.sum(outside_jet, axis=1))}")
    
    # Write histograms
    for hist in histograms.values():
        hist.Write()
    
    out_file.Close()
    print(f"[INFO] Created {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Create histograms from SV analysis ROOT files')
    parser.add_argument('input_file', help='Input ROOT file path')
    parser.add_argument('--dirname', required=True, help='Directory name for cross section lookup')
    parser.add_argument('--output-number', type=int, required=True, help='Sequential output file number')
    parser.add_argument('--output-file', default=None, help='Override output file path')
    
    args = parser.parse_args()
    input_path = args.input_file
    
    # Allow override of output path (for condor jobs)
    if args.output_file:
        output_path = args.output_file
    else:
        current_dir = os.getcwd()
        output_dir = os.path.join(current_dir, "histo", args.dirname)
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"histo_{args.output_number}.root")
    
    try:
        process_file(input_path, output_path, args.dirname)
        return 0
    except Exception as e:
        print(f"[ERROR] Failed to process {input_path}: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
