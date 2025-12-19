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
#https://xsecdb-xsdb-official.app.cern.ch/xsdb/?columns=67108863&currentPage=0&pageSize=10&searchQuery=process_name%3DQCD_Bin-PT-30to50_TuneCP5_13p6TeV_pythia8

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

def get_jet_flavor_mask(partFlav, flavor_type):
    """Get mask for jet flavor"""
    if flavor_type == 'inclusive':
        return ak.ones_like(partFlav, dtype=bool)
    elif flavor_type == 'b':
        return abs(partFlav) == 5
    elif flavor_type == 'c':
        return abs(partFlav) == 4
    elif flavor_type == 'light':
        return (abs(partFlav) < 4) | (partFlav == 21)
    else:
        raise ValueError(f"Unknown flavor type: {flavor_type}")

def process_file(input_file, output_file, dirname):
    """Process a single ROOT file and create histograms"""
    print(f"[INFO] Processing {input_file}")
    
    xsec = get_xsec(dirname)
    
    # Open input file
    with uproot.open(input_file, timeout=300) as f:
        tree = f["analyzer/tree"]
        n_entries = tree.num_entries
        weight = calculate_weight(n_entries, xsec, LUMINOSITY)
        print(f"[INFO] Events: {n_entries}, XSec: {xsec} pb, Weight: {weight}")
        
        # Read data - try reading branches separately to avoid format issues
        try:
            data = tree.arrays(
                ["lxy_SV", "dR_ptl", "invm_ptl", "mass_SV", "pt_trk", "pt_SV", 
                 "eta_trk", "normChi2_SV", "normChi2_trk", "dxySig_trk", 
                 "hits_trk", "mom_ptl", "ip2d_mu", "vtxProb_SV", 
                 "dotProd_ptl", "dotProd_SV", "partFlav_jet"],
                library="ak"
            )
        except KeyError as e:
            print(f"[WARNING] KeyError reading branches, trying alternative method: {e}")
            # Read branches individually
            data = {}
            branch_names = ["lxy_SV", "dR_ptl", "invm_ptl", "mass_SV", "pt_trk", "pt_SV", 
                           "eta_trk", "normChi2_SV", "normChi2_trk", "dxySig_trk", 
                           "hits_trk", "mom_ptl", "ip2d_mu", "vtxProb_SV", 
                           "dotProd_ptl", "dotProd_SV", "partFlav_jet"]
            
            for branch_name in branch_names:
                try:
                    data[branch_name] = tree[branch_name].array(library="ak")
                except Exception as e:
                    print(f"[ERROR] Failed to read branch {branch_name}: {e}")
                    raise
            
            # Convert dict to ak.Record
            data = ak.zip(data, depth_limit=1)
    
    # Quality cuts
    mask_pt = ak.all(data["pt_trk"] >= 5, axis=2)
    mask_pt_evt = ak.all(mask_pt, axis=1)
    mask_nonempty = ak.count(data["pt_SV"], axis=1) > 0
    data = data[mask_nonempty & mask_pt_evt]
    
    n_events_pass = len(data)
    print(f"[INFO] Events passing cuts: {n_events_pass}/{n_entries}")
    
    if n_events_pass == 0:
        print("[WARNING] No events passed cuts")
    
    # Check if mother IDs are same
    same_in_sv = ak.all(data["mom_ptl"] == ak.firsts(data["mom_ptl"], axis=2), axis=2)
    
    # Create output ROOT file
    out_file = ROOT.TFile(output_file, "RECREATE")
    
    # Histogram configuration
    hist_config = {
        "vtxProb_SV": (100, 0, 1, "Vertex Probability", "vtxProb"),
        "dR_ptl": (50, 0, 0.5, "#DeltaR(#mu_{1}, #mu_{2})", "dR_ptl"),
        "invm_ptl": (100, 0, 10, "Invariant Mass [GeV]", "invm_ptl"),
        "pt_trk": (80, 0, 200, "Track p_{T} [GeV]", "pt_trk"),
        "eta_trk": (60, -3, 3, "Track #eta", "eta_trk"),
        "lxy_SV": (50, 0, 10, "L_{xy} [cm]", "lxy_SV"),
    }
    
    flavor_types = ['inclusive', 'b', 'c', 'light']
    mother_types = ['same', 'diff']
    histograms = {}
    
    # Create histograms
    for flavor in flavor_types:
        for mother in mother_types:
            for var, (nbins, xmin, xmax, title, short) in hist_config.items():
                hname = f"h_{short}_{flavor}_{mother}"
                histograms[hname] = ROOT.TH1F(hname, f"{title};{title};Events", 
                                               nbins, xmin, xmax)
                histograms[hname].Sumw2()
    
    # Histograms with vtxProb > 0.2 cut
    for flavor in flavor_types:
        for mother in mother_types:
            for var, (nbins, xmin, xmax, title, short) in hist_config.items():
                if var == "vtxProb_SV":
                    continue
                hname = f"h_{short}_{flavor}_{mother}_vtxcut"
                histograms[hname] = ROOT.TH1F(hname, 
                                               f"{title} (vtxProb>0.2);{title};Events",
                                               nbins, xmin, xmax)
                histograms[hname].Sumw2()
    
    # Invm vs lxy histograms
    lxy_bins = [(0, 0.2), (0.2, 1), (1, 2.4), (2.4, 3.1), (3.1, 7), (7, 10)]
    for flavor in flavor_types:
        for mother in mother_types:
            for i, (lxy_min, lxy_max) in enumerate(lxy_bins):
                hname = f"h_invm_lxy{i}_{flavor}_{mother}"
                histograms[hname] = ROOT.TH1F(hname,
                    f"Invariant Mass (L_{{xy}}=[{lxy_min},{lxy_max}] cm);M(#mu#mu) [GeV];Events",
                    100, 0, 10)
                histograms[hname].Sumw2()
    
    # Mother PDG ID histograms
    for flavor in flavor_types:
        hname = f"h_mom_pdgid_{flavor}_same"
        histograms[hname] = ROOT.TH1F(hname, 
                                       "Mother PDG ID;PDG ID;Events",
                                       600, 0, 600)
        histograms[hname].Sumw2()

    # lxy histo for photon conversion
    hname = "h_lxy_SV_light_same_photon"
    histograms[hname] = ROOT.TH1F(hname,"L_{xy} (Light, Mother=Photon);L_{xy} [cm];Events",50, 0, 10)
    histograms[hname].Sumw2()
    
    # Fill histograms
    for flavor in flavor_types:
        flavor_mask_per_jet = get_jet_flavor_mask(data["partFlav_jet"], flavor)
        flavor_mask = ak.any(flavor_mask_per_jet, axis=1)
        
        for mother in mother_types:
            if mother == 'same':
                mother_mask = ak.all(same_in_sv, axis=1)
            else:
                mother_mask = ~ak.all(same_in_sv, axis=1)
            
            mask = flavor_mask & mother_mask
            if ak.sum(mask) == 0:
                continue
            
            subset = data[mask]
            
            # Fill basic histograms
            for var, (_, _, _, _, short) in hist_config.items():
                hname = f"h_{short}_{flavor}_{mother}"
                values = flat(subset[var])
                for val in values:
                    histograms[hname].Fill(val, weight)
            
            # Fill with vtxProb cut
            vtx_cut_mask = ak.all(subset["vtxProb_SV"] > 0.2, axis=1)
            subset_vtxcut = subset[vtx_cut_mask]
            
            for var, (_, _, _, _, short) in hist_config.items():
                if var == "vtxProb_SV":
                    continue
                hname = f"h_{short}_{flavor}_{mother}_vtxcut"
                values = flat(subset_vtxcut[var])
                for val in values:
                    histograms[hname].Fill(val, weight)
            
            # Fill invm vs lxy
            for i, (lxy_min, lxy_max) in enumerate(lxy_bins):
                lxy_mask = (subset["lxy_SV"] >= lxy_min) & (subset["lxy_SV"] < lxy_max)
                lxy_mask_evt = ak.any(lxy_mask, axis=1)
                subset_lxy = subset[lxy_mask_evt]
                
                invm_vals = flat(subset_lxy["invm_ptl"])
                hname = f"h_invm_lxy{i}_{flavor}_{mother}"
                for val in invm_vals:
                    histograms[hname].Fill(val, weight)
            
            # Fill mother PDG ID
            if mother == 'same':
                mom_vals = flat(subset["mom_ptl"])
                hname = f"h_mom_pdgid_{flavor}_same"
                for val in mom_vals:
                    histograms[hname].Fill(abs(val), weight)

            if flavor == 'light' and mother == 'same':
                photon_mask = ak.all(abs(subset["mom_ptl"]) == 22, axis=2)
                photon_mask_evt = ak.all(photon_mask, axis=1)
                subset_photon = subset[photon_mask_evt]

                lxy_vals = flat(subset_photon["lxy_SV"])
                hname = "h_lxy_SV_light_same_photon"
                for val in lxy_vals:
                    histograms[hname].Fill(val, weight)

            if flavor == 'light' and mother == 'same':
                    photon_mask = ak.all(abs(subset["mom_ptl"]) == 22, axis=2)
                    photon_mask_evt = ak.all(photon_mask, axis=1)
                    subset_photon = subset[photon_mask_evt]

                    lxy_vals = flat(subset_photon["lxy_SV"])
                    hname = "h_lxy_SV_light_same_photon"
                    for val in lxy_vals:
                            histograms[hname].Fill(val, weight)
    
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
