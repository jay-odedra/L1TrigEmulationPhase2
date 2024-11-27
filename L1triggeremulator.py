import ROOT
import itertools
import array

class L1Electron:
    def __init__(self):
        self.electronindex = None
        self.charge = None
        self.hwIso = None
        self.hwQual = None
        self.hwZ0 = None
        self.eta = None
        self.iso = None
        self.phi = None
        self.pt = None
        self.relIso = None
        self.z0 = None

    def __str__(self):
        return (f"electronindex: {self.electronindex}, "
                f"charge: {self.charge}, "
                f"hwIso: {self.hwIso}, "
                f"hwQual: {self.hwQual}, "
                f"hwZ0: {self.hwZ0}, "
                f"eta: {self.eta}, "
                f"iso: {self.iso}, "
                f"phi: {self.phi}, "
                f"pt: {self.pt}, "
                f"relIso: {self.relIso}, "
                f"z0: {self.z0}")

    def set_attributes_from_tree(self, tree, index):
            self.electronindex = index
            self.charge = [x for x in tree.L1GTtkElectron_charge][index]
            self.hwIso = [x for x in tree.L1GTtkElectron_hwIso][index]
            self.hwQual = [x for x in tree.L1GTtkElectron_hwQual][index]
            self.hwZ0 = [x for x in tree.L1GTtkElectron_hwZ0][index]
            self.eta = [x for x in tree.L1GTtkElectron_eta][index]
            self.iso = [x for x in tree.L1GTtkElectron_iso][index]
            self.phi = [x for x in tree.L1GTtkElectron_phi][index]
            self.pt = [x for x in tree.L1GTtkElectron_pt][index]
            self.relIso = [x for x in tree.L1GTtkElectron_relIso][index]
            self.z0 = [x for x in tree.L1GTtkElectron_z0][index]

    
def check_charge(electrons):
    charges = [electron.charge for electron in electrons]
    return 1 in charges and -1 in charges

def create_branches(tree, prefix, attributes):
    branches = {}
    for attr in attributes:
        branches[attr] = array.array("f", [0])
        tree.Branch(f"{prefix}_{attr}", branches[attr], f"{prefix}_{attr}/F")
    return branches

def fill_branches(branches, electron):
    for attr in branches:
        branches[attr][0] = getattr(electron, attr)

def DeltaR(ele1, ele2):
    deta = ele1.eta - ele2.eta
    dphi = ROOT.TVector2.Phi_mpi_pi(ele1.phi - ele2.phi)
    return (deta**2 + dphi**2)**0.5




ele_pt_cut = 0
ele_eta_cut = 2500
ele12_DR_cut_max = 1000
ele12_DR_cut_min = 0
MLL_mass_cut = 6





















# Open the ROOT file
file_path = "/eos/user/j/jodedra/Phase2Analysiswork/L1TrigEmulationPhase2/output_Phase2_L1T_kee.root"
file = ROOT.TFile.Open(file_path)

# Get the tree from the file
tree_name = "Events"  # Replace with the actual name of your tree
tree = file.Get(tree_name)

# Create an output file and tree
output_file = ROOT.TFile("output_with_dielectrons_passing_trig_withnewcuts.root", "RECREATE")
output_tree = ROOT.TTree("Events", "Tree with dielectron pairs")


# Define branches for the output tree
dielectron_mass = array.array("f", [0])
dielectron_pt = array.array("f", [0])
dielectron_eta = array.array("f", [0])
dielectron_phi = array.array("f", [0])
output_tree.Branch("dielectron_mass", dielectron_mass, "dielectron_mass/F")
output_tree.Branch("dielectron_pt", dielectron_pt, "dielectron_pt/F")
output_tree.Branch("dielectron_eta", dielectron_eta, "dielectron_eta/F")
output_tree.Branch("dielectron_phi", dielectron_phi, "dielectron_phi/F")

# Define attributes for electrons
attributes = ["electronindex", "charge", "hwIso", "hwQual", "hwZ0", "eta", "iso", "phi", "pt", "relIso", "z0"]
# Create branches for the first and second electrons in the dielectron pair
dielectron_ele1_branches = create_branches(output_tree, "dielectron_ele1", attributes)
dielectron_ele2_branches = create_branches(output_tree, "dielectron_ele2", attributes)

# Loop over all entries in the tree
for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    l1electronobjects = []
    
    # Skip if there are not more than one electron
    if not tree.nL1GTtkElectron > 1: continue
    for electronindex in range(tree.nL1GTtkElectron):
        # Create L1Electron objects for each electron in the event
        l1electron = L1Electron()
        l1electron.set_attributes_from_tree(tree, electronindex)
        l1electronobjects.append(l1electron)
    
    # Check if there is at least one positive and one negative charge in the list
    if not check_charge(l1electronobjects): continue
    
    
    z0valuemin = 10000
    bestpairid =[None, None]
    dielectron = None
    # Loop over all combinations of electron pairs
    for ele1,ele2 in itertools.combinations(l1electronobjects,2):    
        if ele1.charge + ele2.charge != 0 or ele1.charge * ele2.charge == 0: continue
        
        
        
        tlv1 = ROOT.TLorentzVector()
        tlv1.SetPtEtaPhiM(ele1.pt, ele1.eta, ele1.phi, 0.000511)
            
        tlv2 = ROOT.TLorentzVector()
        tlv2.SetPtEtaPhiM(ele2.pt, ele2.eta, ele2.phi, 0.000511)
        mll_pair = tlv1 + tlv2
        # Skip if the invariant mass is not less than 6

        if not mll_pair.M() < MLL_mass_cut: continue
        # Find the best pair based on the smallest z0 difference
        if abs(ele1.z0 - ele2.z0) < z0valuemin and abs(ele1.z0 - ele2.z0) < 1.0:
            z0valuemin = abs(ele1.z0 - ele2.z0)
            bestpair = sorted([ele1, ele2], key=lambda x: x.pt, reverse=True)
            dielectron = mll_pair
        else:
            continue
    if not dielectron: continue
    
    # Skip if the DR between the two electrons is too big
    if DeltaR(bestpair[0], bestpair[1]) > ele12_DR_cut_max: continue
    
    # Skip if the DR between the two electrons is to small
    if DeltaR(bestpair[0], bestpair[1]) < ele12_DR_cut_min: continue
    
        # Skip if the charges do not sum to zero or if either charge is zero
    if bestpair[0].pt < ele_pt_cut or bestpair[1].pt < ele_pt_cut: continue
    
    # apply pt and eta cuts on the electrons
    if abs(bestpair[0].eta) > ele_eta_cut or abs(bestpair[1].eta) > ele_eta_cut: continue    
    
    # Fill the output tree with the dielectron mass and attributes of the best pair
    dielectron_mass[0] = dielectron.M()
    dielectron_pt[0] = dielectron.Pt()
    dielectron_eta[0] = dielectron.Eta()
    dielectron_phi[0] = dielectron.Phi()
    fill_branches(dielectron_ele1_branches, bestpair[0])
    fill_branches(dielectron_ele2_branches, bestpair[1])
    output_tree.Fill()
# Write and close the output file
output_file.Write()
output_file.Close()
print("Output file created successfully.")



