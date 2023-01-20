import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import numpy as np

def get_mol_descriptors(smiles):
    s = pd.Series()
    try:
        molecule = Chem.MolFromSmiles(smiles)
        descriptors_dic = {}
        descriptors_dic['smiles'] = smiles
        descriptors_dic['molecular_weight'] = Descriptors.ExactMolWt(molecule)
        descriptors_dic['logp'] = Descriptors.MolLogP(molecule)
        descriptors_dic['h_bond_donor'] = Descriptors.NumHDonors(molecule)
        descriptors_dic['h_bond_acceptors'] = Descriptors.NumHAcceptors(molecule)
        descriptors_dic['rotatable_bonds'] = Descriptors.NumRotatableBonds(molecule)
        descriptors_dic['number_of_atoms'] = Chem.rdchem.Mol.GetNumAtoms(molecule)
        descriptors_dic['molar_refractivity'] = Chem.Crippen.MolMR(molecule)
        descriptors_dic['topological_surface_area_mapping'] = Chem.QED.properties(molecule).PSA
        descriptors_dic['formal_charge'] = Chem.rdmolops.GetFormalCharge(molecule)
        descriptors_dic['heavy_atoms'] = Chem.rdchem.Mol.GetNumHeavyAtoms(molecule)
        descriptors_dic['num_of_rings'] = Chem.rdMolDescriptors.CalcNumRings(molecule)
#         descriptors_dic['MolWt']=Descriptors.MolWt(molecule)
#         descriptors_dic['HeavyAtomMolWt']=Descriptors.HeavyAtomMolWt(molecule)
# #         descriptors_dic['ExactMolWt']=Descriptors.ExactMolWt(molecule)
#         descriptors_dic['NumValenceElectrons']=Descriptors.NumValenceElectrons(molecule)
#         descriptors_dic['NumRadicalElectrons']=Descriptors.NumRadicalElectrons(molecule)
#         descriptors_dic['MaxEStateIndex']=Descriptors.MaxEStateIndex(molecule)
#         descriptors_dic['MinEStateIndex']=Descriptors.MinEStateIndex(molecule)
#         descriptors_dic['MaxAbsEStateIndex']=Descriptors.MaxAbsEStateIndex(molecule)
#         descriptors_dic['MinAbsEStateIndex']=Descriptors.MinAbsEStateIndex(molecule)
#         descriptors_dic['BalabanJ']=Descriptors.BalabanJ(molecule)
#         descriptors_dic['BertzCT']=Descriptors.BertzCT(molecule)
#         descriptors_dic['Chi0']=Descriptors.Chi0(molecule)
#         descriptors_dic['Chi0n']=Descriptors.Chi0n(molecule)
#         descriptors_dic['Chi0v']=Descriptors.Chi0v(molecule)
#         descriptors_dic['Chi1']=Descriptors.Chi1(molecule)
#         descriptors_dic['Chi1n']=Descriptors.Chi1n(molecule)
#         descriptors_dic['Chi1v']=Descriptors.Chi1v(molecule)
#         descriptors_dic['Chi2n']=Descriptors.Chi2n(molecule)
#         descriptors_dic['Chi2v']=Descriptors.Chi2v(molecule)
#         descriptors_dic['Chi3n']=Descriptors.Chi3n(molecule)
#         descriptors_dic['Chi3v']=Descriptors.Chi3v(molecule)
#         descriptors_dic['Chi4n']=Descriptors.Chi4n(molecule)
#         descriptors_dic['Chi4v']=Descriptors.Chi4v(molecule)
#         descriptors_dic['EState_VSA1']=Descriptors.EState_VSA1(molecule)
#         descriptors_dic['EState_VSA10']=Descriptors.EState_VSA10(molecule)
#         descriptors_dic['EState_VSA11']=Descriptors.EState_VSA11(molecule)
#         descriptors_dic['EState_VSA2']=Descriptors.EState_VSA2(molecule)
#         descriptors_dic['EState_VSA3']=Descriptors.EState_VSA3(molecule)
#         descriptors_dic['EState_VSA4']=Descriptors.EState_VSA4(molecule)
#         descriptors_dic['EState_VSA5']=Descriptors.EState_VSA5(molecule)
#         descriptors_dic['EState_VSA6']=Descriptors.EState_VSA6(molecule)
#         descriptors_dic['EState_VSA7']=Descriptors.EState_VSA7(molecule)
#         descriptors_dic['EState_VSA8']=Descriptors.EState_VSA8(molecule)
#         descriptors_dic['EState_VSA9']=Descriptors.EState_VSA9(molecule)
#         descriptors_dic['FractionCSP3']=Descriptors.FractionCSP3(molecule)
#         descriptors_dic['HallKierAlpha']=Descriptors.HallKierAlpha(molecule)
#         descriptors_dic['HeavyAtomCount']=Descriptors.HeavyAtomCount(molecule)
#         descriptors_dic['Ipc']=Descriptors.Ipc(molecule)
#         descriptors_dic['Kappa1']=Descriptors.Kappa1(molecule)
#         descriptors_dic['Kappa2']=Descriptors.Kappa2(molecule)
#         descriptors_dic['Kappa3']=Descriptors.Kappa3(molecule)
#         descriptors_dic['LabuteASA']=Descriptors.LabuteASA(molecule)
# #         descriptors_dic['MolLogP']=Descriptors.MolLogP(molecule)
#         descriptors_dic['MolMR']=Descriptors.MolMR(molecule)
#         descriptors_dic['NHOHCount']=Descriptors.NHOHCount(molecule)
#         descriptors_dic['NOCount']=Descriptors.NOCount(molecule)
#         descriptors_dic['NumAliphaticCarbocycles']=Descriptors.NumAliphaticCarbocycles(molecule)
#         descriptors_dic['NumAliphaticHeterocycles']=Descriptors.NumAliphaticHeterocycles(molecule)
#         descriptors_dic['NumAliphaticRings']=Descriptors.NumAliphaticRings(molecule)
#         descriptors_dic['NumAromaticCarbocycles']=Descriptors.NumAromaticCarbocycles(molecule)
#         descriptors_dic['NumAromaticHeterocycles']=Descriptors.NumAromaticHeterocycles(molecule)
#         descriptors_dic['NumAromaticRings']=Descriptors.NumAromaticRings(molecule)
# #         descriptors_dic['NumHAcceptors']=Descriptors.NumHAcceptors(molecule)
# #         descriptors_dic['NumHDonors']=Descriptors.NumHDonors(molecule)
#         descriptors_dic['NumHeteroatoms']=Descriptors.NumHeteroatoms(molecule)
# #         descriptors_dic['NumRotatableBonds']=Descriptors.NumRotatableBonds(molecule)
#         descriptors_dic['NumSaturatedCarbocycles']=Descriptors.NumSaturatedCarbocycles(molecule)
#         descriptors_dic['NumSaturatedHeterocycles']=Descriptors.NumSaturatedHeterocycles(molecule)
#         descriptors_dic['NumSaturatedRings']=Descriptors.NumSaturatedRings(molecule)
#         descriptors_dic['PEOE_VSA1']=Descriptors.PEOE_VSA1(molecule)
#         descriptors_dic['PEOE_VSA10']=Descriptors.PEOE_VSA10(molecule)
#         descriptors_dic['PEOE_VSA11']=Descriptors.PEOE_VSA11(molecule)
#         descriptors_dic['PEOE_VSA12']=Descriptors.PEOE_VSA12(molecule)
#         descriptors_dic['PEOE_VSA13']=Descriptors.PEOE_VSA13(molecule)
#         descriptors_dic['PEOE_VSA14']=Descriptors.PEOE_VSA14(molecule)
#         descriptors_dic['PEOE_VSA2']=Descriptors.PEOE_VSA2(molecule)
#         descriptors_dic['PEOE_VSA3']=Descriptors.PEOE_VSA3(molecule)
#         descriptors_dic['PEOE_VSA4']=Descriptors.PEOE_VSA4(molecule)
#         descriptors_dic['PEOE_VSA5']=Descriptors.PEOE_VSA5(molecule)
#         descriptors_dic['PEOE_VSA6']=Descriptors.PEOE_VSA6(molecule)
#         descriptors_dic['PEOE_VSA7']=Descriptors.PEOE_VSA7(molecule)
#         descriptors_dic['PEOE_VSA8']=Descriptors.PEOE_VSA8(molecule)
#         descriptors_dic['PEOE_VSA9']=Descriptors.PEOE_VSA9(molecule)
#         descriptors_dic['RingCount']=Descriptors.RingCount(molecule)
#         descriptors_dic['SMR_VSA1']=Descriptors.SMR_VSA1(molecule)
#         descriptors_dic['SMR_VSA10']=Descriptors.SMR_VSA10(molecule)
#         descriptors_dic['SMR_VSA2']=Descriptors.SMR_VSA2(molecule)
#         descriptors_dic['SMR_VSA3']=Descriptors.SMR_VSA3(molecule)
#         descriptors_dic['SMR_VSA4']=Descriptors.SMR_VSA4(molecule)
#         descriptors_dic['SMR_VSA5']=Descriptors.SMR_VSA5(molecule)
#         descriptors_dic['SMR_VSA6']=Descriptors.SMR_VSA6(molecule)
#         descriptors_dic['SMR_VSA7']=Descriptors.SMR_VSA7(molecule)
#         descriptors_dic['SMR_VSA8']=Descriptors.SMR_VSA8(molecule)
#         descriptors_dic['SMR_VSA9']=Descriptors.SMR_VSA9(molecule)
#         descriptors_dic['SlogP_VSA1']=Descriptors.SlogP_VSA1(molecule)
#         descriptors_dic['SlogP_VSA10']=Descriptors.SlogP_VSA10(molecule)
#         descriptors_dic['SlogP_VSA11']=Descriptors.SlogP_VSA11(molecule)
#         descriptors_dic['SlogP_VSA12']=Descriptors.SlogP_VSA12(molecule)
#         descriptors_dic['SlogP_VSA2']=Descriptors.SlogP_VSA2(molecule)
#         descriptors_dic['SlogP_VSA3']=Descriptors.SlogP_VSA3(molecule)
#         descriptors_dic['SlogP_VSA4']=Descriptors.SlogP_VSA4(molecule)
#         descriptors_dic['SlogP_VSA5']=Descriptors.SlogP_VSA5(molecule)
#         descriptors_dic['SlogP_VSA6']=Descriptors.SlogP_VSA6(molecule)
#         descriptors_dic['SlogP_VSA7']=Descriptors.SlogP_VSA7(molecule)
#         descriptors_dic['SlogP_VSA8']=Descriptors.SlogP_VSA8(molecule)
#         descriptors_dic['SlogP_VSA9']=Descriptors.SlogP_VSA9(molecule)
#         descriptors_dic['TPSA']=Descriptors.TPSA(molecule)
#         descriptors_dic['VSA_EState1']=Descriptors.VSA_EState1(molecule)
#         descriptors_dic['VSA_EState10']=Descriptors.VSA_EState10(molecule)
#         descriptors_dic['VSA_EState2']=Descriptors.VSA_EState2(molecule)
#         descriptors_dic['VSA_EState3']=Descriptors.VSA_EState3(molecule)
#         descriptors_dic['VSA_EState4']=Descriptors.VSA_EState4(molecule)
#         descriptors_dic['VSA_EState5']=Descriptors.VSA_EState5(molecule)
#         descriptors_dic['VSA_EState6']=Descriptors.VSA_EState6(molecule)
#         descriptors_dic['VSA_EState7']=Descriptors.VSA_EState7(molecule)
#         descriptors_dic['VSA_EState8']=Descriptors.VSA_EState8(molecule)
#         descriptors_dic['VSA_EState9']=Descriptors.VSA_EState9(molecule)
#         descriptors_dic['fr_Al_COO']=Descriptors.fr_Al_COO(molecule)
#         descriptors_dic['fr_Al_OH']=Descriptors.fr_Al_OH(molecule)
#         descriptors_dic['fr_Al_OH_noTert']=Descriptors.fr_Al_OH_noTert(molecule)
#         descriptors_dic['fr_ArN']=Descriptors.fr_ArN(molecule)
#         descriptors_dic['fr_Ar_COO']=Descriptors.fr_Ar_COO(molecule)
#         descriptors_dic['fr_Ar_N']=Descriptors.fr_Ar_N(molecule)
#         descriptors_dic['fr_Ar_NH']=Descriptors.fr_Ar_NH(molecule)
#         descriptors_dic['fr_Ar_OH']=Descriptors.fr_Ar_OH(molecule)
#         descriptors_dic['fr_COO']=Descriptors.fr_COO(molecule)
#         descriptors_dic['fr_COO2']=Descriptors.fr_COO2(molecule)
#         descriptors_dic['fr_C_O']=Descriptors.fr_C_O(molecule)
#         descriptors_dic['fr_C_O_noCOO']=Descriptors.fr_C_O_noCOO(molecule)
#         descriptors_dic['fr_C_S']=Descriptors.fr_C_S(molecule)
#         descriptors_dic['fr_HOCCN']=Descriptors.fr_HOCCN(molecule)
#         descriptors_dic['fr_Imine']=Descriptors.fr_Imine(molecule)
#         descriptors_dic['fr_NH0']=Descriptors.fr_NH0(molecule)
#         descriptors_dic['fr_NH1']=Descriptors.fr_NH1(molecule)
#         descriptors_dic['fr_NH2']=Descriptors.fr_NH2(molecule)
#         descriptors_dic['fr_N_O']=Descriptors.fr_N_O(molecule)
#         descriptors_dic['fr_Ndealkylation1']=Descriptors.fr_Ndealkylation1(molecule)
#         descriptors_dic['fr_Ndealkylation2']=Descriptors.fr_Ndealkylation2(molecule)
#         descriptors_dic['fr_Nhpyrrole']=Descriptors.fr_Nhpyrrole(molecule)
#         descriptors_dic['fr_SH']=Descriptors.fr_SH(molecule)
#         descriptors_dic['fr_aldehyde']=Descriptors.fr_aldehyde(molecule)
#         descriptors_dic['fr_alkyl_carbamate']=Descriptors.fr_alkyl_carbamate(molecule)
#         descriptors_dic['fr_alkyl_halide']=Descriptors.fr_alkyl_halide(molecule)
#         descriptors_dic['fr_allylic_oxid']=Descriptors.fr_allylic_oxid(molecule)
#         descriptors_dic['fr_amide']=Descriptors.fr_amide(molecule)
#         descriptors_dic['fr_amidine']=Descriptors.fr_amidine(molecule)
#         descriptors_dic['fr_aniline']=Descriptors.fr_aniline(molecule)
#         descriptors_dic['fr_aryl_methyl']=Descriptors.fr_aryl_methyl(molecule)
#         descriptors_dic['fr_azide']=Descriptors.fr_azide(molecule)
#         descriptors_dic['fr_azo']=Descriptors.fr_azo(molecule)
#         descriptors_dic['fr_barbitur']=Descriptors.fr_barbitur(molecule)
#         descriptors_dic['fr_benzene']=Descriptors.fr_benzene(molecule)
#         descriptors_dic['fr_benzodiazepine']=Descriptors.fr_benzodiazepine(molecule)
#         descriptors_dic['fr_bicyclic']=Descriptors.fr_bicyclic(molecule)
#         descriptors_dic['fr_diazo']=Descriptors.fr_diazo(molecule)
#         descriptors_dic['fr_dihydropyridine']=Descriptors.fr_dihydropyridine(molecule)
#         descriptors_dic['fr_epoxide']=Descriptors.fr_epoxide(molecule)
#         descriptors_dic['fr_ester']=Descriptors.fr_ester(molecule)
#         descriptors_dic['fr_ether']=Descriptors.fr_ether(molecule)
#         descriptors_dic['fr_furan']=Descriptors.fr_furan(molecule)
#         descriptors_dic['fr_guanido']=Descriptors.fr_guanido(molecule)
#         descriptors_dic['fr_halogen']=Descriptors.fr_halogen(molecule)
#         descriptors_dic['fr_hdrzine']=Descriptors.fr_hdrzine(molecule)
#         descriptors_dic['fr_hdrzone']=Descriptors.fr_hdrzone(molecule)
#         descriptors_dic['fr_imidazole']=Descriptors.fr_imidazole(molecule)
#         descriptors_dic['fr_imide']=Descriptors.fr_imide(molecule)
#         descriptors_dic['fr_isocyan']=Descriptors.fr_isocyan(molecule)
#         descriptors_dic['fr_isothiocyan']=Descriptors.fr_isothiocyan(molecule)
#         descriptors_dic['fr_ketone']=Descriptors.fr_ketone(molecule)
#         descriptors_dic['fr_ketone_Topliss']=Descriptors.fr_ketone_Topliss(molecule)
#         descriptors_dic['fr_lactam']=Descriptors.fr_lactam(molecule)
#         descriptors_dic['fr_lactone']=Descriptors.fr_lactone(molecule)
#         descriptors_dic['fr_methoxy']=Descriptors.fr_methoxy(molecule)
#         descriptors_dic['fr_morpholine']=Descriptors.fr_morpholine(molecule)
#         descriptors_dic['fr_nitrile']=Descriptors.fr_nitrile(molecule)
#         descriptors_dic['fr_nitro']=Descriptors.fr_nitro(molecule)
#         descriptors_dic['fr_nitro_arom']=Descriptors.fr_nitro_arom(molecule)
#         descriptors_dic['fr_nitro_arom_nonortho']=Descriptors.fr_nitro_arom_nonortho(molecule)
#         descriptors_dic['fr_nitroso']=Descriptors.fr_nitroso(molecule)
#         descriptors_dic['fr_oxazole']=Descriptors.fr_oxazole(molecule)
#         descriptors_dic['fr_oxime']=Descriptors.fr_oxime(molecule)
#         descriptors_dic['fr_para_hydroxylation']=Descriptors.fr_para_hydroxylation(molecule)
#         descriptors_dic['fr_phenol']=Descriptors.fr_phenol(molecule)
#         descriptors_dic['fr_phenol_noOrthoHbond']=Descriptors.fr_phenol_noOrthoHbond(molecule)
#         descriptors_dic['fr_phos_acid']=Descriptors.fr_phos_acid(molecule)
#         descriptors_dic['fr_phos_ester']=Descriptors.fr_phos_ester(molecule)
#         descriptors_dic['fr_piperdine']=Descriptors.fr_piperdine(molecule)
#         descriptors_dic['fr_piperzine']=Descriptors.fr_piperzine(molecule)
#         descriptors_dic['fr_priamide']=Descriptors.fr_priamide(molecule)
#         descriptors_dic['fr_prisulfonamd']=Descriptors.fr_prisulfonamd(molecule)
#         descriptors_dic['fr_pyridine']=Descriptors.fr_pyridine(molecule)
#         descriptors_dic['fr_quatN']=Descriptors.fr_quatN(molecule)
#         descriptors_dic['fr_sulfide']=Descriptors.fr_sulfide(molecule)
#         descriptors_dic['fr_sulfonamd']=Descriptors.fr_sulfonamd(molecule)
#         descriptors_dic['fr_sulfone']=Descriptors.fr_sulfone(molecule)
#         descriptors_dic['fr_term_acetylene']=Descriptors.fr_term_acetylene(molecule)
#         descriptors_dic['fr_tetrazole']=Descriptors.fr_tetrazole(molecule)
#         descriptors_dic['fr_thiazole']=Descriptors.fr_thiazole(molecule)
#         descriptors_dic['fr_thiocyan']=Descriptors.fr_thiocyan(molecule)
#         descriptors_dic['fr_thiophene']=Descriptors.fr_thiophene(molecule)
#         descriptors_dic['fr_unbrch_alkane']=Descriptors.fr_unbrch_alkane(molecule)
#         descriptors_dic['fr_urea']=Descriptors.fr_urea(molecule)
        s = pd.Series(descriptors_dic)
    except Exception as ex:
        print(ex)
        print(smiles)
    return s



def get_lipinski(row):
    lipinski = False
    molecular_weight = row['molecular_weight']
    h_bond_donor = row['h_bond_donor']
    h_bond_acceptors = row['h_bond_acceptors']
    logp = row['logp']
    rotatable_bonds = row['rotatable_bonds']
    if molecular_weight <= 500 and logp <= 5 and h_bond_donor <= 5 and h_bond_acceptors <= 5 and rotatable_bonds <= 5:
        lipinski = True
    return lipinski

def get_ghose_filter(row):
    ghose_filter = False
    molecular_weight = row['molecular_weight']
    logp = row['logp']
    molar_refractivity = row['molar_refractivity']
    number_of_atoms = row['number_of_atoms']
    if molecular_weight >= 160 and molecular_weight <= 480 and logp >= 0.4 and logp <= 5.6 and number_of_atoms >= 20 and number_of_atoms <= 70 and molar_refractivity >= 40 and molar_refractivity <= 130:
        ghose_filter = True
    return ghose_filter


def get_veber_filter(row):
    veber_filter = False
    rotatable_bonds = row['rotatable_bonds']
    topological_surface_area_mapping = row['topological_surface_area_mapping']
    if rotatable_bonds <= 10 and topological_surface_area_mapping <= 140:
        veber_filter = True
    return veber_filter


def get_ro3(row):
    rule_of_3 = False
    molecular_weight = row['molecular_weight']
    logp = row['logp']
    h_bond_donor = row['h_bond_donor']
    h_bond_acceptors = row['h_bond_acceptors']
    rotatable_bonds = row['rotatable_bonds']
    if molecular_weight <= 300 and logp <= 3 and h_bond_donor <= 3 and h_bond_acceptors <= 3 and rotatable_bonds <= 3:
        rule_of_3 = True
    return rule_of_3

def get_REOS_filter(row):
    reos_filter = False
    molecular_weight = row['molecular_weight']
    logp = row['logp']
    h_bond_donor = row['h_bond_donor']
    h_bond_acceptors = row['h_bond_acceptors']
    formal_charge = row['formal_charge']
    rotatable_bonds = row['rotatable_bonds']
    heavy_atoms = row['heavy_atoms']
    if molecular_weight >= 200 and molecular_weight <= 500 and logp >= int(0 - 5) and logp <= 5 and h_bond_donor >= 0 and h_bond_donor <= 5 and h_bond_acceptors >= 0 and h_bond_acceptors <= 10 and formal_charge >= int(0-2) and formal_charge <= 2 and rotatable_bonds >= 0 and rotatable_bonds <= 8 and heavy_atoms >= 15 and heavy_atoms <= 50:
        reos_filter = True
    return reos_filter


def get_druglike_filter(row):
    drug_like_filter = False
    molecular_weight = row['molecular_weight']
    num_of_rings = row['num_of_rings']
    rotatable_bonds = row['rotatable_bonds']
    h_bond_donor = row['h_bond_donor']
    h_bond_acceptors = row['h_bond_acceptors']
    logp = row['logp']
    if molecular_weight < 400 and num_of_rings > 0 and rotatable_bonds < 5 and h_bond_donor <= 5 and h_bond_acceptors <= 10 and logp < 5:
        drug_like_filter = True
    return drug_like_filter
        

def get_num_filter(row):
    count = 0
    lipinski = row['lipinski']
    ghose_filter = row['ghose_filter']
    veber_filter = row['veber_filter']
    rule_of_3 = row['rule_of_3']
    reos_filter = row['reos_filter']
    drug_like_filter = row['drug_like_filter']
    
    if lipinski:
        count+=1
    if ghose_filter:
        count+=1
    if veber_filter:
        count+=1
    if rule_of_3:
        count+=1
    if reos_filter:
        count+=1
    if drug_like_filter:
        count+=1
    return count

def get_morgan_fp(smiles, radius, vector_size):
    molecule = Chem.MolFromSmiles(smiles)
    return np.array(AllChem.GetMorganFingerprintAsBitVect(molecule, radius, nBits=vector_size, useFeatures=True)).reshape(-1,1)

def get_hashed_atom_pair_fp(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    return np.array(AllChem.GetHashedAtomPairFingerprintAsBitVect(molecule, nBits=256)).reshape(-1,1)

def get_hashed_topological_torision_fp(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    return np.array(AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(molecule , nBits=256)).reshape(-1,1)