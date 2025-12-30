import json
import math
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
# Suppress RDKit warnings
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

def add_dicts(dict1, dict2):
    """
    Combines two dictionaries by adding the values of common keys.

    Returns
    -------
    dict
        A new dictionary containing the summed values. Keys with a final value of 0 are removed.
    """
    # Start with a copy of the first dictionary
    result_dict = dict1.copy()
    
    # Iterate over the second dictionary
    for key, value in dict2.items():
        # Add the value from the second dictionary to the result,
        # using .get(key, 0) to handle keys that are not in the first dictionary.
        result_dict[key] = result_dict.get(key, 0) + value

    # Create a new dictionary, excluding any keys where the summed value is 0
    return {key: value for key, value in result_dict.items() if value != 0}

def subtract_dicts(dict1, dict2):
    """
    Subtracts the values of one dictionary from another based on common keys.

    Returns
    -------
    dict
        A new dictionary containing the results of the subtraction. Keys with a final value of 0 are removed.
    """
    # Start with a copy of the first dictionary (the minuend)
    result_dict = dict1.copy()

    # Iterate over the second dictionary (the subtrahend)
    for key, value in dict2.items():
        # Subtract the value from the second dictionary from the result,
        # using .get(key, 0) to handle keys not present in the first.
        result_dict[key] = result_dict.get(key, 0) - value

    # Create a new dictionary, excluding any keys where the final value is 0
    return {key: value for key, value in result_dict.items() if value != 0}

def reaction_string_to_moiety_change_dict(reaction_string, radius):
    """
    Parses a reaction SMILES string and calculates the net change in moieties.

    Parameters
    ----------
    reaction_string : str
        A reaction SMILES string (e.g., "CCO.O>>CC(=O)O").
    radius : int or str
        The radius for substructure definition, or 'MAX' to treat entire molecules as moieties.

    Returns
    -------
    dict
        A dictionary representing the net change of moieties {moiety_smiles: count}.
        Positive counts are products, negative counts are reactants.
    """
    # Split the reaction string into reactant and product sides
    initial_state_str, final_state_str = reaction_string.split(">>")
    
    # Split each side into individual molecule SMILES
    initial_state_smiles = initial_state_str.split(".")
    final_state_smiles = final_state_str.split(".")
    
    # Check if we are using substructure moieties or whole-molecule moieties
    if radius != 'MAX':
        # --- Substructure (moiety) mode ---
        state_i = {} # Dictionary to hold reactant moiety counts
        for smile in initial_state_smiles:
            # Count substructures for each reactant molecule and aggregate them
            state_i = add_dicts(state_i, MolSmilesToMoietyDict(smile,radius))
            
        state_f = {} # Dictionary to hold product moiety counts
        for smile in final_state_smiles:
            # Count substructures for each product molecule and aggregate them
            state_f = add_dicts(state_f, MolSmilesToMoietyDict(smile,radius))
        
        # Calculate the net change by subtracting reactant counts from product counts
        return subtract_dicts(state_f, state_i)
    else:
        # --- Whole-molecule mode ---
        state_i = {} # Dictionary to hold reactant molecule counts
        for smile in initial_state_smiles:
            state_i = add_dicts(state_i, {smile: 1})

        state_f = {} # Dictionary to hold product molecule counts
        for smile in final_state_smiles:
            state_f = add_dicts(state_f, {smile: 1})
            
        # Calculate the net change by subtracting reactant counts from product counts
        return subtract_dicts(state_f, state_i)

def moiety_dict_to_reaction_smiles(d):
    """
    Constructs a reaction SMILES string from a net change dictionary.

    Parameters
    ----------
    d : dict
        A dictionary of {moiety_smiles: count}.

    Returns
    -------
    str
        A reaction SMILES string in the format "reactants>>products".
    """
    final = []    # List to hold product SMILES
    initial = []  # List to hold reactant SMILES
    
    # Iterate through the net change dictionary
    for key, value in d.items():
        if value > 0:
            # Positive values correspond to products; add them to the 'final' list
            final.extend([key] * abs(value))
        elif value < 0:
            # Negative values correspond to reactants; add them to the 'initial' list
            initial.extend([key] * abs(value))
            
    # Join the lists of SMILES with '.' to form the reactant and product strings
    Initial_state = '.'.join(initial)
    Final_state = '.'.join(final)
    
    # Combine them into a full reaction SMILES string
    return Initial_state + '>>' + Final_state

# creates moiety count vector for various radii
def MolSmilesToMoietyDict(smiles,radius):
    
    if radius != math.ceil(radius):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    else:
        mol = Chem.MolFromSmiles(smiles)

    
    # Initialize the dictionary to store counts of each substructure's SMILES
    smi_count = {}
    # Iterate over each atom in the molecule to treat it as a potential substructure center
    for i in range(len([atom for atom in mol.GetAtoms()])):
        if mol.GetAtomWithIdx(i).GetAtomicNum() != 1: # find moieties for all non-hydrogen atoms
            
            mod_mol = Chem.RWMol(mol)
    
            # For non integer radii info about the bond type is kept 
            # while masking all other features of the atom not in ceil(r)-1 bond radius away
            # The atomic number 0 corresponds to a wildcard atom ('*').
            # this clears its label, charge, and number of radical electrons 
            # Use a set to store the indices of atoms in the current search to make sure they are not masked 
            if radius != math.ceil(radius):
                
                # Iteratively find the largest possible chemical environment up to the given radius.
                # This loop starts at the max radius and decreases, finding the first valid environment.
                # this caps the largest moiety to be only the molecule itself
                for r in range(math.floor(radius), -1, -1):
                    env = Chem.FindAtomEnvironmentOfRadiusN(mod_mol, r, i,useHs=True)
                    if env:
                        break
    
                atoms_to_keep = set()
                atoms_to_keep.add(i)
                # ...collect all atom indices connected by those bonds
                for bidx in env:
                    atoms_to_keep.add(mod_mol.GetBondWithIdx(bidx).GetBeginAtomIdx())
                    atoms_to_keep.add(mod_mol.GetBondWithIdx(bidx).GetEndAtomIdx())
    
                for j in range(mod_mol.GetNumAtoms()):
                    if j not in atoms_to_keep :
                        mod_mol.GetAtomWithIdx(j).SetAtomicNum(0) # atom *
                        mod_mol.GetAtomWithIdx(j).SetAtomMapNum(0) # label
                        mod_mol.GetAtomWithIdx(j).SetFormalCharge(0) # charge
                        mod_mol.GetAtomWithIdx(j).SetNumRadicalElectrons(0) # number of radical electrons  
           
            for r in range(math.ceil(radius), -1, -1):
                env = Chem.FindAtomEnvironmentOfRadiusN(mod_mol, r, i,useHs=True)
                if env:
                    break # Stop once the largest possible environment is found
    
            # Use a set to store the indices of atoms in the current substructure to avoid duplicates
            atoms = set()
            atoms.add(i)
    
            # ...collect all atom indices connected by those bonds
            for bidx in env:
                atoms.add(mod_mol.GetBondWithIdx(bidx).GetBeginAtomIdx())
                atoms.add(mod_mol.GetBondWithIdx(bidx).GetEndAtomIdx())
    
            # Generate a canonical SMILES string for the identified substructure
            substructure = Chem.MolFragmentToSmiles(mod_mol, atomsToUse=list(atoms), bondsToUse=env, canonical=True, isomericSmiles=True)

            #mark a moiety from a labeled atom with an !
            label = mod_mol.GetAtomWithIdx(i).GetAtomMapNum()
            if label != 0:
                substructure = substructure + '!'
                
            # Increment the count for this substructure in the dictionary
            if substructure in smi_count:
                smi_count[substructure] += 1
            else:
                smi_count[substructure] = 1
                 
    return smi_count

# these are protonation reactions of various substrates in M-CSA 
solvent_rxns_names = {
'[OH2:1]>>[OH-:1].[H+]' : 'Water',
'[OH3+:1]>>[OH2:1].[H+]' : 'Oxonium',
'[NH4+].O>>[OH3+].N' : 'Ammonium',
'[CH3:2][CH2:2][c:2]1[cH:2][nH:2][cH:2][nH+:2]1.[OH2:1]>>[OH3+:1].[CH3:2][CH2:2][c:2]1[cH:2][nH:2][cH:2][n:2]1' : 'HisN1',
'[CH3:2][CH2:2][c:2]1[cH:2][nH+:2][cH:2][nH:2]1.[OH2:1]>>[OH3+:1].[CH3:2][CH2:2][c:2]1[cH:2][n:2][cH:2][nH:2]1' : 'HisN3',
'[CH3:3][CH2:3][CH2:3][C:3](=[O:3])[OH:3].[OH2:1]>>[OH3+:1].[CH3:3][CH2:3][CH2:3][C:3](=[O:3])[O-:3]' : 'Asp/Glu',
'[CH3:4][CH2:4][CH2:4][CH2:4][CH2:4][NH3+:4].[OH2:1]>>[OH3+:1].[CH3:4][CH2:4][CH2:4][CH2:4][CH2:4][NH2:4]' : 'Lys',
'[CH3:5][CH2:5][SH:5].[OH2:1]>>[OH3+:1].[CH3:5][CH2:5][S-:5]' : 'Cys',
'[CH3:6][CH2:6][c:6]1[cH:6][cH:6][c:6]([OH:6])[cH:6][cH:6]1.[OH2:1]>>[OH3+:1].[CH3:6][CH2:6][c:6]1[cH:6][cH:6][c:6]([O-:6])[cH:6][cH:6]1' : 'Tys',
'[CH3:7][CH2:7][OH:7].[OH2:1]>>[OH3+:1].[CH3:7][CH2:7][O-:7]' : 'Ser',
'[CH3:8][CH2:8][CH2:8][CH2:8][NH:8][C:8]([NH2:8])=[NH2+:8].[OH2:1]>>[OH3+:1].[CH3:8][CH2:8][CH2:8][CH2:8][NH:8][C:8](=[NH:8])[NH2:8]' : 'Arg',
'CC([NH3+])C(=O)[O-].[OH2:1]>>[OH3+:1].CC(N)C(=O)[O-]' : 'Amine_secondary',
'CC(N)C(=O)O.[OH2:1]>>[OH3+:1].CC(N)C(=O)[O-]' : 'Carboxylic',
'O=S(=O)([O-])S.[OH2:1]>>[OH3+:1].O=S(=O)([O-])[S-]' : 'Thiosulfate_sulfur',
'O=P(O)(O)O.[OH2:1]>>[OH3+:1].O=P([O-])(O)O' : 'Phosphate',
'O=P([O-])(O)O.[OH2:1]>>[OH3+:1].O=P([O-])([O-])O' : 'Phosphate_-1',
'O=P([O-])([O-])O.[OH2:1]>>[OH3+:1].O=P([O-])([O-])[O-]' : 'Phosphate_-2',
'O=P([O-])([O-])O.[OH2:1]>>[OH3+:1].O=[P-](=O)([O-])[O-]' : 'Phosphate_-3',
'[NH3+]CC(=O)COP(=O)([O-])[O-].[OH2:1]>>[OH3+:1].NCC(=O)COP(=O)([O-])[O-]' : 'amine_primary',
'Cc1[nH+]cc(COP(=O)([O-])[O-])c(CO)c1O.[OH2:1]>>[OH3+:1].Cc1ncc(COP(=O)([O-])[O-])c(CO)c1O' : 'PLP_N',
'C[n+]1cnc(N)c2nc[nH]c21.[OH2:1]>>[OH3+:1].Cn1cnc(N)c2ncnc1-2' : 'Adenine_methylN',
'CP(=O)([O-])OCC1OC(n2cnc3c(N)ncnc32)C(O)C1O.[OH2:1]>>[OH3+:1].CP(=O)([O-])OCC1OC(n2cnc3c(N)ncnc32)C(O)C1[O-]' : 'alcohol_secondary*',
'O=C([O-])O.[OH2:1]>>[OH3+:1].O=C([O-])[O-]' : 'Bicarbonate_-1',
'C[NH2+]CC(=O)[O-].[OH2:1]>>[OH3+:1].CNCC(=O)[O-]' : 'Methylamine',
'O=C([O-])C(=O)CP(=O)([O-])O.[OH2:1]>>[OH3+:1].O=C([O-])C(=O)CP(=O)([O-])[O-]' : '3-Phosphonopyruvate',
'O=NO.[OH2:1]>>[OH3+:1].O=N[O-]' : 'Nitrite',
'Nc1nc2c(c(=O)[nH]1)[NH+]=C(C(O)C(O)COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])CN2.[OH2:1]>>[OH3+:1].Nc1nc2c(c(=O)[nH]1)N=C(C(O)C(O)COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])CN2' : 'THFA',
'[NH3+]C1OC(COP(=O)([O-])[O-])C(O)C1O.[OH2:1]>>[OH3+:1].NC1OC(COP(=O)([O-])[O-])C(O)C1O' : 'AmineC1_sugar',
'Cl.[OH2:1]>>[OH3+:1].[Cl-]' : 'HCl',
'O=C(C=Cc1ccc(O)cc1)c1c(O)cc(O)cc1O.[OH2:1]>>[OH3+:1].O=C(C=Cc1ccc(O)cc1)c1c([O-])cc(O)cc1O' : 'Phenol',
'CC(C)(CO)C(O)C(=O)[NH2+]CCC(=O)[O-].[OH2:1]>>[OH3+:1].CC(C)(CO)C(O)C(=O)NCCC(=O)[O-]' : 'Methylamide',
'[NH3+]c1ncnc2[nH]cnc12.[OH2:1]>>[OH3+:1].Nc1ncnc2[nH]cnc12' : 'Adenine_primary',
'O=[As]([O-])(O)O.[OH2:1]>>[OH3+:1].O=[As]([O-])([O-])O' : 'Arsenate',
'C#N.[OH2:1]>>[OH3+:1].[C-]#N' : 'Cyanide',
'O=S(=O)(O)O.[OH2:1]>>[OH3+:1].O=S(=O)([O-])O' : 'Sulfate',
'O=S(=O)([O-])O.[OH2:1]>>[OH3+:1].O=S(=O)([O-])[O-]' : 'Sulfate_-1',
'Nc1c(NCC(O)C(O)C(O)CO)[nH]c(=O)[nH]c1=O.[OH2:1]>>[OH3+:1].Nc1c(NCC(O)C(O)C(O)CO)[nH]c(=O)[n-]c1=O' : 'Uracil',
'O=S([O-])O.[OH2:1]>>[OH3+:1].O=S([O-])[O-]' : 'Sulfite_-1',
'Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(CC(O)C(O)C(O)COP(=O)([O-])OP(=O)([O-])OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)c2cc1C.[OH2:1]>>[OH3+:1].Cc1cc2nc3c(=O)[n-]c(=O)nc-3n(CC(O)C(O)C(O)COP(=O)([O-])OP(=O)([O-])OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)c2cc1C' : 'FAD_N3',
'CC(=O)C(C)(O)C(=O)[O-].[OH2:1]>>[OH3+:1].CC(=O)C(C)([O-])C(=O)[O-]' : '2-Acetyllactate*',
'OCC1OC([OH+]CC2OC(O)(CO)C(O)C2O)C(O)C(O)C1O.[OH2:1]>>[OH3+:1].OCC1OC(OCC2OC(O)(CO)C(O)C2O)C(O)C(O)C1O' : 'EsterC1C5_sugar',
'O=C[NH2+]CC(=O)NC1OC(COP(=O)([O-])[O-])C(O)C1O.[OH2:1]>>[OH3+:1].O=CNCC(=O)NC1OC(COP(=O)([O-])[O-])C(O)C1O' : 'Formyl-amine',
'O=CO.[OH2:1]>>[OH3+:1].O=C[O-]' : 'Formate',
'CNC(CCCC[NH2+]CCCC[NH3+])C(C)=O.[OH2:1]>>[OH3+:1].CNC(CCCCNCCCC[NH3+])C(C)=O' : 'Diethylamine',
'CN(CC(=O)[O-])C(N)=[NH2+].[OH2:1]>>[OH3+:1].CN(CC(=O)[O-])C(=N)N' : 'Creatine_N',
'Cc1cc2c(cc1C)N(CC(O)C(O)C(O)COP(=O)([O-])[O-])c1[nH]c(=O)[nH]c(=O)c1N2.[OH2:1]>>[OH3+:1].Cc1cc2c(cc1C)N(CC(O)C(O)C(O)COP(=O)([O-])[O-])c1[n-]c(=O)[nH]c(=O)c1N2' : 'FAD_N1',
'C[NH3+].[OH2:1]>>[OH3+:1].CN' : 'Methylamine',
'O=C(O)c1cccnc1.[OH2:1]>>[OH3+:1].O=C([O-])c1cccnc1' : 'Benzoate',
'C[NH+](C)C.[OH2:1]>>[OH3+:1].CN(C)C' : 'Trimethylamine',
}

## mechanisms from the M-CSA database and 
with open('MCSA_rules_rMAX_1005.json', "r") as file:
    MCSA_rules = json.load(file)
with open('MCSA_mols_rMAX_1005.json', "r") as file:
    MCSA_mols = json.load(file)

labeled_entry_strings = []
for j in range(len(MCSA_rules)):
    labeled_mech_strings = []
    for i in range(len(MCSA_rules[j])):
        labeled_step_string = []
        for s in range(len(MCSA_rules[j][i])):
            labeled_step_string.append(moiety_dict_to_reaction_smiles(MCSA_rules[j][i][s]))
        labeled_mech_strings.append(labeled_step_string)
    labeled_entry_strings.append(labeled_mech_strings)

labeled_rxn_strings = []
for j in range(len(MCSA_mols)):
    labeled_rxn_strings.append(moiety_dict_to_reaction_smiles(MCSA_mols[j]))

for radius in [0,1,2]:
    
    full_name_solvent_rxn_dict = {}
    for rxn in solvent_rxns_names:
        reactants = rxn.split('>>')[0]
        products = rxn.split('>>')[1]
        protonation_rxn = products + '>>' + reactants
        full_name_solvent_rxn_dict[solvent_rxns_names[rxn] + '_depro'] = reaction_string_to_moiety_change_dict(rxn,radius)
        full_name_solvent_rxn_dict[solvent_rxns_names[rxn] + '_pro'] = reaction_string_to_moiety_change_dict(protonation_rxn,radius)
    
    Rules = pd.DataFrame()
    for j in range(len(labeled_entry_strings)):
        if labeled_entry_strings[j] != []:
            for i in range(len(labeled_entry_strings[j])):
                if labeled_entry_strings[j][i] != []:
                    Overall_Rxn = reaction_string_to_moiety_change_dict(labeled_rxn_strings[j],radius)
                    Overall_Mech = {}
                    step_dicts = []
                    for step in labeled_entry_strings[j][i]:
                        step_dict = reaction_string_to_moiety_change_dict(step,radius)
                        step_dicts.append(step_dict)
                        Overall_Mech = add_dicts(step_dict,Overall_Mech)
                    
                    if Overall_Rxn == Overall_Mech:
                        for s in range(len(step_dicts)):
                            if labeled_entry_strings[j][i][s] != '>>':
                                df = pd.DataFrame(step_dicts[s], index=[str(j+1) + "_" + str(i+1) + "_" + str(s+1)])
                                df = df.T
                                Rules = pd.concat([Rules, df], axis=1)
                                df = pd.DataFrame(subtract_dicts({},step_dicts[s]), index=["(" + str(j+1) + "_" + str(i+1) + "_" + str(s+1) + ")"])
                                df = df.T
                                Rules = pd.concat([Rules, df], axis=1)

    for name in full_name_solvent_rxn_dict:
        df = pd.DataFrame(full_name_solvent_rxn_dict[name], index=[name])
        df = df.T
        Rules = pd.concat([Rules, df], axis=1)
                   
    # Fill missing values with 0
    Rules= Rules.fillna(0).astype(int)
    
    
    # Dictionary to store unique columns
    unique_columns = {}
    original_to_unique = {}
    # Identify and store unique columns
    for col in Rules:
        col_tuple = tuple(Rules[col])
        if col_tuple not in unique_columns.values():
            unique_columns[col] = col_tuple
            original_to_unique[col] = col  # Map original to itself initially
        else:
            # Find the unique column name that matches the current column
            for unique_col, unique_tuple in unique_columns.items():
                if unique_tuple == col_tuple:
                    original_to_unique[col] = unique_col
                    break
    
    # Create a new DataFrame from the unique columns
    Unique_Rules = pd.DataFrame(unique_columns, index=Rules.index)
    
    # Create a dictionary to map unique column values to concatenated names
    unique_to_concat = {unique_col: [] for unique_col in unique_columns.keys()}
    
    # Populate the dictionary with original column names
    for original_col, unique_col in original_to_unique.items():
        unique_to_concat[unique_col].append(original_col)
    
    # Concatenate the names
    unique_to_concat = {k: '&'.join(v) for k, v in unique_to_concat.items()}
    
    # Rename the columns in the new DataFrame
    Unique_Rules.rename(columns=unique_to_concat, inplace=True)
    
    Unique_Rules.to_csv('Unique_Rules/Unique_Rules_'+str(int(radius))+'.csv', index=True, header=True)
    
