import json
import pandas as pd
import time
import math
import pulp

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


with open('M-CSA_arrow_rules_r0.json', 'r') as f:
    MCSA_arrow_rules = json.load(f)

MCSA_mechanism_arrows = {}
for j in range(len(MCSA_arrow_rules)):
    for i in range(len(MCSA_arrow_rules[j])):
        MCSA_mechanism_arrows[str(j+1)+'_'+str(i+1)] = set()
        for s in range(len(MCSA_arrow_rules[j][i])):
            for dir in range(len(MCSA_arrow_rules[j][i][s])):
                for a in range(len(MCSA_arrow_rules[j][i][s][dir])):
                    arr = MCSA_arrow_rules[j][i][s][dir][a]
                    MCSA_mechanism_arrows[str(j+1)+'_'+str(i+1)].add(arr)    


# A hardcoded dictionary mapping common protonation/deprotonation rule names
# to a set of generic "arrow environment" reaction SMILES. This is used
# in the similarity scoring to represent these common steps.
Protonation_arrow_rules = {
     'Phenol_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},'Phenol_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'Asp/Glu_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},'Asp/Glu_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'HisN3_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},'HisN3_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},
     'Lys_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},'Lys_depro': {'[OH]>>[O]', '[O]>>[H].[O]' '[NH+]>>[N+]', '[N]>>[H].[N]'},
     'Arg_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},'Arg_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},
     'HisN1_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},'HisN1_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},
     'THFA_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},'THFA_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},
     'Tys_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'Formate_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},'Formate_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'Ser_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'Cys_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[S-]>>[H].[S-]', '[SH]>>[S]'},
     'alcohol_secondary_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     '2-Acetyllactate_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'Amine_secondary_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},'Amine_secondary_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},
     'Water_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},'Water_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'Oxonium_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O]>>[H].[O]', '[OH+]>>[O+]'},'Oxonium_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O]>>[H].[O]', '[OH+]>>[O+]'},
     'Ammonium_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[N]>>[H].[N]', '[NH+]>>[N+]'},'Ammonium_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[N]>>[H].[N]', '[NH+]>>[N+]'},
     'Carboxylic_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},'Carboxylic_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'Thiosulfate_sulfur_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[S-]>>[H].[S-]', '[SH]>>[S]'},'Thiosulfate_sulfur_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[S-]>>[H].[S-]', '[SH]>>[S]'},
     'Phosphate_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},'Phosphate_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'Phosphate_-1_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},'Phosphate_-1_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'Phosphate_-2_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},'Phosphate_-2_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'Phosphate_-3_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},'Phosphate_-3_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'amine_primary_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},'amine_primary_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},
     'PLP_N_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},'PLP_N_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},
     'Adenine_methylN_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH]>>[N]', '[N]>>[H].[N]', '[N]=[C]>>[N][C]'},'Adenine_methylN_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH]>>[N]', '[N]>>[H].[N]', '[N]=[C]>>[N][C]'},
     'Bicarbonate_-1_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},'Bicarbonate_-1_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'Methylamine_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},'Methylamine_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},
     '3-Phosphonopyruvate_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},'3-Phosphonopyruvate_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'Nitrite_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},'Nitrite_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'AmineC1_sugar_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},'AmineC1_sugar_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},
     'HCl_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[Cl-]>>[H].[Cl-]', '[ClH]>>[Cl]'},'HCl_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[Cl-]>>[H].[Cl-]', '[ClH]>>[Cl]'},
     'Methylamide_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},'Methylamide_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},
     'Adenine_primary_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},'Adenine_primary_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},
     'Arsenate_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},'Arsenate_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'Cyanide_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[C-]>>[H].[C-]', '[CH]>>[C]'},'Cyanide_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[C-]>>[H].[C-]', '[CH]>>[C]'},
     'Sulfate_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},'Sulfate_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'Sulfate_-1_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},'Sulfate_-1_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'Uracil_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[N-]>>[H].[N-]', '[NH]>>[N]'},'Uracil_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[N-]>>[H].[N-]', '[NH]>>[N]'},
     'Sulfite_-1_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},'Sulfite_-1_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'FAD_N3_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[N-]>>[H].[N-]', '[NH]>>[N]'},'FAD_N3_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[N-]>>[H].[N-]', '[NH]>>[N]'},
     'EsterC1C5_sugar_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O]>>[H].[O]', '[OH+]>>[O+]'},'EsterC1C5_sugar_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O]>>[H].[O]', '[OH+]>>[O+]'},
     'Formyl-amine_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},'Formyl-amine_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},
     'Diethylamine_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},'Diethylamine_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},
     'Creatine_N_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},'Creatine_N_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},
     'FAD_N1_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[N-]>>[H].[N-]', '[NH]>>[N]'},'FAD_N1_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[N-]>>[H].[N-]', '[NH]>>[N]'},
     'Benzoate_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},'Benzoate_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[O-]>>[H].[O-]', '[OH]>>[O]'},
     'Trimethylamine_depro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'},'Trimethylamine_pro': {'[OH]>>[O]', '[O]>>[H].[O]', '[NH+]>>[N+]', '[N]>>[H].[N]'}}

def similar_mechs(solution):
    """
    Calculates the similarity of a candidate mechanism (`solution`) to a global
    database of known M-CSA mechanisms.
    NOTE: Depends on global variables `MCSA_arrow_rules` and `MCSA_mechanism_arrows`.

    Returns
    -------
    dict
        A dictionary of {mcsa_id: similarity_score}, sorted by score descending.
    """
    # Use a set to store the unique arrow environments for the candidate solution
    solution_arrows = set()
    for rule in solution:
        steps = rule.split('&')
        for step in steps:
            # Handle special protonation rules by looking them up in the dictionary
            if 'pro' in step or 'depro' in step:
                for a in Protonation_arrow_rules.get(step, []):
                    solution_arrows.add(a)
            # Handle standard rule notation (e.g., 'mcsaID_step_substep')
            else:
                try:
                    # Parse the rule ID components
                    if '(' in step: # Handle alternative formatting
                        j, i, s = step.split('(')[1].split(')')[0].split('_')
                    else:
                        j, i, s = step.split('_')
                    # Look up the arrow environments for this rule in the global M-CSA data
                    for arr in MCSA_arrow_rules[int(j)-1][int(i)-1][int(s)-1]:
                        solution_arrows.update(arr)
                except (ValueError, IndexError):
                    # If parsing or lookup fails, just skip this step
                    continue
    
    solution_scores = {}
    # Iterate over every known mechanism in the global M-CSA arrow database
    for m in MCSA_mechanism_arrows:
        if MCSA_mechanism_arrows[m]: # Ensure the known mechanism is not empty
            # The similarity score is the Jaccard index: |A ) B| / |A * B|
            # This is calculated here using a mathematical identity:
            # |A ) B| = |A| + |B| - |A * B|
            # Here, added_lists represents |A| + |B|
            # And set(added_lists) represents |A * B|
            known_arrows_set = MCSA_mechanism_arrows[m]
            intersection_size = len(solution_arrows.intersection(known_arrows_set))
            union_size = len(solution_arrows.union(known_arrows_set))
            
            # Calculate the score, avoiding division by zero
            score = intersection_size / union_size if union_size > 0 else 0.0
            solution_scores[m] = score
    
    # Return the scores sorted from highest to lowest similarity
    return dict(sorted(solution_scores.items(), key=lambda item: item[1], reverse=True))


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

def MechFind(desired_reaction,radius,max_steps,iterations,time_limit):
    start_time = time.time()
    Unique_Rules= pd.read_csv('Unique_Rules/Unique_Rules_'+str(radius)+'.csv', index_col=0)
    
    old_names = Unique_Rules.columns.tolist()
    new_names = list(range(len(old_names)))
    column_mapping = dict(zip(old_names, new_names))
    Unique_Rules_renamed = Unique_Rules.rename(columns=column_mapping)
    
    ###### Sets ######
    moiety_index = Unique_Rules_renamed.index.tolist()  # moiety set
    rules_index = Unique_Rules_renamed.columns.values.tolist() # rule set
    #print("Number of unique rules used in this search:", len(rules_index))
    
    ###### parameters ######
        # T(m,r) contains atom stoichiometry for each rule
    T = Unique_Rules_renamed.to_dict(orient="index")
    
    
    # overall reaction input
    T_o = reaction_string_to_moiety_change_dict(desired_reaction,radius)
    for index in moiety_index:
        if index not in T_o:
            T_o[index] = 0.0
    
    # reactant input
    C_R = MolSmilesToMoietyDict(desired_reaction.split('>>')[0],radius)

    # makes sure if a reactant moiety is labeled that moiety is still included in constraint 3
    labeled_moieties = [] # M*
    for moiety in moiety_index:
        if '!' in moiety:
            # water is a fundametal part of life, it is assumed to always be available
            # sometimes water behaves as a cofactor and a reactant in the same mechanism
            # this makes sure even if it a reactant it is always available as a cofactor
            #  radius = 0            radius = 0.5            radius >= 1
            if moiety == '[O:1]!' or moiety == '*[O:1]*!' or moiety == '[H][O:1][H]!' : 
                labeled_moieties.append(moiety)
            else:
                # some reactants are labeled this specifies the only moieties of the reactants 
                if moiety in C_R.keys():
                    continue
                else:
                    labeled_moieties.append(moiety)
                
    # adds the rest of the moieties in moiety_index to C_R as zeros        
    for index in moiety_index:
        if index not in C_R:
            C_R[index] = 0.0
            
    
    ###### variables ######
    y = pulp.LpVariable.dicts("y", rules_index, lowBound=0, cat="Integer")
    
    # create minRules MILP problem
    minRules = pulp.LpProblem("minRules", pulp.LpMinimize)
    
    ####### objective function ####
    minRules += pulp.lpSum([y[r] for r in rules_index])
    
    
    ####### constraints ####
    
    # constraint 1: moiety change balance
    for m in moiety_index:
        minRules += (pulp.lpSum([T[m][r] * y[r] for r in rules_index if T[m][r] != 0]) == [T_o[m]]
                    , "moiety_balance_" + str(moiety_index.index(m)))
        
    # constraint 2: customized constraints
    # the number of steps of the pathway
    minRules += pulp.lpSum([y[r] for r in rules_index]) <= max_steps
    
    solutions = []
    sol_num = 0
    while len(solutions) < iterations:
        
        elapsed_time = time.time() - start_time
        if time_limit != 'none':
            if elapsed_time > time_limit:
                solutions.append('took too long')
                break
        
        try:
            minRules.solve(pulp.PULP_CBC_CMD(msg = 0, timeLimit = time_limit - elapsed_time))
        except:
            solutions.append('error in solver')
        
        if pulp.LpStatus[minRules.status] == 'Infeasible':

            solutions.append('infeasible')
            break
        elif pulp.LpStatus[minRules.status] == 'Not Solved':

            solutions.append('took too long')
            break
        
        # constraint 3: integer cuts
        integer_cut_rules = []
        solution = []
        for r in rules_index:
            for n in range(int(y[r].varValue)):
                integer_cut_rules.append(r)
                solution.append(r)
    
        sol_num += 1        
        length = len(integer_cut_rules) - 1
        minRules += (pulp.lpSum([y[r] for r in integer_cut_rules]) <= length,
                    "integer_cut_" + str(sol_num),)
    
        ## Set K
        K_steps = list(range(1, len(solution) + 1))
    
        z = pulp.LpVariable.dicts("z", (K_steps,rules_index), lowBound=0, upBound=1, cat="Binary")
    
        # create OrderRules MILP problem
        OrderRules = pulp.LpProblem("OrderRules", pulp.LpMinimize)
        
        ####### objective function ####
        OrderRules += 0
        
        # constraint 4: the cumlitive sum of all rules and reactant to be greater than or equal to zero
        for m in moiety_index:
            if m not in labeled_moieties:
                for current_k in K_steps:
                    OrderRules += (C_R[m] + pulp.lpSum([[T[m][r] * z[k][r] for r in set(solution) if T[m][r] != 0] for k in range(1,current_k+1) ]) >= 0
                                ,"cumulative_sum_" + str(moiety_index.index(m)) + '_' + str(current_k))
        
    
        # constraint 5: fixes one rule to be used per step
        for k in K_steps:
            OrderRules += pulp.lpSum([z[k][r] for r in set(solution)]) == 1 
    
        # constraint 6: fixes z_kr == y_r for all rules in solution
        for r in set(solution):
            OrderRules += pulp.lpSum([z[k][r] for k in K_steps] ) == int(y[r].varValue)
        
        OrderRules.solve(pulp.PULP_CBC_CMD(msg = 0))
        
       
        if pulp.LpStatus[OrderRules.status] != 'Optimal':
            continue
            
        solution = [] 
        for k in K_steps:
            for r in rules_index:
                if z[k][r].varValue == 1:
                    solution.append(Unique_Rules.columns.tolist()[r])
        solutions.append(solution)
    
    scores = {}
    for n in range(len(solutions)):
        solution = solutions[n]
        if solution == 'took too long' or solution == 'infeasible' or solution == 'error in solver':
            scores[n] = 0.0
        else:
            mechs = similar_mechs(solution)
            scores[n] = mechs[list(mechs.keys())[0]]

    sorted_scores = dict(sorted(scores.items(), key=lambda item: item[1], reverse=True))
    
    sorted_solutions = []
    for idx in sorted_scores:
        sorted_solutions.append(solutions[idx])
        
    elapsed_time = time.time() - start_time
    print('Final time: '+ str(round(elapsed_time))+' seconds')
    print('Number of solutions from minRules: '+ str(sol_num))

    return sorted_solutions

