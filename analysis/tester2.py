
from rdkit import Chem
from rdkit.Chem import rdChemReactions

# Define a reaction with two reactants: carboxylic acid and amine → amide
rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')

# Define two substrate molecules
reactant1 = Chem.MolFromSmiles('CC(=O)O')  # acetic acid
reactant2 = Chem.MolFromSmiles('CN')       # methylamine

# Run the reaction
products = rxn.RunReactants((reactant1, reactant2))

# Display the result
if products:
    for product_set in products:
        for product in product_set:
            print(Chem.MolToSmiles(product))
else:
    print("No products generated.")


# Run the reaction
products = rxn.RunReactants((reactant2, reactant1))

# Display the result
if products:
    for product_set in products:
        for product in product_set:
            print(Chem.MolToSmiles(product))
else:
    print("No products generated.")

