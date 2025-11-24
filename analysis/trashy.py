'''
Description:  Inputs the substrate and product names for a reaction and returns an object with:

{
substrate_name: "...",
product_name: "...",
substrate_smiles: "...",
product_smiles: "...",
matches: [], # The list of EVODEX IDs and their links in EVODEX
conclusion: "No match found to any known enzymatic reactions"
}

'''


def evaluation_reaction(substrate_name: str, product_name: str) -> object:
    """
    Evaluate a reaction based on substrate and product names,
    returning details and matching results from EVODEX.
    """
    # Look up smiles of substrate (if no hit, return "Could not resolve substrate to SMILES")

    # Look up smiles of product (if no hit, return "Could not resolve product to SMILES")

    # Run it through EVODEX at each level of abstraction
    # Check if there is an F match.  If none, return "No match found to any known enzymatic reactions"
    # Check if there is an C match.  If none, return "The reaction has precedent with other reactions with a similar formula difference, but the specific mechanism is unprecedented"
    # Check if there is an N match.  If none, return "The reaction has precedent with other reactions sharing similar reactive groups, but the specific mechanism is unprecedented"
    # Check if there is an E match.  If none, return "The reaction matches known reaction mechanisms partially but not the entire electronic manifold is present."

    # If got this far, it has an E match, so return that
    pass
