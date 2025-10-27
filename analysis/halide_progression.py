from evodex.operators import extract_operator

if __name__ == "__main__":
    smirks = "[C:4][C:1][Cl]>>[C:4][C:1][O][C]"
    operator = extract_operator(smirks)
    print("ethanol",operator)

    smirks = "[C:4]=[C:1][Cl]>>[C:4]=[C:1][O][C]"
    operator = extract_operator(smirks)
    print("vinyl alcohol",operator)

    smirks = "[C:4]=[C:5][C:1][Cl]>>[C:4]=[C:5][C:1][O][C]"
    operator = extract_operator(smirks)
    print("allyl alcohol", operator)

    smirks = "[C:6]=[C:5][C:4]=[C:1][Cl]>>[C:6]=[C:5][C:4]=[C:1][O][C]"
    operator = extract_operator(smirks)
    print("dibutenyl alcohol",operator)

    smirks = "[C:7]1=[C:6]-[C:5]=[C:4]-[C:3]=[C:2]1[Cl]>>[C:7]1=[C:6]-[C:5]=[C:4]-[C:3]=[C:2]1[O][C]"
    operator = extract_operator(smirks)
    print("phenol",operator)

    smirks = "[C:7]1=[C:6]-[C:5]=[C:4]-[C:3]=[C:2]1[C:8][Cl]>>[C:7]1=[C:6]-[C:5]=[C:4]-[C:3]=[C:2]1[C:8][O][C]"
    operator = extract_operator(smirks)
    print("benzyl alcohol",operator)

    smirks = "[C:7]1=[C:6]-[C:5]=[C:4]-[C:3]=[C:2]1[C:8][C:9][Cl]>>[C:7]1=[C:6]-[C:5]=[C:4]-[C:3]=[C:2]1[C:8][C:9][O][C]"
    operator = extract_operator(smirks)
    print("phenethyl alcohol",operator)