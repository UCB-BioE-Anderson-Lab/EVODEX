

from evodex.operators import extract_operator

if __name__ == "__main__":
    smirks = "[C:4][C:1][O:2]>>[C:4][C:1][O:2][C]"
    operator = extract_operator(smirks)
    print("ethanol",operator)

    smirks = "[C:4]=[C:1][O:2]>>[C:4]=[C:1][O:2][C]"
    operator = extract_operator(smirks)
    print("vinyl alcohol",operator)

    smirks = "[C:4]=[C:5][C:1][O:2]>>[C:4]=[C:5][C:1][O:2][C]"
    operator = extract_operator(smirks)
    print("allyl alcohol", operator)

    smirks = "[C:6]=[C:5][C:4]=[C:1][O:2]>>[C:6]=[C:5][C:4]=[C:1][O:2][C]"
    operator = extract_operator(smirks)
    print("dibutenyl alcohol",operator)

    smirks = "[c:7]1[c:6][c:5][c:4][c:3][c:2]1[O:1]>>[c:7]1[c:6][c:5][c:4][c:3][c:2]1[O:1][C]"
    operator = extract_operator(smirks)
    print("phenol",operator)

    smirks = "[c:7]1[c:6][c:5][c:4][c:3][c:2]1[C:8][O:1]>>[c:7]1[c:6][c:5][c:4][c:3][c:2]1[C:8][O:1][C]"
    operator = extract_operator(smirks)
    print("benzyl alcohol",operator)

    smirks = "[c:7]1[c:6][c:5][c:4][c:3][c:2]1[C:8][C:9][O:1]>>[c:7]1[c:6][c:5][c:4][c:3][c:2]1[C:8][C:9][O:1][C]"
    operator = extract_operator(smirks)
    print("phenethyl alcohol",operator)