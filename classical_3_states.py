def equality(a, b, c, letter):
    if a == b and b == c and c == letter:
        return 1
    return 0


def equality_two(a, b, c, letter):
    if (a == b and a == letter) or (a == c and a == letter) or (b == c and b == letter):
        return 1
    return 0


def one_equality(a, b, c, letter):
    if (a == letter and b != letter and c != letter) or (b == letter and a != letter and c != letter) or (
            c == letter and a != letter and b != letter):
        return 1
    return 0


def classical_3_states(dssp_struct, stride_struct, promotif_struct, len):
    classical_3_states_tab = []

    for i in range(len):
        a = dssp_struct[i]
        b = stride_struct[i]
        c = promotif_struct[i]

        if equality(a, b, c, "H") == 1:
            classical_3_states_tab.append("H")
        elif equality(a, b, c, "E") == 1:
            classical_3_states_tab.append("E")
        elif equality(a, b, c, "C") == 1:
            classical_3_states_tab.append("C")

        elif equality_two(a, b, c, "H") == 1:
            classical_3_states_tab.append("h")
        elif equality_two(a, b, c, "E") == 1:
            classical_3_states_tab.append("e")
        elif equality_two(a, b, c, "C") == 1:
            classical_3_states_tab.append("c")

        elif one_equality(a, b, c, "H") == 1:
            classical_3_states_tab.append("*")
        elif one_equality(a, b, c, "E") == 1:
            classical_3_states_tab.append("*")
        elif one_equality(a, b, c, "C") == 1:
            classical_3_states_tab.append("*")

        else:
            classical_3_states_tab.append("-")

    return classical_3_states_tab
