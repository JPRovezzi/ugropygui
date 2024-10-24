# =============================================================================
# Sulfur test cases
# =============================================================================
from .case import Case


sulfur_cases = [
    Case(
        identifier="C(=S)=S",
        identifier_type="smiles",
        cases_module="sulfur",
        unifac_result={"CS2": 1},
        psrk_result={"CS2": 1},
        joback_result={},
    ),
    Case(
        identifier="CCCC1=CC=C(C=C1)C(C)S",
        identifier_type="smiles",
        cases_module="sulfur",
        unifac_result={},
        psrk_result={
            "CH3": 2,
            "CH2": 1,
            "ACH": 4,
            "ACCH2": 1,
            "CHSH": 1,
            "AC": 1,
        },
        joback_result={
            "-CH3": 2,
            "-CH2-": 2,
            ">CH-": 1,
            "ring=CH-": 4,
            "ring=C<": 2,
            "-SH": 1,
        },
    ),
    Case(
        identifier="CCCC1=CC=C(CS)C=C1",
        identifier_type="smiles",
        cases_module="sulfur",
        unifac_result={
            "CH3": 1,
            "ACH": 4,
            "ACCH2": 1,
            "CH2SH": 1,
            "CH2": 1,
            "AC": 1,
        },
        psrk_result={
            "CH3": 1,
            "ACH": 4,
            "ACCH2": 1,
            "CH2SH": 1,
            "CH2": 1,
            "AC": 1,
        },
        joback_result={
            "-CH3": 1,
            "-CH2-": 3,
            "ring=CH-": 4,
            "ring=C<": 2,
            "-SH": 1,
        },
    ),
    Case(
        identifier="CC1=CC=C(CS)C=C1",
        identifier_type="smiles",
        cases_module="sulfur",
        unifac_result={"ACH": 4, "ACCH3": 1, "CH2SH": 1, "AC": 1},
        psrk_result={"ACH": 4, "ACCH3": 1, "CH2SH": 1, "AC": 1},
        joback_result={
            "-CH3": 1,
            "-CH2-": 1,
            "ring=CH-": 4,
            "ring=C<": 2,
            "-SH": 1,
        },
    ),
    Case(
        identifier="SCC1=CC=NC=C1",
        identifier_type="smiles",
        cases_module="sulfur",
        unifac_result={"C5H4N": 1, "CH2SH": 1},
        psrk_result={"C5H4N": 1, "CH2SH": 1},
        joback_result={
            "-CH2-": 1,
            "ring=CH-": 4,
            "ring=C<": 1,
            "-N= (ring)": 1,
            "-SH": 1,
        },
    ),
    Case(
        identifier="CS",
        identifier_type="smiles",
        cases_module="sulfur",
        unifac_result={"CH3SH": 1},
        psrk_result={"CH3SH": 1},
        joback_result={"-CH3": 1, "-SH": 1},
    ),
    Case(
        identifier="CCS",
        identifier_type="smiles",
        cases_module="sulfur",
        unifac_result={"CH3": 1, "CH2SH": 1},
        psrk_result={"CH3": 1, "CH2SH": 1},
        joback_result={"-CH3": 1, "-CH2-": 1, "-SH": 1},
    ),
    Case(
        identifier="CSC",
        identifier_type="smiles",
        cases_module="sulfur",
        unifac_result={"CH3": 1, "CH3S": 1},
        psrk_result={"CH3": 1, "CH3S": 1},
        joback_result={"-CH3": 2, "-S- (non-ring)": 1},
    ),
    Case(
        identifier="CCSCC",
        identifier_type="smiles",
        cases_module="sulfur",
        unifac_result={"CH3": 2, "CH2": 1, "CH2S": 1},
        psrk_result={"CH3": 2, "CH2": 1, "CH2S": 1},
        joback_result={"-CH3": 2, "-CH2-": 2, "-S- (non-ring)": 1},
    ),
    Case(
        identifier="CC(C)SC(C)C",
        identifier_type="smiles",
        cases_module="sulfur",
        unifac_result={"CH3": 4, "CH": 1, "CHS": 1},
        psrk_result={"CH3": 4, "CH": 1, "CHS": 1},
        joback_result={"-CH3": 4, ">CH-": 2, "-S- (non-ring)": 1},
    ),
    Case(
        identifier="C1CCS(=O)(=O)C1",
        identifier_type="smiles",
        cases_module="sulfur",
        unifac_result={"CH2": 2, "(CH2)2SU": 1},
        psrk_result={},
        joback_result={},
    ),
    Case(
        identifier="CC1CC(S(=O)(=O)C1)C",
        identifier_type="smiles",
        cases_module="sulfur",
        unifac_result={"CH3": 2, "CH2": 1, "CH": 1, "CH2CHSU": 1},
        psrk_result={},
        joback_result={},
    ),
    Case(
        identifier="CC(C)S",
        identifier_type="smiles",
        cases_module="sulfur",
        unifac_result={},
        psrk_result={"CH3": 2, "CHSH": 1},
        joback_result={"-CH3": 2, ">CH-": 1, "-SH": 1},
    ),
    Case(
        identifier="CC(S)C1=CC=CC=C1",
        identifier_type="smiles",
        cases_module="sulfur",
        unifac_result={},
        psrk_result={"CH3": 1, "ACH": 5, "CHSH": 1, "AC": 1},
        joback_result={
            "-CH3": 1,
            ">CH-": 1,
            "ring=CH-": 5,
            "ring=C<": 1,
            "-SH": 1,
        },
    ),
    Case(
        identifier="CC(C)(C)S",
        identifier_type="smiles",
        cases_module="sulfur",
        unifac_result={},
        psrk_result={"CH3": 3, "CSH": 1},
        joback_result={"-CH3": 3, ">C<": 1, "-SH": 1},
    ),
]
