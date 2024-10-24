# =============================================================================
# Description: Test cases for aldehydes
# =============================================================================
from .case import Case


aldehydes_cases = [
    Case(
        identifier="C(=O)C=O",
        identifier_type="smiles",
        cases_module="aldhydes",
        commentary="",
        unifac_result={"HCO": 2},
        psrk_result={"HCO": 2},
        joback_result={"O=CH- (aldehyde)": 2},
    ),
    Case(
        identifier="C1=CC=C(C(=C1)C=O)O",
        identifier_type="smiles",
        cases_module="aldhydes",
        commentary="salicylaldehyde",
        unifac_result={"ACH": 4, "ACOH": 1, "AC": 1, "HCO": 1},
        psrk_result={"ACH": 4, "ACOH": 1, "AC": 1, "HCO": 1},
        joback_result={
            "ring=CH-": 4,
            "ring=C<": 2,
            "-OH (phenol)": 1,
            "O=CH- (aldehyde)": 1,
        },
    ),
    Case(
        identifier="CC(C=C)C=O",
        identifier_type="smiles",
        cases_module="aldhydes",
        commentary="2-Methyl-3-butenal",
        unifac_result={"CH3": 1, "CH": 1, "CH2=CH": 1, "HCO": 1},
        psrk_result={"CH3": 1, "CH": 1, "CH2=CH": 1, "HCO": 1},
        joback_result={
            "-CH3": 1,
            ">CH-": 1,
            "=CH2": 1,
            "=CH-": 1,
            "O=CH- (aldehyde)": 1,
        },
    ),
    Case(
        identifier="C1=CC=C(C=C1)C=CC=O",
        identifier_type="smiles",
        cases_module="aldhydes",
        commentary="Cinnamaldehyde",
        unifac_result={"ACH": 5, "AC": 1, "CH=CH": 1, "HCO": 1},
        psrk_result={"ACH": 5, "AC": 1, "CH=CH": 1, "HCO": 1},
        joback_result={
            "=CH-": 2,
            "ring=CH-": 5,
            "ring=C<": 1,
            "O=CH- (aldehyde)": 1,
        },
    ),
    Case(
        identifier="C1=CC=C(C=C1)C=O",
        identifier_type="smiles",
        cases_module="aldhydes",
        commentary="benzaldehyde",
        unifac_result={"ACH": 5, "AC": 1, "HCO": 1},
        psrk_result={"ACH": 5, "AC": 1, "HCO": 1},
        joback_result={"ring=CH-": 5, "ring=C<": 1, "O=CH- (aldehyde)": 1},
    ),
    Case(
        identifier="C1CCC(CC1)C=O",
        identifier_type="smiles",
        cases_module="aldhydes",
        commentary="cyclohexanecarbaldehyde",
        unifac_result={"CH2": 5, "CH": 1, "HCO": 1},
        psrk_result={"CH2": 5, "CH": 1, "HCO": 1},
        joback_result={"ring-CH2-": 5, "ring>CH-": 1, "O=CH- (aldehyde)": 1},
    ),
    Case(
        identifier="CCCCC=O",
        identifier_type="smiles",
        cases_module="aldhydes",
        commentary="pentanal",
        unifac_result={"CH3": 1, "CH2": 3, "HCO": 1},
        psrk_result={"CH3": 1, "CH2": 3, "HCO": 1},
        joback_result={"-CH3": 1, "-CH2-": 3, "O=CH- (aldehyde)": 1},
    ),
    Case(
        identifier="CC(C)CC=O",
        identifier_type="smiles",
        cases_module="aldhydes",
        commentary="3-methylbutanal",
        unifac_result={"CH3": 2, "CH2": 1, "CH": 1, "HCO": 1},
        psrk_result={"CH3": 2, "CH2": 1, "CH": 1, "HCO": 1},
        joback_result={
            "-CH3": 2,
            "-CH2-": 1,
            ">CH-": 1,
            "O=CH- (aldehyde)": 1,
        },
    ),
    Case(
        identifier="CC=O",
        identifier_type="smiles",
        cases_module="aldhydes",
        commentary="acetaldehyde",
        unifac_result={"CH3": 1, "HCO": 1},
        psrk_result={"CH3": 1, "HCO": 1},
        joback_result={"-CH3": 1, "O=CH- (aldehyde)": 1},
    ),
    Case(
        identifier=r"CCCCCC\C(C=O)=C/C1=CC=CC=C1",
        identifier_type="smiles",
        cases_module="aldhydes",
        commentary="2-Hexyl-3-Phenyl-2-Propenal",
        unifac_result={
            "ACH": 5,
            "AC": 1,
            "CH=C": 1,
            "CH2": 5,
            "CH3": 1,
            "HCO": 1,
        },
        psrk_result={
            "ACH": 5,
            "AC": 1,
            "CH=C": 1,
            "CH2": 5,
            "CH3": 1,
            "HCO": 1,
        },
        joback_result={
            "-CH3": 1,
            "-CH2-": 5,
            "=CH-": 1,
            "=C<": 1,
            "ring=CH-": 5,
            "ring=C<": 1,
            "O=CH- (aldehyde)": 1,
        },
    ),
]
