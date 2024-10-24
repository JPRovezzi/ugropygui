# =============================================================================
# Description: Test cases for hydrocarbons module.
#
# This module contains the test cases for hydrocarbons.
#
# Only aliphatic and saturated hydrocarbons are must be considered here,
# linear o cyclic.
# =============================================================================
from .case import Case


hydrocarbons_cases = [
    Case(  # Ethane
        "CC",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH3": 2},
        psrk_result={"CH3": 2},
        joback_result={"-CH3": 2},
    ),
    Case(  # 2,2-dimethylpropane
        "CC(C)(C)C",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH3": 4, "C": 1},
        psrk_result={"CH3": 4, "C": 1},
        joback_result={"-CH3": 4, ">C<": 1},
    ),
    Case(  # Cyclohexane
        "C1CCCCC1",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH2": 6},
        psrk_result={"CH2": 6},
        joback_result={"ring-CH2-": 6},
    ),
    Case(  # 2-methylpropane
        "CC(C)C",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH3": 3, "CH": 1},
        psrk_result={"CH3": 3, "CH": 1},
        joback_result={"-CH3": 3, ">CH-": 1},
    ),
    Case(  # Hexane
        "CCCCCC",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH3": 2, "CH2": 4},
        psrk_result={"CH3": 2, "CH2": 4},
        joback_result={"-CH3": 2, "-CH2-": 4},
    ),
    Case(
        "C1CC2CCCC3CCCC1C23",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH2": 8, "CH": 4},
        psrk_result={"CH2": 8, "CH": 4},
        joback_result={"ring-CH2-": 8, "ring>CH-": 4},
    ),
    Case(
        "C1C2CCCCC2C2CCCCC12",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH2": 9, "CH": 4},
        psrk_result={"CH2": 9, "CH": 4},
        joback_result={"ring-CH2-": 9, "ring>CH-": 4},
    ),
    Case(
        "C1C2CC1CCCC2",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH2": 6, "CH": 2},
        psrk_result={"CH2": 6, "CH": 2},
        joback_result={"ring-CH2-": 6, "ring>CH-": 2},
    ),
    Case(
        "C1CCCCCCCC1",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH2": 9},
        psrk_result={"CH2": 9},
        joback_result={"ring-CH2-": 9},
    ),
    Case(
        "C1C2CC3CC1CC(C2)C3",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH2": 6, "CH": 4},
        psrk_result={"CH2": 6, "CH": 4},
        joback_result={"ring-CH2-": 6, "ring>CH-": 4},
    ),
    Case(
        "C12C3C1C1C2C31",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH": 6},
        psrk_result={"CH": 6},
        joback_result={"ring>CH-": 6},
    ),
    Case(
        "C1CC2CC1CCC2",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH2": 6, "CH": 2},
        psrk_result={"CH2": 6, "CH": 2},
        joback_result={"ring-CH2-": 6, "ring>CH-": 2},
    ),
    Case(
        "C1CC2CC3CCC2CC13",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH2": 6, "CH": 4},
        psrk_result={"CH2": 6, "CH": 4},
        joback_result={"ring-CH2-": 6, "ring>CH-": 4},
    ),
    Case(
        "C12C3C4C1C1C2C3C41",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH": 8},
        psrk_result={"CH": 8},
        joback_result={"ring>CH-": 8},
    ),
    Case(
        "C1CC1",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH2": 3},
        psrk_result={"CH2": 3},
        joback_result={"ring-CH2-": 3},
    ),
    Case(
        "C1CCC1",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH2": 4},
        psrk_result={"CH2": 4},
        joback_result={"ring-CH2-": 4},
    ),
    Case(
        "CC12C3CCC4CCC1C234",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH3": 1, "CH2": 4, "CH": 3, "C": 2},
        psrk_result={"CH3": 1, "CH2": 4, "CH": 3, "C": 2},
        joback_result={"-CH3": 1, "ring-CH2-": 4, "ring>CH-": 3, "ring>C<": 2},
    ),
    Case(
        "CCC(CC)C(C)(C)C",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH3": 5, "CH2": 2, "CH": 1, "C": 1},
        psrk_result={"CH3": 5, "CH2": 2, "CH": 1, "C": 1},
        joback_result={"-CH3": 5, "-CH2-": 2, ">CH-": 1, ">C<": 1},
    ),
    Case(
        "C1CCC2CCCCC2C1",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH2": 8, "CH": 2},
        psrk_result={"CH2": 8, "CH": 2},
        joback_result={"ring-CH2-": 8, "ring>CH-": 2},
    ),
    Case(
        "C1CCC(CC1)CC2CCCCC2",
        "smiles",
        "hydrocarbons",
        unifac_result={"CH2": 11, "CH": 2},
        psrk_result={"CH2": 11, "CH": 2},
        joback_result={"-CH2-": 1, "ring-CH2-": 10, "ring>CH-": 2},
    ),
]
