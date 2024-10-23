# =============================================================================
# Description: This file contains the test cases for the molecules that
# contains nitrogen groups.
#
# TODO: add tests for amines concatenations.
# =============================================================================
from .case import Case


nitrogen_cases = [
    Case(
        identifier="CC1(CC(CC(C1)(C)CN=C=O)N=C=O)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 3, "CH2": 4, "CH": 1, "C": 2, "NCO": 2},
        psrk_result={},
        joback_result={
            "-CH3": 3,
            "-CH2-": 1,
            "=C=": 2,
            "ring-CH2-": 3,
            "ring>CH-": 1,
            "ring>C<": 2,
            "=O (other than above)": 2,
            "-N= (non-ring)": 2,
        },
    ),
    Case(
        identifier="CC(=O)N",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 1, "AMH2": 1},
        psrk_result={"CH3": 1, "AMH2": 1},
        joback_result={"-CH3": 1, ">C=O (non-ring)": 1, "-NH2": 1},
    ),
    Case(
        identifier="CC(=O)NC",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 1, "AMHCH3": 1},
        psrk_result={"CH3": 1, "AMHCH3": 1},
        joback_result={"-CH3": 2, ">C=O (non-ring)": 1, ">NH (non-ring)": 1},
    ),
    Case(
        identifier="CCNC(=O)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 2, "AMHCH2": 1},
        psrk_result={"CH3": 2, "AMHCH2": 1},
        joback_result={
            "-CH3": 2,
            "-CH2-": 1,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 1,
        },
    ),
    Case(
        identifier="CC(=O)N(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 1, "AM(CH3)2": 1},
        psrk_result={"CH3": 1, "AM(CH3)2": 1},
        joback_result={"-CH3": 3, ">C=O (non-ring)": 1, ">N- (non-ring)": 1},
    ),
    Case(
        identifier="CCN(C)C(=O)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 2, "AMCH3CH2": 1},
        psrk_result={"CH3": 2, "AMCH3CH2": 1},
        joback_result={
            "-CH3": 3,
            "-CH2-": 1,
            ">C=O (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CCN(CC)C(=O)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 3, "AM(CH2)2": 1},
        psrk_result={"CH3": 3, "AM(CH2)2": 1},
        joback_result={
            "-CH3": 3,
            "-CH2-": 2,
            ">C=O (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CCN(C(C)C)C(=O)NC(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={},
        psrk_result={},
        joback_result={
            "-CH3": 5,
            "-CH2-": 1,
            ">CH-": 2,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CC(C)NC(=O)N(C)C(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={},
        psrk_result={},
        joback_result={
            "-CH3": 5,
            ">CH-": 2,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CCN(CC)C(=O)NC(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 4, "CHNH": 1, "AM(CH2)2": 1},
        psrk_result={"CH3": 4, "CHNH": 1, "AM(CH2)2": 1},
        joback_result={
            "-CH3": 4,
            "-CH2-": 2,
            ">CH-": 1,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CCN(C)C(=O)NC(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 3, "CHNH": 1, "AMCH3CH2": 1},
        psrk_result={"CH3": 3, "CHNH": 1, "AMCH3CH2": 1},
        joback_result={
            "-CH3": 4,
            "-CH2-": 1,
            ">CH-": 1,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CC(C)NC(=O)N(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 2, "CHNH": 1, "AM(CH3)2": 1},
        psrk_result={"CH3": 2, "CHNH": 1, "AM(CH3)2": 1},
        joback_result={
            "-CH3": 4,
            ">CH-": 1,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CCNC(=O)N(CC)C(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 4, "CH": 1, "CH2N": 1, "AMHCH2": 1},
        psrk_result={"CH3": 4, "CH": 1, "CH2N": 1, "AMHCH2": 1},
        joback_result={
            "-CH3": 4,
            "-CH2-": 2,
            ">CH-": 1,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CCNC(=O)N(C)C(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 3, "CH": 1, "CH3N": 1, "AMHCH2": 1},
        psrk_result={"CH3": 3, "CH": 1, "CH3N": 1, "AMHCH2": 1},
        joback_result={
            "-CH3": 4,
            "-CH2-": 1,
            ">CH-": 1,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CCNC(=O)N(CC)CC",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 3, "CH2NH": 1, "AM(CH2)2": 1},
        psrk_result={"CH3": 3, "CH2NH": 1, "AM(CH2)2": 1},
        joback_result={
            "-CH3": 3,
            "-CH2-": 3,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CCNC(=O)N(C)CC",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 2, "CH2NH": 1, "AMCH3CH2": 1},
        psrk_result={"CH3": 2, "CH2NH": 1, "AMCH3CH2": 1},
        joback_result={
            "-CH3": 3,
            "-CH2-": 2,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CCNC(=O)N(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH2NH": 1, "AM(CH3)2": 1, "CH3": 1},
        psrk_result={"CH2NH": 1, "AM(CH3)2": 1, "CH3": 1},
        joback_result={
            "-CH3": 3,
            "-CH2-": 1,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CCN(C(C)C)C(=O)NC",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 3, "CH": 1, "CH2N": 1, "AMHCH3": 1},
        psrk_result={"CH3": 3, "CH": 1, "CH2N": 1, "AMHCH3": 1},
        joback_result={
            "-CH3": 4,
            "-CH2-": 1,
            ">CH-": 1,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CNC(=O)N(C)C(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 2, "CH": 1, "CH3N": 1, "AMHCH3": 1},
        psrk_result={"CH3": 2, "CH": 1, "CH3N": 1, "AMHCH3": 1},
        joback_result={
            "-CH3": 4,
            ">CH-": 1,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CCN(CC)C(=O)NC",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 2, "CH3NH": 1, "AM(CH2)2": 1},
        psrk_result={"CH3": 2, "CH3NH": 1, "AM(CH2)2": 1},
        joback_result={
            "-CH3": 3,
            "-CH2-": 2,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CCN(C)C(=O)NC",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 1, "CH3NH": 1, "AMCH3CH2": 1},
        psrk_result={"CH3": 1, "CH3NH": 1, "AMCH3CH2": 1},
        joback_result={
            "-CH3": 3,
            "-CH2-": 1,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CNC(=O)N(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3NH": 1, "AM(CH3)2": 1},
        psrk_result={"CH3NH": 1, "AM(CH3)2": 1},
        joback_result={
            "-CH3": 3,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CC(C)NC(=O)NC(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={},
        psrk_result={},
        joback_result={
            "-CH3": 4,
            ">CH-": 2,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 2,
        },
    ),
    Case(
        identifier="CCNC(=O)NC(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 3, "CHNH": 1, "AMHCH2": 1},
        psrk_result={"CH3": 3, "CHNH": 1, "AMHCH2": 1},
        joback_result={
            "-CH3": 3,
            "-CH2-": 1,
            ">CH-": 1,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 2,
        },
    ),
    Case(
        identifier="CCNC(=O)NCC",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 2, "CH2NH": 1, "AMHCH2": 1},
        psrk_result={"CH3": 2, "CH2NH": 1, "AMHCH2": 1},
        joback_result={
            "-CH3": 2,
            "-CH2-": 2,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 2,
        },
    ),
    Case(
        identifier="CNC(=O)NC(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 2, "CHNH": 1, "AMHCH3": 1},
        psrk_result={"CH3": 2, "CHNH": 1, "AMHCH3": 1},
        joback_result={
            "-CH3": 3,
            ">CH-": 1,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 2,
        },
    ),
    Case(
        identifier="CCNC(=O)NC",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result=[
            {"CH2NH": 1, "AMHCH3": 1, "CH3": 1},
            {"CH3": 1, "CH3NH": 1, "AMHCH2": 1},
        ],
        psrk_result=[
            {"CH2NH": 1, "AMHCH3": 1, "CH3": 1},
            {"CH3": 1, "CH3NH": 1, "AMHCH2": 1},
        ],
        joback_result={
            "-CH3": 2,
            "-CH2-": 1,
            ">C=O (non-ring)": 1,
            ">NH (non-ring)": 2,
        },
    ),
    Case(
        identifier="CNC(=O)NC",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3NH": 1, "AMHCH3": 1},
        psrk_result={"CH3NH": 1, "AMHCH3": 1},
        joback_result={"-CH3": 2, ">C=O (non-ring)": 1, ">NH (non-ring)": 2},
    ),
    Case(
        identifier="CN",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3NH2": 1},
        psrk_result={"CH3NH2": 1},
        joback_result={"-CH3": 1, "-NH2": 1},
    ),
    Case(
        identifier="CC(C)N",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 2, "CHNH2": 1},
        psrk_result={"CH3": 2, "CHNH2": 1},
        joback_result={"-CH3": 2, ">CH-": 1, "-NH2": 1},
    ),
    Case(
        identifier="CCCN",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 1, "CH2": 1, "CH2NH2": 1},
        psrk_result={"CH3": 1, "CH2": 1, "CH2NH2": 1},
        joback_result={"-CH3": 1, "-CH2-": 2, "-NH2": 1},
    ),
    Case(
        identifier="COC1=C(OC)C2=C3C(CC22CCC(=O)C=C2)NCCC3=C1",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={
            "CH2": 2,
            "C": 1,
            "CH=CH": 1,
            "ACH": 1,
            "AC": 3,
            "ACCH2": 1,
            "ACCH": 1,
            "CH2CO": 1,
            "CH3O": 2,
            "CH2NH": 1,
        },
        psrk_result={
            "CH2": 2,
            "C": 1,
            "CH=CH": 1,
            "ACH": 1,
            "AC": 3,
            "ACCH2": 1,
            "ACCH": 1,
            "CH2CO": 1,
            "CH3O": 2,
            "CH2NH": 1,
        },
        joback_result={
            "-CH3": 2,
            "ring-CH2-": 5,
            "ring>CH-": 1,
            "ring>C<": 1,
            "ring=CH-": 3,
            "ring=C<": 5,
            "-O- (non-ring)": 2,
            ">C=O (ring)": 1,
            ">NH (ring)": 1,
        },
    ),
    Case(
        identifier="C1C2=C(C=CN1)C3=CC=CC=C3N2",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={},
        psrk_result={},
        joback_result={
            "ring-CH2-": 1,
            "ring=CH-": 6,
            "ring=C<": 4,
            ">NH (ring)": 2,
        },
    ),
    Case(
        identifier="CC(C)N(CN)C(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={},
        psrk_result={},
        joback_result={
            "-CH3": 4,
            "-CH2-": 1,
            ">CH-": 2,
            "-NH2": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CC(C)N(C)CN",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 2, "CH": 1, "CH2NH2": 1, "CH3N": 1},
        psrk_result={"CH3": 2, "CH": 1, "CH2NH2": 1, "CH3N": 1},
        joback_result={
            "-CH3": 3,
            "-CH2-": 1,
            ">CH-": 1,
            "-NH2": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CC(C)NC(C)NC(C)(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 6, "C": 1, "CHNH": 2},
        psrk_result={"CH3": 6, "C": 1, "CHNH": 2},
        joback_result={"-CH3": 6, ">CH-": 2, ">C<": 1, ">NH (non-ring)": 2},
    ),
    Case(
        identifier="CC(C)NC(C)N",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 3, "CHNH2": 1, "CHNH": 1},
        psrk_result={"CH3": 3, "CHNH2": 1, "CHNH": 1},
        joback_result={"-CH3": 3, ">CH-": 2, "-NH2": 1, ">NH (non-ring)": 1},
    ),
    Case(
        identifier="CCC(C)(C)NC(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 5, "CH2": 1, "C": 1, "CHNH": 1},
        psrk_result={"CH3": 5, "CH2": 1, "C": 1, "CHNH": 1},
        joback_result={
            "-CH3": 5,
            "-CH2-": 1,
            ">CH-": 1,
            ">C<": 1,
            ">NH (non-ring)": 1,
        },
    ),
    Case(
        identifier="CCNC(C)CC",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result=[
            {"CH3": 3, "CH2": 2, "CHNH": 1},
            {"CH3": 3, "CH2": 1, "CH": 1, "CH2NH": 1},
        ],
        psrk_result=[
            {"CH3": 3, "CH2": 2, "CHNH": 1},
            {"CH3": 3, "CH2": 1, "CH": 1, "CH2NH": 1},
        ],
        joback_result={"-CH3": 3, "-CH2-": 2, ">CH-": 1, ">NH (non-ring)": 1},
    ),
    Case(
        identifier="CCCNC",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result=[
            {"CH3": 2, "CH2": 1, "CH2NH": 1},
            {"CH3": 1, "CH2": 2, "CH3NH": 1},
        ],
        psrk_result=[
            {"CH3": 2, "CH2": 1, "CH2NH": 1},
            {"CH3": 1, "CH2": 2, "CH3NH": 1},
        ],
        joback_result={"-CH3": 2, "-CH2-": 2, ">NH (non-ring)": 1},
    ),
    Case(
        identifier="CN1CCCC1CC(=O)CC1CCCN1C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result=[
            {"CH3": 1, "CH2": 6, "CH": 2, "CH2CO": 1, "CH3N": 1, "CH2N": 1},
            {"CH2": 7, "CH": 2, "CH2CO": 1, "CH3N": 2},
            {"CH3": 2, "CH2": 5, "CH": 2, "CH2CO": 1, "CH2N": 2},
        ],
        psrk_result=[
            {"CH3": 1, "CH2": 6, "CH": 2, "CH2CO": 1, "CH3N": 1, "CH2N": 1},
            {"CH2": 7, "CH": 2, "CH2CO": 1, "CH3N": 2},
            {"CH3": 2, "CH2": 5, "CH": 2, "CH2CO": 1, "CH2N": 2},
        ],
        joback_result={},
    ),
    Case(
        identifier="CNC",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 1, "CH3NH": 1},
        psrk_result={"CH3": 1, "CH3NH": 1},
        joback_result={"-CH3": 2, ">NH (non-ring)": 1},
    ),
    Case(
        identifier="CCNCC",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 2, "CH2": 1, "CH2NH": 1},
        psrk_result={"CH3": 2, "CH2": 1, "CH2NH": 1},
        joback_result={"-CH3": 2, "-CH2-": 2, ">NH (non-ring)": 1},
    ),
    Case(
        identifier="CC(C)NC(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 4, "CH": 1, "CHNH": 1},
        psrk_result={"CH3": 4, "CH": 1, "CHNH": 1},
        joback_result={"-CH3": 4, ">CH-": 2, ">NH (non-ring)": 1},
    ),
    Case(
        identifier="CC(C)NCN",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 2, "CH2NH2": 1, "CHNH": 1},
        psrk_result={"CH3": 2, "CH2NH2": 1, "CHNH": 1},
        joback_result={
            "-CH3": 2,
            "-CH2-": 1,
            ">CH-": 1,
            "-NH2": 1,
            ">NH (non-ring)": 1,
        },
    ),
    Case(
        identifier="C1CN2CCC1CC2",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH2": 5, "CH": 1, "CH2N": 1},
        psrk_result={"CH2": 5, "CH": 1, "CH2N": 1},
        joback_result={},
    ),
    Case(
        identifier="CCN(CC(=O)CC)C1=CC=CC=C1",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={
            "CH3": 2,
            "CH2": 1,
            "ACH": 5,
            "AC": 1,
            "CH2CO": 1,
            "CH2N": 1,
        },
        psrk_result={
            "CH3": 2,
            "CH2": 1,
            "ACH": 5,
            "AC": 1,
            "CH2CO": 1,
            "CH2N": 1,
        },
        joback_result={
            "-CH3": 2,
            "-CH2-": 3,
            "ring=CH-": 5,
            "ring=C<": 1,
            ">C=O (non-ring)": 1,
            ">N- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CCN(C(C)C)C(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 5, "CH": 2, "CH2N": 1},
        psrk_result={"CH3": 5, "CH": 2, "CH2N": 1},
        joback_result={"-CH3": 5, "-CH2-": 1, ">CH-": 2, ">N- (non-ring)": 1},
    ),
    Case(
        identifier="CCN(C)CC",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result=[
            {"CH3": 2, "CH2": 2, "CH3N": 1},
            {"CH3": 3, "CH2": 1, "CH2N": 1},
        ],
        psrk_result=[
            {"CH3": 2, "CH2": 2, "CH3N": 1},
            {"CH3": 3, "CH2": 1, "CH2N": 1},
        ],
        joback_result={"-CH3": 3, "-CH2-": 2, ">N- (non-ring)": 1},
    ),
    Case(
        identifier="CN(C)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 2, "CH3N": 1},
        psrk_result={"CH3": 2, "CH3N": 1},
        joback_result={"-CH3": 3, ">N- (non-ring)": 1},
    ),
    Case(
        identifier="CCN(CC)CC",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 3, "CH2": 2, "CH2N": 1},
        psrk_result={"CH3": 3, "CH2": 2, "CH2N": 1},
        joback_result={"-CH3": 3, "-CH2-": 3, ">N- (non-ring)": 1},
    ),
    Case(
        identifier="C1=CC=C(C=C1)N",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"ACH": 5, "ACNH2": 1},
        psrk_result={"ACH": 5, "ACNH2": 1},
        joback_result={"ring=CH-": 5, "ring=C<": 1, "-NH2": 1},
    ),
    Case(
        identifier="C1=CC2=C(C=C1)C1=C(C=NC=C1)C1=C2C=CC=C1",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"ACH": 8, "AC": 4, "C5H3N": 1},
        psrk_result={"ACH": 8, "AC": 4, "C5H3N": 1},
        joback_result={"ring=CH-": 11, "ring=C<": 6, "-N= (ring)": 1},
    ),
    Case(
        identifier="CC1=CC(C)=CN=C1",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 2, "C5H3N": 1},
        psrk_result={"CH3": 2, "C5H3N": 1},
        joback_result={
            "-CH3": 2,
            "ring=CH-": 3,
            "ring=C<": 2,
            "-N= (ring)": 1,
        },
    ),
    Case(
        identifier="O=C1CCC2=C(O1)C=CN=C2",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH2COO": 1, "C5H3N": 1, "CH2": 1},
        psrk_result={"CH2COO": 1, "C5H3N": 1, "CH2": 1},
        joback_result={
            "ring-CH2-": 2,
            "ring=CH-": 3,
            "ring=C<": 2,
            "-COO- (ester)": 1,
            "-N= (ring)": 1,
        },
    ),
    Case(
        identifier="CCCC1=CC=C(CC2=CC=NC=C2)C=C1",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 1, "ACH": 4, "ACCH2": 2, "C5H4N": 1, "CH2": 1},
        psrk_result={"CH3": 1, "ACH": 4, "ACCH2": 2, "C5H4N": 1, "CH2": 1},
        joback_result={
            "-CH3": 1,
            "-CH2-": 3,
            "ring=CH-": 8,
            "ring=C<": 3,
            "-N= (ring)": 1,
        },
    ),
    Case(
        identifier="CN1CCCC1C2=CN=CC=C2",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result=[
            {"CH2": 3, "CH3N": 1, "C5H4N": 1, "CH": 1},
            {"CH3": 1, "CH2": 2, "CH": 1, "CH2N": 1, "C5H4N": 1},
        ],
        psrk_result=[
            {"CH2": 3, "CH3N": 1, "C5H4N": 1, "CH": 1},
            {"CH3": 1, "CH2": 2, "CH": 1, "CH2N": 1, "C5H4N": 1},
        ],
        joback_result={},
    ),
    Case(
        identifier="CC1=CC(C)=C(C)C=N1",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={},
        psrk_result={},
        joback_result={
            "-CH3": 3,
            "ring=CH-": 2,
            "ring=C<": 3,
            "-N= (ring)": 1,
        },
    ),
    Case(
        identifier="C1=CC=C(C=C1)C1=CC=NC=C1",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"ACH": 5, "AC": 1, "C5H4N": 1},
        psrk_result={"ACH": 5, "AC": 1, "C5H4N": 1},
        joback_result={"ring=CH-": 9, "ring=C<": 2, "-N= (ring)": 1},
    ),
    Case(
        identifier="CC(=C)C1=CC=NC=C1",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 1, "CH2=C": 1, "C5H4N": 1},
        psrk_result={"CH3": 1, "CH2=C": 1, "C5H4N": 1},
        joback_result={
            "-CH3": 1,
            "=CH2": 1,
            "=C<": 1,
            "ring=CH-": 4,
            "ring=C<": 1,
            "-N= (ring)": 1,
        },
    ),
    Case(
        identifier="CC1=NC=CC(O)=C1",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"C5H3N": 1, "CH3": 1, "OH": 1},
        psrk_result={"C5H3N": 1, "CH3": 1, "OH": 1},
        joback_result={
            "-CH3": 1,
            "ring=CH-": 3,
            "ring=C<": 2,
            "-OH (phenol)": 1,
            "-N= (ring)": 1,
        },
    ),
    Case(
        identifier="C1=CC=NC=C1",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"C5H5N": 1},
        psrk_result={"C5H5N": 1},
        joback_result={"ring=CH-": 5, "-N= (ring)": 1},
    ),
    Case(
        identifier="CC1=CN=CC=C1",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"C5H4N": 1, "CH3": 1},
        psrk_result={"C5H4N": 1, "CH3": 1},
        joback_result={
            "-CH3": 1,
            "ring=CH-": 4,
            "ring=C<": 1,
            "-N= (ring)": 1,
        },
    ),
    Case(
        identifier="CC1=C(N=CC=C1)C",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"C5H3N": 1, "CH3": 2},
        psrk_result={"C5H3N": 1, "CH3": 2},
        joback_result={
            "-CH3": 2,
            "ring=CH-": 3,
            "ring=C<": 2,
            "-N= (ring)": 1,
        },
    ),
    Case(
        identifier="CC#N",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3CN": 1},
        psrk_result={"CH3CN": 1},
        joback_result={"-CH3": 1, "-CN": 1},
    ),
    Case(
        identifier="CCC#N",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 1, "CH2CN": 1},
        psrk_result={"CH3": 1, "CH2CN": 1},
        joback_result={"-CH3": 1, "-CH2-": 1, "-CN": 1},
    ),
    Case(
        identifier="CCCC1=CC=C(C[N+]([O-])=O)C=C1",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={
            "CH3": 1,
            "ACH": 4,
            "ACCH2": 1,
            "CH2NO2": 1,
            "CH2": 1,
            "AC": 1,
        },
        psrk_result={
            "CH3": 1,
            "ACH": 4,
            "ACCH2": 1,
            "CH2NO2": 1,
            "CH2": 1,
            "AC": 1,
        },
        joback_result={
            "-CH3": 1,
            "-CH2-": 3,
            "ring=CH-": 4,
            "ring=C<": 2,
            "-NO2": 1,
        },
    ),
    Case(
        identifier="[O-][N+](=O)CC1=CC=CC=C1",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"ACH": 5, "CH2NO2": 1, "AC": 1},
        psrk_result={"ACH": 5, "CH2NO2": 1, "AC": 1},
        joback_result={"-CH2-": 1, "ring=CH-": 5, "ring=C<": 1, "-NO2": 1},
    ),
    Case(
        identifier="C[N+](=O)[O-]",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3NO2": 1},
        psrk_result={"CH3NO2": 1},
        joback_result={"-CH3": 1, "-NO2": 1},
    ),
    Case(
        identifier="CCC[N+](=O)[O-]",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 1, "CH2": 1, "CH2NO2": 1},
        psrk_result={"CH3": 1, "CH2": 1, "CH2NO2": 1},
        joback_result={"-CH3": 1, "-CH2-": 2, "-NO2": 1},
    ),
    Case(
        identifier="CC(C)[N+](=O)[O-]",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"CH3": 2, "CHNO2": 1},
        psrk_result={"CH3": 2, "CHNO2": 1},
        joback_result={"-CH3": 2, ">CH-": 1, "-NO2": 1},
    ),
    Case(
        identifier="C1=CC=C(C=C1)[N+](=O)[O-]",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={"ACH": 5, "ACNO2": 1},
        psrk_result={"ACH": 5, "ACNO2": 1},
        joback_result={"ring=CH-": 5, "ring=C<": 1, "-NO2": 1},
    ),
    Case(
        identifier="[O-][N+](=O)C1=CC=NC=C1",
        identifier_type="smiles",
        cases_module="nitrogen",
        unifac_result={},
        psrk_result={},
        joback_result={
            "ring=CH-": 4,
            "ring=C<": 1,
            "-N= (ring)": 1,
            "-NO2": 1,
        },
    ),
]
