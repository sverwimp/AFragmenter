import pytest
from afragmenter import sequence_reader as sr
from pathlib import Path


def test_translate_three2one_iter():
    # Test valid sequences
    assert sr.translate_three2one_iter(["Ala"]) == "A"
    assert sr.translate_three2one_iter(["Ala", "Cys"]) == "AC"
    assert sr.translate_three2one_iter(["Ala", "Cys", "Asp"]) == "ACD"
    assert sr.translate_three2one_iter(["Ala", "Cys", "Asp", "Glu"]) == "ACDE"
    assert sr.translate_three2one_iter(["Ala", "Cys", "Asp", "Glu", "Phe"]) == "ACDEF"
    # Test case insensitivity
    assert sr.translate_three2one_iter(["ala"]) == "A"
    assert sr.translate_three2one_iter(["ALA"]) == "A"
    assert sr.translate_three2one_iter(["aLa"]) == "A"
    with pytest.raises(KeyError):
        # Invalid residue
        assert sr.translate_three2one_iter(["Cy"])
        assert sr.translate_three2one_iter(["C"])
        assert sr.translate_three2one_iter(["alacysasp"])
        # String passed instead of list
        assert sr.translate_three2one_iter("Ala") == "A"


def test_read_mmcif_sequence():
# TODO: find a really small protein structure file to test this and other functions
    P0A9R4_seq = "MPKIVILPHQDLCPDGAVLEANSGETILDAALRNGIEIEHACEKSCACTTCHCIVR"\
                 "EGFDSLPESSEQEDDMLDKAWGLEPESRLSCQARVTDEDLVVEIPRYTINHAREH"
    mmcif_file = Path(__file__).parent / "data" / "P0A9R4" / "P0A9R4.cif"
    # Test valid sequence
    assert sr._read_mmcif_sequence(mmcif_file) == P0A9R4_seq
    assert sr.read_mmcif_sequence(mmcif_file, len(P0A9R4_seq)) == P0A9R4_seq

    # Test invalid sequence length
    with pytest.raises(ValueError):
        sr.read_mmcif_sequence(mmcif_file, len(P0A9R4_seq) - 1)
        sr.read_mmcif_sequence(mmcif_file, len(P0A9R4_seq) + 1)
    
    # Test invalid file path
    with pytest.raises(FileNotFoundError):
        sr._read_mmcif_sequence("invalid.cif")
        sr.read_mmcif_sequence("invalid.cif", len(P0A9R4_seq))
    
    # Test invalid file type
    pdb_file = Path(__file__).parent / "data" / "P0A9R4" / "P0A9R4.pdb"
    empty_file = Path(__file__).parent / "data" / "empty.txt"
    with pytest.raises(ValueError):
        sr._read_mmcif_sequence(pdb_file)
        sr._read_mmcif_sequence(empty_file)
        sr.read_mmcif_sequence(pdb_file, len(P0A9R4_seq))
        sr.read_mmcif_sequence(empty_file, len(P0A9R4_seq))
    

@pytest.mark.filterwarnings("ignore")
def test_read_pdb_sequence():
    P0A9R4_seq = "MPKIVILPHQDLCPDGAVLEANSGETILDAALRNGIEIEHACEKSCACTTCHCIVR"\
                 "EGFDSLPESSEQEDDMLDKAWGLEPESRLSCQARVTDEDLVVEIPRYTINHAREH"
    
    pdb_file = Path(__file__).parent / "data" / "P0A9R4" / "P0A9R4.pdb"
    # Test valid sequence
    assert sr._read_pdb_sequence(pdb_file, "A") == P0A9R4_seq
    assert sr.read_pdb_sequence(pdb_file, "A", len(P0A9R4_seq)) == P0A9R4_seq

    # Test invalid sequence length
    with pytest.raises(ValueError):
        sr.read_pdb_sequence(pdb_file, "A", len(P0A9R4_seq) - 1)
        sr.read_pdb_sequence(pdb_file, "A", len(P0A9R4_seq) + 1)
    
    # Test invalid file path
    with pytest.raises(FileNotFoundError):
        sr._read_pdb_sequence("invalid.pdb", "A")
        sr.read_pdb_sequence("invalid.pdb", "A", len(P0A9R4_seq))
    
    # Test invalid file type
    # Reading the cif file with PDBParser will raise a warning:
    # PDBConstructionWarning: Ignoring unrecognized record '#' at line 1657
    mmcif_file = Path(__file__).parent / "data" / "P0A9R4" / "P0A9R4.cif"
    empty_file = Path(__file__).parent / "data" / "empty.txt"
    with pytest.raises(ValueError):
        sr._read_pdb_sequence(mmcif_file, "A")
        sr._read_pdb_sequence(empty_file, "A")
        sr.read_pdb_sequence(mmcif_file, "A", len(P0A9R4_seq))
        sr.read_pdb_sequence(empty_file, "A", len(P0A9R4_seq))
