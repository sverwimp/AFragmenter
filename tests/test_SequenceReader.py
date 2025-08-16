import pytest
from afragmenter import sequence_reader
from afragmenter.sequence_reader import SequenceReader

def test_translate_three2one_iter():
    """Test translating a list of three-letter amimo acids to one-letter amimo acids."""
    assert sequence_reader._translate_three2one_iter(['ALA', 'CYS', 'ASP']) == 'ACD'
    assert sequence_reader._translate_three2one_iter(['GLY', 'HIS', 'ILE']) == 'GHI'
    with pytest.raises(KeyError):
        sequence_reader._translate_three2one_iter(['XYZ'])


def test_get_content(tmp_path):
    """Test getting the content from a file."""
    content = "test content"
    file = tmp_path / "test.txt"
    with open(file, 'w') as f:
        f.write(content)

    assert sequence_reader._get_content(file) == content
    assert sequence_reader._get_content(content) == content


def test_determine_file_format():
    """Test determining the file format from the content."""
    fasta_content = ">test\nACDGH"
    pdb_content = "HEADER    TEST PDB\nATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00  0.00           N\n"
    mmcif_content = "data_test\n_entity_poly.pdbx_seq_one_letter_code\nA\n_entity_poly_seq.mon_id\nALA CYS ASP GLY HIS ILE\n"
    assert sequence_reader.determine_file_format(fasta_content) == 'fasta'
    assert sequence_reader.determine_file_format(pdb_content) == 'pdb'
    assert sequence_reader.determine_file_format(mmcif_content) == 'mmcif'
    with pytest.raises(ValueError, match="Unsupported file format, please provide a FASTA, PDB or mmCIF file"):
        sequence_reader.determine_file_format("invalid content")


def test_read_first_valid_line():
    """Test reading the first valid line from a file."""
    content = "# Comment\n\nValid line\n"
    assert sequence_reader._read_first_valid_line(content) == 'Valid line'
    assert sequence_reader._read_first_valid_line("# Only comments\n\n") is None


def test_read_fasta_sequence_content():
    """Test reading a FASTA file content."""
    fasta_content = ">test\nACDGH"
    reader = SequenceReader(fasta_content, format='FASTA')
    assert reader.sequence == 'ACDGH'
    assert reader.name == 'test'
    assert reader.format == 'fasta'
    assert reader.seq_length == 5


def test_read_fasta_sequence_file(tmp_path):
    """Test reading a FASTA file."""
    fasta_content = ">test\nACDGH"
    fasta_file = tmp_path / "test.fasta"
    with open(fasta_file, 'w') as f:
        f.write(fasta_content)

    reader = SequenceReader(fasta_file, format='FASTA')
    assert reader.sequence == 'ACDGH'
    assert reader.name == 'test'
    assert reader.format == 'fasta'
    assert reader.seq_length == 5


def test_read_fasta_sequence_seq_length():
    """Test reading a FASTA file with correct sequence length passed."""
    fasta_content = ">test\nACDGH"
    reader = SequenceReader(fasta_content, format='FASTA', seq_length=5)
    assert reader.sequence == 'ACDGH'
    assert reader.name == 'test'
    assert reader.format == 'fasta'
    assert reader.seq_length == 5


def test_read_fasta_sequence_wrong_seq_length():
    """Test reading a FASTA file with wrong sequence length passed."""
    fasta_content = ">test\nACDGH"
    with pytest.raises(ValueError) as excinfo:
        SequenceReader(fasta_content, format='FASTA', seq_length=4)
    assert "No sequence of length" in str(excinfo.value)
    assert "found in the FASTA content" in str(excinfo.value)


def test_read_fasta_sequence_empty():
    """Test reading an empty FASTA file."""
    fasta_content = ""
    with pytest.raises(ValueError, match="No sequence found in FASTA content"):
        SequenceReader(fasta_content, format='FASTA')


def test_read_pdb_sequence():
    """Test reading sequence from PDB file content."""
    pdb_content = "HEADER    TEST PDB\n" \
                  "COMPND   2 MOLECULE: TEST PDB;\n" \
                  "COMPND   3 CHAIN: A\n" \
                  "SOURCE    MOL_ID: 1;\n" \
                  "ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00  0.00           N\n" \
                  "ATOM      2  CA  ALA A   1      12.104  14.207  11.000  1.00  0.00           C\n" \
                  "ATOM      3  C   ALA A   1      13.104  15.207  12.000  1.00  0.00           C\n" \
                  "ATOM      4  O   ALA A   1      14.104  16.207  13.000  1.00  0.00           O\n" \
                  "ATOM      5  CB  ALA A   1      15.104  17.207  14.000  1.00  0.00           C\n" \
                  "ATOM      6  N   CYS A   2      16.104  18.207  15.000  1.00  0.00           N\n"
    reader = SequenceReader(pdb_content, format='PDB')
    assert reader.sequence == 'AC'
    assert reader.name == 'test_pdb'


def test_read_pdb_sequence_file(tmp_path):
    """Test reading sequence from PDB file."""
    pdb_content = "HEADER    TEST PDB\n" \
                  "COMPND   2 MOLECULE: TEST PDB;\n" \
                  "COMPND   3 CHAIN: A\n" \
                  "SOURCE    MOL_ID: 1;\n" \
                  "ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00  0.00           N\n" \
                  "ATOM      2  CA  ALA A   1      12.104  14.207  11.000  1.00  0.00           C\n" \
                  "ATOM      3  C   ALA A   1      13.104  15.207  12.000  1.00  0.00           C\n" \
                  "ATOM      4  O   ALA A   1      14.104  16.207  13.000  1.00  0.00           O\n" \
                  "ATOM      5  CB  ALA A   1      15.104  17.207  14.000  1.00  0.00           C\n" \
                  "ATOM      6  N   CYS A   2      16.104  18.207  15.000  1.00  0.00           N\n"
    pdb_file = tmp_path / "test.pdb"
    with open(pdb_file, 'w') as f:
        f.write(pdb_content)

    reader = SequenceReader(pdb_file, format='PDB')
    assert reader.sequence == 'AC'
    assert reader.name == 'test_pdb'
    

def test_read_pdb_sequence_seq_length():
    """Test reading sequence from PDB file with correct sequence length passed."""
    pdb_content = "HEADER    TEST PDB\n" \
                  "COMPND   2 MOLECULE: TEST PDB;\n" \
                  "COMPND   3 CHAIN: A\n" \
                  "SOURCE    MOL_ID: 1;\n" \
                  "ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00  0.00           N\n" \
                  "ATOM      2  CA  ALA A   1      12.104  14.207  11.000  1.00  0.00           C\n" \
                  "ATOM      3  C   ALA A   1      13.104  15.207  12.000  1.00  0.00           C\n" \
                  "ATOM      4  O   ALA A   1      14.104  16.207  13.000  1.00  0.00           O\n" \
                  "ATOM      5  CB  ALA A   1      15.104  17.207  14.000  1.00  0.00           C\n" \
                  "ATOM      6  N   CYS A   2      16.104  18.207  15.000  1.00  0.00           N\n"
    reader = SequenceReader(pdb_content, format='PDB', seq_length=2)
    assert reader.sequence == 'AC'
    assert reader.name == 'test_pdb'


def test_read_pdb_sequence_seq_wrong_length():
    """Test reading a PDB file with wrong sequence length passed."""
    pdb_content = "HEADER    TEST PDB\n" \
                  "COMPND   2 MOLECULE: TEST PDB;\n" \
                  "COMPND   3 CHAIN: A\n" \
                  "SOURCE    MOL_ID: 1;\n" \
                  "ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00  0.00           N\n" \
                  "ATOM      2  CA  ALA A   1      12.104  14.207  11.000  1.00  0.00           C\n" \
                  "ATOM      3  C   ALA A   1      13.104  15.207  12.000  1.00  0.00           C\n" \
                  "ATOM      4  O   ALA A   1      14.104  16.207  13.000  1.00  0.00           O\n" \
                  "ATOM      5  CB  ALA A   1      15.104  17.207  14.000  1.00  0.00           C\n" \
                  "ATOM      6  N   CYS A   2      16.104  18.207  15.000  1.00  0.00           N\n"
    with pytest.raises(ValueError) as excinfo:
        SequenceReader(pdb_content, format='PDB', seq_length=12)
    assert "Sequence length mismatch" in str(excinfo.value)


def test_read_pdb_sequence_empty(tmp_path):
    """Test reading an empty PDB file."""
    pdb_content = ""
    pdb_file = tmp_path / "test.pdb"
    with open(pdb_file, 'w') as f:
        f.write(pdb_content)
    with pytest.raises(ValueError, match="Unable to infer file format, file is empty or contains only comments"):
        SequenceReader(pdb_file, format='auto')


def test_read_pdb_sequence_no_chain_b():
    """Test reading a PDB file withouth the asked chain."""
    pdb_content = "HEADER    TEST PDB\n" \
                  "COMPND   2 MOLECULE: TEST PDB;\n" \
                  "COMPND   3 CHAIN: A\n" \
                  "SOURCE    MOL_ID: 1;\n" \
                  "ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00  0.00           N\n" \
                  "ATOM      2  CA  ALA A   1      12.104  14.207  11.000  1.00  0.00           C\n" \
                  "ATOM      3  C   ALA A   1      13.104  15.207  12.000  1.00  0.00           C\n" \
                  "ATOM      4  O   ALA A   1      14.104  16.207  13.000  1.00  0.00           O\n" \
                  "ATOM      5  CB  ALA A   1      15.104  17.207  14.000  1.00  0.00           C\n" \
                  "ATOM      6  N   CYS A   2      16.104  18.207  15.000  1.00  0.00           N\n"
    with pytest.raises(ValueError, match="Chain B not found in PDB file."):
        SequenceReader(pdb_content, format='PDB', query_chain='B')


def test_read_mmcif_sequence_content():
    """Test reading sequence from mmCIF file content. using _entity_poly.pdbx_seq_one_letter_code"""
    mmcif_content = "data_testcase\n" \
                    "_entry.id testcase\n" \
                    "_entity_poly.pdbx_seq_one_letter_code     \n" \
                    ";ACDGHI\n" \
                    ";\n" \
                    "#\n" \
                    "loop_\n" \
                    "_entity_poly_seq.entity_id\n" \
                    "_entity_poly_seq.hetero\n" \
                    "_entity_poly_seq.mon_id\n" \
                    "_entity_poly_seq.num\n" \
                    "1 n ALA 1   \n" \
                    "1 n CYS 2   \n" \
                    "1 n ASP 3   \n" \
                    "1 n GLY 4   \n" \
                    "1 n HIS 5   \n" \
                    "#"  # If entity_poly_sea.mon_id is used, the sequence is incomplete
    reader = SequenceReader(mmcif_content, format='mmCIF')
    assert reader.sequence == 'ACDGHI'
    assert reader.name == 'testcase'


def test_read_mmcif_sequence_file(tmp_path):
    """Test reading sequence from mmCIF file. using _entity_poly_seq.mon_id"""
    mmcif_content = "data_testcase\n" \
                    "_entry.id testcase\n" \
                    "_entity_poly.pdbx_seq_one_letter_code     \n" \
                    ";ACDGHI\n" \
                    ";\n" \
                    "#\n" \
                    "loop_\n" \
                    "_entity_poly_seq.entity_id\n" \
                    "_entity_poly_seq.hetero\n" \
                    "_entity_poly_seq.mon_id\n" \
                    "_entity_poly_seq.num\n" \
                    "1 n ALA 1   \n" \
                    "1 n CYS 2   \n" \
                    "1 n ASP 3   \n" \
                    "1 n GLY 4   \n" \
                    "1 n HIS 5   \n" \
                    "#"  # If entity_poly_sea.mon_id is used, the sequence is incomplete
    
    mmcif_file = tmp_path / "testcase.cif"
    with open(mmcif_file, 'w') as f:
        f.write(mmcif_content)
    
    reader = SequenceReader(mmcif_file, format='mmCIF')
    assert reader.sequence == 'ACDGHI'
    assert reader.name == 'testcase'


def test_read_mmcif_sequence_mon_id_content():
    """Test reading sequence from mmCIF file content. using _entity_poly_seq.mon_id"""
    mmcif_content = "data_testcase2\n" \
                    "_entry.id testcase2\n" \
                    "loop_\n" \
                    "_entity_poly_seq.entity_id\n" \
                    "_entity_poly_seq.hetero\n" \
                    "_entity_poly_seq.mon_id\n" \
                    "_entity_poly_seq.num\n" \
                    "1 n ALA 1   \n" \
                    "1 n CYS 2   \n" \
                    "1 n ASP 3   \n" \
                    "1 n GLY 4   \n" \
                    "1 n HIS 5   \n" \
                    "1 n ILE 6   \n" \
                    "#\n"
    reader = SequenceReader(mmcif_content, format='mmCIF')
    assert reader.sequence == 'ACDGHI'
    assert reader.name == 'testcase2'


def test_read_mmcif_sequence_mon_id_file(tmp_path):
    """Test reading sequence from mmCIF file content. using _entity_poly_seq.mon_id"""
    mmcif_content = "data_testcase2\n" \
                    "_entry.id testcase2\n" \
                    "loop_\n" \
                    "_entity_poly_seq.entity_id\n" \
                    "_entity_poly_seq.hetero\n" \
                    "_entity_poly_seq.mon_id\n" \
                    "_entity_poly_seq.num\n" \
                    "1 n ALA 1   \n" \
                    "1 n CYS 2   \n" \
                    "1 n ASP 3   \n" \
                    "1 n GLY 4   \n" \
                    "1 n HIS 5   \n" \
                    "1 n ILE 6   \n" \
                    "#\n"
    mmcif_file = tmp_path / "testcase2.cif"
    with open(mmcif_file, 'w') as f:
        f.write(mmcif_content)
    
    reader = SequenceReader(mmcif_file, format='mmCIF')
    assert reader.sequence == 'ACDGHI'
    assert reader.name == 'testcase2'   


def test_read_mmcif_sequence_seq_length():
    """Test reading sequence from mmCIF file with correct sequence length passed."""
    mmcif_content = "data_testcase\n" \
                    "_entry.id testcase\n" \
                    "_entity_poly.pdbx_seq_one_letter_code     \n" \
                    ";ACDGHI\n" \
                    ";\n" \
                    "#\n" \
                    "loop_\n" \
                    "_entity_poly_seq.entity_id\n" \
                    "_entity_poly_seq.hetero\n" \
                    "_entity_poly_seq.mon_id\n" \
                    "_entity_poly_seq.num\n" \
                    "1 n ALA 1   \n" \
                    "1 n CYS 2   \n" \
                    "1 n ASP 3   \n" \
                    "1 n GLY 4   \n" \
                    "1 n HIS 5   \n" \
                    "1 n ILE 6   \n"  
    reader = SequenceReader(mmcif_content, format='mmCIF', seq_length=6)
    assert reader.sequence == 'ACDGHI'
    assert reader.name == 'testcase'


def test_read_mmcif_sequence_seq_wrong_length():
    """Test reading sequence from mmCIF file with wrong sequence length passed."""
    mmcif_content = "data_testcase\n" \
                    "_entry.id testcase\n" \
                    "_entity_poly.pdbx_seq_one_letter_code     \n" \
                    ";ACDGHI\n" \
                    ";\n" \
                    "#\n" \
                    "loop_\n" \
                    "_entity_poly_seq.entity_id\n" \
                    "_entity_poly_seq.hetero\n" \
                    "_entity_poly_seq.mon_id\n" \
                    "_entity_poly_seq.num\n" \
                    "1 n ALA 1   \n" \
                    "1 n CYS 2   \n" \
                    "1 n ASP 3   \n" \
                    "1 n GLY 4   \n" \
                    "1 n HIS 5   \n" \
                    "1 n ILE 6   \n"  
    with pytest.raises(ValueError) as excinfo:
        SequenceReader(mmcif_content, format='mmCIF', seq_length=7)
    assert "Sequence length mismatch" in str(excinfo.value)


def test_read_mmcif_sequence_empty(tmp_path):
    """Test reading an empty mmCIF file."""
    mmcif_content = ""
    mmcif_file = tmp_path / "test.cif"
    with open(mmcif_file, 'w') as f:
        f.write(mmcif_content)
    with pytest.raises(ValueError, match="Unable to infer file format, file is empty or contains only comments"):
        SequenceReader(mmcif_file, format='auto')


def test_read_mmcif_sequence_no_sequence():
    """Test reading mmCIF file with not sequence."""
    mmcif_content = "data_testcase\n" \
                    "_entry.id testcase\n"
    with pytest.raises(ValueError, match="Could not find the sequence in the mmCIF dictionary."):
        SequenceReader(mmcif_content, format='mmCIF')
                    
     
def test_clusters_to_fasta():
    """Test parsing a sequence with cluster intervals."""
    fasta_content = ">test\nABCDEFG"
    reader = SequenceReader(fasta_content, format='FASTA')
    cluster_intervals = {0: [(0, 2), (5, 6)], 1: [(3, 4)]}
    parsed_sequences = reader.clusters_to_fasta(cluster_intervals)
    assert parsed_sequences == {0: 'ABC--FG', 1: 'DE'}


def test_clusters_to_fasta_empty():
    """Test an empty cluster interval."""
    fasta_content = ">test\nABCDEFG"
    reader = SequenceReader(fasta_content, format='FASTA')
    cluster_intervals = {}
    parsed_sequences = reader.clusters_to_fasta(cluster_intervals)
    assert parsed_sequences == {}


def test_clusters_to_fasta_wrong_interval():
    """Test an interval outside the sequence length.
    By default, python does not raise an error when slicing is out of range."""
    fasta_content = ">test\nABCDEFG"
    reader = SequenceReader(fasta_content, format='FASTA')
    cluster_intervals = {0: [(0, 2), (5, 16)], 1: [(3, 4)]}
    parsed_sequences = reader.clusters_to_fasta(cluster_intervals)
    assert parsed_sequences == {0: 'ABC--FG', 1: 'DE'}
    


def test_get_name_cif():
    """Test getting the name from a mmCIF dictionary."""
    mmcif_dict = {'_entry.id': ['test']}
    assert sequence_reader._get_name_cif(mmcif_dict) == 'test'


def test_get_name_pdb():
    """Test getting the name from a PDB structure."""
    class MockStructure:
        header = {'compound': {'1': {'molecule': 'test_molecule'}}, 'name': 'test_name'}
    assert sequence_reader._get_name_pdb(MockStructure()) == 'test_molecule'
    MockStructure.header = {'name': 'test_name'}
    assert sequence_reader._get_name_pdb(MockStructure()) == 'test_name'


def test_get_name_fasta():
    """Test getting the name from a FASTA record."""
    class MockRecord:
        id = 'test_id'
    assert sequence_reader._get_name_fasta(MockRecord()) == 'test_id'


def test_clean_name():
    """Test cleaning the name."""
    assert sequence_reader._clean_name('AF-test-F1') == 'test'
    assert sequence_reader._clean_name('test name') == 'test_name'
    assert sequence_reader._clean_name('') == ''