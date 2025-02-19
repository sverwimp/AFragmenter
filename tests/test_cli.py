from afragmenter import cli
from click.testing import CliRunner
from pathlib import Path
import pytest

@pytest.fixture
def json_file():
    return Path(__file__).parent / "data" / "B1LFV6" / "B1LFV6.json"

@pytest.fixture
def cif_file():
    return Path(__file__).parent / "data" / "B1LFV6" / "B1LFV6.cif"

@pytest.fixture
def pdb_file():
    return Path(__file__).parent / "data" / "B1LFV6" / "B1LFV6.pdb"

def test_help():
    runner = CliRunner()
    result = runner.invoke(cli.main, ["--help"])
    assert result.exit_code == 0

def test_missing_required_args(cif_file):
    runner = CliRunner()
    result = runner.invoke(cli.main, ["--structure", cif_file])
    assert result.exit_code != 0
    
    either_required = ["--json", "--afdb"]
    either_required_message = "Either {} is required.".format(" or ".join(either_required))
    assert either_required_message in result.output

def test_invalid_json_path():
    runner = CliRunner()
    result = runner.invoke(cli.main, ["--json", "invalid.json"])
    assert result.exit_code != 0
    assert "Invalid value for '--json'" in result.output

def test_invalid_structure_path(json_file):
    runner = CliRunner()
    result = runner.invoke(cli.main, ["--json", json_file, "--structure", "invalid.pdb"])
    assert result.exit_code != 0
    assert "Invalid value for '--structure' / '-s'" in result.output

def test_valid_arguments(json_file):
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(cli.main, ["--json", json_file, "--n-iterations", "2"])
        assert result.exit_code == 0

def test_save_outputs_and_name(json_file, cif_file):
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(cli.main, ["--json", json_file, "--structure", cif_file,
                                                 "--save-result", "results.csv", "--plot-result", "results.png",
                                                 "--save-fasta", "output.fasta", "--n-iterations", "2",
                                                 "--name", "test_name"])
        assert result.exit_code == 0
        assert Path("results.csv").exists()
        assert Path("results.png").exists()
        assert Path("output.fasta").exists()
        with open("output.fasta") as f:
            assert f.read().startswith(">test_name")

def test_emtpy_name(json_file, cif_file):
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(cli.main, ["--json", json_file, "--structure", cif_file,
                                                 "--save-fasta", "output.fasta", "--n-iterations", "2",
                                                 "--name", ""])
        assert result.exit_code == 0
        assert Path("output.fasta").exists()
        with open("output.fasta") as f:
            assert f.read().startswith(">1")

def test_all_arguments(json_file, cif_file):
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(cli.main, [
            "--json", json_file,
            "--structure", cif_file,
            "--afdb", "B1LFV6",
            "--resolution", "0.6",
            "--objective-function", "CPM",
            "--n-iterations", "2",
            "--threshold", "5",
            "--min-size", "10",
            "--no-merge",
            "--name", "B1LFV6",
            "--save-result", "results.csv",
            "--plot-result", "results.png",
            "--save-fasta", "output.fasta"
        ])
        assert result.exit_code == 0
        assert Path("results.csv").exists()
        assert Path("results.png").exists()
        assert Path("output.fasta").exists()

def test_invalid_argument():
    runner = CliRunner()
    json_file = Path(__file__).parent / "data" / "B1LFV6" / "B1LFV6.json"
    result = runner.invoke(cli.main, ["--json", json_file, "--invalid"])
    assert result.exit_code != 0
    assert "No such option" in result.output
    
