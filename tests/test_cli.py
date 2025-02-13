from afragmenter import cli
from click.testing import CliRunner
from pathlib import Path

def test_help():
    runner = CliRunner()
    result = runner.invoke(cli.afragmenter, ["--help"])
    assert result.exit_code == 0
    assert "Usage: afragmenter [OPTIONS]" in result.output

def test_version():
    runner = CliRunner()
    result = runner.invoke(cli.afragmenter, ["--version"])
    assert result.exit_code == 0
    assert "afragmenter, version" in result.output

def test_missing_required_args():
    runner = CliRunner()
    cif_file = Path(__file__).parent / "data" / "B1LFV6" / "B1LFV6.cif"
    result = runner.invoke(cli.afragmenter, ["--structure", cif_file])
    assert result.exit_code != 0
    
    either_required = ["--json", "--afdb"]
    either_required_message = "Either {} is required.".format(" or ".join(either_required))
    assert either_required_message in result.output

def test_invalid_json_path():
    runner = CliRunner()
    result = runner.invoke(cli.afragmenter, ["--json", "invalid.json"])
    assert result.exit_code != 0
    assert "Invalid value for '--json'" in result.output

def test_invalid_structure_path():
    runner = CliRunner()
    json_file = Path(__file__).parent / "data" / "B1LFV6" / "B1LFV6.json"
    result = runner.invoke(cli.afragmenter, ["--json", json_file, "--structure", "invalid.pdb"])
    assert result.exit_code != 0
    assert "Invalid value for '--structure' / '-s'" in result.output

def test_valid_arguments():
    runner = CliRunner()
    json_file = Path(__file__).parent / "data" / "B1LFV6" / "B1LFV6.json"
    with runner.isolated_filesystem():
        result = runner.invoke(cli.afragmenter, ["--json", json_file, "--n-iterations", "10"])
        assert result.exit_code == 0

def test_all_arguments():
    runner = CliRunner()
    with runner.isolated_filesystem():
        json_file = Path(__file__).parent / "data" / "B1LFV6" / "B1LFV6.json"
        structure_file = Path(__file__).parent / "data" / "B1LFV6" / "B1LFV6.cif"
        result = runner.invoke(cli.afragmenter, [
            "--json", json_file,
            "--structure", structure_file,
            "--afdb", "afdb",
            "--resolution", "0.6",
            "--objective-function", "CPM",
            "--n-iterations", "10",
            "--threshold", "5",
            "--min-size", "10",
            "--no-merge",
            "--plot-results", "results.png",
            "--output-fasta", "output.fasta"
        ])
        assert result.exit_code == 0
        assert Path("results.png").exists()
        assert Path("output.fasta").exists()

def test_invalid_argument():
    runner = CliRunner()
    json_file = Path(__file__).parent / "data" / "B1LFV6" / "B1LFV6.json"
    result = runner.invoke(cli.afragmenter, ["--json", json_file, "--invalid"])
    assert result.exit_code != 0
    assert "No such option" in result.output
    
