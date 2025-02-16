import pytest
from afragmenter.format_result import format_csv, is_notebook, _build_rich_table



@pytest.fixture
def sample_intervals():
    return {
        0: [(0, 5), (10, 15)],
        1: [(20, 25), (30, 35)]
    }

def test_format_csv(sample_intervals):
    result = format_csv(sample_intervals)
    expected = "domain,nres,chopping\n1,12,1-6_11-16\n2,12,21-26_31-36\n"
    for r, e in zip(result.replace('\r', '').split('\n'), expected.split('\n')):
        assert r == e

def test_format_csv_delimiter(sample_intervals):
    result = format_csv(sample_intervals, delimiter=';')
    expected = "domain;nres;chopping\n1;12;1-6_11-16\n2;12;21-26_31-36\n"
    for r, e in zip(result.replace('\r', '').split('\n'), expected.split('\n')):
        assert r == e

def test_is_notebook():
    assert is_notebook() == False

def test_format_as_rich_table(sample_intervals):
    table = _build_rich_table(sample_intervals)
    assert table.row_count == 2

# I'm not sure how to test rich tables... so I'm going to skip it for now