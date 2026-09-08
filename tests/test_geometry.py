from wolf_elispot import ELISpotAnalyzer

def test_96_well_centers():
    analyzer = ELISpotAnalyzer()
    centers = analyzer.well_centers()
    assert len(centers) == 96
    assert "A1" in centers
    assert "H12" in centers
