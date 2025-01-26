import pytest
from specdata import SpecDB


@pytest.fixture
def oh_ax():
    return SpecDB("OHAX.db")


def test_get_spectrum(oh_ax):
    spec = oh_ax.get_spectrum(Trot=3000, Tvib=3000, wmin=310, wmax=320)
    assert spec.y.sum() == pytest.approx(734360, rel=1)
    assert spec.y[:10].sum() == pytest.approx(11599, rel=1)
    assert spec.y[-10:].sum() == pytest.approx(1170, rel=1)
    assert spec.y[100:200].sum() == pytest.approx(198190, rel=1)
