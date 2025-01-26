import pytest

from data_models.enums import SpectralSystemEnum
from data_models.models import SimulatedSpectrumParams, SpectralSystemParams
from specdata import SpecDB, generate_spectrum


@pytest.fixture
def oh_ax():
    return SpecDB("OHAX.db")


def test_get_spectrum(oh_ax):
    spec = oh_ax.get_spectrum(Trot=3000, Tvib=3000, wmin=310, wmax=320)
    assert spec.y.sum() == pytest.approx(734360, rel=1)
    assert spec.y[:10].sum() == pytest.approx(11599, rel=1)
    assert spec.y[-10:].sum() == pytest.approx(1170, rel=1)
    assert spec.y[100:200].sum() == pytest.approx(198190, rel=1)


def test_generate_spectrum(oh_ax):
    params = SimulatedSpectrumParams(
        slitf_gauss=0.0396,
        slitf_lorentz=0.0138,
        baseline=0.55,
        baseline_slope=0,
        species_params={
            SpectralSystemEnum.OHAX: SpectralSystemParams(
                Tvib=3000,
                Trot=3000,
                intensity=1,
            )
        },
    )

    spec = generate_spectrum(
        simspec_params=params,
        step=0.0318,
        wmin=310,
        wmax=320,
        simulations={
            "OHAX": SpecDB("OHAX.db"),
        },
    )

    assert spec.y.sum() == pytest.approx(23244173, rel=1)
    assert spec.y[:10].sum() == pytest.approx(11599, rel=1)
    assert spec.y[-10:].sum() == pytest.approx(1170, rel=1)
    assert spec.y[100:200].sum() == pytest.approx(198190, rel=1)
