from typing import Dict

from pydantic import BaseModel

from data_models.enums import SpectralSystemEnum


class SpectralSystemParams(BaseModel):
    Tvib: float
    Trot: float
    intensity: float


class SpectralSystemParamFitSettings(BaseModel):
    Tvib_min: float | None
    Tvib_max: float | None
    Tvib_fixed: bool

    Trot_min: float | None
    Trot_max: float | None
    Trot_fixed: bool

    intensity_min: float | None
    intensity_max: float | None
    intensity_fixed: bool


class SimulatedSpectrumParams(BaseModel):
    slitf_gauss: float
    slitf_lorentz: float
    baseline: float
    baseline_slope: float
    species_params: Dict[SpectralSystemEnum, SpectralSystemParams]


class SimulatedSpectrumParamFitSettings(BaseModel):
    slitf_gauss_min: float | None
    slitf_gauss_max: float | None
    slitf_gauss_fixed: bool

    slitf_lorentz_min: float | None
    slitf_lorentz_max: float | None
    slitf_lorentz_fixed: bool

    baseline_min: float | None
    baseline_max: float | None
    baseline_fixed: bool

    baseline_slope_min: float | None
    baseline_slope_max: float | None
    baseline_slope_fixed: bool


class MeasuredSpectrumParams(BaseModel):
    wav_start: float
    wav_step: float
    wav_2nd: float


class MeasuredSpectrumParamFitSettings(BaseModel):
    wav_start_min: float | None
    wav_start_max: float | None
    wav_start_fixed: bool

    wav_step_min: float | None
    wav_step_max: float | None
    wav_step_fixed: bool

    wav_2nd_min: float | None
    wav_2nd_max: float | None
    wav_2nd_fixed: bool
