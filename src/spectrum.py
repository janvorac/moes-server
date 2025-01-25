import logging
import warnings

import numpy as np
from scipy.signal import fftconvolve  # type: ignore [import-untyped]
from scipy.special import wofz  # type: ignore [import-untyped]


class Spectrum:
    """An object holding x and y axis and implememnting some methods
    often used in processing of spectroscopic data.
    """

    def __init__(
        self, x: np.typing.NDArray[np.float64], y: np.typing.NDArray[np.float64]
    ):
        """

        kwargs:
        -------
        x: 1D numpy array of wavelengths (or wavenumbers, if you prefffer)

        y: 1D numpy array of intensities (or fluxes or signal strengths, whatever)

        x and y should be equally long (though it is not strictly
        forbidden, most of the methods will fail otherwise)

        """
        self.x = x
        self.y = y

        if len(self.x) != len(self.y):
            raise ValueError("Spectrum: length of x and y mismatch!")

    def __len__(self):
        return len(self.x)

    def convolve_with_slit_function(
        self,
        gauss: float = 0.1,
        lorentz: float = 1e-9,
        instrumental_step: float | None = None,
    ):
        """
        Broaden the peaks in the spectrum by voigt profile and by a rectangle of given width.
        Changes state of the instance, returns nothing.

        **kwargs:
        ---------
        gauss: *float* gaussian HWHM, defaults to 0.1
        lorentz: *float* lorentzian HWHM, defaults to 1e-9
        step: *float* distance between pixels in nm

        return:
        -------
        None, modifies the spectrum in place
        """

        slit = voigt(self.x, gauss, lorentz, np.mean(self.x), 1.0)  # type: ignore [arg-type]
        slit /= np.sum(slit)

        numpoints = len(self.y)
        if instrumental_step is None:
            convolution_profile = slit[slit > max(slit) / 1000.0]
        else:
            # covolve with a rectangle of width equal to the instrumental step
            # to avoid losing thin lines
            simulated_step = self.x[1] - self.x[0]
            if instrumental_step / simulated_step < 1:
                msg = "Your simulated spectra resolution is more rough than experimental data."
                warnings.warn(msg, UserWarning)
                convolution_profile = slit[slit > max(slit) / 1000.0]
            else:
                instrumental_step_profile = np.ones(
                    int(instrumental_step / simulated_step) + 1
                )
                if len(slit) >= len(instrumental_step_profile):
                    convolution_profile_uncut = fftconvolve(
                        slit, instrumental_step_profile, mode="same"
                    )
                else:
                    convolution_profile_uncut = fftconvolve(
                        instrumental_step_profile, slit, mode="same"
                    )
                convolution_profile = convolution_profile_uncut[
                    convolution_profile_uncut > max(convolution_profile_uncut) / 1000
                ]

        self.y = fftconvolve(self.y, convolution_profile, mode="same")

        if len(self.y) == 0:
            self.y = (
                np.ones(numpoints) * 1e100
            )  # if the array gets destroyed by fftconvolve,
            # set array to ridiculously huge values
        if any(np.isnan(self.y)):
            self.y[:] = 1e100
        return

    def refine_mesh(self, points_per_nm: int = 3000):
        """
        adds artificial zeros in between lines. Usually used after creating
        a simulated spectrum before convolution with slit function.
        Necessary for later comparing simulation with measurement.

        return:
        Spectrum objects with pretty many points (or fine mesh, if you prefer)
        """

        start_spec = np.min(self.x) - 2  # prevent lines from falling to edges
        end_spec = np.max(self.x) + 2

        no_of_points = int(np.abs(end_spec - start_spec) * points_per_nm)

        spec = np.zeros((no_of_points, 2))

        spec[:, 0] = np.linspace(start_spec, end_spec, no_of_points)

        for i in range(len(self.x)):
            index = int((self.x[i] - start_spec) * points_per_nm + 0.5)
            spec[index, 1] += self.y[i]
        self.x = spec[:, 0]
        self.y = spec[:, 1]
        return spec


def match_spectra(sim_spec: Spectrum, exp_spec: Spectrum) -> tuple[Spectrum, Spectrum]:
    """
    Take two Spectrum objects with different x-axes
    (the ranges must partially overlap)
    and return a tuple of Spectrum objects defined at the same x-points.
    This enables comparing.
    Args:
    ----
    exp_spec: *Spectrum object*, experimental spectrum
    sim_spec: *Spectrum object*, simulated spectrum.

    Returns:
    (Spectrum_simulated, Spectrum_experimental): a tuple of Spectrum objects, with identical x-axes.
    """
    if len(sim_spec.x) == 0:
        return (Spectrum(x=exp_spec.x, y=np.zeros_like(exp_spec.x)), exp_spec)

    if np.min(exp_spec.x) < np.min(sim_spec.x):
        sim_spec.x = np.concatenate(
            [[min(exp_spec.x)], [min(sim_spec.x) - 1e-3], sim_spec.x]
        )
        sim_spec.y = np.concatenate([[0], [0], sim_spec.y])
    if np.max(exp_spec.x) > np.max(sim_spec.x):
        sim_spec.x = np.concatenate(
            [sim_spec.x, [max(sim_spec.x) + 1e-3], [max(exp_spec.x)]]
        )
        sim_spec.y = np.concatenate([sim_spec.y, [0], [0]])

    y_interp = np.interp(exp_spec.x, sim_spec.x, sim_spec.y)

    return (Spectrum(x=exp_spec.x, y=y_interp), exp_spec)


def compare_spectra(spectrum_exp, spectrum_sim):
    """
    spectrum_exp: Spectrum object
    spectrum_sim: Spectrum object
    The wavelength axes are expected to differ, but overlap

    returns: sqrt of (sum of squares of differences of the spectra divided
             by (the number of points)**2 )
    """
    matched = match_spectra(spectrum_sim, spectrum_exp)
    dif = matched[0].y - matched[1].y
    logging.debug("len(dif) = %i", len(dif))
    logging.debug("dif = ")
    logging.debug("%s", dif)
    return dif


def _voigt(x: np.typing.NDArray[np.float64], y: np.typing.NDArray[np.float64]):
    """
    Taken from `astro.rug.nl <http://www.astro.rug.nl/software/kapteyn-beta/kmpfittutorial.html?highlight=voigt#voigt-profiles/>`_

    The voigt function is also the real part of
    `w(z) = exp(-z^2) erfc(iz)`, the complex probability function,
    which is also known as the Faddeeva function. Scipy has
    implemented this function under the name `wofz()`
    """
    z = x + 1j * y
    I = wofz(z).real
    return I


def voigt(
    nu: np.typing.NDArray[np.float64],
    alphaD: float,
    alphaL: float,
    nu_0: float,
    A: float,
):
    """
    Taken from `astro.rug.nl <http://www.astro.rug.nl/software/kapteyn-beta/kmpfittutorial.html?highlight=voigt#voigt-profiles/>`_

    The voigt line shape in terms of its physical parameters

    Args:
      **nu**:  light frequency axis

      **alphaD**:  Doppler broadening HWHM

      **alphaL**:  Lorentzian broadening HWHM

      **nu_0**:  center of the line

      **A**:  integral under the line

    Returns:
      **V**: The voigt profile on the nu axis
    """
    if alphaD == 0:
        alphaD = 1e-10
    if alphaL == 0:
        alphaL = 1e-10
    f = np.sqrt(np.log(2))
    x = (nu - nu_0) / alphaD * f
    y = alphaL / alphaD * f
    V = A * f / (alphaD * np.sqrt(np.pi)) * _voigt(x, y)
    return V
