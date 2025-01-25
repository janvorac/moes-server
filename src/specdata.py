import pathlib
import sqlite3 as sqlite
import warnings
from copy import copy
from typing import TYPE_CHECKING, Any, Literal, cast

import numpy
import pandas as pd

import spectrum
from data_models.enums import SpectralSystemEnum
from data_models.models import SimulatedSpectrumParams, SpectralSystemParams

DATA_DIR = pathlib.Path(__file__).parent / "data"


BOLTZMANN_CONSTANT = 69.50348004861274  # 1/(mK)
kB = BOLTZMANN_CONSTANT / 100  # inverse cm


WAV_RESERVE = 2


class DatabaseError(Exception):
    pass


class SpecDB:
    """
    Class for working with spectral databases, using pandas
    """

    def __init__(self, filename: str):
        """
        args:
        -----
        filename: name of the database file. Just the file, not the path. The files
        will be ALWAYS searched for exclusively in lighteroes/data
        directory. Providing full path will result in error.

        """

        self.specie_name = filename.replace(".db", "")
        self.filename = filename
        self.uorl = "upper"  # default, fuck off lower

        to_open = DATA_DIR / filename

        if not SpecDB.isSQLite3(to_open):
            raise DatabaseError(f"{to_open} is not a valid sqllite database!")

        self.conn = sqlite.connect(to_open)

        self.states = pd.read_sql_query("select J,E_J,E_v from upper_states", self.conn)
        self.last_Trot: float | None = None
        self.last_Tvib: float | None = None
        self.norm: float | None = None
        self.spec: spectrum.Spectrum | numpy.ndarray | None = None
        self.last_wmin: float = 0
        self.last_wmax: float = numpy.inf
        self.table: pd.DataFrame | None = None

    @staticmethod
    def isSQLite3(filename: pathlib.Path) -> bool:
        from os.path import getsize, isfile

        if not isfile(filename):
            return False
        if getsize(filename) < 100:  # SQLite database file header is 100 bytes
            return False

        with open(filename, "rb") as fd:
            header = fd.read(100)
        return header[:16] == b"SQLite format 3\x00"

    def __getstate__(self) -> str:
        return self.filename

    def __setstate__(self, filename: str) -> "SpecDB":
        return self.__init__(filename)  # type: ignore[misc]

    def calculate_norm(
        self,
        Trot: float,
        Tvib: float,
    ) -> float:
        parts = (2 * self.states.J + 1) * numpy.exp(
            -self.states.E_J / (kB * Trot) - self.states.E_v / (kB * Tvib)
        )
        return numpy.sum(parts)

    def get_spectrum(
        self,
        Trot: float,
        Tvib: float,
        wmin: float,
        wmax: float,
        y_scaling: Literal["intensity", "photon_flux"] = "photon_flux",
        refractive_index: Literal["vacuum", "air"] = "air",
    ) -> spectrum.Spectrum:
        """
        kwargs:
           wavelength: either 'vacuum' or 'air', returns the wavelength in vacuum or air (default is 'air')

           as_spectrum: bool, if True object Spectrum is returned, otherwise a numpy array is returned

            y: either 'photon_flux' (default), i.e. the y*axis is population * (emission coefficient)
              or 'intensity', i.e. the y-axis is population * (emission coefficient) * wavenumber
        """

        wav = cast(
            Literal["vacuum_wavelength", "air_wavelength"],
            refractive_index + "_wavelength",
        )
        recalculate_pops = False

        if wmin < self.last_wmin or wmax > self.last_wmax or self.table is None:

            self.last_wmin = wmin - WAV_RESERVE
            self.last_wmax = wmax + WAV_RESERVE
            self.table = self.get_table_from_DB(self.last_wmin, self.last_wmax, wav=wav)
            recalculate_pops = True

        if self.last_Trot != Trot or self.last_Tvib != Tvib or recalculate_pops:
            self.norm = self.calculate_norm(Trot, Tvib)
            self.last_Trot = Trot
            self.last_Tvib = Tvib
            self.table["pops"] = (
                (2 * self.table["J"] + 1)
                * numpy.exp(
                    -self.table["E_v"] / (kB * Tvib) - self.table["E_J"] / (kB * Trot)  # type: ignore[operator]
                )
                / self.norm
            )

        self.table["y"] = self.table["pops"] * self.table["A"]

        if y_scaling == "intensity":
            self.table["y"] *= self.table["wavenumber"]
        self.spec = spectrum.Spectrum(x=self.table[wav], y=self.table["y"])
        return copy(self.spec)

    def get_table_from_DB(
        self,
        wmin=None,
        wmax=None,
        wav: Literal["air_wavelength", "vacuum_wavelength"] = "air_wavelength",
    ) -> pd.DataFrame:
        """ """
        q = "SELECT air_wavelength, vacuum_wavelength, A, J, E_J, E_v, wavenumber"
        q += " FROM "
        if wmin != 0 and wmax != numpy.inf:
            q += " (SELECT * FROM lines WHERE " + wav + " BETWEEN ? AND ?)"
            params = (wmin, wmax)
        else:
            q += " lines "
            params = ()  # type: ignore[assignment]
        q += " INNER JOIN "
        q += " upper_states on upper_state=upper_states.id"
        q += " ORDER BY " + wav
        table = pd.read_sql_query(q, self.conn, params=params)
        return table

    def get_lines_by_states(
        self,
        wmin: float,
        wmax: float,
        minlines: int = 1,
        refractive_index: Literal["air", "vacuum"] = "air",
        max_J: float | None = None,
        max_v: float | None = None,
        singlet_like: bool = False,
    ) -> dict[str, Any]:
        wav = cast(
            Literal["vacuum_wavelength", "air_wavelength"],
            refractive_index + "_wavelength",
        )

        q = "select "
        q += "lines.id, A, " + wav + ", upper_state, branch, wavenumber, lower_state, "
        q += "E_J, J, component, E_v, v"
        q += (
            " from lines inner join "
            + self.uorl
            + "_states on "
            + self.uorl
            + "_state="
            + self.uorl
            + "_states.id"
        )
        q += " where lines." + wav + " between ? and ?"
        params = [wmin, wmax]
        if max_v is not None:
            q += " and v <= ?"
            params.append(max_v)
        if max_J is not None:
            q += " and J <= ?"
            params.append(max_J)

        big_table = pd.read_sql_query(q, self.conn, params=params)

        if singlet_like:
            gr = big_table.groupby(["v", "J"])
            states = big_table.loc[:, ["E_J", "J", "E_v", "v"]].drop_duplicates()
            states = states.groupby(["v", "J"]).mean()
        else:
            gr = big_table.groupby([self.uorl + "_state"])
            states = big_table.loc[
                :, [self.uorl + "_state", "E_J", "J", "component", "E_v", "v"]
            ].drop_duplicates()
            states.set_index(self.uorl + "_state", inplace=True)

        specs = []
        numlines = []
        state_names = []
        for name, group in gr:
            if len(group) >= minlines:
                numlines.append(len(group))
                specs.append(group.loc[:, [wav, "A"]])
                state_names.append(name)
        if singlet_like:
            states_ordered = states.loc[state_names].reset_index()
        else:
            states_ordered = states.loc[state_names, :]  # reorder
        print(len(numlines), len(states_ordered))
        states_ordered["numlines"] = numlines

        return {
            "specs": specs,
            "states": states_ordered.loc[states_ordered["numlines"] >= minlines],
        }


def generate_spectrum(
    simspec_params: SimulatedSpectrumParams,
    step: float,
    wmin: float,
    wmax: float,
    simulations: dict[SpectralSystemEnum, SpecDB],
    points_per_nm: int = 1000,
) -> spectrum.Spectrum:

    if (len(simulations) == 0) or (len(simspec_params.species_params) == 0):
        warnings.warn("No simulation files given, returning empty spectrum!", Warning)
        return spectrum.Spectrum(x=[], y=[])

    spectra_xs = []
    spectra_ys = []
    for species, species_params in simspec_params.species_params.items():
        sim = simulations[species]
        temp_spec = sim.get_spectrum(
            Trot=species_params.Trot,
            Tvib=species_params.Tvib,
            wmin=wmin,
            wmax=wmax,
        )
        temp_spec.y *= species_params.intensity
        spectra_xs.append(temp_spec.x)
        spectra_ys.append(temp_spec.y)

    spec_x = numpy.concatenate(spectra_xs)
    spec_y = numpy.concatenate(spectra_ys)
    spec_y = spec_y[spec_x.argsort()]
    spec_x = numpy.sort(spec_x)
    spec = spectrum.Spectrum(x=spec_x, y=spec_y)
    spec.refine_mesh(points_per_nm=points_per_nm)
    spec.convolve_with_slit_function(
        gauss=simspec_params.slitf_gauss,
        lorentz=simspec_params.slitf_lorentz,
        instrumental_step=step,
    )
    if len(spec.y) > 0:
        spec.y += simspec_params.baseline
        spec.y += simspec_params.baseline_slope * (spec.x - wmin)

    return spec


if __name__ == "__main__":
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

    test_spectrum = generate_spectrum(
        simspec_params=params,
        step=0.0318,
        wmin=300,
        wmax=320,
        simulations={
            "OHAX": SpecDB("OHAX.db"),
        },
    )
    with open("test_spectrum.txt", "w") as f:
        for x, y in zip(test_spectrum.x, test_spectrum.y):
            f.write(f"{x} {y}\n")
