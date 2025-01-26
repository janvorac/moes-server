from enum import Enum


class SpectralSystemEnum(str, Enum):
    C2SWAN: str = "C2_swan"
    N2CB: str = "N2CB"
    N2PlusBX: str = "N2PlusBX"
    NHAX: str = "NHAX"
    NOBX: str = "NOBX"
    OHAX: str = "OHAX"
