"""
Object description of CSG surfaces

"""

from __future__ import annotations

"""
Copyright 2024 F4E | European Joint Undertaking for ITER and the Development of
Fusion Energy (‘Fusion for Energy’). Licensed under the EUPL, Version 1.2 or - 
as soon they will be approved by the European Commission - subsequent versions
of the EUPL (the “Licence”). You may not use this work except in compliance
with the Licence. You may obtain a copy of the Licence at: 
    https://eupl.eu/1.2/en/ 
Unless required by applicable law or agreed to in writing, software distributed
under the Licence is distributed on an “AS IS” basis, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the Licence permissions
and limitations under the Licence.
"""

import re
import math

from f4enix.input.auxiliary import _get_comment

pat_space = re.compile(r"\s+")
pat_geom = re.compile(r"\d+\s+[A-Z][A-Z]")
pat_newline = re.compile(r"\s*\n")


class Torus:
    def __init__(
        self,
        id_num: int,
        geom: str,
        origin: list[float],
        major_radius: float,
        minor_radius_A: float,
        minor_radius_B: float,
        comment: str = None,
    ):
        """Represent a CSG definition of a torus

        Parameters
        ----------
        id_num : int
            number ID of the surface
        geom : str
            one of the MCNP torus allowed geomtries (e.g. TZ)
        origin : list[float]
            x, y, z coordinates of the origin of the torus in cm
        major_radius : float
            major radius of the torus in cm
        minor_radius_A : float
            minor A radius of the torus in cm
        minor_radius_B : float
            minor B radius of the torus in cm
        comment : str, optional
            dollar comment associated with the geometry, by default None

        Attributes
        ----------
        id_num : int
            number ID of the surface
        geom : str
            one of the MCNP torus allowed geomtries (e.g. TZ)
        origin : list[float]
            x, y, z coordinates of the origin of the torus in cm
        major_radius : float
            major radius of the torus in cm
        minor_radius_A : float
            minor A radius of the torus in cm
        minor_radius_B : float
            minor B radius of the torus in cm
        comment : str, optional
            dollar comment associated with the geometry, by default None

        """
        geom = geom.upper()
        self.id_num = id_num
        self.geom = geom
        self.origin = origin
        self.major_radius = major_radius
        self.minor_radius_A = minor_radius_A
        self.minor_radius_B = minor_radius_B
        self.comment = comment

    @classmethod
    def from_text(cls, text: str) -> Torus:
        """Create a Torus object from an MCNP CSG text description

        Parameters
        ----------
        text : str
            MCNP CSG text description

        Returns
        -------
        Torus
            CSG object representing the torus
        """
        pieces = pat_space.split(text)
        id_num = int(pieces[0])
        geom = pieces[1]
        origin = [float(x) for x in pieces[2:5]]
        major_radius = float(pieces[5])
        minor_radius_A = float(pieces[6])
        minor_radius_B = float(pieces[7])
        comment = _get_comment(text)

        return Torus(
            id_num, geom, origin, major_radius, minor_radius_A, minor_radius_B, comment
        )

    def __repr__(self) -> str:
        return self.to_text()

    def to_text(self) -> str:
        """print the CSG MCNP text format for the torus

        Returns
        -------
        str
            MCNP CSG text description of the torus
        """
        return f"{self.id_num} {self.geom} {self.origin[0]} {self.origin[1]} {self.origin[2]} {self.major_radius} {self.minor_radius_A} {self.minor_radius_B} {self.comment}"

    def offset(self, value: float) -> None:
        """Create an offset of the torus surface. This means varying the minor radii

        Parameters
        ----------
        value : float
            offset value in cm
        """
        # for a torus, an offset means changing the minor radius
        self.minor_radius_A += value
        self.minor_radius_B += value


class Cone:
    def __init__(
        self,
        id_num: int,
        geom: str,
        vertex: float,
        t_square: float,
        inclination: float,
        comment: str = None,
    ) -> None:
        """Represent a CSG definition of a cone

        Parameters
        ----------
        id_num : int
            number ID of the surface
        geom : str
            one of the MCNP cone allowed geometries (e.g. CZ). Only cones
            parallel to main axes are supported
        vertex : float
            the offset from origin of the cone vertex in cm
        t_square : float
            the square of the tangent of the half angle of the cone
        inclination : float
            inclination of the cone (see MCNP manual)
        comment : str, optional
            dollar comment associated with the surface, by default None
        """
        geom = geom.upper()
        if geom not in ["KX", "KY", "KZ"]:
            raise ValueError(
                f"Only cones parallel to the main axes are supported, not {geom}"
            )
        self.id_num = id_num
        self.geom = geom
        self.vertex = vertex
        self.t_square = t_square
        self.inclination = inclination
        self.comment = comment
        self.half_teta = math.atan(math.sqrt(t_square))

    @classmethod
    def from_text(cls, text: str) -> Cone:
        """Create a Cone object from an MCNP CSG text description

        Parameters
        ----------
        text : str
            MCNP CSG text description of a cone

        Returns
        -------
        Cone
            CSG object representing the cone

        Raises
        ------
        ValueError
            if the cone is not parallel to the main axes
        """
        # only parallel cones are supported
        pieces = pat_space.split(text)
        id_num = int(pieces[0])
        geom = pieces[1].upper()
        quota = float(pieces[2])
        t_square = float(pieces[3])
        inclination = float(pieces[4])
        comment = _get_comment(text)
        return cls(id_num, geom, quota, t_square, inclination, comment)

    def to_text(self) -> str:
        """print the CSG MCNP text format for the cone"""
        return f"{self.id_num} {self.geom} {self.vertex} {self.t_square} {self.inclination} {self.comment}"

    def __repr__(self) -> str:
        return self.to_text()

    def offset(self, value: float) -> None:
        """Create an offset of the cone surface. This means varying the cone vertex
        position on the axis

        Parameters
        ----------
        value : float
            offset value in cm
        """
        self.vertex += value / math.sin(self.half_teta)


class Cylinder:
    def __init__(
        self, id_num: int, geom: str, radius: float, comment: str = None
    ) -> None:
        """Represent a CSG definition of a cylinder

        Parameters
        ----------
        id_num : int
            number ID of the surface
        geom : str
            one of the MCNP cylinder allowed geometries. Only CZ are supported
            for the moment being
        radius : float
            radius of the cylinder in cm
        comment : str | None, optional
            dollar comment associated with the geometry, by default None
        """
        geom = geom.upper()
        if geom != "CZ":
            raise ValueError(f"Only CZ cylinders are supported, not {geom}")
        self.id_num = id_num
        self.geom = geom
        self.radius = radius
        self.comment = comment

    @classmethod
    def from_text(cls, text: str) -> Cylinder:
        """Create a Cylinder object from an MCNP CSG text description

        Parameters
        ----------
        text : str
            MCNP CSG text description of a cylinder

        Returns
        -------
        Cylinder
            CSG object representing the cylinder
        """
        # only CZ supported
        pieces = pat_space.split(text)
        id_num = int(pieces[0])
        geom = pieces[1].upper()
        radius = float(pieces[2])
        comment = _get_comment(text)
        return cls(id_num, geom, radius, comment)

    def to_text(self) -> str:
        """text representation of the cylinder in MCNP CSG format

        Returns
        -------
        str
            MCNP CSG text description of the cylinder
        """
        return f"{self.id_num} {self.geom} {self.radius} {self.comment}"

    def __repr__(self) -> str:
        return self.to_text()

    def offset(self, value: float) -> None:
        """Create an offset of the cylinder surface. This means varying the radius

        Parameters
        ----------
        value : float
            offset value in cm
        """
        self.radius += value


class Plane:
    def __init__(
        self, id_num: int, geom: str, quota: float, comment: str | None
    ) -> None:
        """Represent a CSG definition of a plane

        Parameters
        ----------
        id_num : int
            number ID of the surface
        geom : str
            one of the MCNP plane allowed geometries. Only PZ, PX, PY are supported
            for the moment being
        quota : float
            offset of the plane in cm with respect to the origin
        comment : str | None, optional
            dollar comment associated with the geometry, by default None

        Attributes
        ----------
        id_num : int
            number ID of the surface
        geom : str
            one of the MCNP plane allowed geometries. Only PZ are supported
            for the moment being
        quota : float
            offset of the plane in cm with respect to the origin
        comment : str | None, optional
            dollar comment associated with the geometry, by default None

        """
        geom = geom.upper()
        if geom not in ["PZ", "PX", "PY"]:
            raise ValueError(f"Only PZ planes are supported, not {geom}")
        self.id_num = id_num
        self.geom = geom
        self.quota = quota
        self.comment = comment

    @classmethod
    def from_text(cls, text: str) -> Plane:
        """Create a Plane object from an MCNP CSG text description

        Parameters
        ----------
        text : str
            MCNP CSG text description of a plane

        Returns
        -------
        Plane
            CSG object representing the plane
        """
        # Only planes parallel to the main axes are supported
        pieces = pat_space.split(text)
        id_num = int(pieces[0])
        geom = pieces[1]
        offset = float(pieces[2])
        comment = _get_comment(text)
        return cls(id_num, geom, offset, comment)

    def to_text(self) -> str:
        """text representation of the plane in MCNP CSG format"""
        return f"{self.id_num} {self.geom} {self.quota} {self.comment}"

    def __repr__(self) -> str:
        return self.to_text()

    def offset(self, value: float) -> None:
        """Create an offset of the plane surface. This means varying the quota

        Parameters
        ----------
        value : float
            offset value in cm
        """
        self.quota += value
