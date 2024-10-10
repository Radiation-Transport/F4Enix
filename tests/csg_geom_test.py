"""
C --- LOWER SECTION ---
1 TZ 0 0 -352.5 514.5 158 158 $ inner torus (positive offset)
2 TZ 0 0 -261 512.5 249.5 249.5 $ outer torus (positive offset)
3 KZ 23185.875 0.00047776881961240957 -1 $ splitter cone (no offset)
C
C --- UPPER SECTION ---
11 TZ 0 0 356.3 516.0 159.5 159.5 $ inner torus (positive offset)
12 TZ 0 0 19.5   398.5 503.4 503.4 $ outer torus (positive offset)
13 KZ 906.5 2.096588 -1 $ (positive offset)
14 KZ -557.512037762645 0.47696543145338977 1 $ outer splitter cone (no offset)
15 KZ -390.8473312058339 0.47696543145338977 1 $ inner splitter cone (no offset)
C
C --- MID SECTION
21 TZ 0 0 79.5 588.30 304.4 304.4 $ eq1 (positive offset)
22 TZ 0 0 49.1 340.2 554.4 554.4 $ eq2 (positive offset)
23 TZ 0 0 15.2 688.45 204.4 204.4 $ eq3 (positive offset)
24 KZ -1980.9 0.209324 1 $ (negative offset)
25 CZ 356.5 $ inboard cylinder (negative offset)
C
C --- horizontal planes ---
101 PZ -352.5 $ lower plane (no offset)
102 PZ 356.3 $ (no offset)
C cones for mid-section
103 KZ -106.47471022128559 10.006677777777762 1 $ (no offset)
104 KZ 7.414873035066506 66.60492771814405 1 $ (no offset)
105 KZ 82.21638190954774 105.53168045875019 -1 $ (no offset)
106 KZ 330.17922989446475 4.777283063576085 -1 $ (no offset)
"""

from f4enix.input.csg_geom import Torus, Cone, Plane, Cylinder


class TestTorus:
    text = "1 TZ 0.0 0.0 -352.5 514.5 158.0 155.0 $ inner torus (positive offset)\n"

    def test_from_text(self):
        torus = Torus.from_text(self.text)
        assert torus.id_num == 1
        assert torus.geom == "TZ"
        assert torus.origin == [0, 0, -352.5]
        assert torus.major_radius == 514.5
        assert torus.minor_radius_A == 158
        assert torus.minor_radius_B == 155
        assert torus.comment == "$ inner torus (positive offset)\n"

    def test_to_text(self):
        torus = Torus.from_text(self.text)
        assert torus.to_text() == self.text

    def test_offset(self):
        torus = Torus.from_text(self.text)
        torus.offset(10)
        assert torus.minor_radius_A == 168
        assert torus.minor_radius_B == 165


class TestCone:
    text = "3 kz 23185.875 0.00047776881961240957 -1 $ splitter cone (no offset)\n"

    def test_from_text(self):
        cone = Cone.from_text(self.text)
        assert cone.id_num == 3
        assert cone.geom == "KZ"
        assert cone.vertex == 23185.875
        assert cone.t_square == 0.00047776881961240957
        assert cone.inclination == -1
        assert cone.comment == "$ splitter cone (no offset)\n"

    def test_to_text(self):
        cone = Cone.from_text(self.text)
        assert (
            cone.to_text()
            == "3 KZ 23185.875 0.00047776881961240957 -1.0 $ splitter cone (no offset)\n"
        )

    def test_offset(self):
        cone = Cone.from_text(self.text)
        cone.offset(10)
        assert cone.vertex == 23643.48427656681


class TestCylinder:
    text = "25 CZ 356.5 $ inboard cylinder (negative offset)\n"

    def test_from_text(self):
        cylinder = Cylinder.from_text(self.text)
        assert cylinder.id_num == 25
        assert cylinder.geom == "CZ"
        assert cylinder.radius == 356.5
        assert cylinder.comment == "$ inboard cylinder (negative offset)\n"

    def test_to_text(self):
        cylinder = Cylinder.from_text(self.text)
        assert cylinder.to_text() == self.text

    def test_offset(self):
        cylinder = Cylinder.from_text(self.text)
        cylinder.offset(-10)
        assert cylinder.radius == 346.5


class TestPlane:
    text = "101 PZ -352.5 $ lower plane (no offset)\n"

    def test_from_text(self):
        plane = Plane.from_text(self.text)
        assert plane.id_num == 101
        assert plane.geom == "PZ"
        assert plane.quota == -352.5
        assert plane.comment == "$ lower plane (no offset)\n"

    def test_to_text(self):
        plane = Plane.from_text(self.text)
        assert plane.to_text() == self.text

    def test_offset(self):
        plane = Plane.from_text(self.text)
        plane.offset(10)
        assert plane.quota == -342.5
