#!/usr/bin/env python3
# Copyright (C) 2021, Miklos Maroti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from typing import Any
import csv
import math
import random
import re
import sys
from freecad_scripts.libs import FreeCAD, Units, GmshTools, FemToolsCcx


class PressureVessel(object):
    """
    The base class to work with parametric pressure vessel models.
    """

    def __init__(self, filename: str, debug=True):
        """
        Creates a pressure vessel analysis class that can be used to run
        multiple simulations for the given design template by changing its
        parameters.
        """
        self.filename = filename
        self.debug = debug

        print("Opening:", filename)
        self.doc = FreeCAD.open(filename)

        if self.doc.getObject('FEMMeshGmsh').FemMesh.NodeCount:
            print("WARNING: clean the mesh in the model to save space")
        if self.doc.getObject('CCX_Results'):
            print("WARNING: remove the CCX results from the model to save space")

        self.sketch_params = []
        obj = self.doc.getObject('Sketch')
        for c in obj.Constraints:
            if c.Name:
                self.sketch_params.append(str(c.Name))

    def print_info(self):
        """
        Prints out all relevant information from the design template
        and the output of design analysis.
        """
        names = [obj.Name for obj in self.doc.Objects]
        print("Object names:", ", ".join(names))

        obj = self.doc.getObject('Sketch')
        print("Sketch parameters:")
        for c in obj.Constraints:
            if not c.Name:
                continue
            elif 'angle' in c.Name:
                print("  {} = {} rad".format(
                    c.Name, self.sketch_get_angle(c.Name)))
            else:
                print("  {} = {} m".format(
                    c.Name, self.sketch_get_length(c.Name)))

        print("Body properties:")
        print("  body_area = {:.6f} m^2".format(self.get_body_area()))
        print("  body_volume = {:.9f} m^3".format(self.get_body_volume()))
        print("  body_mass = {:.3f} kg".format(self.get_body_mass()))
        print("  outer_length = {:.9f} m".format(self.get_outer_length()))
        print("  outer_diameter = {:.9f} m".format(self.get_outer_diameter()))
        print("  outer_area = {:.6f} m^2".format(self.get_outer_area()))
        print("  outer_volume = {:.9f} m^3".format(self.get_outer_volume()))
        print("  inner_area = {:.6f} m^2".format(self.get_inner_area()))
        print("  inner_volume = {:.9f} m^3".format(self.get_inner_volume()))

        print("FEM parameters:")
        print("  test_pressure =", self.get_test_pressure(), "MPa")
        print("  mesh_length =", self.get_mesh_length(), "m")

        print("Material parameters:")
        print("  youngs_modulus =", self.get_youngs_modulus(), "MPa")
        print("  poisson_ratio =", self.get_poisson_ratio())
        print("  tensile_strength =", self.get_tensile_strength(), "MPa")
        print("  density =", self.get_density(), "kg/m^3")

        if self.has_mesh_properties():
            print("Mesh properties:")
            print("  node_count =", self.get_node_count())
            print("  edge_count =", self.get_edge_count())
            print("  face_count =", self.get_face_count())
            print("  volume_count =", self.get_volume_count())
        else:
            print("Mesh properties: none")

        if self.has_fem_properties():
            print("FEM results:")
            print("  vonmises_stress = {} MPa".format(
                self.get_vonmises_stress()))
            print("  tresca_stress = {} MPa".format(
                self.get_tresca_stress()))
            print("  max_displacement = {} m".format(
                self.get_max_displacement()))
            print("  structural_failure =",
                  self.get_structural_failure())
            print("  maximum_pressure = {} MPa".format(
                self.get_maximum_pressure()))
        else:
            print("FEM Results: none")

    @staticmethod
    def print_properties(obj):
        print(obj.Name, "properties:")
        for name in obj.PropertiesList:
            print(" ", name, "=", getattr(obj, name))

    def sketch_set_length(self, param: str, value: float):
        """
        Sets a length constraint of the sketch object in meters. This method
        may throw an exception if an internal constraint cannot be satisfied.
        """
        obj = self.doc.getObject('Sketch')
        obj.setDatum(param, Units.Quantity(
            value * 1e3, Units.Unit('mm')))

    def sketch_set_angle(self, param: str, value: float):
        """
        Sets a length constraint of the sketch object in radians. This method
        may throw an exception if an internal constraint cannot be satisfied.
        """
        obj = self.doc.getObject('Sketch')
        obj.setDatum(param, Units.Quantity(value, Units.Unit('rad')))

    def sketch_get_length(self, param: str) -> float:
        obj = self.doc.getObject('Sketch')
        return float(obj.getDatum(param).getValueAs('mm')) * 1e-3

    def sketch_get_angle(self, param: str) -> float:
        obj = self.doc.getObject('Sketch')
        return float(obj.getDatum(param).getValueAs('rad'))

    def set_test_pressure(self, value: float):
        """
        Sets the outside pressure acting on the vessel in mega pascals.
        """
        obj = self.doc.getObject('ConstraintPressure')
        obj.Pressure = float(value)

    def get_test_pressure(self) -> float:
        obj = self.doc.getObject('ConstraintPressure')
        return float(obj.Pressure)

    def set_mesh_length(self, value: float):
        """
        Sets the maximum edge length for the meshing algorithm in meters.
        """
        obj = self.doc.getObject('FEMMeshGmsh')
        obj.CharacteristicLengthMax = Units.Quantity(
            value * 1e3, Units.Unit('mm'))

    def get_mesh_length(self) -> float:
        obj = self.doc.getObject('FEMMeshGmsh')
        return obj.CharacteristicLengthMax.getValueAs('mm') * 1e-3

    def set_youngs_modulus(self, value: float):
        """
        Sets the Youngs modulus of the material in mega pascals.
        """
        obj = self.doc.getObject('MaterialSolid')
        mat = dict(obj.Material)
        mat['YoungsModulus'] = str(value) + ' MPa'
        obj.Material = mat

    def get_youngs_modulus(self) -> float:
        obj = self.doc.getObject('MaterialSolid')
        return Units.Quantity(obj.Material['YoungsModulus']).getValueAs('MPa')

    def set_poisson_ratio(self, value: float):
        """
        Sets the poisson ratio of the material.
        """
        obj = self.doc.getObject('MaterialSolid')
        mat = dict(obj.Material)
        mat['PoissonRatio'] = str(value)
        obj.Material = mat

    def get_poisson_ratio(self):
        obj = self.doc.getObject('MaterialSolid')
        return float(obj.Material['PoissonRatio'])

    def set_tensile_strength(self, value: float):
        """
        Sets the ultimate tensile strength of the material in mega pascals.
        """
        obj = self.doc.getObject('MaterialSolid')
        mat = dict(obj.Material)
        mat['UltimateTensileStrength'] = str(value) + ' MPa'
        obj.Material = mat

    def get_tensile_strength(self) -> float:
        obj = self.doc.getObject('MaterialSolid')
        return Units.Quantity(obj.Material['UltimateTensileStrength']).getValueAs('MPa')

    def set_density(self, value: float):
        """
        Sets the density of the material in kg/m3.
        """
        obj = self.doc.getObject('MaterialSolid')
        mat = dict(obj.Material)
        mat['Density'] = str(value) + ' kg/m^3'
        obj.Material = mat

    def get_density(self):
        obj = self.doc.getObject('MaterialSolid')
        return float(Units.Quantity(obj.Material['Density']).getValueAs('kg/m^3'))

    def recompute(self):
        self.doc.recompute()

    def get_body_area(self):
        """
        Returns the body volume in square meters.
        """
        obj = self.doc.getObject('Body')
        return obj.Shape.Area * 1e-6

    def get_body_volume(self):
        """
        Returns the body volume in cubic meters.
        """
        obj = self.doc.getObject('Body')
        return obj.Shape.Volume * 1e-9

    def get_body_mass(self):
        """
        Returns the body mass in kilogramms.
        """
        return self.get_body_volume() * self.get_density()

    def get_outer_area(self):
        """
        Returns the outer area in square meters.
        """
        obj = self.doc.getObject('Body')
        return obj.Shape.OuterShell.Area * 1e-6

    def get_outer_volume(self):
        """
        Returns the outer area in cubic meters.
        """
        obj = self.doc.getObject('Body')
        return obj.Shape.OuterShell.Volume * 1e-9

    def get_outer_length(self):
        """
        Returns the outer length (along the x-axis) in meters.
        """
        obj = self.doc.getObject('Body')
        return obj.Shape.BoundBox.XLength * 1e-3

    def get_outer_diameter(self):
        """
        Returns the outer diameter (along maximum of the y and z-axis diameter)
        in meters.
        """
        obj = self.doc.getObject('Body')
        return max(obj.Shape.BoundBox.YLength, obj.Shape.BoundBox.ZLength) * 1e-3

    def get_inner_area(self):
        obj = self.doc.getObject('Body')
        return self.get_body_area() - self.get_outer_area()

    def get_inner_volume(self):
        obj = self.doc.getObject('Body')
        return self.get_outer_volume() - self.get_body_volume()

    def clean(self):
        """
        Removes all temporary artifacts from the model.
        """
        if self.doc.getObject('CCX_Results'):
            self.doc.removeObject('CCX_Results')
        if self.doc.getObject('ResultMesh'):
            self.doc.removeObject('ResultMesh')
        if self.doc.getObject('ccx_dat_file'):
            self.doc.removeObject('ccx_dat_file')

    RE_MESH_WARNINGS = re.compile(
        "^Warning :\s*(\d+) warnings?$", re.MULTILINE)
    RE_MESH_ERRORS = re.compile(
        "^Warning :\s*(\d+) errors?$", re.MULTILINE)

    def run_analysis(self):
        """
        Set the various parameters, then call this method and query the results.
        """
        self.clean()
        self.doc.recompute()

        if self.debug:
            print("Running GMSH mesher ...", end=' ', flush=True)
        self.mesh_warnings = 0
        self.mesh_errors = 0
        mesher = GmshTools(self.doc.getObject('FEMMeshGmsh'))
        err = mesher.create_mesh()
        if isinstance(err, str):
            match = PressureVessel.RE_MESH_WARNINGS.search(err)
            self.mesh_warnings = int(match.group(1)) if match else 1e99
            if not match:
                print(err)
            match = PressureVessel.RE_MESH_ERRORS.search(err)
            self.mesh_errors = int(match.group(1)) if match else 1e99
        elif err:
            raise err
        obj = self.doc.getObject('FEMMeshGmsh').FemMesh
        if self.debug:
            print(obj.NodeCount, "nodes,",
                  obj.EdgeCount, "edges,",
                  obj.FaceCount, "faces,",
                  obj.VolumeCount, "volumes,",
                  self.mesh_errors, "mesh errors,",
                  self.mesh_warnings, "mesh warnings")

        if self.debug:
            print("Running FEM analysis ...", end=' ', flush=True)
        fea = FemToolsCcx(
            self.doc.getObject('Analysis'),
            self.doc.getObject('SolverCcxTools'))
        fea.purge_results()
        fea.update_objects()
        fea.setup_working_dir()
        fea.setup_ccx()
        err = fea.check_prerequisites()
        if err:
            raise ValueError("FEM error: " + err)
        fea.write_inp_file()
        fea.ccx_run()
        fea.load_results()
        obj = self.doc.getObject('CCX_Results')
        assert obj.ResultType == 'Fem::ResultMechanical'
        if self.debug:
            print("vonMises stress: {:.2f} MPa".format(
                self.get_vonmises_stress()))

    def has_mesh_properties(self):
        obj = self.doc.getObject('FEMMeshGmsh').FemMesh
        return True if obj and obj.NodeCount else False

    def get_node_count(self):
        obj = self.doc.getObject('FEMMeshGmsh').FemMesh
        return obj.NodeCount

    def get_edge_count(self):
        obj = self.doc.getObject('FEMMeshGmsh').FemMesh
        return obj.EdgeCount

    def get_face_count(self):
        obj = self.doc.getObject('FEMMeshGmsh').FemMesh
        return obj.FaceCount

    def get_volume_count(self):
        obj = self.doc.getObject('FEMMeshGmsh').FemMesh
        return obj.VolumeCount

    def has_fem_properties(self):
        obj = self.doc.getObject('CCX_Results')
        return True if obj else False

    def get_vonmises_stress(self) -> float:
        """
        Returns the maximum vonMises stress in mega pascals.
        """
        obj = self.doc.getObject('CCX_Results')
        return max(obj.vonMises)

    def get_tresca_stress(self) -> float:
        """
        Returns the maximum tresca (shear) stress in mega pascals.
        """
        obj = self.doc.getObject('CCX_Results')
        return max(obj.MaxShear)

    def get_max_displacement(self) -> float:
        """
        Returns the maximum displacement in meters.
        """
        obj = self.doc.getObject('CCX_Results')
        return max(obj.DisplacementLengths) * 1e-3

    def get_structural_failure(self) -> bool:
        """
        Returns if the maximum vonMises stress is larger than the tensile strength.
        """
        return self.get_vonmises_stress() >= self.get_tensile_strength()

    def get_maximum_pressure(self) -> float:
        """
        Returns the maximum pressure the vessel can withstand without deformation.
        """
        return self.get_test_pressure() * self.get_tensile_strength() / self.get_vonmises_stress()

    def set(self, name: str, value: Any):
        """
        Allows to set any parameter value of the model. Do not name methods starting
        with set_ if you do not want it to be a parameter.
        """
        if name in self.sketch_params:
            self.sketch_set_length(name, value)
        elif hasattr(self, 'set_' + name):
            getattr(self, 'set_' + name)(value)
        else:
            raise ValueError("Unknown parameter: " + name)

    def get(self, name: str) -> Any:
        """
        Allows to get any parameter value of the model. Do not name methods starting
        with set_ if you do not want it to be a parameter.
        """
        if name in self.sketch_params:
            if 'angle' in name:
                return self.sketch_get_angle(name)
            else:
                return self.sketch_get_length(name)
        elif hasattr(self, name):
            return getattr(self, name)
        elif hasattr(self, 'get_' + name):
            return getattr(self, 'get_' + name)()
        else:
            raise ValueError("Unknown parameter: " + name)

    def csv_open_output(self, filename: str) -> csv.DictWriter:
        if filename == '-':
            file = sys.stdout
        else:
            file = open(filename, 'w', newline='', encoding='utf-8')

        fieldnames = list(self.sketch_params)
        fieldnames.extend([
            'test_pressure', 'mesh_length',
            'youngs_modulus', 'poisson_ratio', 'tensile_strength', 'density',
            'body_area', 'body_volume', 'body_mass',
            'outer_length', 'outer_diameter', 'outer_area', 'outer_volume',
            'inner_area', 'inner_volume',
            'node_count', 'edge_count', 'face_count', 'volume_count',
            'mesh_warnings', 'mesh_errors',
            'vonmises_stress', 'tresca_stress', 'max_displacement',
            'structural_failure', 'maximum_pressure'
        ])

        writer = csv.DictWriter(file, fieldnames)
        if filename != '-':
            writer.file = file
        writer.writeheader()
        return writer

    def csv_write_row(self, writer: csv.DictWriter):
        row = {name: self.get(name) for name in writer.fieldnames}
        writer.writerow(row)

    def csv_flush_output(self, writer: csv.DictWriter):
        if hasattr(writer, 'file'):
            writer.file.flush()

    def csv_close_output(self, writer: csv.DictWriter):
        if hasattr(writer, 'file'):
            writer.file.close()

    def study_random(self, count: int, output: str):
        writer = self.csv_open_output(output)
        while count > 0:
            try:
                # set sketch lengths
                for name in self.sketch_params:
                    if 'thickness' in name:
                        value = random.uniform(0.005, 0.02)
                        print("setting", name, "=", value)
                        self.sketch_set_length(name, value)
                    elif 'length' in name:
                        value = random.uniform(0.0, 2.0)
                        print("setting", name, "=", value)
                        self.sketch_set_length(name, value)
                    elif 'radius' in name:
                        value = random.uniform(0.1, 1.0)
                        print("setting", name, "=", value)
                        self.sketch_set_length(name, value)
                    elif 'angle' in name:
                        value = random.uniform(0.0, 0.5 * math.pi)
                        print("setting", name, "=", value)
                        self.sketch_set_angle(name, value)
                    else:
                        value = random.uniform(0.0, 1.0)
                        print("setting", name, "=", value)
                        self.sketch_set_length(name, value)

                self.recompute()
                mesh_len = 0.1 * (self.get_body_volume() ** 0.333)
                print("setting mesh_len =", mesh_len)
                self.set_mesh_length(mesh_len)
            except ValueError:
                # we just skip bad combinations for now
                print("failed constraints")
                continue

            try:
                self.run_analysis()
            except RuntimeError:
                print("failed analysis")
                continue

            self.csv_write_row(writer)
            count -= 1
        self.csv_close_output(writer)


def run(args=None):
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('model', type=str, metavar='FILE',
                        help="a parametric FreeCAD model of the pressure vessel model")
    parser.add_argument('--study', type=str, choices=['random'], default=None,
                        help="generate a CSV file for the given study")
    parser.add_argument('--output', type=str, metavar='FILE', default='output.csv',
                        help="output CSV filename")
    parser.add_argument('--count', type=int, metavar='NUM', default=1000,
                        help="generate this many random samples")
    args = parser.parse_args(args)

    vessel = PressureVessel(args.model)

    if args.study == 'random':
        vessel.study_random(args.count, args.output)
    else:
        vessel.run_analysis()
        vessel.print_info()


if __name__ == '__main__':
    run()
