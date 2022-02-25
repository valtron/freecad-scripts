from pathlib import Path
from time import perf_counter as time
import pickle

from freecad_scripts import pressure_vessel as pv
import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

def main():
    t = Analysis()
    
    rng = np.random.default_rng(0)
    inputs = [
        random_input(rng)
        for _ in range(300)
    ]
    mesh_scales = [0.2, 0.15, 0.1, 0.05]
    
    field = 'vonmises_stress'
    
    data = []
    for i, inp in tqdm(list(enumerate(inputs)), "input"):
        l = []
        for mesh_scale in mesh_scales:
            out = run(i, t, inp, mesh_scale)
            out_fem = out.get('out_fem')
            if out_fem is None:
                out = None
            else:
                if out_fem.get(field) is None:
                    out = None
            l.append(out)
        
        if l[-1] is None:
            continue
        
        for j in range(len(mesh_scales) - 1):
            if l[j] is None:
                continue
            err = l[j]['out_fem'][field][-5] / l[-1]['out_fem'][field][-5]
            if err < 1:
                err = 1/err
            vol = l[j]['out_meshing']['body_volume']
            mesh_len = mesh_scales[j] * vol**(1/3)
            thick1 = inp['thick1']
            thick2 = inp['thick2']
            aspect_ratio = inp['L'] / inp['D']
            theta = inp['theta']
            R = inp['D'] / 2
            k = 1/np.hypot(R * np.tan(theta), R)
            data.append([
                i, j, err, mesh_len, thick1, thick2, aspect_ratio, theta, k, vol,
            ])
    data = np.array(data)
    print(data.shape)
    
    err_max = 2
    E = np.clip(data[:,2], -np.inf, err_max)
    fig, ax = plt.subplots(constrained_layout = True, figsize = (8, 5))
    ax.set_title("Relative error distribution vs. mesh scale")
    ax.set_xlabel("Relative error")
    ax.set_ylabel("Density")
    ax.axvline(err_max, linestyle = '--', color = 'gray')
    bins = np.linspace(1, err_max + 0.1, 30)
    for j, mesh_scale in enumerate(mesh_scales[:-1]):
        ax.hist(
            E[data[:,1] == j], bins = bins, density = True, alpha = 0.5,
            label = "mesh_scale = {:.2f}".format(mesh_scale),
        )
    ax.legend()
    plt.savefig("fig/0 - Histogram.pdf")
    
    xs = [
        (3, "Mesh length"),
        (4, "Thickness1"),
        (5, "Thickness2"),
        (6, "Aspect ratio"),
        (7, "Theta"),
        (8, "Curvature"),
        (9, "Body volume"),
    ]
    
    for x_idx, x_label in xs:
        fig, ax = plt.subplots(constrained_layout = True, figsize = (8, 5))
        ax.set_title("Relative error vs. {}".format(x_label))
        ax.set_xlabel(x_label)
        ax.set_ylabel("Relative error")
        ax.axhline(err_max, linestyle = '--', color = 'gray')
        for j, mesh_scale in enumerate(mesh_scales[:-1]):
            m = (data[:,1] == j)
            ax.scatter(data[m,x_idx], E[m], marker = '.', alpha = 0.3)
        if x_label in ("Body volume", "Aspect ratio", "Curvature"):
            ax.set_xscale('log')
        plt.savefig("fig/{} - {}.pdf".format(x_idx, x_label))

def run(i, t, inp, mesh_scale):
    cache_file = Path('.cache/{:03}-{:.2f}.pickle'.format(i, mesh_scale))
    if not cache_file.exists():
        out = t.do(inp, mesh_scale)
        with cache_file.open('wb') as fh:
            pickle.dump(out, fh)
    with cache_file.open('rb') as fh:
        return pickle.load(fh)

class Analysis:
    def __init__(self):
        self.impl = pv.PressureVessel('freecad_scripts/models/pressure_vessel0.FCStd', debug = False)
        self.impl.sketch_params = None
        self._docs = { 0: self.impl.doc }
    
    def do(self, inp, mesh_scale):
        out = {}
        out['input'] = inp
        out['mesh_scale'] = mesh_scale
        out['error'] = None
        
        try:
            sketch, params = convert_inputs(**inp)
        except AssertionError as ex:
            out['error'] = ex
            return out
        
        out['sketch'] = sketch
        out['params'] = params
        
        self._set_sketch(sketch)
        try:
            self._set_params(params)
        except ValueError as ex:
            out['error'] = ex
            return out
        
        mesh_len = mesh_scale * (self.impl.get_body_volume()**(1/3))
        self.impl.set_mesh_length(mesh_len)
        
        t_meshing = time()
        out_meshing = self._run_meshing()
        t_meshing = time() - t_meshing
        out['out_meshing'] = out_meshing
        out['time_meshing'] = t_meshing
        
        t_fem = time()
        out_fem = self._run_fem()
        t_fem = time() - t_fem
        out['out_fem'] = out_fem
        out['time_fem'] = t_fem
        
        return out
    
    def _run_meshing(self):
        from freecad_scripts.libs import GmshTools
        
        self.impl.clean()
        self.impl.recompute()
        mesher = GmshTools(self.impl.doc.getObject('FEMMeshGmsh'))
        err = mesher.create_mesh()
        
        out = {
            'node_count': self.impl.get_node_count(),
            'edge_count': self.impl.get_edge_count(),
            'face_count': self.impl.get_face_count(),
            'err': err,
        }
        fields = [
            'test_pressure', 'poisson_ratio',
            'density', 'body_area', 'body_volume', 'body_mass',
            'outer_length', 'outer_diameter', 'outer_area', 'outer_volume',
            'inner_area', 'inner_volume', 'volume_count',
        ]
        for f in fields:
            out[f] = getattr(self.impl, 'get_{}'.format(f))()
        out['youngs_modulus'] = self.impl.get_youngs_modulus().Value
        out['tensile_strength'] = self.impl.get_tensile_strength().Value
        
        return out
    
    def _run_fem(self):
        from freecad_scripts.libs import FemToolsCcx
        
        out = {}
        out['quantiles'] = QUANTILES
        
        doc = self.impl.doc
        obj = doc.getObject('FEMMeshGmsh').FemMesh
        fea = FemToolsCcx(doc.getObject('Analysis'), doc.getObject('SolverCcxTools'))
        fea.purge_results()
        fea.update_objects()
        fea.setup_working_dir()
        fea.setup_ccx()
        err = fea.check_prerequisites()
        if err:
            out['error'] = err
            return out
        fea.write_inp_file()
        fea.ccx_run()
        fea.load_results()
        ccx = doc.getObject('CCX_Results')
        
        if ccx is None:
            out['error'] = ValueError("missing CCX_Results")
            return out
        
        try:
            out['vonmises_stress'] = np.quantile(ccx.vonMises, QUANTILES)
            out['tresca_stress'] = np.quantile(ccx.MaxShear, QUANTILES)
            out['max_displacement'] = np.quantile(ccx.DisplacementLengths, QUANTILES) * 1e-3
        except RuntimeError as ex:
            out['error'] = ex
        
        return out
    
    def _set_sketch(self, sketch):
        if sketch not in self._docs:
            import FreeCAD
            self._docs[sketch] = FreeCAD.open('freecad_scripts/models/pressure_vessel{}.FCStd'.format(sketch))
        self.impl.doc = self._docs[sketch]
    
    def _set_params(self, params):
        self.impl.disable_recompute()
        sk = self.impl.doc.getObject('Sketch')
        for i in range(len(sk.Constraints)):
            if sk.Constraints[i].IsActive:
                sk.toggleActive(i)
        for k, v in params.items():
            if 'angle' in k:
                self.impl.sketch_set_angle(k, v)
            else:
                self.impl.sketch_set_length(k, v)
        for i in range(len(sk.Constraints)):
            if not sk.Constraints[i].IsActive:
                sk.toggleActive(i)
        self.impl.recompute()

QUANTILES = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 0.999, 1]

def convert_inputs(L, D, thick1, thick2, theta, nw, nh):
    R = D/2
    
    assert 0 <= theta <= np.pi/2
    assert thick1 <= thick2 <= R
    
    x = R * np.tan(theta)
    r = np.hypot(x, R)
    
    l1 = L - 2 * (r - x)
    assert l1 >= 0, l1
    
    radius = R - thick2
    
    if nw >= 1e-3 and nh >= 1e-3:
        assert np.isclose(theta, 0)
        assert np.isclose(thick1, thick2)
        assert nh <= radius
        assert nw <= l1
        return 3, {
            'thickness1': thick1,
            'length': max(l1, 1e-3),
            'thickness2': nw,
            'thickness3': nh,
            'radius': radius,
        }
    
    l2 = (r - thick1) * (1 - np.sqrt(1 - ((R - thick2)/(r - thick1))**2))
    
    if theta > 0:
        return 4, {
            'thickness1': thick1,
            'length1': max(l1, 1e-3),
            'thickness2': thick2,
            'radius': radius,
            'length2': l2,
        }
    
    if thick2 > thick1:
        return 2, {
            'thickness1': thick1,
            'length': max(l1, 1e-3),
            'thickness2': thick2,
            'radius': radius,
        }
    
    if l1 >= 1e-3:
        return 1, {
            'radius': radius,
            'thickness': thick1,
            'length': l1,
        }
    
    return 0, {
        'radius': radius,
        'thickness': thick1,
    }

def random_input(rng):
    L = rng.uniform(0.3, 2)
    D = rng.uniform(0.3, 1)
    aspect_ratio = L/D
    theta_min = max(0, np.arctan((1 - aspect_ratio**2)/(2*aspect_ratio)))
    theta = rng.uniform(theta_min, np.pi/2)
    thick1 = rng.uniform(0.005, 0.02)
    thick2 = rng.uniform(thick1, 2 * thick1)
    nw = 0
    nh = 0
    
    return {
        'L': L,
        'D': D,
        'thick1': thick1,
        'thick2': thick2,
        'theta': theta,
        'nw': 0,
        'nh': 0,
    }

if __name__ == '__main__':
    main()
