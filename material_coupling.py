"""
Adapter coupling heatConduction → MaterialStateSolver (material_engine).

This module is opt-in and lazy-imported: heatConduction never imports
material_engine at module load time. If MaterialStateSolver is not
installed, heatConduction remains fully functional on its static-JSON path.

Typical usage:

    from material_coupling import make_layered_coupler
    coupler = make_layered_coupler(para, yaml_paths=['pica', 'steel_304'])
    TProfile, cache = hc.solve(para, material_hook=coupler.hook)

One MaterialEngine + MaterialState is created per TPS layer. The hook
is called each converged timestep with cache['T0']=T_old, cache['T']=T_new;
it calls material_update per layer and writes the updated per-node
k/rho/cp/Q arrays back into para. Shared-boundary nodes use
last-write-wins (layer l+1 overwrites the interface node), matching
parameter.normalize_conductivity.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np


# ------------------------------------------------------------------------ #
#  Lazy import
# ------------------------------------------------------------------------ #

def _import_material_engine():
    try:
        from material_engine import MaterialEngine  # noqa: F401
    except ImportError as e:
        raise ImportError(
            "material_engine not found. Install MaterialStateSolver "
            "(pip install -e /path/to/MaterialStateSolver) to use "
            "evolving material properties."
        ) from e
    return MaterialEngine


# ------------------------------------------------------------------------ #
#  YAML path resolution
# ------------------------------------------------------------------------ #

_DEFAULT_DB_REL = Path(__file__).parent / '..' / '..' / 'MaterialStateSolver' / 'tps_material_db' / 'models'


def _candidate_roots():
    roots = []
    env = os.environ.get('MATERIAL_STATE_SOLVER_PATH')
    if env:
        roots.append(Path(env) / 'tps_material_db' / 'models')
    roots.append(_DEFAULT_DB_REL.resolve())
    return roots


def resolve_material_yaml(name_or_path):
    """Resolve a material name or path to an absolute YAML path.

    Accepts:
        - absolute or relative path to a .yaml file
        - bare name e.g. 'pica' → searches under DB roots
        - category-qualified name e.g. 'ablators/pica'
        - any of the above with or without '.yaml' suffix

    Search order:
        1. The literal path as given (cwd-relative or absolute).
        2. $MATERIAL_STATE_SOLVER_PATH/tps_material_db/models/**/<name>.yaml
        3. Default DB at ../../MaterialStateSolver/tps_material_db/models/
    """
    p = Path(name_or_path)
    # direct-path check
    if p.suffix == '.yaml' and p.exists():
        return p.resolve()

    # try as <root>/<name_or_path>.yaml and <root>/<name_or_path>
    stem = str(name_or_path)
    if not stem.endswith('.yaml'):
        stem_yaml = stem + '.yaml'
    else:
        stem_yaml = stem

    searched = []
    for root in _candidate_roots():
        searched.append(str(root))
        if not root.exists():
            continue
        # category-qualified: 'ablators/pica' → root/ablators/pica.yaml
        direct = root / stem_yaml
        if direct.exists():
            return direct.resolve()
        # recursive glob by basename
        basename = Path(stem_yaml).name
        matches = list(root.rglob(basename))
        if matches:
            return matches[0].resolve()

    raise FileNotFoundError(
        f"Material YAML '{name_or_path}' not found. Searched roots: {searched}. "
        "Set MATERIAL_STATE_SOLVER_PATH env var or pass an absolute path."
    )


# ------------------------------------------------------------------------ #
#  Coupler
# ------------------------------------------------------------------------ #

class LayeredMaterialCoupler:
    """One MaterialEngine + MaterialState per TPS layer.

    The coupler owns one engine per layer; each layer's MaterialState
    spans npl nodes (matching parameter.normalize_conductivity's layer
    → node mapping). On __init__ it seeds initial per-node k/rho/cp/Q
    into para via evaluate_properties(). Each timestep hook() calls
    material_update() per layer and writes back to para.

    Shared boundary nodes are written last-write-wins so layer l+1's
    properties win at the interface, matching the static pipeline.
    """

    def __init__(self, para, yaml_paths, *, zero_Q=False):
        MaterialEngine = _import_material_engine()

        k_layers = np.atleast_1d(np.asarray(para['layerConductivities']))
        n_layers = len(k_layers)
        if len(yaml_paths) != n_layers:
            raise ValueError(
                f"yaml_paths has {len(yaml_paths)} entries but para has "
                f"{n_layers} layers"
            )

        self.yaml_paths = [resolve_material_yaml(p) for p in yaml_paths]
        self.npl = int(para['nodesPerLayer'])
        self.n_layers = n_layers
        self.layer_ranges = [
            (l * (self.npl - 1), l * (self.npl - 1) + (self.npl - 1))
            for l in range(n_layers)
        ]
        self.zero_Q = bool(zero_Q)

        self.engines = [MaterialEngine(str(p)) for p in self.yaml_paths]
        self.states = [e.create_initial_state(n_points=self.npl) for e in self.engines]

        self._initialize_properties(para)

    # -------------------------------------------------------------- #
    def _num_nodes(self):
        return self.n_layers * (self.npl - 1) + 1

    def _initialize_properties(self, para):
        """Evaluate properties at IC temperature and write per-node arrays
        into para['conductivity' / 'density' / 'heatCapacity' / 'volumetricHeatSource'].
        """
        N = self._num_nodes()
        T_ic = float(para['IC value'])
        k_arr = np.zeros(N)
        rho_arr = np.zeros(N)
        cp_arr = np.zeros(N)
        Q_arr = np.zeros(N)
        T_vec = np.full(self.npl, T_ic, dtype=np.float64)
        for l, (s, e) in enumerate(self.layer_ranges):
            props = self.engines[l].evaluate_properties(T_vec, self.states[l])
            # last-write-wins at shared node (s..e inclusive; e == (l+1)*(npl-1))
            k_arr[s:e + 1] = np.asarray(props.k)
            rho_arr[s:e + 1] = np.asarray(props.rho)
            cp_arr[s:e + 1] = np.asarray(props.cp)
            Q_arr[s:e + 1] = 0.0 if self.zero_Q else np.asarray(props.Q)

        para['conductivity'] = k_arr
        para['density'] = rho_arr
        para['heatCapacity'] = cp_arr
        para['volumetricHeatSource'] = Q_arr

    # -------------------------------------------------------------- #
    def hook(self, para, cache, timeStep):
        """Called by hc.solve() after each converged Newton step.

        At entry: cache['T0'] is T_old (previous-step), cache['T'] is
        T_new (just converged). This method MUST be called BEFORE
        storeUpdateResult (which would overwrite T0 with T_new).
        """
        if timeStep == 0:
            return  # initial properties seeded in constructor
        dt = float(para['deltaTime'])
        T_old = np.asarray(cache['T0']).ravel()
        T_new = np.asarray(cache['T']).ravel()

        k_arr = para['conductivity']
        rho_arr = para['density']
        cp_arr = para['heatCapacity']
        Q_arr = para.get('volumetricHeatSource', None)
        if Q_arr is None:
            Q_arr = np.zeros_like(k_arr)
            para['volumetricHeatSource'] = Q_arr

        for l, (s, e) in enumerate(self.layer_ranges):
            T_n = T_old[s:e + 1]
            T_np1 = T_new[s:e + 1]
            result = self.engines[l].material_update(T_n, T_np1, self.states[l], dt)
            self.states[l] = result.state
            k_arr[s:e + 1] = np.asarray(result.properties.k)
            rho_arr[s:e + 1] = np.asarray(result.properties.rho)
            cp_arr[s:e + 1] = np.asarray(result.properties.cp)
            Q_arr[s:e + 1] = 0.0 if self.zero_Q else np.asarray(result.properties.Q)

    # -------------------------------------------------------------- #
    def reset(self, para):
        """Reset states to virgin and re-seed initial properties (for
        optimizer FD perturbations that start each solve from IC)."""
        self.states = [e.create_initial_state(n_points=self.npl) for e in self.engines]
        self._initialize_properties(para)


# ------------------------------------------------------------------------ #
#  Factory
# ------------------------------------------------------------------------ #

def make_layered_coupler(para, yaml_paths=None, *, zero_Q=False):
    """Build a LayeredMaterialCoupler.

    If yaml_paths is None, resolve from para['materials'] (which may be
    material names, category/name, or paths).
    """
    if yaml_paths is None:
        yaml_paths = list(para['materials'])
    return LayeredMaterialCoupler(para, yaml_paths, zero_Q=zero_Q)
