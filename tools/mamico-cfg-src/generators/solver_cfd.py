import shutil
from pathlib import Path

from generators.cell_size import get_timesteps_per_coupling_cycle
from generators.domain_size import get_domain_sizes, get_md_offset_z
from generators.simulation import get_coupling_cycles
from generators.solver_md import get_wall_velocity
from utils import (
    DYNAMIC_VISCOSITY,
    MD_DENSITY,
    MD_TIMESTEP,
    check_if_replacing,
    get_asset_text,
    get_domain_size,
)
from xml_templating import PartialXml

# preCICE participant/mesh/data names, must match coupling/scenario/precice/CouetteScenario.h
_precice_participant = "mamico"
_precice_foam_participant = "Fluid"
_precice_foam_mesh_outer = "Fluid-Mesh-Outer"
_precice_foam_mesh_inner = "Fluid-Mesh-Inner"
_precice_macro2md_mesh = "mamico-macro2MD-mesh"
_precice_md2macro_mesh = "mamico-MD2macro-mesh"
_precice_macro2md_data = "VelocityMacro"
_precice_md2macro_data = "VelocityMicro"
# CFD cells surrounding the MD hole; they supply CFD -> MD values
_foam_cell_set_shell = "mdShell"
# The six faces of the hole; they receive MD -> CFD values
_foam_inner_patches = ("lefti", "righti", "fronti", "backi", "buttomi", "topi")
# Overlap layers MaMiCo keeps outside its md2macro region
# (CouettePreciceInterface::getOuterRegion in coupling/scenario/precice/CouetteScenario.h)
_precice_outer_region = 3
# Extra cell layers by which the OpenFOAM hole is pulled INSIDE the md2macro
# region. With 0 the hole face coincides with the md2macro boundary, so a
# nearest-neighbor read mapping hands the CFD one cell half a cell away. With 1
# the face lies strictly between two genuinely sampled MaMiCo cells, which is
# what makes linear-cell-interpolation a real (inner+outer)/2 average over a
# volume centred on the boundary -- the operator Nie 2004 section 2 specifies
# and the one the direct IcoFoamInterface coupling used. getOuterRegion() stays
# at 3, so the Nie momentum layers are untouched.
_precice_hole_inset = 1  # default; overridden by the hole_inset option
_precice_config_name = "precice-config.xml"
_precice_foam_case_name = "FoamCase"


def _create_foam_setup(get_config_value) -> Path:
    src_path = Path(__file__).parent.parent / "assets" / f"FoamSetup.template"
    dst_path = Path(get_config_value("output_dir")) / "FoamSetup"
    check_if_replacing(dst_path, get_config_value)
    shutil.copytree(src_path, dst_path, dirs_exist_ok=True)
    return dst_path


def _get_geometry(get_config_value) -> dict:
    """Channel, MD-box and hole geometry for the OpenFOAM participant.

    The hole is the md2macro region, i.e. the MD box inset by MaMiCo's overlap
    layers -- not the MD box itself. The overlap ring therefore keeps its CFD
    cells and both solvers cover it, which is what makes this a Schwarz coupling.
    At cell_size 5 on the medium domain this reproduces FoamSetup.template exactly.
    """
    cell_size = get_config_value("cell_size")
    md_size, cfd_size = get_domain_sizes(get_config_value)
    offset_xy = (cfd_size - md_size) / 2
    offset_z = get_md_offset_z(get_config_value)
    md_min = (offset_xy, offset_xy, offset_z)
    md_max = (offset_xy + md_size, offset_xy + md_size, offset_z + md_size)

    overlap = _precice_outer_region * cell_size
    hole_overlap = (_precice_outer_region + int(get_config_value("hole_inset"))) * cell_size
    hole_min = tuple(v + hole_overlap for v in md_min)
    hole_max = tuple(v - hole_overlap for v in md_max)

    def counts(lo, hi):
        # cells below / inside / above the hole along one axis
        return (max(1, round(lo / cell_size)),
                max(1, round((hi - lo) / cell_size)),
                max(1, round((cfd_size - hi) / cell_size)))

    return dict(
        cell_size=cell_size,
        md_size=md_size,
        cfd_size=cfd_size,
        md_min=md_min,
        md_max=md_max,
        hole_min=hole_min,
        hole_max=hole_max,
        nx=counts(hole_min[0], hole_max[0]),
        ny=counts(hole_min[1], hole_max[1]),
        nz=counts(hole_min[2], hole_max[2]),
        # the overlap ring: a box over the whole MD domain picks up exactly the
        # cells between the hole surface and the MD box edge, since the hole is empty
        shell_min=md_min,
        shell_max=md_max,
        # Cells between the hole surface and the md2macro boundary exist in the
        # CFD mesh but must not join mdShell, which has to stay exactly the
        # _precice_outer_region macro2MD layers getOuterRegion() defines.
        # With _precice_hole_inset = 0 this box holds no cells and the
        # corresponding topoSet delete action is a no-op.
        shell_cut_min=tuple(v + overlap for v in md_min),
        shell_cut_max=tuple(v - overlap for v in md_max),
    )




def _get_side_conditions(get_config_value) -> dict:
    """Lateral CFD boundaries: open, with Neumann pressure everywhere.

    This is the arrangement of the only published two-way MaMiCo/OpenFOAM study,
    Wittenberg & Neumann, Computation 9 (2021) 128 appendix A.3: "A Neumann boundary
    condition is applied at every boundary for the pressure and at the outer boundary
    for the velocity."

    The cyclic alternative is faithful to Neumann & Bian, CPC 220 (2017) section 2.3,
    but that reference is one-way, where the hole injects nothing the CFD must absorb.
    Two-way it deadlocks the continuity treatment: cyclic patches are coupled, so
    adjustPhi skips them, and pinning p at the far wall to compensate makes
    p.needReference() false, so adjustPhi is never called at all. The divergence error
    in the MD interface data is then left for the pressure solve to redistribute, which
    drives exactly the streamwise blow-up seen in sept20PreciceFoamTest.

    With p Neumann on every patch, p.needReference() is true and adjustPhi runs; the
    four open sides have fixesValue() == false, so they can absorb the hole's net flux.
    Pressure level is set by pRefCell/pRefValue in system/fvSolution.
    """
    return {
        "side-mesh-xmin": "type patch;",
        "side-mesh-xmax": "type patch;",
        "side-mesh-ymin": "type patch;",
        "side-mesh-ymax": "type patch;",
        "side-u-bc": "type            zeroGradient;",
        "side-p-bc": "type            zeroGradient;",
        "fixedwall-p-bc": "type            zeroGradient;",
    }


def _get_hole_condition(get_config_value) -> dict:
    """Whether MD velocities actually reach the CFD solution.

    Two-way: Dirichlet on the hole faces; the adapter overwrites them each window.
    One-way: zeroGradient, so the BC re-derives itself from the interior and the values
    the adapter writes are discarded. MaMiCo and preCICE still exchange both directions
    -- the published arrangement is "technically bidirectional, physically
    unidirectional" (Neumann & Bian section 4), and PreciceAdapter only registers
    MD2macro vertices under twoWayCoupling(), which preCICE needs non-empty.
    """
    if get_config_value("coupling_2way"):
        return {"hole-u-bc": "type            fixedValue;\n        value           uniform (0 0 0);"}
    return {"hole-u-bc": "type            zeroGradient;"}


def _create_precice_config(get_config_value) -> Path:
    """Write the preCICE XML. MaMiCo provides both coupling meshes; OpenFOAM maps onto them."""
    xml = PartialXml(get_asset_text("precice-config.xml.template"))
    xml.substitute("mamico-participant", _precice_participant)
    xml.substitute("foam-participant", _precice_foam_participant)
    xml.substitute("foam-mesh-outer", _precice_foam_mesh_outer)
    xml.substitute("foam-mesh-inner", _precice_foam_mesh_inner)
    xml.substitute("macro2md-mesh", _precice_macro2md_mesh)
    xml.substitute("md2macro-mesh", _precice_md2macro_mesh)
    xml.substitute("macro2md-data", _precice_macro2md_data)
    xml.substitute("md2macro-data", _precice_md2macro_data)
    xml.substitute("max-time-windows", get_coupling_cycles(get_config_value))
    xml.substitute(
        "time-window-size",
        get_timesteps_per_coupling_cycle(get_config_value) * MD_TIMESTEP,
    )
    path = Path(get_config_value("output_dir")) / _precice_config_name
    check_if_replacing(path, get_config_value)
    with open(path, "w") as file:
        file.write(xml.get_raw())
    return path


def _create_precice_foam_case(get_config_value) -> Path:
    """Copy the OpenFOAM participant template and substitute its geometry/timing placeholders."""
    src_path = Path(__file__).parent.parent / "assets" / "PreciceFoamCase.template"
    dst_path = Path(get_config_value("output_dir")) / _precice_foam_case_name
    check_if_replacing(dst_path, get_config_value)
    shutil.copytree(src_path, dst_path, dirs_exist_ok=True)

    geometry = _get_geometry(get_config_value)
    delta_t = min(get_timesteps_per_coupling_cycle(get_config_value) * MD_TIMESTEP, 0.25)
    substitutions = {
        "chan": geometry["cfd_size"],
        "hole-x0": geometry["hole_min"][0], "hole-x1": geometry["hole_max"][0],
        "hole-y0": geometry["hole_min"][1], "hole-y1": geometry["hole_max"][1],
        "hole-z0": geometry["hole_min"][2], "hole-z1": geometry["hole_max"][2],
        "shell-min": " ".join(str(v) for v in geometry["shell_min"]),
        "shell-max": " ".join(str(v) for v in geometry["shell_max"]),
        "shell-cut-min": " ".join(str(v) for v in geometry["shell_cut_min"]),
        "shell-cut-max": " ".join(str(v) for v in geometry["shell_cut_max"]),
        "delta-t": delta_t,
        "end-time": get_coupling_cycles(get_config_value) * get_timesteps_per_coupling_cycle(get_config_value) * MD_TIMESTEP,
        "kinematic-viscosity": DYNAMIC_VISCOSITY / MD_DENSITY,
        "wall-velocity": get_wall_velocity(get_config_value),
        # Both participants are started from the output directory
        "precice-config": _precice_config_name,
        "foam-case-name": _precice_foam_case_name,
        "foam-participant": _precice_foam_participant,
        "foam-mesh-outer": _precice_foam_mesh_outer,
        "foam-mesh-inner": _precice_foam_mesh_inner,
        "cell-set-shell": _foam_cell_set_shell,
        "inner-patches": " ".join(_foam_inner_patches),
        "macro2md-data": _precice_macro2md_data,
        "md2macro-data": _precice_md2macro_data,
    }
    substitutions.update(_get_side_conditions(get_config_value))
    substitutions.update(_get_hole_condition(get_config_value))
    for axis, ns in (("nx", geometry["nx"]), ("ny", geometry["ny"]), ("nz", geometry["nz"])):
        for idx, n in enumerate(ns):
            substitutions["%s%d" % (axis, idx)] = n

    for file_path in sorted(dst_path.rglob("*")):
        if not file_path.is_file():
            continue
        text = file_path.read_text()
        for placeholder, value in substitutions.items():
            text = text.replace("{" + placeholder + "}", str(value))
        file_path.write_text(text)
    for script in ("Allrun", "Allclean"):
        (dst_path / script).chmod(0o755)
    return dst_path


def validate(get_config_value) -> str:
    """MaMiCo config validation:
    The analytical CFD solver does not support two-way coupling.
    Only the OpenFOAM solver should be used with two-way coupling.
    """
    key = __name__.split(".")[-1]
    solver = get_config_value(key)
    use_2way_coupling = get_config_value("coupling_2way")
    if solver == "foam" and get_domain_size(get_config_value) != 2:
        return f"OpenFOAM can only be used with the medium domain size."
    # FoamSetup.template is a fixed mesh built for medium + cell size 5
    # (hole 35..65 / 20..50, cell counts 7,6,7 and 4,6,10); any other cell
    # size silently mismatches the MD domain it is supposed to overlap.
    if solver == "foam" and get_config_value("cell_size") != 5:
        return f"The OpenFOAM setup template only matches cell size 5."
    if solver == "analytical" and use_2way_coupling:
        return f"Cannot use two-way coupling with analytical CFD solver."
    if (
        solver == "lb"
        and get_config_value("simulation") != "test"
        and use_2way_coupling
    ):
        return f"Two-way coupling is not stable with LBCouette solver for more than a quick test."
    if (
        solver == "fd"
        and get_config_value("simulation") != "test"
        and use_2way_coupling
    ):
        return f"Two-way coupling is not stable with FDCouette solver for more than a quick test."
    return None


def apply(partial_xml, get_config_value) -> None:
    key = __name__.split(".")[-1]
    solver = get_config_value(key)
    partial_xml.substitute("cfd-solver", solver)
    partial_xml.substitute("viscosity", DYNAMIC_VISCOSITY)
    print("Substituted CFD solver")
    macro_solver_extra = ""
    if solver == "foam":
        path = _create_foam_setup(get_config_value)
        print(f"Also created files under {path}")
        output_dir = get_config_value("output_dir")
        macro_solver_extra = (
            f'foam-setup-directory="{output_dir}"\n' + 'foam-setup-folder="FoamSetup"'
        )
    elif solver == "precice":
        macro_solver_extra = f'precice-config-file="{_precice_config_name}"'
        print(f"Also created file {_create_precice_config(get_config_value)}")
        print(f"Also created files under {_create_precice_foam_case(get_config_value)}")
    partial_xml.substitute("macro-solver-extra", macro_solver_extra)
