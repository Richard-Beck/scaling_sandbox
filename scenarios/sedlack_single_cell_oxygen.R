# Parameters for a Sedlack et al.-inspired, single-cell oxygen model.

sedlack_single_cell_oxygen_scenario <- function() {
  solubility_aqueous <- 1.29e-3 # mol / m^3 / mmHg
  list(
    name = "Homogenized single-cell oxygen transport",
    source = "Sedlack et al. (2022), adapted to a uniform intracellular sink",
    geometry = list(
      cell_radii_um = c(x = 30, y = 30, z = 12),
      # Equal-volume radius for the initial radial screening implementation.
      equivalent_spherical_radius_um = (30 * 30 * 12) ^ (1 / 3),
      medium_height_um = 150,
      lateral_half_width_um = 50,
      calibration_spacing_um = 100,
      grid_dx_um = 10,
      grid_dz_um = 3,
      substrate = "glass_no_flux"
    ),
    transport = list(
      diffusion_medium_um2_s = 1900,
      diffusion_cytoplasm_um2_s = 900,
      solubility_medium_mol_m3_mmhg = solubility_aqueous,
      solubility_cytoplasm_mol_m3_mmhg = solubility_aqueous,
      membrane_permeability_um_s = 42 * 1e4
    ),
    consumption = list(
      # Sedlack's mitochondrial-domain value: an intentionally extreme upper bound.
      vmax_mito_mol_m3_s = 10,
      # Baseline calibrated in the 3D Cartesian model to the reported 5% O2,
      # 100-um-spacing mean cytoplasmic pO2 of 21.4 mmHg.
      vmax_effective_mol_m3_s = 0.09869059,
      calibration_target_cytoplasm_po2_mmhg = 21.4,
      radial_screening_vmax_mol_m3_s = 0.194,
      # Former 0.10 x 10 setting: retain only as an extreme sensitivity case.
      vmax_effective_extreme_mol_m3_s = 1.0,
      km_mol_m3 = 1e-4
    ),
    boundary_conditions = list(
      applied_po2_mmhg = c(`0.5% O2` = 3.8, `5% O2` = 38, `10% O2` = 76, `20% O2` = 152),
      screening_po2_mmhg = 38,
      top = "Dirichlet: c = S_medium * pO2_applied",
      lateral = "No flux (single-cell unit-cell symmetry)",
      bottom = "No flux (glass)",
      membrane = "Permeability coupling in pO2, not concentration continuity"
    )
  )
}
