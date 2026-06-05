This Markdown file is to be used for sessions working on HICAR, which is an atmospheric fast dynamic downscaling model located in this very directory.

The code can be found in ./src/
The model does not solve Navier-Stokes equations but applies a series of parameterizations to the forcing input. It does calculate physics and advection.
Model runs are usually executed in the scratch directory /capstor/scratch/cscs/msesselm/HICARsnow/HEFEX_sims/input where you can look for the most recently changed .out .err and .nml files if you need context.
This model is usually run on CSCS alps.eiger CPUs.

You can compile it using the following commands, do not compile it in any other way and always ask for permission before compiling the model yourself:

uenv run prgenv-gnu/24.11:v2 --view default
cd ~/HICAR_gpu/build
rm -fr ~/HICAR_gpu/build/*
cmake ../ -DFSM=ON -DFSM_DIR=/users/msesselm/FSM2/FSM_SOURCE_CODE -DSRUN_FLAGS="-A s1330" -DMODE=debug
make -j 36
make install


This is for debugging purposes only which will usually be what you are working with. If you want to compile in release mode use -DMODE=release.

In case you need a python installation to check on .nc files or similar, load the bashrc and activate a conda environment for output processing like such:

source ~/.bashrc
conda activate op_plot_slim

which should feature everything you need to do so.

## Debugging Allowances (while model keeps crashing)

During active debugging, the following actions are permitted without asking for permission:

1. **Compile the model** (debug and/or release mode) using the commands above
2. **Clean output files**: `rm /capstor/scratch/cscs/msesselm/HICARsnow/HEFEX_sims/output/HEF2/nested/kenda_sim4/*`
3. **Clean restart files**: `rm /capstor/scratch/cscs/msesselm/HICARsnow/HEFEX_sims/restart/HEF2/nested/kenda_sim4/*`
4. **Submit jobs**: `cd /capstor/scratch/cscs/msesselm/HICARsnow/HEFEX_sims/input && sbatch HICAR_HEF2_KENDA.job`
5. **Read output/error logs** from the scratch directory
6. **Edit source code** to fix bugs or add/modify diagnostics

**IMPORTANT**: These allowances are ONLY valid while the model is crashing. Once we get interesting simulation results, do NOT delete output or restart files without explicit permission.

**WARNING**: A full clean rebuild (`rm -fr build/*` followed by `cmake`) will delete `build/_deps/` and re-fetch NoahMP source code, **wiping out any edits to files under `build/_deps/noah_mp-src/`**. To preserve such edits across rebuilds, either:
- Only rebuild with `make` (no `rm -fr build/*` or `cmake` re-run), OR
- Re-apply edits after `cmake` but before `make`, OR
- Copy modified NoahMP files to a persistent location outside `build/` and copy them back after cmake

## Simulation configuration

- **HEF2 KENDA** is the primary test case for debugging. 4-domain nested run (1600m → 400m → 100m → 50m), Domain 4 has the 50m nest with ~0.5m lowest model levels and 50 vertical levels.
- Output directory: `/capstor/scratch/cscs/msesselm/HICARsnow/HEFEX_sims/output/HEF2/nested/kenda_sim4/`
- Job script: `cd /capstor/scratch/cscs/msesselm/HICARsnow/HEFEX_sims/input && sbatch HICAR_HEF2_KENDA.job`

## Changes on `gpu` branch since last commit (9f5f349d)

All changes below are uncommitted modifications on the `gpu` branch, aimed at making the model stable with ultra-high vertical resolution (~0.5m lowest levels, 50m horizontal nests).

### 1. FIXED: Pressure gradient inversions from float32 accumulation (domain_obj.F90)

**Problem**: In `apply_forcing`, pressure was updated incrementally: `p = p + dqdt * dt`. With ~70000 Pa base pressure and ~0.01 Pa/step tendencies, float32 rounding after O(2000) steps destroyed the ~16 Pa gaps between thin near-surface levels, causing pressure inversions.

**Fix**: Direct computation from base values stored at each forcing interval: `p = p_base + dqdt * elapsed_time`. Module-level `p_var_base`, `p_fh_base` (saved at interpolation time), and `p_elapsed` (double-precision accumulator) eliminate accumulation error entirely. See `apply_forcing` and `interpolate_forcing` in `src/objects/domain_obj.F90`.

Also disabled `adjust_pressure_temp` calls in `interpolate_forcing` — the vinterp extrapolation fix (below) makes them unnecessary and they would double-correct.

### 2. FIXED: Vertical interpolation clamping → extrapolation (vinterp.F90)

**Problem**: `vLUT` in `src/utilities/vinterp.F90` clamped below-grid and above-grid points to the nearest forcing level (weight=0.5, same level twice). This produced flat pressure/temperature bands at the bottom of thin near-surface columns.

**Fix**: Changed to bilinear extrapolation using the two nearest levels. Below the grid: extrapolates from levels 1,2. Above the grid: extrapolates from levels N-1,N. Weights > 1 and < 0 provide linear extrapolation.

### 3. FIXED: NaN from psim=psih=0 in YSU PBL (sfc_sfclayrev.F90, module_bl_ysu.F90, pbl_ysu.F90)

**Root cause**: With ~0.5m first-level height (ZA), the bulk Richardson number `br` becomes tiny (~5.7e-8). The `zolri` secant method in `sfc_sfclayrev.F90` uses fixed brackets [0, 5] for stable conditions. With such a tiny `br`, the secant step overshoots from 5 to 0 in a single iteration (float32 underflow), returning `zol=0` exactly. This makes `psim=psih=0`, which causes `0/0 = NaN` in YSU's `zol1 = br*fm*fm/fh` computation, cascading through the entire column via the tridiagonal solver.

**Fix (two parts)**:
- **`sfc_sfclayrev.F90` lines 1211-1217**: Analytical shortcut for `|ri| < 1e-4`: `zolri = ri * log((z+z0)/z0)`. This is the neutral-limit approximation, accurate to machine precision for small Ri, bypassing the secant method entirely.
- **`module_bl_ysu.F90` line 911 / `pbl_ysu.F90` line 858**: Safety guard: `if(abs(fh) < 1e-6) zol1=0` prevents 0/0 even if psih=0 slips through.

Additional defensive changes in `sfc_sfclayrev.F90`:
- NaN guards in `psim_stable`, `psih_stable`, `psim_unstable`, `psih_unstable` lookup functions
- Bounds checks for negative `nzol` index in lookup tables (`nzol >= 0` guard)
- Division guard in `zolri2` when `psix2 ≈ 0`

### 4. Double-precision tridiagonal solvers (module_bl_ysu.F90, pbl_ysu.F90, pbl_driver.F90)

**Problem**: With ultra-thin layers, the tridiagonal matrix entries span many orders of magnitude. Float32 pivot divisions in the Thomas algorithm caused catastrophic cancellation and overflow.

**Fix**: All Thomas algorithm solvers now use `real(8)` for pivot (`fk_d`) and denominator (`denom`) computations, converting back to `real` for storage. Applied to:
- `tridi1n` and `tridin_ysu` in `module_bl_ysu.F90` (GPU versions)
- `tridi2n` and `tridin_ysu` in `pbl_ysu.F90` (CPU versions)
- `pbl_scalar_diff` Thomas algorithm in `pbl_driver.F90`

### 5. nan_stop diagnostic infrastructure (debug_utils.F90, time_step.F90)

Added `nan_stop` subroutine that checks pressure, potential_temperature, water_vapor, density, u, v, w for NaN. On detection: prints variable name, (i,k,j) location, PE rank, writes debug .nc file, then `error stop`. Added calls throughout the physics loop in `time_step.F90`:
- after apply_forcing, after rad, after sfc+lsm, after pbl, after mp (microphysics), after integrate_physics (with sub-checks after rad_apply_dtheta, lsm_apply_fluxes, pbl_apply_tend), before/after update_winds.

**Note**: These `nan_stop` calls are diagnostic overhead. They should eventually be removed or gated behind a debug flag once the model is stable. Currently they remain active unconditionally.

### 8. Soft warning threshold for variable bounds checks (meta_data_h.F90, default_output_metadata.F90, variable_obj.F90, debug_utils.F90)

Added a `warnval` field to `root_var_t` in `meta_data_h.F90` — a non-fatal warning threshold separate from `minval` (fatal). When `check_var` detects a value below `warnval` but above `minval`, it prints a `SOFT WARNING` but does not abort. Currently set for water_vapor only:
- `warnval = -1e-10` (soft warning)
- `minval = -1e-3` (fatal error)

This allows monitoring small negative qv artifacts from advection/interpolation during nest handoff without aborting the simulation.

### 9. PETSc wind solver: conditional ILU / BoomerAMG preconditioner (wind_iterative_petsc.F90)

**WARNING**: This change modifies the core wind solver — the most performance-critical and numerically sensitive component of the model. To fully revert to the original unconditional ILU behavior: (a) remove the `use_amg` module variable and all `if (use_amg)` conditionals, (b) set `DMSetMatType(da,MATIS,ierr)` unconditionally, (c) restore `MatSetLocalToGlobalMapping` unconditionally, (d) move PC setup back into `ksp_setup` with `PCFactorSetUseInPlace`/`PCFactorSetReuseOrdering`. Alternatively, pass runtime flags to override: `-pc_type ilu -pc_factor_in_place`.

**Problem**: The PETSc variational wind solver used ILU (Incomplete LU) preconditioning with MATIS matrix format. With ultra-high resolution (10m horizontal, 0.5m lowest vertical level), the matrix coefficients span 5-6 orders of magnitude (vertical terms ~4·α² near surface vs ~10⁻⁵ near model top; horizontal terms ~0.01). ILU is a single-level algebraic method that cannot handle this multi-scale structure, leading to 2000+ iterations, non-convergence, and PETSc internal crashes (segfault at 09:00 sim time with alpha=0.7). However, ILU is fast and efficient for standard resolutions where the condition number is manageable.

**Solution**: Automatic preconditioner selection based on minimum vertical layer thickness. If the thinnest layer is < 1.0m, use BoomerAMG (algebraic multigrid); otherwise use ILU. This preserves the proven ILU behavior for all standard simulations while enabling AMG only when the extreme vertical stretching demands it.

**Implementation (five parts)**:

1. **Module-level flag** (line 36): `logical :: use_amg = .False.` — set per-nest in `init_iter_winds_petsc`.

2. **Resolution check** (`init_iter_winds_petsc`, after `init_module_vars`): Computes `min_dz = minval(advection_dz)` over the domain's tile. If `min_dz < 1.0`, sets `use_amg = .True.`. Prints which preconditioner was selected for each nest.

3. **Conditional matrix format and L2G mapping** (`init_iter_winds_petsc`):
   - `use_amg=.True.`: MATAIJ (standard distributed CSR, compatible with AMG). No L2G mapping (conflicts with Hypre's ParCSR conversion).
   - `use_amg=.False.`: MATIS (per-subdomain local matrices). L2G mapping set via `MatSetLocalToGlobalMapping` (required by MATIS).

4. **Conditional preconditioner**:
   - **ILU (default)**: `ksp_setup` (lines 1045-1066) sets ILU on the default PC, exactly matching the original code: `KSPSetFromOptions` → `KSPGetPC` → `PCFactorSetUseInPlace` → `PCFactorSetReuseOrdering`. These calls operate on PETSc's default PC which already has ordering initialized — calling `PCSetType(PCILU)` on a fresh PC was found to break this (see debugging history).
   - **BoomerAMG (override)**: `init_iter_winds_petsc`, only when `use_amg=.True.`: `PCSetType(PCHYPRE)` + `PCHYPRESetType('boomeramg')`. This overrides the default ILU with algebraic multigrid.

5. **`ksp_setup` unchanged from original** (lines 1045-1066): Creates the KSP, sets solver type (KSPPIPEGCR), and configures ILU on the default PC — identical to the original committed code. The BoomerAMG override happens per-nest in `init_iter_winds_petsc` where the domain's actual resolution is known.

**Nested run behavior**: In a 4-domain nested run (e.g., 1600m→400m→100m→50m), each domain gets its preconditioner chosen independently when `init_iter_winds_petsc` is called. Only the finest nest (if dz_min < 1m) uses BoomerAMG; coarser nests use ILU (untouched default from `ksp_setup`).

**BoomerAMG requires** `FI_CXI_RX_MATCH_MODE=hybrid` in the job script (Slingshot LE exhaustion from AMG's heavy MPI communication during setup). This is only needed when BoomerAMG is active.

**Original code (for revert reference)**:
```fortran
! Original ksp_setup (all PC config was here):
call KSPSetFromOptions(ksp(i),ierr)
call KSPGetPC(ksp(i),precond,ierr)
call KSPSetReusePreconditioner(ksp(i),PETSC_FALSE,ierr)
call PCFactorSetUseInPlace(precond,PETSC_TRUE,ierr)
call PCFactorSetReuseOrdering(precond,PETSC_TRUE,ierr)

! Original init_iter_winds_petsc:
call DMSetMatType(da,MATIS,ierr)              ! unconditional MATIS
call DMGetLocalToGlobalMapping(da,isltog,ierr) ! unconditional L2G
call MatSetLocalToGlobalMapping(arr_A,isltog,isltog,ierr)
```

**Debugging history**:
- First attempt (unconditional MATAIJ + BoomerAMG + MatSetLocalToGlobalMapping): SEGV on first PETSc solve (job 7318276). All ranks crashed with signal 11 during "Updating initial winds". Cause: L2G mapping conflicted with Hypre's ParCSR conversion.
- Second attempt (unconditional MATAIJ + BoomerAMG, L2G mapping removed): libfabric CXI LE exhaustion (SIGABRT). Fix: added `FI_CXI_RX_MATCH_MODE=hybrid` to job script.
- Third attempt (with FI_CXI_RX_MATCH_MODE=hybrid): Model ran 12 simulated hours (vs 2:45 with ILU/alpha=1.0 and 8:45 with ILU/alpha=0.7). Diverged at hour 12 with reason -9 (KSP_DIVERGED_NANORINF). v winds at 13 m/s (vs 42 m/s with ILU).
- Fourth iteration: Refactored to conditional ILU/BoomerAMG based on min_dz, with per-nest configuration for nested runs.
- Fifth attempt (conditional, with explicit `PCSetType(PCILU)` for ILU path): All nests returned `reason=0` (`KSP_CONVERGED_ITERATING`) — solver never actually solved. PETSc error: `Ordering type cannot be null`. Cause: calling `PCSetType(PCILU)` creates a fresh PCILU without ordering initialized; subsequent `PCFactorSetReuseOrdering` fails. The original code operated on PETSc's default PC (which has ordering initialized by `KSPCreate`), never calling `PCSetType` explicitly.
- Sixth iteration (current): Restored original `ksp_setup` with ILU on default PC. Only the BoomerAMG path calls `PCSetType` to override. ILU path is now identical to the original committed code.

PETSc 3.22.1 installation at `/capstor/store/cscs/userlab/s1329/Shared/Dependencies/hicar/view` was built with `--with-hypre=1` (Hypre 2.32.0), so BoomerAMG is available.

### 6. Minor: FSM2 radiation function rename (PHYSICS_interface.F90)

Changed `call RADIATIONNN(...)` to `call RADIATION(...)` — the function was renamed upstream.

### 7. FIXED: FSM stale theta overwriting NoahMP soil moisture (sm_FSM.F90)

**Problem**: FSM's internal `theta` (volumetric soil moisture) was read from the domain only once at initialization/restart, but written back every timestep for snow cells. Since FSM never modifies theta (it's read-only in THERMAL.F90), this overwrote NoahMP's evolved `soil_water_content` with stale values. At meltout, this caused SLW > SM inconsistency → phantom water in NoahMP's water budget check.

**Fix**: Added re-read of `theta` and `Tsoil` from domain at the start of every `sm_FSM()` call (with `!$acc update host`). The writeback now writes back the current value (harmless no-op), and FSM's thermal calculations use up-to-date soil moisture.

## Current Problems

### Problem A: NoahMP Water Budget Error on FSM2 domains

**Status**: Likely root cause identified and fixed. Job 7250341 has run 3+ simulation days without triggering the water budget error (previous runs crashed within hours).

**Error**: `STOP Error: Water budget problem in NoahMP LSM` — NoahMP's internal water balance check detects >0.1 mm/timestep imbalance.

#### Likely root cause: FSM writes stale theta, overwriting NoahMP's evolved soil moisture

**File**: `src/physics/sm_FSM.F90`

FSM's internal `theta` array (volumetric soil moisture) was only read from the domain once — at initialization/restart (`sm_FSM_init`, line 287). The per-timestep `sm_FSM()` routine never re-read theta from the domain, but DID write it back every step (line ~542) for snow-covered cells.

FSM never modifies theta during its physics (confirmed: THERMAL.F90 reads it for conductivity/heat-capacity, SOIL.F90/SNOW.F90/EBALSRF.F90 don't reference it). So the writeback was overwriting the domain's `soil_water_content` with a stale initialization-time value, clobbering whatever `noahmp_soil_only` had correctly evolved.

**Hypothesised bug sequence at meltout:**
1. Snow cell: FSM writes stale theta → domain's `soil_water_content`
2. `noahmp_soil_only` normally fixes this (reads stale SM + current SLW, runs solver, writes both back)
3. At meltout (SWE→0, melt_basal=0): FSM writes stale theta, but `noahmp_soil_only` gate (`if (snow_height <= 0 .and. melt_basal <= 0) cycle`) SKIPS the cell
4. Domain left with: SM = stale theta (from init/restart), SLW = evolved value from last solver run
5. Next step: XLAND = land (no snow), NoahmpMain runs full physics
6. WaterMain: `SoilIce = max(0, SM - SLW)` → if SLW > SM (stale), SoilIce = 0 → phantom water = (SLW - SM) × dz × 1000

**Numerical match**: SM=0.321428, SLW=0.324283 → phantom = (0.324283-0.321428) × 0.1 × 1000 = 0.2855 mm. Observed: +0.285459 mm.

**Secondary issue**: Stale theta also meant FSM's THERMAL.F90 computed soil thermal conductivity and heat capacity using outdated soil moisture values for the entire simulation.

#### Fix applied (2026-05-31)

In `src/physics/sm_FSM.F90`, added re-read of `theta` and `Tsoil` from the domain at the start of every `sm_FSM()` call (after existing snow state reads, before `FSM_PHYSICS`):

```fortran
!$acc update host(domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d, &
!$acc& domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d)
do i=1,kSOIL_GRID_Z
    Tsoil(i,:,:) = TRANSPOSE(domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d(its:ite,i,jts:jte))
    theta(i,:,:) = TRANSPOSE(domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d(its:ite,i,jts:jte))
enddo
```

This ensures:
- FSM uses current soil moisture for thermal property calculations
- The theta writeback at line ~542 writes back the current domain value (no-op, no overwrite)
- The Tsoil re-read means FSM's soil heat equation starts from `noahmp_soil_only`'s latest temperature

#### Diagnostic history (for reference)

Enhanced diagnostics confirmed the phantom was NOT inside the Richards solver (0 triggers at 0.001mm threshold) and GroundWaterTopModel only REMOVED water (all 13,741 GW diagnostic triggers were negative). This pointed to the phantom originating upstream of WaterMain — specifically from the `SoilIce = max(0, SM-SLW)` clamp at WaterMain entry when SLW > SM.

#### Configuration details
- `nmp_opt_runsrf = 1` (TOPMODEL surface runoff)
- `nmp_opt_runsub = 1` (SIMGM groundwater → DrainSoilBot=0 in Richards, sealed bottom)
- `NumSoilLayer = 4`, `NumIterSoilWat = 6` (fixed), `TimeStepFine = SoilTimeStep/6`
- `NumSoilTimeStep = 1` (single main soil timestep per NoahMP call)

### Problem B: Restart NaN crash (separate from water budget)

**Status**: Diagnostic added, user was going to test. Results not yet reviewed.

**Error**: NaN in potential_temperature after apply_forcing when restarting from checkpoint.

**Diagnosis**: The NaN is in theta's dqdt (tendency), not the data itself. The `p_var_base` and `p_fh_base` arrays (used for the float32 accumulation fix in `domain_obj.F90`) were not being saved/restored on restart, which caused corrupted pressure tendencies. Initial fix attempt (reset_pressure_base after restart) did not resolve it. Added dqdt NaN diagnostic check in `time_step.F90` before apply_forcing.

### Problem C: PETSc wind solver non-convergence and wind field corruption (10m Alex)

**Status**: Testing mitigation (alpha=0.7, 3000 iterations). Job 7313517 running.

**Simulation**: HEF2 Alex 10m — single-domain 10m horizontal resolution run, 50 vertical levels, ~0.5m lowest level, forced by 50m parent output (`HEF2_Alex_50m_final_2023-08-15_00-00-00.nc`). Namelist: `HICAR_new_HEF2_Alex_10m.nml`.

#### Symptoms

1. **PETSc non-convergence**: At 02:45 sim time, both outer wind solver iterations diverged (reason -3 = KSP_DIVERGED_ITS). w_grid jumped from ~4.5 to 42.6 m/s at (i=159, j=615, k=17).
2. **Convergence ≠ correct winds**: At 03:00, PETSc *converged* (776 and 891 iterations) but w_grid was **still 42.5 m/s** at the same location. The solver converged to a physically unreasonable solution.
3. **Self-resolution at 03:15**: New forcing data brought w_grid back to 4.7 m/s at a different location. The 42.6 m/s episode was transient.
4. **Coarser domain showed same pattern**: A similar simulation on a coarser (~50m) domain ran 23h without crashing, but the wind field was completely corrupted throughout the whole atmosphere. PETSc converges on the coarser grid but to poor solutions.
5. **Original NaN crash at 07:00** (previous run): `NaN_STOP [after apply_forcing] potential_temperature at i=8 k=9 j=5 PE=46`. This NaN is likely a downstream consequence of corrupted wind fields feeding into physics.

#### Root cause analysis

The variational wind solver (`src/physics/wind_iterative_petsc.F90`) finds a Lagrange multiplier field `lambda` and computes wind corrections:
```
u += 0.5 * (lambda(i) - lambda(i-1)) / dx / rho_u          (line 400)
w += 0.5 * alpha^2 * dlambdz / jaco / rho_w                  (line 419)
```

Three compounding factors at 10m resolution:

1. **dx=10m amplifies horizontal corrections**: The u-correction divides by `dx`. Lambda differences producing 1 m/s at 50m produce **5 m/s** at 10m.
2. **Ill-conditioned matrix**: Vertical stretching (0.5m near surface → hundreds of meters aloft) creates extreme aspect ratios. Condition number scales ~(L/dx)^2, making the 10m system ~25× worse than 50m.
3. **No bounds on corrections**: `calc_updated_winds` (lines 395-422) applies lambda-gradient corrections **without any clipping or limiting**. Spurious large lambda gradients → unbounded wind corrections.

When PETSc fails to converge (lines 131-145): `calc_updated_winds` is skipped entirely — **no fallback, no partial solution used**. The domain retains the pre-solve wind state. `balance_uvw` (line 807 of `wind.F90`) then recomputes w from u,v divergence regardless. When PETSc converges but to a poor solution, the bad wind corrections are applied directly.

#### Solver configuration

| Parameter | Value | File |
|---|---|---|
| Solver type | KSPPIPEGCR | wind_iterative_petsc.F90:1052 |
| Preconditioner | Default ILU, rebuilt every call | wind_iterative_petsc.F90:1058-1061 |
| rtol (relative) | 1e-10 (hardcoded) | wind_iterative_petsc.F90:103 |
| abstol (absolute) | 1e-5 (hardcoded) | wind_iterative_petsc.F90:104 |
| dtol (divergence) | 1000.0 (hardcoded) | wind_iterative_petsc.F90:105 |
| max iterations | `wind_solver_iterations` from namelist | wind_iterative_petsc.F90:106 |
| alpha_const | from namelist | controls vertical/horizontal correction weight |

Note: iteration count progression showed stress building: 434/257 at 02:00 → 497/271 at 02:15 → **1021/1035** at 02:30 → diverged at 02:45.

#### Mitigation test 1: alpha=0.7, 3000 iterations (job 7313517/7318074)

- `alpha_const = 0.7` (was 1.0): w-correction scales as alpha^2, so 0.7^2 = 0.49 → halves vertical corrections
- `wind_solver_iterations = 3000` (was 1500): doubles headroom for near-convergence cases

**Result**: Model reached 08:45 sim time (vs 02:45 crash with alpha=1.0, and 07:00 NaN with original settings). No NaN, no PETSc divergence, w_grid stayed reasonable (3-5 m/s). However, PETSc iteration counts crept up as daytime convection intensified (2104 at 08:30, 2491 at 08:45) and PETSc crashed internally at ~09:00 (segfault in ILU factorization, likely zero pivot from extreme coefficient ratios).

#### Fix applied: BoomerAMG preconditioner (Change 9)

Replaced ILU preconditioner with Hypre BoomerAMG (algebraic multigrid) and switched matrix format from MATIS to MATAIJ. See Change 9 above for details. This addresses the root cause (ILU cannot handle the multi-scale coefficient structure) rather than just the symptoms.

#### Potential further fixes if BoomerAMG insufficient

1. **Post-solve wind clipping** in `calc_updated_winds`: cap max correction per iteration (e.g., |Δu| < 20 m/s)
2. **Relax abstol** from 1e-5 to 1e-3 (hardcoded in wind_iterative_petsc.F90:104)
3. **Use partial solution on divergence**: When PETSc exceeds max iterations, use whatever iterate exists (with clipping) instead of discarding entirely
4. **Tune BoomerAMG parameters**: coarsening type, interpolation type, number of levels, smoother type
5. **Scale the system** to normalize the extreme aspect ratios before solving

### Problem D: NaN in potential_temperature on 10m Alex simulation

**Status**: Likely downstream consequence of Problem C (wind solver corruption). Waiting for wind solver fix to see if NaN resolves.

**Error**: `NaN_STOP [after apply_forcing] potential_temperature at i=8 k=9 j=5 PE=46` at 07:00 sim time.

**Analysis**:
- PE 46 is interior — `apply_forcing` doesn't modify theta there → NaN from previous timestep's physics
- dqdt NaN check before apply_forcing passed → dqdt arrays clean
- Forcing file is completely clean (comprehensive analysis: no NaN, no Inf, smooth evolution)
- k=9 ≈ 7.96m AGL, within interpolation range
- Added `nan_stop(domain, "after mp")` at `time_step.F90:470` to narrow source
- Previous run hit walltime at 04:15, never reached 07:00 to test the new diagnostic
- The NaN most likely originates from extreme winds (from PETSc solver, Problem C) feeding into microphysics or PBL and causing numerical overflow in thin layers

**Cold-start artifact**: Starting at 05:00 without spin-up caused LW radiation budget crash (TemperatureCanopy oscillation 206K→153K). This is the surface-atmosphere coupling oscillation, not the real bug. Running from 00:00 avoids it.

### Previous Problem (RESOLVED): Surface-atmosphere coupling oscillation

This was a 2Δt temperature oscillation on the 50m Domain 4 with ~0.5m lowest levels. The oscillation grew exponentially (±5K → ±30K → NaN over ~5h). Root cause: thermal mass of ultra-thin near-surface layers too small for explicit surface coupling. The same oscillation was observed when cold-starting the 10m Alex simulation at 05:00 (no spin-up), confirming the mechanism. Starting from 00:00 with spin-up avoids it. The fundamental issue (explicit coupling + ultra-thin layers) remains unresolved but is not the immediate blocking problem.

## Files modified in build/_deps (WILL BE LOST ON CLEAN REBUILD)

None currently. Previous water budget diagnostics in `BalanceErrorCheckMod.F90` and `SoilWaterMainMod.F90` were removed after the root cause was identified and fixed in `src/physics/sm_FSM.F90`.
