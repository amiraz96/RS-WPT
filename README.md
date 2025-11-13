# On the Deployment of Multiple Radio Stripes for Large-Scale Near-Field RF Wireless Power Transfer

Reference MATLAB implementation for the algorithms introduced in the paper "On the Deployment of Multiple Radio Stripes for Large-Scale Near-Field RF Wireless Power Transfer" available at [arXiv:2508.21640](https://arxiv.org/abs/2508.21640). The code explores how reconfigurable surfaces (RSs) can be used to improve multi-user wireless power transfer (WPT) by optimizing the aperture geometry, deployment strategy, and transmit beamforming.

This README summarizes the repository structure, external dependencies, and the workflow for reproducing the numerical studies reported in the paper.

## Key capabilities

- **Successive convex/geometric optimization of RS layouts.** Core routines such as `Do_SCA.m`, `Do_SGP.m`, `Do_SGP_Line.m`, and `Do_SGP_Polygon.m` jointly optimize RS element coordinates and power splitting across users under minimum-spacing and length constraints using CVX-based solvers.
- **Beamforming for optimized apertures.** `Do_Beamforming.m` implements both semidefinite programming (SDP) and maximum-ratio transmission (MRT) strategies to evaluate the harvested power once the RS placement is fixed.  Channel coefficients are generated through near-field spherical propagation models in `Do_Channels.m`.
- **Flexible geometry generation.** Helper scripts create rectangular (`Create_Rect_RS.m`), fully-digital (`Create_FD.m`), line, and polygonal apertures that are reused across simulations and deployment studies.
- **Monte Carlo evaluation pipelines.** Simulation drivers such as `SGP_simulations.m`, `SCA_simulations.m`, and `MonteCarlo_Shape_Line.m` iterate over user drops, aperture sizes, and carrier frequencies to generate the datasets that underpin the figures in the paper.

## Repository layout

```
RS-WPT/
├── Create_FD.m                     # Grid generator for fully-digital baselines
├── Create_Rect_RS.m                # Rectangular RS initialization
├── Do_*.m                          # Optimization, channel, and beamforming primitives
├── *_simulations.m                 # Entry points for reproducing SGP/SCA case studies
├── *_MonteCarlo.m                  # Monte Carlo evaluation scripts across scenarios
├── plot_*.m                        # Post-processing and visualization helpers
└── README.md                       # Project documentation (this file)
```

A number of scripts expect supporting `.mat` files (for example `final_locs_defined.mat`, `final_locs_defined_montecarlo_100.mat`, and pre-computed optimization results). These datasets are generated during earlier simulation stages and must be placed in the repository root before running the corresponding experiments.

## Requirements

- MATLAB R2021b or later (earlier versions may work but are untested).
- The [CVX toolbox](http://cvxr.com/cvx/) with a solver that supports geometric programming and semidefinite programming; install and run `cvx_setup` before launching any of the optimization scripts.
- MATLAB Parallel Computing Toolbox for the `parfor` loops used in several simulation drivers (for example `SGP_simulations.m`). You can switch to `for` loops if the toolbox is unavailable.

## Reproducing the main experiments

The numerical results in the paper are organized around three stages: generating user deployments, optimizing RS placements, and evaluating the resulting beamforming strategies.

1. **Prepare user deployment datasets.** Scripts such as `clustering_simul.m` and `cluster_power_monte_opt.m` create Poisson-based user layouts and cluster assignments saved to `.mat` files. Ensure the generated files (e.g., `final_locs_defined.mat`, `final_locs_defined_montecarlo_100.mat`) are available before running the optimization routines.
2. **Optimize RS geometries.**
   - Run `SGP_simulations.m` to compute successive geometric programming layouts over different carrier frequencies for a fixed RS aperture length. Results are stored in `SGP_OPT_Results_overF_*.mat` files.
   - Run `SCA_simulations.m` for the successive convex approximation baseline using the same deployment assumptions. Outputs are saved as `SCA_OPT_Results_overF_*.mat`.
   - Geometry-specific routines (`Polygon_simulations.m`, `Line_simulations.m`) call `Do_SGP_Polygon.m` or `Do_SGP_Line.m` to explore alternative shapes.
3. **Monte Carlo evaluation.**
   - Use `SGP_SCA_MonteCarlo.m` to average the harvested power of SGP and SCA designs across many random user drops, reporting both SDP and MRT beamforming metrics.
   - `MonteCarlo_Shape_Line.m` and `Angle_Line_MonteCarlo.m` compare polygonal or angular line layouts against baselines over varying aperture lengths and frequencies.
   - The `Do_Cluster_BF.m` and `Do_Deploy.m` utilities evaluate per-cluster deployment strategies by contrasting RS and fully-digital apertures under identical user sets.

Each simulation script saves intermediate `.mat` files that the subsequent analysis and plotting functions consume. Inspect the individual scripts to adjust carrier frequencies, aperture sizes, and the number of RS elements according to the scenarios described in the paper.


## Tips for extending the code

- The optimization functions accept arbitrary 3D user coordinates; adjust the `M_locs` tensors in the `.mat` files to evaluate new deployment patterns.
- Constraints on inter-element spacing and RS aperture length are explicitly parameterized (`inter_dist`, `max_L`), making it straightforward to study different hardware capabilities.【F:RS-WPT/Do_SCA.m†L21-L46】【F:RS-WPT/Do_SGP_Line_Angle.m†L1-L55】
- Additional beamforming policies can be integrated by reusing the channel matrices returned by `Do_Channels.m` and adding new processing blocks after the optimization stage.

## Citing

If you use this code, please cite the associated article. A BibTeX entry is provided below.

```bibtex
@misc{azarbahram2025deploymentmultipleradiostripes,
      title={On the Deployment of Multiple Radio Stripes for Large-Scale Near-Field RF Wireless Power Transfer}, 
      author={Amirhossein Azarbahram and Onel L. A. López and Petar Popovski and Matti Latva-aho},
      year={2025},
      eprint={2508.21640},
      archivePrefix={arXiv},
      primaryClass={eess.SP},
      url={https://arxiv.org/abs/2508.21640}, 
}
```

## License

The repository inherits the usage rights granted by the original authors. Refer to the paper or contact the authors for explicit licensing terms.
