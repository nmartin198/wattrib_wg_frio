# Incorporating Weather Attribution to Future Water Budget Projections – Frio Basin

This repository provides the source code and other files for the journal article [**Incorporating Weather Attribution to Future Water Budget Projections**](https://doi.org/10.3390/hydrology10120219). The graphical abstract for this article is displayed below.
<br/>
<figure>
    <img src="/assets/Graphical-Abstract.png"
         width="900"
         alt="Graphical Abstract">
</figure>
<br/>

In this study, a custom weather generator (WG) is created and calibrated that includes extreme events. A schematic for the WG framework is displayed below. The purpose of this WG framework is to enable calibration to drought magnitudes and likelihoods, under human induced climate change. After calibration, the attribution constrained WG will provide a climate description that provides increased probability of future occurrence for what was historically severe and extreme drought relative to Coupled Model Intercomparison Project (CMIP), global climate model (GCM) simulation results. 
<br/>
<figure>
    <img src="/assets/WG_Framework.png"
         width="1000"
         alt="Weather Generator with Events Schematic">
</figure>


## Source Code

The WG source code is available at [**src**](https://github.com/nmartin198/wattrib_wg_frio/tree/main/src).


## Calibration Example

An example calibration to a weather attribution study findings is presented in [**calibration**](https://github.com/nmartin198/wattrib_wg_frio/tree/main/calibration).


## Implementation Example

Implementation of the attribution constrained WG to produce an ensemble of stochastic day-to-day weather realizations is provided in [**final**](https://github.com/nmartin198/wattrib_wg_frio/tree/main/final).


## Author

* **Nick Martin** nick.martin@alumni.stanford.edu
<br/>

## License

This project is licensed under the GNU Affero General Public License v.3.0 - see the [LICENSE](LICENSE) file for details.
<br/>
