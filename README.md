# PV System Analysis (Without Economic Analysis)

**Author:** Derek / Team  
**Location:** Sacramento International Airport (38.7° N, –121.6 E)  
**Last Updated:** June 17, 2025  

---

## Overview

This MATLAB script performs a detailed performance and environmental analysis of a photovoltaic (PV) array at Sacramento International Airport. It covers:

1. **Irradiance Calculations**  
   - Hourly extraterrestrial, global, diffuse, and beam irradiance  
   - Tilted‐surface irradiance components (beam, diffuse, ground-reflected)  

2. **PV Array Performance Modeling**  
   - Cell temperature estimates  
   - Hourly and daily power output  
   - System losses (shading, soiling, inverter, wiring)  
   - Degradation over a 25-year lifetime  

3. **Environmental Impact Assessment**  
   - Annual and lifetime CO₂ emissions reduction  
   - “Trees planted” and “cars removed” equivalents  

4. **Additional Metrics**  
   - Performance Ratio (PR)  
   - Capacity Factor (CF)  
   - Specific yield and energy yield  
   - Loss breakdown and pie-chart visualization  

---

## Prerequisites

- **MATLAB R2021a** or later  
- No specialized toolboxes required (relies on core MATLAB functions)  

---

## Getting Started

1. **Place the script**  
   - Save the code into a file named `PVSystemAnalysis.m` (or your preferred name).  
2. **Adjust site and system parameters**  
   - At the top of the file, modify latitude, longitude, tilt, panel specs, installed capacity, albedo, etc.  
3. **Run the script**  
   - In the MATLAB Command Window:
     ```matlab
     >> PVSystemAnalysis
     ```

---

## File Structure

- **Sections marked by comments:**  
  1. **Parameters Initialization** – define constants, panel specs, site data  
  2. **Solar Irradiance Calculations** – compute hourly `GHI`, `DHI`, `DNI`, and tilted components  
  3. **PV Array Performance Modeling** – temperature correction, power output, loss factors  
  4. **Environmental Impact Assessment** – CO₂ reduction and equivalencies  
  5. **Additional Performance Metrics** – PR, CF, yields, loss analysis  
  6. **Lifetime Projections** – degradation and annual energy over 25 years  
  7. **Summary & Conclusions** – printed results and key figures  

- **Plots Generated:**  
  - Hourly irradiance components  
  - Zenith and incidence angles  
  - Efficiency degradation curve  
  - Annual energy production over time  
  - Annual CO₂ reduction over time  
  - Loss analysis pie chart  

- **Tables Displayed:**  
  - Hourly irradiance & angles  
  - Efficiency and energy production per year  

---

## Key Outputs

- **Daily energy yield** (kWh/day)  
- **Annual energy production** (MWh)  
- **Total CO₂ reduction** (tons over lifetime)  
- **Performance Ratio (PR)** and **Capacity Factor (CF)**  
- **Equivalent trees planted** and **cars removed**  

---

## Customization Tips

- **Change the analysis day** by updating `day = xxx;` (Julian day).  
- **Use real measured irradiance** by replacing `H_actual` with site data.  
- **Extend lifetime** by adjusting `years`.  
- **Include economic analysis** by adding cost and revenue calculations in a new section.  

---

## License & Contribution

Feel free to adapt and extend this script for your own PV performance studies. If you improve or spot issues, please share your updates back to the team.

---

**Enjoy your PV analysis!**
