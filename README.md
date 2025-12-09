# Plantarflexor-StaticOptimization

Extends OpenSim static optimization to analyze how muscle morphology affects plantarflexor power during gait. Tests ±1 SD variations in tendon compliance and fiber length using Rajagopal 2016 data. Built for ME EN 6960 - Simulation of Human Movement.

## Purpose

This project investigates how muscle modeling assumptions affect predicted plantarflexor fiber power during walking gait. Specifically, it examines the sensitivity of soleus and gastrocnemius medialis fiber power to variations in:
- Tendon compliance (tendon slack length)
- Optimal fiber length
- Combined morphological parameters (Fmax, Lopt, Lts, pennation angle)

Understanding these sensitivities is critical for:
- Accurate biomechanical analysis of push-off mechanics
- Informing prosthesis and assistive device design
- Quantifying uncertainty in musculoskeletal simulations

This work extends the static optimization framework from Lab 4 of ME EN 6960 (Simulation of Human Movement) using physiologically-grounded parameter variations based on population standard deviations from Rajagopal et al. (2016).

---

## Software Requirements

### Required Software
- **MATLAB** (tested on R2020a or later)
- **OpenSim 4.x** with MATLAB API configured
  - Download: https://simtk.org/projects/opensim
  - Setup guide: https://simtk-confluence.stanford.edu/display/OpenSim/Scripting+with+Matlab

### Required Files (Included in Repository)
- `CustomStaticOptimization_Sensitivity.m` - Main analysis script
- `loadFilterCropArray.m` - Helper function for loading and filtering data
- `subject_walk_adjusted.osim` - Rajagopal2016 musculoskeletal model
- `coordinates.mot` - Joint angle kinematics from walking trial
- `grf_walk.xml` - Ground reaction force data and configuration

---

## How to Use

### Setup
1. **Install OpenSim 4.x** and configure the MATLAB API
```matlab
   % Test OpenSim installation in MATLAB
   import org.opensim.modeling.*;
   Model();  % Should run without error
```

2. **Clone this repository**
```bash
   git clone https://github.com/yourusername/plantarflexor-sensitivity.git
   cd plantarflexor-sensitivity
```

3. **Set MATLAB path** to the repository directory

### Running the Analysis

1. **Open MATLAB** and navigate to the repository directory

2. **Run the main script**
```matlab
   CustomStaticOptimization_Sensitivity
```

3. **What the script does:**
   - Computes kinematics using OpenSim AnalyzeTool (if not already computed)
   - Runs inverse dynamics to calculate joint moments
   - Performs static optimization for 5 conditions:
     - Baseline (model defaults)
     - +1 SD Morphology (all parameters increased)
     - -1 SD Morphology (all parameters decreased)
     - +1 SD Tendon Compliance (longer/more compliant tendon)
     - -1 SD Tendon Compliance (shorter/stiffer tendon)
   - Extracts fiber force, velocity, length, and activation
   - Computes fiber power (force × velocity)
   - Generates comprehensive visualizations and statistics

4. **Expected runtime:** 5-15 minutes depending on your computer

### Output Files

The script generates:
- `sensitivity_results.mat` - All results data structure
- Multiple figures showing:
  - Fiber power time series for each muscle and condition
  - Activation, fiber length, and velocity comparisons
  - Bar charts of peak power and total work
  - Percent change analyses
- Console output with detailed numerical results tables

---

## Results

### Key Findings

**Soleus**
- Baseline peak power: 27.495 W
- +1 SD Morphology: -98.1474% change
- -1 SD Morphology: 1102.42% change  
- +1 SD Tendon Compliance: -82.7879% change
- -1 SD Tendon Compliance: 252.982% change

**Gastrocnemius Medialis**
- Baseline peak power: 36.3863 W
- +1 SD Morphology: 35.7323% change
- -1 SD Morphology: +1086.25% change
- +1 SD Tendon Compliance: 36.6241% change
- -1 SD Tendon Compliance: 606.298% change

## How to Reproduce Results

### Exact Reproduction

1. **Use the same software versions:**
   - MATLAB R2020a or later
   - OpenSim 4.1 or 4.2

2. **Use provided data files** (unchanged):
   - `subject_walk_adjusted.osim`
   - `coordinates.mot`
   - `grf_walk.xml`

3. **Run the script** as described above

4. **Verify outputs:**
   - Check console output tables match reported values
   - Compare generated figures to those in `/results` folder

### Parameter Values Used

From Rajagopal et al. (2016):

**Soleus:**
- Fmax: 6195 ± 1606 N
- Lopt: 4.4 ± 1.0 cm  
- Lts: 27.7 ± 1.0 cm
- Pennation: 21.9 ± 8.0°

**Gastrocnemius Medialis:**
- Fmax: 3116 ± 727 N
- Lopt: 5.1 ± 1.0 cm
- Lts: 39.9 ± 1.1 cm  
- Pennation: 9.5 ± 4.3°


---

## Methods Summary

### Static Optimization Formulation

**Objective:** Minimize muscle activation squared
```
min: Σ(a_i²) + w_reserves × Σ(τ_reserve²)
```

**Subject to:** Moment balance constraint
```
Σ[r_i × F_max,i × a_i × f_l,i × f_v,i × cos(α_i)] + τ_reserve = τ_ID
```

Where:
- `a_i` = muscle activation (0-1)
- `r_i` = muscle moment arm
- `f_l,i` = force-length multiplier
- `f_v,i` = force-velocity multiplier
- `α_i` = pennation angle
- `τ_ID` = inverse dynamics joint moment

### Fiber Power Calculation
```
Fiber Force = a × f_l × f_v × F_max + f_p × F_max
Fiber Power = Fiber Force × Fiber Velocity
```

Where `f_p` is the passive force multiplier.

### Model Properties
- **Muscle model:** Millard2012EquilibriumMuscle
- **Tendon compliance:** Disabled (`ignore_tendon_compliance = true`)
- **Activation dynamics:** Disabled (`ignore_activation_dynamics = true`)
- **Note:** Tendon slack length variations still affect muscle operating point on force-length curve

---

## Limitations and Future Work

### Current Limitations
1. **Static optimization** - Does not capture muscle-tendon dynamics or series elastic effects
2. **Generic model** - Not subject-specific; uses scaled Rajagopal2016 model
3. **Limited validation** - No experimental fiber velocity/force measurements for direct comparison

### Future Directions
1. Enable tendon compliance dynamics to capture realistic muscle-tendon interaction
2. Validate against subject-specific ultrasound measurements
3. Extend to running and other activities
4. Apply individualized morphology for patient-specific analysis
5. Investigate implications for prosthetic foot design

---

## References

1. **Rajagopal, A., et al.** (2016). Full-body musculoskeletal model for muscle-driven simulation of human gait. *IEEE Trans Biomed Eng*, 63(10), 2068-2079.

2. **Neptune, R. R., et al.** (2001). Muscle contributions to whole-body sagittal plane angular momentum during walking. *J Biomech*, 34(6), 783-790.

3. **Dorn, T. W., et al.** (2012). Muscular strategy shift in human running. *J Exp Biol*, 215(11), 1944-1956.

4. **Arnold, E. M., et al.** (2013). Fibre operating lengths of human lower limb muscles during walking. *Phil Trans R Soc B*, 366, 1530-1539.

5. **Zajac, F. E.** (1989). Muscle and tendon: properties, models, scaling, and application to biomechanics and motor control. *Crit Rev Biomed Eng*, 17(4), 359-411.

6. **Uchida, T. K., et al.** (2016). Simulating ideal assistive devices to reduce the metabolic cost of running. *PLOS ONE*, 11(9), e0163417.

7. **Krishnaswamy, P., et al.** (2011). Human leg model predicts ankle muscle-tendon morphology, state, roles and energetics in walking. *PLOS Comput Biol*, 7(3), e1001107.

8. **Farris, D. J., & Sawicki, G. S.** (2012). Human medial gastrocnemius force-velocity behavior shifts with locomotion speed and gait. *Proc Natl Acad Sci USA*, 109(3), 977-982.

9. **Lichtwark, G. A., & Wilson, A. M.** (2006). Interactions between the human gastrocnemius muscle and the Achilles tendon during incline, level and decline locomotion. *J Exp Biol*, 209(21), 4379-4388.

10. **Fukunaga, T., et al.** (2001). In vivo behaviour of human muscle tendon during walking. *Proc R Soc Lond B*, 268(1464), 229-233.

---

## Author

Jack Donohoe  
ME EN 6960 - Simulation of Human Movement  
University of Utah  
Final Project, Fall 2025

---

## Acknowledgments

- OpenSim development team at Stanford University
- Dr. Uhlrich for course guidance and Lab 4 framework
- Rajagopal et al. for providing the musculoskeletal model and population data
