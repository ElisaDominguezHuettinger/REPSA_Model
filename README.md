# REPSA_Model
A data-based mathematical model predicts dynamic trajectories of a protected plant community under different conservation regimes


•	Notebook with the full symbolic analysis of the 4 models, analytical derivation of cases in which we observe different stability regimes 
  (extension, local extension, coexistence; mono / metastability) for each of the models: Analisis_simbolico_estabilidad.mlx (in Spanish)
  
  
•	Code to run example simulations of the 4 models for the different cases (as explained above):  Simulations_4_models_different_cases.m 


•	Code to analyse and plot phase space with example trajectories for model 1A (Native only, with negative interactions):  Example trajectories:  DosSP_c_int_n.m 


•	Code in which by manually tuning parameters we find a reasonable fit to data of model 1B: DosSP_c_int_p_OPTIMIZACION.m 


•	Code to analyse and plot phase space with example trajectories for model 1B (Native only, with positive interactions); 
  also a visualization of a bifurcation diagram (with an arbitrarily chosen bifurcation parameter) is given:  DosSP_c_int_p.m
  
  
•	Code with the odes of model 2A (Native and exotic, with negative interactions (Eq. 2A)  m3SPintNEG.m Note, code contains no further análisis of this model, only the equations


•	Symbolic analysis of model 2B (Native and exotic, with positive interactions) TresSP_c_int_p_nat_n_exo.m Note: no plotting of phase space


•	Simulates 3D model (two native with one exotic, with positive interactions between natives) with empirical data; 
  for a parameter set that shows reasonable fit to data, depending on initial conditions we observe convergence to either coexistence or extension of local species with, 
  respectively, extension or dominance of exotic species: TresSP_c_int_p_nat_n_exo_OPTIMIZACION.m
  
  
•	Global optimization of model 2B (aka 4) i.e. positive interactions between native and negative with exotic: OptimizationModel4_TresSP_con_int_pos_nat_n_exo.m 


•	Runs several trajectories (initial conditions) of the model 2B (native species with positive interactions + exotic species) with given parameters, 
  and colours the trajectories depending on where they converge to, to visualize in 3D the regions of attraction of the two stable states: AtractingRegions.m
  --> to be substituted by Pablos code (visualize / calculate / estimate / show the sepparatrix in a more elegant form)
  
 * Stability analysis and numerical integration of model II-A.   m3SPintNEG.m


* generation of matrices with qualitative behaviours: 
m3SPintNEG_Proportion_of_behaviours.m (model IIB)
TresSP_c_int_p_nat_n_exo_Proportion_of_behaviours.m (Model IIA)

