within Kidney;

model NephronSF
  parameter Real L_CortexOuterMedula = 0.04 "cortex and outer medula thicknes";
  parameter Real L_ProximalTubulus = 0.03 "proximal tubulus length";
  DomainLineSegment1D cortexOuterMedula(L = L_CortexOuterMedula, N = 100);
  field Real qD_H2O(domain = cortexOuterMedula) "wather flux in descending part of LOH";
  field Real qA_H2O(domain = cortexOuterMedula) "wather flux in  ascending part of LOH";
  parameter Real reabsFrac = 5.0/6.0 "rebsorption fraction w.r.t. q_JG_in";
  field Real qD_out_H2O(domain = cortexOuterMedula) "water reabsbsorption per unit length";
  parameter Real q_total_in = 0.180 /24/360 "total volume of filtrate produced per second (180L/day)";
  parameter Real q_JG_in = q_total_in*0.85 "flux of filtrate in juxtaglomerular nephrons, is 85% of total";
equation
//descending
 qD_out_H2O = if (cortexOuterMedula.x < L_ProximalTubulus) then q_JG_in*reabsFrac/L_ProximalTubulus else 0 indomain cortexOuterMedula "water reabsorption";
 qD_H2O = q_JG_in indomain cortexOuterMedula.left "input boundary condition";
 pder(qD_H2O,x) = - qD_out_H2O indomain cortexOuterMedula "water flow decrease due to reabsorption";
 qD_H2O = extrapolateField(qD_H2O) indomain cortexOuterMedula.right "output bc extrapolation";
//ascending
 qD_H2O = - qA_H2O indomain cortexOuterMedula.right "bc between descending and ascending part of LOH";
 pder(qA_H2O,x) = 0 indomain cortexOuterMedula;
 qA_H2O = extrapolateField(qA_H2O) indomain cortexOuterMedula.left;
end NephronSF;