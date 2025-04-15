---
status: started
---
Many of these were easy to do as there was a direct tie to the TS1 manual conversion of user reactions.

Some are a little more complex.

To-do: 
- [x] pan_b
- [x] mpan_b
- [ ] so2h
- [x] oh_dms
- [x] nh3h
- [x] strat38
- [x] strat 72-80


Notes: 

~~The pans just aren't finished.~~

We checked with Larry over the issue with MPAN and there being 2 forward rates. Change to using the forward from \*.inp and ignore the forward in mo_usrrxt.F90.

To just feed in a rate constant, use type "USER_DEFINED", and give the reactants/products. 

--- 
nh3h 
has a gamma of 0, so ? Going to put it in and see what happens

---
oh_dms 
possibly using a constant for one of the m(:,k). Becky is checking.

2/11: it wouldn't be that. I worked out how they got there in the TS1 version and then worked it out for this implementation. Still missing a m(:,k).

Some time later: 

---
so2h
Structure of the user rxn is that there are several options of how to calculate gamma. Horowitz et al. 2020 and the namelist for UFS-AM4 configuration indicate that the Zheng et al. 2015 calculation should be used. 

- In gfdl_am3chem_mod.F90, within the gfdl_am4chem_init subroutine/lines 1433-1441 appear to fill trop_option%gSO2_dynamic to -1, 1, 2 depending on what text was given in the namelist. 1 = wang2014, 2= zheng2015.
---
strat38 
~~Did not finish reviewing this one.~~
Complete. Followed logic/formatting from PBZNIT_M and HO2NO2_M
- followup from 2/11: need to recheck because I messed up some of these. 
---
strat 72->80

These rxns follow a similar format in the user rxns file in that they all follow from strat_chem_get_gamma and strat_chem_get_hetrates. 

Within each of these routines there is a lot that happens. 

In several points, there is a dependence of the calculation on values contained in "psc". The fill of psc looks like it's in the subroutine strat_chem_get_psc

Gamma calculation: 
indexes 2,6,7 set to 0 (strat73, strat77, strat78)
index 8 set to 0.25 (strat79)
index 9 (strat80) is given as $ 1.0/(1.2422 +1.0/(0.114 + EXP(29.24-0.396x psc%wh2so4))) $
indexes 1, 3, 4, 5 ae calculated from a bunch of things starting with temperature. 

hetrates calculation: 
Gammas are part of the input. Couple of stages of case selection. 
Ultimately leads to grouping of treatment based on if reacting with HCl, H2O, or HBr. Within each of these, there is a calculation based on NAT, ICE, and liquid aerosol. 

