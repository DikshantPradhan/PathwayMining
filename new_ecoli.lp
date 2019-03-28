\* Problem: FBA_Ecoli_core_model *\

Maximize
 obj: + Biomass_Ecoli_core_w_GAM

Subject To
 r13dpg(c): + PGK + GAPD = 0
 r2pg(c): - PGM - ENO = 0
 r3pg(c): + PGM - PGK - 1.496 Biomass_Ecoli_core_w_GAM = 0
 r6pgc(c): + PGL - GND = 0
 r6pgl(c): - PGL + G6PDH2r = 0
 ac(c): + ACt2r - ACKr = 0
 ac(e): - EX_ac(e) - ACt2r = 0
 acald(c): + ALCD2x + ACALDt - ACALD = 0
 acald(e): - EX_acald(e) - ACALDt = 0
 accoa(c): - PTAr + PFL + PDH - MALS - CS
 - 3.7478 Biomass_Ecoli_core_w_GAM + ACALD = 0
 acon_C(c): - ACONTb + ACONTa = 0
 actp(c): + PTAr + ACKr = 0
 adp(c): + SUCOAS - PYK + PPCK + PGK + PFK + GLNabc + GLNS
 + 59.81 Biomass_Ecoli_core_w_GAM - ATPS4r + ATPM + 2 ADK1 + ACKr = 0
 akg(c): + ICDHyr - GLUSy + GLUDy + 4.1182 Biomass_Ecoli_core_w_GAM
 + AKGt2r - AKGDH = 0
 akg(e): - EX_akg(e) - AKGt2r = 0
 amp(c): + PPS - ADK1 = 0
 atp(c): - SUCOAS + PYK - PPS - PPCK - PGK - PFK - GLNabc - GLNS
 - 59.81 Biomass_Ecoli_core_w_GAM + ATPS4r - ATPM - ADK1 - ACKr = 0
 cit(c): + CS - ACONTa = 0
 co2(c): + PPCK - PPC + PDH + ME2 + ME1 + ICDHyr + GND + CO2t + AKGDH
 = 0
 co2(e): - EX_co2(e) - CO2t = 0
 coa(c): - SUCOAS + PTAr - PFL - PDH + MALS + CS
 + 3.7478 Biomass_Ecoli_core_w_GAM - AKGDH - ACALD = 0
 dhap(c): - TPI + FBA = 0
 e4p(c): - TKT2 + TALA - 0.361 Biomass_Ecoli_core_w_GAM = 0
 etoh(c): + ETOHt2r - ALCD2x = 0
 etoh(e): - EX_etoh(e) - ETOHt2r = 0
 f6p(c): + TKT2 + TALA + PGI - PFK + FRUpts2 + FBP
 - 0.0709 Biomass_Ecoli_core_w_GAM = 0
 fdp(c): + PFK - FBP - FBA = 0
 for(c): + PFL - FORti + FORt2 = 0
 for(e): + FORti - FORt2 - EX_for(e) = 0
 fru(e): - FRUpts2 - EX_fru(e) = 0
 fum(c): + SUCDi + FUMt2_2 - FUM - FRD7 = 0
 fum(e): - FUMt2_2 - EX_fum(e) = 0
 g3p(c): + TPI + TKT2 + TKT1 - TALA - GAPD + FBA
 - 0.129 Biomass_Ecoli_core_w_GAM = 0
 g6p(c): - PGI + GLCpts - G6PDH2r - 0.205 Biomass_Ecoli_core_w_GAM = 0
 glc_D(e): - GLCpts - EX_glc(e) = 0
 gln_L(c): - GLUSy - GLUN + GLNabc + GLNS
 - 0.2557 Biomass_Ecoli_core_w_GAM = 0
 gln_L(e): - GLNabc - EX_gln_L(e) = 0
 glu_L(c): + GLUt2r + 2 GLUSy + GLUN - GLUDy - GLNS
 - 4.9414 Biomass_Ecoli_core_w_GAM = 0
 glu_L(e): - GLUt2r - EX_glu_L(e) = 0
 glx(c): - MALS + ICL = 0
 h2o(c): - PPS - PPC - PGL - MALS + H2Ot - GLUN - GLUDy - GLNabc - FUM
 - FBP + ENO + CYTBD - CS - 59.81 Biomass_Ecoli_core_w_GAM + ATPS4r
 - ATPM - ACONTb + ACONTa = 0
 h2o(e): - H2Ot - EX_h2o(e) = 0
 h(c): + 2 THD2 + SUCCt3 + 2 SUCCt2_2 + PYRt2r - PYK + 2 PPS + PPC
 + PIt2r + PGL + PFK - 4 NADH16 + MDH + 2 MALt2_2 + MALS + LDH_D
 + GLUt2r - GLUSy + GLUDy + GLNabc + GLNS + GAPD + G6PDH2r + 2 FUMt2_2
 + FORt2 + ETOHt2r + D_LACt2 - 2 CYTBD + CS
 + 59.81 Biomass_Ecoli_core_w_GAM + 3 ATPS4r + ATPM + ALCD2x + AKGt2r
 + ACt2r + ACALD = 0
 h(e): - 2 THD2 - SUCCt3 - 2 SUCCt2_2 - PYRt2r - PIt2r + 3 NADH16
 - 2 MALt2_2 - GLUt2r - 2 FUMt2_2 - FORt2 - EX_h(e) - ETOHt2r - D_LACt2
 + 2 CYTBD - 4 ATPS4r - AKGt2r - ACt2r = 0
 icit(c): - ICL - ICDHyr + ACONTb = 0
 lac_D(c): - LDH_D + D_LACt2 = 0
 lac_D(e): - EX_lac_D(e) - D_LACt2 = 0
 mal_L(c): - ME2 - ME1 - MDH + MALt2_2 + MALS + FUM = 0
 mal_L(e): - MALt2_2 - EX_mal_L(e) = 0
 nad(c): + THD2 - PDH - NADTRHD + NADH16 - ME1 - MDH - LDH_D - GAPD
 - 3.547 Biomass_Ecoli_core_w_GAM - ALCD2x - AKGDH - ACALD = 0
 nadh(c): - THD2 + PDH + NADTRHD - NADH16 + ME1 + MDH + LDH_D + GAPD
 + 3.547 Biomass_Ecoli_core_w_GAM + ALCD2x + AKGDH + ACALD = 0
 nadp(c): - THD2 + NADTRHD - ME2 - ICDHyr - GND + GLUSy - GLUDy
 - G6PDH2r + 13.0279 Biomass_Ecoli_core_w_GAM = 0
 nadph(c): + THD2 - NADTRHD + ME2 + ICDHyr + GND - GLUSy + GLUDy
 + G6PDH2r - 13.0279 Biomass_Ecoli_core_w_GAM = 0
 nh4(c): + NH4t + GLUN + GLUDy - GLNS = 0
 nh4(e): - NH4t - EX_nh4(e) = 0
 o2(c): + O2t - 0.5 CYTBD = 0
 o2(e): - O2t - EX_o2(e) = 0
 oaa(c): - PPCK + PPC + MDH - CS - 1.7867 Biomass_Ecoli_core_w_GAM = 0
 pep(c): - PYK + PPS + PPCK - PPC - GLCpts - FRUpts2 + ENO
 - 0.5191 Biomass_Ecoli_core_w_GAM = 0
 pi(c): + SUCOAS - PTAr + PPS + PPC + PIt2r + GLNabc + GLNS - GAPD + FBP
 + 59.81 Biomass_Ecoli_core_w_GAM - ATPS4r + ATPM = 0
 pi(e): - PIt2r - EX_pi(e) = 0
 pyr(c): + PYRt2r + PYK - PPS - PFL - PDH + ME2 + ME1 + LDH_D + GLCpts
 + FRUpts2 - 2.8328 Biomass_Ecoli_core_w_GAM = 0
 pyr(e): - PYRt2r - EX_pyr(e) = 0
 q8(c): - SUCDi - NADH16 + FRD7 + CYTBD = 0
 q8h2(c): + SUCDi + NADH16 - FRD7 - CYTBD = 0
 r5p(c): - TKT1 - RPI - 0.8977 Biomass_Ecoli_core_w_GAM = 0
 ru5p_D(c): + RPI - RPE + GND = 0
 s7p(c): + TKT1 - TALA = 0
 succ(c): - SUCOAS - SUCDi - SUCCt3 + SUCCt2_2 + ICL + FRD7 = 0
 succ(e): + SUCCt3 - SUCCt2_2 - EX_succ(e) = 0
 succoa(c): + SUCOAS + AKGDH = 0
 xu5p_D(c): - TKT2 - TKT1 + RPE = 0

Bounds
 -1000 <= ACALD <= 1000
 -1000 <= ACALDt <= 1000
 -1000 <= ACKr <= 1000
 -1000 <= ACONTa <= 1000
 -1000 <= ACONTb <= 1000
 -1000 <= ACt2r <= 1000
 -1000 <= ADK1 <= 1000
 0 <= AKGDH <= 1000
 -1000 <= AKGt2r <= 1000
 -1000 <= ALCD2x <= 1000
 8.39 <= ATPM <= 1000
 -1000 <= ATPS4r <= 1000
 0 <= Biomass_Ecoli_core_w_GAM <= 1000
 -1000 <= CO2t <= 1000
 0 <= CS <= 1000
 0 <= CYTBD <= 1000
 -1000 <= D_LACt2 <= 1000
 -1000 <= ENO <= 1000
 -1000 <= ETOHt2r <= 1000
 0 <= EX_ac(e) <= 1000
 0 <= EX_acald(e) <= 1000
 0 <= EX_akg(e) <= 1000
 -1000 <= EX_co2(e) <= 1000
 0 <= EX_etoh(e) <= 1000
 0 <= EX_for(e) <= 1000
 0 <= EX_fru(e) <= 1000
 0 <= EX_fum(e) <= 1000
 -10 <= EX_glc(e) <= 1000
 0 <= EX_gln_L(e) <= 1000
 0 <= EX_glu_L(e) <= 1000
 -1000 <= EX_h(e) <= 1000
 -1000 <= EX_h2o(e) <= 1000
 0 <= EX_lac_D(e) <= 1000
 0 <= EX_mal_L(e) <= 1000
 -1000 <= EX_nh4(e) <= 1000
 -1000 <= EX_o2(e) <= 1000
 -1000 <= EX_pi(e) <= 1000
 0 <= EX_pyr(e) <= 1000
 0 <= EX_succ(e) <= 1000
 -1000 <= FBA <= 1000
 0 <= FBP <= 1000
 0 <= FORt2 <= 1000
 0 <= FORti <= 1000
 0 <= FRD7 <= 1000
 0 <= FRUpts2 <= 1000
 -1000 <= FUM <= 1000
 0 <= FUMt2_2 <= 1000
 -1000 <= G6PDH2r <= 1000
 -1000 <= GAPD <= 1000
 0 <= GLCpts <= 1000
 0 <= GLNS <= 1000
 0 <= GLNabc <= 1000
 -1000 <= GLUDy <= 1000
 0 <= GLUN <= 1000
 0 <= GLUSy <= 1000
 -1000 <= GLUt2r <= 1000
 0 <= GND <= 1000
 -1000 <= H2Ot <= 1000
 -1000 <= ICDHyr <= 1000
 0 <= ICL <= 1000
 -1000 <= LDH_D <= 1000
 0 <= MALS <= 1000
 0 <= MALt2_2 <= 1000
 -1000 <= MDH <= 1000
 0 <= ME1 <= 1000
 0 <= ME2 <= 1000
 0 <= NADH16 <= 1000
 0 <= NADTRHD <= 1000
 -1000 <= NH4t <= 1000
 -1000 <= O2t <= 1000
 0 <= PDH <= 1000
 0 <= PFK <= 1000
 0 <= PFL <= 1000
 -1000 <= PGI <= 1000
 -1000 <= PGK <= 1000
 0 <= PGL <= 1000
 -1000 <= PGM <= 1000
 -1000 <= PIt2r <= 1000
 0 <= PPC <= 1000
 0 <= PPCK <= 1000
 0 <= PPS <= 1000
 -1000 <= PTAr <= 1000
 0 <= PYK <= 1000
 -1000 <= PYRt2r <= 1000
 -1000 <= RPE <= 1000
 -1000 <= RPI <= 1000
 0 <= SUCCt2_2 <= 1000
 0 <= SUCCt3 <= 1000
 0 <= SUCDi <= 1000
 -1000 <= SUCOAS <= 1000
 -1000 <= TALA <= 1000
 0 <= THD2 <= 1000
 -1000 <= TKT1 <= 1000
 -1000 <= TKT2 <= 1000
 -1000 <= TPI <= 1000

End
