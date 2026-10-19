/*QFracture-2016 10-year major osteoporotic fracture calculation*/

*Identify records with all inputs required for their sex 
generate byte qfracture_calculable = ///
    inlist(lower(strtrim(sex)), "female", "male") & ///
    inrange(age, 30, 99) & qf_bmi > 0 & !missing(qf_bmi) & ///
    inrange(qf_alcohol_cat6, 0, 5) & ///
    inrange(qf_ethrisk, 0, 9) & ///
    inrange(qf_smoke_cat, 0, 4)

*Construct the sex-specific linear predictor 
tempvar linear_predictor dage age_1 age_2 dbmi bmi_1 bmi_2
generate double `linear_predictor' = 0 if qfracture_calculable
generate double `dage' = age / 10 if qfracture_calculable
generate double `dbmi' = qf_bmi / 10 if qfracture_calculable

/* Female QFracture_2016_fracture4 equation. */
generate double `age_1' = (`dage'^2) - 25.463895797729492 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female"
generate double `age_2' = (`dage'^3) - 128.495315551757810 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female"
generate double `bmi_1' = (`dbmi'^(-1)) - 0.382189363241196 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female"

replace `linear_predictor' = `linear_predictor' ///
    + (`age_1' * 0.148823061721650830) ///
    + (`age_2' * -0.0095516624764288762) ///
    + (`bmi_1' * 2.8180291389827810) ///
    if qfracture_calculable & lower(strtrim(sex)) == "female"

replace `linear_predictor' = `linear_predictor' -0.016119659842711558 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female" & ///
    qf_alcohol_cat6 == 1
replace `linear_predictor' = `linear_predictor' +0.018142191954688299 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female" & ///
    qf_alcohol_cat6 == 2
replace `linear_predictor' = `linear_predictor' +0.087039813091311105 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female" & ///
    qf_alcohol_cat6 == 3
replace `linear_predictor' = `linear_predictor' +0.485087668164837120 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female" & ///
    qf_alcohol_cat6 == 4
replace `linear_predictor' = `linear_predictor' +0.452147004572386320 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female" & ///
    qf_alcohol_cat6 == 5

replace `linear_predictor' = `linear_predictor' -0.425660692163662540 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female" & qf_ethrisk == 2
replace `linear_predictor' = `linear_predictor' -0.554320911950214160 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female" & qf_ethrisk == 3
replace `linear_predictor' = `linear_predictor' -0.918260109780693060 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female" & qf_ethrisk == 4
replace `linear_predictor' = `linear_predictor' -0.681936065314830420 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female" & qf_ethrisk == 5
replace `linear_predictor' = `linear_predictor' -1.466848340498807700 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female" & qf_ethrisk == 6
replace `linear_predictor' = `linear_predictor' -0.910123811422844600 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female" & qf_ethrisk == 7
replace `linear_predictor' = `linear_predictor' -0.642178331754473920 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female" & qf_ethrisk == 8
replace `linear_predictor' = `linear_predictor' -0.503682943263451090 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female" & qf_ethrisk == 9

replace `linear_predictor' = `linear_predictor' +0.055735693430561166 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female" & qf_smoke_cat == 1
replace `linear_predictor' = `linear_predictor' +0.163389566170135280 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female" & qf_smoke_cat == 2
replace `linear_predictor' = `linear_predictor' +0.154048833869658700 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female" & qf_smoke_cat == 3
replace `linear_predictor' = `linear_predictor' +0.232977159175790450 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female" & qf_smoke_cat == 4

replace `linear_predictor' = `linear_predictor' ///
    + b_antidepressant  *  0.293589157888128390 ///
    + b_anycancer       *  0.117552273314779310 ///
    + b_asthmacopd      *  0.199719375335289990 ///
    + b_corticosteroids *  0.218702009424677490 ///
    + b_cvd             *  0.141964372550305170 ///
    + b_dementia        *  0.469738726308515210 ///
    + b_endocrine       *  0.110539402821759230 ///
    + b_epilepsy2       *  0.402446009860403300 ///
    + b_falls           *  0.332232162630358370 ///
    + b_hrt_oest        * -0.202786445635823240 ///
    + b_liver           *  0.483116157696158620 ///
    + b_malabsorption   *  0.168747780183557460 ///
    + b_parkinsons      *  0.474223935803918140 ///
    + b_ra_sle          *  0.226705932747190420 ///
    + b_renal           *  0.250864872379400640 ///
    + b_type1           *  0.783288716093229360 ///
    + b_type2           *  0.236386965781406080 ///
    + fh_osteoporosis   *  0.383794975586049470 ///
    if qfracture_calculable & lower(strtrim(sex)) == "female"


/* Male QFracture_2016_fracture4 equation. */
replace `age_1' = (`dage'^0.5) - 2.201908826828003 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male"
replace `age_2' = `dage' - 4.848402023315430 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male"
replace `bmi_1' = (`dbmi'^(-1)) - 0.375702142715454 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male"
generate double `bmi_2' = (`dbmi'^(-0.5)) - 0.612945437431335 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male"

replace `linear_predictor' = `linear_predictor' ///
    + (`age_1' * -8.596211339574217900) ///
    + (`age_2' *  2.353458556759995700) ///
    + (`bmi_1' * 18.353475363901577000) ///
    + (`bmi_2' * -19.090527346773332000) ///
    if qfracture_calculable & lower(strtrim(sex)) == "male"

replace `linear_predictor' = `linear_predictor' -0.112964871341206640 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male" & ///
    qf_alcohol_cat6 == 1
replace `linear_predictor' = `linear_predictor' -0.101993314738360990 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male" & ///
    qf_alcohol_cat6 == 2
replace `linear_predictor' = `linear_predictor' -0.0087172177838055233 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male" & ///
    qf_alcohol_cat6 == 3
replace `linear_predictor' = `linear_predictor' +0.237836729189438680 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male" & ///
    qf_alcohol_cat6 == 4
replace `linear_predictor' = `linear_predictor' +0.590067482832426940 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male" & ///
    qf_alcohol_cat6 == 5

replace `linear_predictor' = `linear_predictor' -0.249775936169492060 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male" & qf_ethrisk == 2
replace `linear_predictor' = `linear_predictor' -0.338872340936377190 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male" & qf_ethrisk == 3
replace `linear_predictor' = `linear_predictor' -1.178065231275621700 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male" & qf_ethrisk == 4
replace `linear_predictor' = `linear_predictor' -0.563753612841744520 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male" & qf_ethrisk == 5
replace `linear_predictor' = `linear_predictor' -0.880267422140305620 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male" & qf_ethrisk == 6
replace `linear_predictor' = `linear_predictor' -0.622083906805602100 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male" & qf_ethrisk == 7
replace `linear_predictor' = `linear_predictor' -0.870024076008438780 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male" & qf_ethrisk == 8
replace `linear_predictor' = `linear_predictor' -0.358693460399084340 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male" & qf_ethrisk == 9

replace `linear_predictor' = `linear_predictor' +0.037462254711292393 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male" & qf_smoke_cat == 1
replace `linear_predictor' = `linear_predictor' +0.253665462428324720 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male" & qf_smoke_cat == 2
replace `linear_predictor' = `linear_predictor' +0.272854168990491700 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male" & qf_smoke_cat == 3
replace `linear_predictor' = `linear_predictor' +0.349825978043770030 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male" & qf_smoke_cat == 4

replace `linear_predictor' = `linear_predictor' ///
    + b_antidepressant  * 0.418500980256445660 ///
    + b_anycancer       * 0.292384966784961790 ///
    + b_asthmacopd      * 0.265824315015200320 ///
    + b_carehome        * 0.094672730961445892 ///
    + b_corticosteroids * 0.260159289395475500 ///
    + b_cvd             * 0.215542097682090070 ///
    + b_dementia        * 0.625306639047599420 ///
    + b_epilepsy2       * 0.669570806430142660 ///
    + b_falls           * 0.447714134262650220 ///
    + b_liver           * 0.797093608102493630 ///
    + b_malabsorption   * 0.250892103549892440 ///
    + b_parkinsons      * 0.825427634021887350 ///
    + b_ra_sle          * 0.384658583273509320 ///
    + b_renal           * 0.447967288006248630 ///
    + b_type1           * 1.012168914506730500 ///
    + b_type2           * 0.243925443249027930 ///
    + fh_osteoporosis   * 0.989871751246659380 ///
    if qfracture_calculable & lower(strtrim(sex)) == "male"


/*Convert the linear predictor to 10-year absolute risk*/
generate double qfracture_mof_10y_pct = .

replace qfracture_mof_10y_pct = ///
    100 * (1 - (0.978452086448669 ^ exp(`linear_predictor'))) ///
    if qfracture_calculable & lower(strtrim(sex)) == "female"

replace qfracture_mof_10y_pct = ///
    100 * (1 - (0.991695225238800 ^ exp(`linear_predictor'))) ///
    if qfracture_calculable & lower(strtrim(sex)) == "male"

label variable qfracture_mof_10y_pct ///
    "QFracture 2016 10-year major osteoporotic fracture risk (%)"

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	