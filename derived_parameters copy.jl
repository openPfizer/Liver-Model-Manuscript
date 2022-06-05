### A Pluto.jl notebook ###
# v0.19.6

using Markdown
using InteractiveUtils

# ╔═╡ 4a97d273-8195-4493-b13d-2e41de123d5d
begin
	using Markdown
	using InteractiveUtils
	using Printf
	using DataFrames
	using XLSX, CSV
	
	md"""
	#### Load necessary packages
	"""
end

# ╔═╡ 2cad0dea-172e-4779-96b1-259d8393612c
	md"""
	# Derived Rate Constants Workbook
	Here we derive the basal fluxes and rate constants of the model. We rely on literature and steady state constraints to determine the typical value (TV) of the parameters of the model. We then allow these parameters to float to a lower bound value (LV) and higher bound value (HV) in our plausible patient searching.
	"""

# ╔═╡ 8f830c81-e3ce-445f-a338-154d4a3f9038
begin
	df = DataFrame(XLSX.readtable("parameters.xlsx","parameters")...)
	#dict1 = Dict(Pair.(df.Parameter, df.ParamNum))
	md"""
	#### Load parameters.xlsx into a DataFrame
	"""
end

# ╔═╡ 68c18624-20a2-4085-85c0-5efb600ee505
begin
	kg_to_g = 1000.0 # [g/kg]
	L_to_mL = 1000.0 # [mL/L]
	rho_tg = 0.9 # [g/cm^3], https://bionumbers.hms.harvard.edu/files/Densities%20of%20triglycerides%20at%2020%C2%B0C%20&%2040%C2%B0C.pdf
	tg_MW = 772 # molecular weight of an average TG, https://laney.edu/bill_trego/wp-content/uploads/sites/242/2015/01/Calculating-FATG-Molar-Mass.pdf

	glycerol_MW = 92.1

	FA_TG = 3 # [mol FA/mol TG], conversion
	fa_MW = (tg_MW - glycerol_MW)/FA_TG
	
	rhoe_tg = 9*tg_MW/1000.0 # kcal/g --> kcal/mmol
	rhoe_fa = 9.6*fa_MW/1000.0 # [kcal/mmol], includes unit conversions

	vol_frac_cytosol = 0.7 # [fraction], 				https://www.proteinatlas.org/humanproteome/cell/cytosol
	vol_ratio_er_to_cyt = 0.16/0.5 # [ratio], http://book.bionumbers.org/how-big-is-the-endoplasmic-reticulum-of-cells/
	vol_hepat = 3.4e-9 # [cm^3/cell], https://en.wikipedia.org/wiki/Hepatocyte
	num_hepat = 139*10^6 # [cells/g liver], https://www.sciencedirect.com/science/article/pii/S088723330600124X?via%3Dihub
	mass_liver = 1500 # [g], https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/liver-weight

	mass_liver_low = (1677 - 2*396)/1677 # [multiplier], Grandmaison et al. FSI. 2001. 2SD
	mass_liver_high = (1677 + 2*396)/1677 # [multiplier], Grandmaison et al. FSI. 2001. 2SD
	
	vol_hepat_tot = vol_hepat*num_hepat*mass_liver # [mL]

	diet_kcals = 2400 # [kcal/day], middle-aged male at normal BMI, moderate exercise (https://www.calculator.net/calorie-calculator.html?ctype=standard&cage=45&csex=m&cheightfeet=5&cheightinch=10&cpound=165&cheightmeter=180&ckg=65&cactivity=1.465&cmop=0&coutunit=c&cformula=m&cfatpct=20&printit=0&x=22&y=15)

	diet_kcals_low = 1596 # [kcal/day], 60kg, 80 year old male with no regular exercise, low-fat diet
	diet_kcals_high = 3230 # [kcal/day], 25 year old, 120 kg male, regular exercise, high-fat diet
	
	diet_frac_fat = 0.35 # [fraction]
	lcfa_frac = 0.85 # [fraction]
	ba_fat = 0.95 # [fraction], bioavailability of fat, https://www.sciencedirect.com/topics/medicine-and-dentistry/fat-intake

	bmr_liver = 300 # [kcal/kg/day], https://www.fao.org/3/m2845e/m2845e00.htm
	frac_liver_fatox = 0.6 # [fraction], portion of liver TEE that is from fat oxidation
	fa_basal = 0.15 # [mM], Holzhütter H-G, Berndt N. Computational Hypothesis: How Intra-Hepatic Functional Heterogeneity May Influence the Cascading Progression of Free Fatty Acid-Induced Non-Alcoholic Fatty Liver Disease (NAFLD). Cells. 2021; 10(3):578. https://doi.org/10.3390/cells10030578

	vol_frac_tg = 0.05 # [fraction]
	rho_liver = 1.07 # [g/mL], https://eurjmedres.biomedcentral.com/articles/10.1186/2047-783X-15-8-345#:~:text=Median%20liver%20density%20was%201.07%20g%2Fml.

	tg_plasma_basal = 1.5385 # [mM], approximately a mean daily concentration
	vd_tg_plasma = 4.515 # [L], https://doi.org/10.1007/s11745-001-0696-6
	vd_tg_plasma_low = 4.515 - 2*0.308*sqrt(10) # [L], mean - 2*SD
	vd_tg_plasma_high = 4.515 + 2*0.308*sqrt(10)# [L], mean + 2*SD

	emax_vldl_prod = 33.62629534 # [mmols/day], Adiels et al. 2005.
	ec50_lit = 2.2879 # [%-volume], Adiels et al. 2005.

	frac_fa_dnl_basal = 0.05 # [fraction], Lambert et al. 2014.
	frac_fa_diet_basal = 0.05 # [fraction], Lambert et al. 2014.

	dnl_low = 0.1 #(3.3 - 2*sqrt(6)*0.8)/3.3 # [multiplier], Schwarz et al. JCI. 1995. 2SD of eucaloric feed value overlaps zero, assume the lower-bound is ~0. Lower bound 0 is not useful for parameter searching, so set at low, but non-zero.
	dnl_high = 8 # [multiplier], Schwarz et al. JCI. 1995. 6 - 10x increase observed from basal with high CHO intake.

	# Useful values for setting LV/HV of parameters
	tg_plasma_low = 30/tg_MW*10.0 # [mM], low plasma TG value, mg/dL --> mM
	tg_plasma_high = 750/tg_MW*10.0 # [mM], high plasma TG value, mg/dL --> mM

	lf_prct_low = 0.3 # [%], low liver fat as vol. %. Functionally, there are many people without observable liver fat, set very low, but non-zero
	lf_prct_high = 67.5 # [%], high liver fast as vol. %

	Q_cardiac = 5*1440.0 # [L/day], cardiac output of adult converted to per day
	nefa_high = 0.900 # [mM], https://doi.org/10.1161/01.ATV.20.6.1588, Table 1, Diabetic group.
	
	md"""
	### Useful constants for calculations
	Various physical and biological constants used in the calculations below. Some could be treated as variables for inter-patient variability (e.g., diet kcals, diet frac fat), if desired.
	"""
end

# ╔═╡ 0bcd76cf-e186-4c73-b701-801c933c14b9
begin
	vd_cyt = vol_frac_cytosol * vol_hepat * num_hepat * mass_liver * 1/L_to_mL # [L]
	vd_er = vd_cyt*vol_ratio_er_to_cyt

	md"""
	### Cytosolic and ER Volumes (L)
	Determine the volume of the cytosol of a hepatocyte. This value can be floated for plausible patients, but reasonably not over orders-of-magnitude.
	"""
end

# ╔═╡ 0feb0564-0a63-48ef-a748-3d3fbfe76d7d
begin
	liver_fat_ox = bmr_liver*frac_liver_fatox*mass_liver/kg_to_g # [kcal/day]
	liver_fat_ox_mmol = liver_fat_ox/rhoe_fa # [mmols-FA/day]
	kbetaox = liver_fat_ox_mmol/(fa_basal*vd_cyt) # [1/day], FA_basal is assumed from a literature source (see reference in Description).

	md"""
	### Rate Constant of Beta Oxidation
	"""
end

# ╔═╡ 7c06328e-6039-4077-b25e-d2285ba3c4ff
begin
	chylo_basal_flux = diet_kcals*diet_frac_fat*lcfa_frac*ba_fat/rhoe_tg # [mmol-TG/day]

	md"""
	### Daily Chylomicron Flux
	Chylomicron flux is essentially a proxy for fat intake over 24 hours, adjusted slightly for bioavailability.
	"""
end

# ╔═╡ 1eec7d02-2720-44b7-9f15-fde740c0443e
begin
	tg_cy_basal = vol_hepat_tot*vol_frac_tg*rho_tg/tg_MW*1e6/vol_hepat_tot # [mM]

	md"""
	### Basal Hepatocyte Triglyceride Concentrations
	Concentration of triglyceride in the cytosol and ER. Straight calculation of amount of TG, density of TG, and volume of cytosol (calculated previously). ER is assumed similar to cytosol.
	"""
end

# ╔═╡ 407c2009-d18e-48b1-89d6-04ff4cafb1bc
begin
	
	ec50_vldl_voltg_mL = ec50_lit/100*vol_hepat_tot # [mL]
	ec50_vldl_mM = ec50_vldl_voltg_mL*rho_tg/tg_MW/vol_hepat_tot*1e6 # [mM], includes unit conversions for mol/mL --> mmol/L
	
	md"""
	### EC50 VLDL Production -- Unit Conversion
	"""
end

# ╔═╡ ea93a794-3847-4828-8b47-bd72e9d3cb44
begin
	vldl_flux_basal = emax_vldl_prod*tg_cy_basal/(ec50_vldl_mM+tg_cy_basal) # [mmols-TG/day], we assumed ER ~ Cytosolic concentration at basal
	
	tg_input_flux = vldl_flux_basal + chylo_basal_flux # [mmols-TG]/day	

	# Total output of FA from the liver/day:
	fa_liver_output_flux = vldl_flux_basal*FA_TG + liver_fat_ox_mmol # [mmol-FA/day]
	fa_liver_input_flux = fa_liver_output_flux # [mmol-FA/day], mass balance at steady state (input == output)

	# Break down/out the various FA input fluxes (i.e., DNL, plasma uptake, NEFAs)
	frac_fa_other = 1 - (frac_fa_dnl_basal + frac_fa_diet_basal) # [fraction]
	dnl_flux_basal = fa_liver_input_flux*frac_fa_dnl_basal # [mmol/day]
	
	# Dietary contribution is trickier because VLDL is originating from multiple sources. _Assume_ that chylos are diet and everything else is "fatty acid" derived (adipose storage)
	frac_tg_plasma_diet = chylo_basal_flux/tg_input_flux # [fraction], approximation
	kuptake_liver_tg = frac_fa_diet_basal*fa_liver_input_flux/(FA_TG*
		tg_plasma_basal*vd_tg_plasma*frac_tg_plasma_diet) # [1/day]

	tot_tg_liver_uptake = kuptake_liver_tg*tg_plasma_basal*vd_tg_plasma # [mmols-TG/day]

	nefa_uptake_flux = fa_liver_input_flux - dnl_flux_basal - FA_TG*tot_tg_liver_uptake # [mmols-FA/day]

	md"""
	### TG Mass Balance
	We work backwards through several of the fluxes in the model to use the steady state constraints to determine the needed flux values and rate constants for the system to be at steady state.
	"""	
end

# ╔═╡ a7080cb5-6afc-4572-9c8f-adf0c70ef071
begin
	
	tg_lipase_clear = tg_input_flux - tot_tg_liver_uptake # [mmols/day]
	klipase_clear = tg_lipase_clear/(tg_plasma_basal*vd_tg_plasma) # [1/day], clearance by lipases in plasma

	md"""
	### TG Plasma Clearance Rate
	Calculate the clearance rate of triglyceride from the plasma via lipoprotein lipase and endothelial lipase
	"""
end

# ╔═╡ 7c1c980e-b248-4932-86b8-f2fb878a7eb5
begin
	fa_er_basal = fa_basal # [mM], assumption
	tg_synth_er_basal = vldl_flux_basal/vd_er # [mM-TG/day]
	ksynth_er_tg = vldl_flux_basal*FA_TG/(fa_er_basal^3*vd_er) # [1/(mM^2*days)]
	kuptake_er = vldl_flux_basal*FA_TG/(fa_basal*vd_cyt) # [1/days]

	#er_uptake_flux = kuptake_er*fa_basal*vd_cyt # [mmols-FA-cyt/day]
	#vd_er

	md"""
	### Remaining Mass Balances
	We back-solved for the majority of the major fluxes at baseline. Now back-solve for the remaining rate constants needed to support the flux.
	"""
end

# ╔═╡ bb78de7c-a91f-42f5-8658-6e0669b108ec
begin
	vldl_frac_ltg_cyt = 0.35 # [fraction], based off Fabrinni 2008 to approximate a residence time of 3 - 4 days
	lipolysis_flux_mmols = vldl_flux_basal*vldl_frac_ltg_cyt # [mmols-TG/day], estimation of a residence time
	klipo_cy_tg = lipolysis_flux_mmols/(tg_cy_basal*vd_cyt) # [1/day]
	ksynth_cy_tg = lipolysis_flux_mmols*FA_TG/(fa_basal^3*vd_cyt) # [1/(mM^2*day]
	ester_cyt_flux = ksynth_cy_tg*fa_basal^3*vd_cyt

	md"""
	### Lipid Droplet Dynamics
	In the basal or baseline state, the lipid droplet should be approximately in equilibrium, neither growing nor shrinking dramatically. Thus, estimating one arm (i.e., lipolysis or esterification) will determine the second arm. We base our estimate on estimates of liver TG being the primary contributor to VLDL, which 	indicates a need for a high cycling rate.
	"""
end

# ╔═╡ f87951a8-b411-498c-ac35-cb3ccb395331
begin
	
	df[df.Parameter .== "klipase_clear", :TV] .= klipase_clear # [1/day]
	df[df.Parameter .== "dnl_basal_flux", :TV] .= dnl_flux_basal # [mmols/day]
	df[df.Parameter .== "kuptake_liver_tg", :TV] .= kuptake_liver_tg # [mmols/day]
	df[df.Parameter .== "nefa_uptake_flux", :TV] .= nefa_uptake_flux # [mmols/day]
	df[df.Parameter .== "ec50_vldl_prod", :TV] 	.= ec50_vldl_mM # [mM-TG]
	df[df.Parameter .== "tg_cy_basal", :TV] .= tg_cy_basal # [mM]
	df[df.Parameter .== "tg_er_basal", :TV] .= tg_cy_basal # [mM], assumption
	df[df.Parameter .== "chylo_basal_flux", :TV] .= chylo_basal_flux # [mmols/day]
	df[df.Parameter .== "kbetaox", :TV] .= kbetaox # [1/day]
	df[df.Parameter .== "vd_cyt", :TV] .= vd_cyt # [L]
	df[df.Parameter .== "vd_er", :TV] .= vd_er # [L]
	df[df.Parameter .== "fa_cy_basal", :TV] .= fa_basal # [mM]
	df[df.Parameter .== "fa_er_basal", :TV] .= fa_basal # [mM]
	df[df.Parameter .== "tg_p_basal", :TV] 	.= tg_plasma_basal # [mM]
	df[df.Parameter .== "vd_tg_p", :TV] 	.= vd_tg_plasma # [L]
	df[df.Parameter .== "emax_vldl_prod", :TV] 	.= emax_vldl_prod # [mmol-TG/day]
	df[df.Parameter .== "kuptake_er", :TV] 	.= kuptake_er # [1/day]
	df[df.Parameter .== "ksynth_er_tg", :TV] 	.= ksynth_er_tg # [1/(mM^2*day)]
	df[df.Parameter .== "ksynth_cy_tg", :TV] 	.= ksynth_cy_tg # [1/(mM^2*day)]
	df[df.Parameter .== "klipo_cy_tg", :TV] 	.= klipo_cy_tg # [1/day]
	
	# Set the bounds for the plausible patient searches
	# LV_M, HV_M (low value multiplier, high value multiplier) should be pre-set (i.e., before the plausible patient search) based on biological constraints.
	df[df.Parameter .== "klipase_clear", :LVM] .= tg_plasma_basal/tg_plasma_high
	df[df.Parameter .== "klipase_clear", :HVM] .= tg_plasma_basal/tg_plasma_low
	df[df.Parameter .== "chylo_basal_flux", :LVM] .= diet_kcals_low/diet_kcals*0.25/0.35
	df[df.Parameter .== "chylo_basal_flux", :HVM] .= diet_kcals_high/diet_kcals*0.45/0.35
	df[df.Parameter .== "dnl_basal_flux", :LVM] .= dnl_low
	df[df.Parameter .== "dnl_basal_flux", :HVM] .= dnl_high
	df[df.Parameter .== "nefa_uptake_flux", :LVM] .= 0.5 # If NEFAs are allowed to go too low, they become a non-contributor to hepatic FAs, which is in contrast to Lambert et al. Narrow the lower window a bit. 
	#df[df.Parameter .== "nefa_uptake_flux", :HVM] .= 4.0 #Alternative (much wider) Q_cardiac*nefa_high/nefa_uptake_flux # Based on limitations only by cardiac output
	df[df.Parameter .== "vd_tg_p", :LVM] .= vd_tg_plasma_low/vd_tg_plasma
	df[df.Parameter .== "vd_tg_p", :HVM] .= vd_tg_plasma_high/vd_tg_plasma
	df[df.Parameter .== "vd_cyt", :LVM] .= mass_liver_low
	df[df.Parameter .== "vd_cyt", :HVM] .= mass_liver_high
	df[df.Parameter .== "vd_er", :LVM] .= mass_liver_low
	df[df.Parameter .== "vd_er", :HVM] .= mass_liver_high
	
	
	
	lf_frac_high = lf_prct_high/100
	lf_frac_low = lf_prct_low/100
	
	#df[df.Parameter .== "fa_cy_basal", :LVM]  No lower bound, keep at 0.25 multiplier
	df[df.Parameter .== "fa_cy_basal", :HVM] .= (1-vol_frac_tg)/(1-lf_frac_high)*lf_frac_high/vol_frac_tg
	df[df.Parameter .== "tg_cy_basal", :LVM] .= (1-vol_frac_tg)/(1-lf_frac_low)*lf_frac_low/vol_frac_tg
	df[df.Parameter .== "tg_cy_basal", :HVM] .= (1-vol_frac_tg)/(1-lf_frac_high)*lf_frac_high/vol_frac_tg
	#df[df.Parameter .== "fa_er_basal", :LVM]  No lower bound, keep at 0.25 multiplier
	df[df.Parameter .== "fa_er_basal", :HVM] .= (1-vol_frac_tg)/(1-lf_frac_high)*lf_frac_high/vol_frac_tg
	df[df.Parameter .== "tg_er_basal", :LVM] .= (1-vol_frac_tg)/(1-lf_frac_low)*lf_frac_low/vol_frac_tg
	df[df.Parameter .== "tg_er_basal", :HVM] .= (1-vol_frac_tg)/(1-lf_frac_high)*lf_frac_high/vol_frac_tg
	df[df.Parameter .== "tg_p_basal", :LVM] .= tg_plasma_low/tg_plasma_basal
	df[df.Parameter .== "tg_p_basal", :HVM] .= tg_plasma_high/tg_plasma_basal
		
	df.LV = df.TV .* df.LVM
	df.HV = df.TV .* df.HVM

	# Write the values back to the original Excel sheet: "parameters_pluto.csv".
	CSV.write("parameters_pluto.csv",df)
	
	md"""
	### Rewrite the DataFrame Values
	Set the low and high values for the plausible patient search based on the TV values determined in this notebook. Save back the values to parameters.xlsx for easier access in the main program scripts.
	"""

end


