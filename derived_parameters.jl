### A Pluto.jl notebook ###
# v0.19.4

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
	dict1 = Dict(Pair.(df.Parameter, df.ParamNum))
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

	lf_prct_low = 0.1 # [%], low liver fat as vol. %. Functionally, there are many people without observable liver fat, set very low, but non-zero
	lf_prct_high = 90 # [%], high liver fast as vol. %

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
	#df[df.Parameter .== "nefa_uptake_flux", :LVM] .= No lower bound, keep at 0.25 multiplier
	df[df.Parameter .== "nefa_uptake_flux", :HVM] .= Q_cardiac*nefa_high/nefa_uptake_flux # Based on limitations only by cardiac output
	df[df.Parameter .== "vd_tg_p", :LVM] .= vd_tg_plasma_low/vd_tg_plasma
	df[df.Parameter .== "vd_tg_p", :HVM] .= vd_tg_plasma_high/vd_tg_plasma
	df[df.Parameter .== "vd_cyt", :LVM] .= mass_liver_low
	df[df.Parameter .== "vd_cyt", :HVM] .= mass_liver_high
	df[df.Parameter .== "vd_er", :LVM] .= mass_liver_low
	df[df.Parameter .== "vd_er", :HVM] .= mass_liver_high
	#df[df.Parameter .== "fa_cy_basal", :LVM]  No lower bound, keep at 0.25 multiplier
	df[df.Parameter .== "fa_cy_basal", :HVM] .= lf_prct_high/(vol_frac_tg*100)
	df[df.Parameter .== "tg_cy_basal", :LVM] .= lf_prct_low/(vol_frac_tg*100)
	df[df.Parameter .== "tg_cy_basal", :HVM] .= lf_prct_high/(vol_frac_tg*100)
	#df[df.Parameter .== "fa_er_basal", :LVM]  No lower bound, keep at 0.25 multiplier
	df[df.Parameter .== "fa_er_basal", :HVM] .= lf_prct_high/(vol_frac_tg*100)
	df[df.Parameter .== "tg_er_basal", :LVM] .= lf_prct_low/(vol_frac_tg*100)
	df[df.Parameter .== "tg_er_basal", :HVM] .= lf_prct_high/(vol_frac_tg*100)
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
InteractiveUtils = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
Markdown = "d6f4376e-aef5-505a-96c1-9c027394607a"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
XLSX = "fdbf4ff8-1666-58a4-91e7-1b58723a45e0"

[compat]
CSV = "~0.10.4"
DataFrames = "~1.3.4"
XLSX = "~0.7.10"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "b153278a25dd42c65abbf4e62344f9d22e59191b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.43.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "daa21eb85147f72e41f6352a57fccea377e310a9"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.4"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "cc1a8e22627f33c789ab60b36a9132ac050bbf75"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.12"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EzXML]]
deps = ["Printf", "XML2_jll"]
git-tree-sha1 = "0fa3b52a04a4e210aeb1626def9c90df3ae65268"
uuid = "8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615"
version = "1.1.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "129b104185df66e408edd6625d480b7f9e9823a0"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.18"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "61feba885fac3a407465726d0c330b3055df897f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "1285416549ccfcdf0c50d4997a94331e88d68413"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "6a2f7d70512d205ca8c7ee31bfa9f142fe74310c"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.12"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.XLSX]]
deps = ["Dates", "EzXML", "Printf", "Tables", "ZipFile"]
git-tree-sha1 = "7fa8618da5c27fdab2ceebdff1da8918c8cd8b5d"
uuid = "fdbf4ff8-1666-58a4-91e7-1b58723a45e0"
version = "0.7.10"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.ZipFile]]
deps = ["Libdl", "Printf", "Zlib_jll"]
git-tree-sha1 = "3593e69e469d2111389a9bd06bac1f3d730ac6de"
uuid = "a5390f91-8eb1-5f08-bee0-b1d1ffed6cea"
version = "0.9.4"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─2cad0dea-172e-4779-96b1-259d8393612c
# ╠═4a97d273-8195-4493-b13d-2e41de123d5d
# ╠═8f830c81-e3ce-445f-a338-154d4a3f9038
# ╠═68c18624-20a2-4085-85c0-5efb600ee505
# ╠═0bcd76cf-e186-4c73-b701-801c933c14b9
# ╠═0feb0564-0a63-48ef-a748-3d3fbfe76d7d
# ╠═7c06328e-6039-4077-b25e-d2285ba3c4ff
# ╠═1eec7d02-2720-44b7-9f15-fde740c0443e
# ╠═407c2009-d18e-48b1-89d6-04ff4cafb1bc
# ╠═ea93a794-3847-4828-8b47-bd72e9d3cb44
# ╠═a7080cb5-6afc-4572-9c8f-adf0c70ef071
# ╠═7c1c980e-b248-4932-86b8-f2fb878a7eb5
# ╠═bb78de7c-a91f-42f5-8658-6e0669b108ec
# ╠═f87951a8-b411-498c-ac35-cb3ccb395331
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
