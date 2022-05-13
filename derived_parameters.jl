### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ 4a97d273-8195-4493-b13d-2e41de123d5d
begin
	using Markdown
	using InteractiveUtils
	using Printf

	md"""
	# Derived Rate Constants
	"""
end

# ╔═╡ 68c18624-20a2-4085-85c0-5efb600ee505
begin
	kg_to_g = 1000.0
	mL_to_L = 1000.0
	rho_tg = 0.9 # [g/cm^3], https://bionumbers.hms.harvard.edu/files/Densities%20of%20triglycerides%20at%2020%C2%B0C%20&%2040%C2%B0C.pdf
	tg_MW = 772 # molecular weight of an average TG, https://laney.edu/bill_trego/wp-content/uploads/sites/242/2015/01/Calculating-FATG-Molar-Mass.pdf

	glycerol_MW = 92.1
	
	fa_MW = (tg_MW - glycerol_MW)/3
	
	rhoe_tg = 9*tg_MW/1000.0 # kcal/g --> kcal/mmol
	rhoe_fa = 9.6*fa_MW/1000.0 # [kcal/mmol], includes unit conversions
	
	md"""
	### Useful Constants
	"""
end

# ╔═╡ fb85dfb0-8e7d-4a1f-89b8-1ffab7549e89
rhoe_fa # [kcal/mmol-fat]

# ╔═╡ 0bcd76cf-e186-4c73-b701-801c933c14b9
begin
	vol_frac_cytosol = 0.7 # [fraction], 				https://www.proteinatlas.org/humanproteome/cell/cytosol
	vol_hepat = 3.4e-9 # [cm^3/cell], https://en.wikipedia.org/wiki/Hepatocyte
	num_hepat = 139*10^6 # [cells/g liver], https://www.sciencedirect.com/science/article/pii/S088723330600124X?via%3Dihub
	mass_liver = 1500 # [g], https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/liver-weight
	vol_hepat_tot = vol_hepat*num_hepat*mass_liver # [mL]
	
	md"""
	### Cytosolic Volume (L)
	"""
end

# ╔═╡ dcabdd7a-f3da-41db-82f7-b0f6aba9ca02
vd_cyt = vol_frac_cytosol * vol_hepat * num_hepat * mass_liver * 1/mL_to_L # [L]

# ╔═╡ 0feb0564-0a63-48ef-a748-3d3fbfe76d7d
begin
	bmr_liver = 300 # [kcal/kg/day], https://www.fao.org/3/m2845e/m2845e00.htm
	frac_liver_fatox = 0.6 # [fraction], portion of liver TEE that is from fat oxidation
	fa_basal = 0.15 # [mM], Holzhütter H-G, Berndt N. Computational Hypothesis: How Intra-Hepatic Functional Heterogeneity May Influence the Cascading Progression of Free Fatty Acid-Induced Non-Alcoholic Fatty Liver Disease (NAFLD). Cells. 2021; 10(3):578. https://doi.org/10.3390/cells10030578
	
	md"""
	### Rate of Beta Oxidation
	"""
end

# ╔═╡ 0244ae8f-b600-4087-82a3-3c104c51342a
begin
	liver_fat_ox = bmr_liver*frac_liver_fatox*mass_liver/kg_to_g # [kcal/day]
	liver_fat_ox_mmol = liver_fat_ox/rhoe_fa # [mmols-FA/day]
	kbetaox = liver_fat_ox_mmol/(fa_basal*vd_cyt) # [1/day]
	
end

# ╔═╡ 61be4a19-5772-4e31-8113-3d9cce2824d9
md"""
### Daily Chylomicron Flux
"""

# ╔═╡ 7c06328e-6039-4077-b25e-d2285ba3c4ff
begin
	diet_kcals = 2400 # [kcal/day]
	diet_frac_fat = 0.35 # [fraction]
	lcfa_frac = 0.85 # [fraction]
	ba_fat = 0.95 # [fraction], bioavailability of fat, https://www.sciencedirect.com/topics/medicine-and-dentistry/fat-intake
end

# ╔═╡ 2c37a8b3-3c8d-43ff-bdd3-4b71e5e1f9e1
chylo_flux_basal = diet_kcals*diet_frac_fat*lcfa_frac*ba_fat/rhoe_tg # [mmol-TG/day]

# ╔═╡ 7fa8770d-f9b0-4534-ba89-476987209ab4
md"""
### Basal Cytosolic Triglyceride Concentration
"""

# ╔═╡ 1eec7d02-2720-44b7-9f15-fde740c0443e
begin
	vol_frac_tg = 0.05 # [fraction]
	rho_liver = 1.07 # [g/mL], https://eurjmedres.biomedcentral.com/articles/10.1186/2047-783X-15-8-345#:~:text=Median%20liver%20density%20was%201.07%20g%2Fml.
end

# ╔═╡ 0bba348a-a829-416f-b2bf-88c588270641
tg_cyt_basal = vol_hepat_tot*vol_frac_tg*rho_tg/tg_MW*1e6/vol_hepat_tot # [mM]

# ╔═╡ f4fdacf8-3c98-45b7-be71-dc2b2fc5e9c5
md"""
### EC50 VLDL Production -- Unit Conversion
"""

# ╔═╡ 407c2009-d18e-48b1-89d6-04ff4cafb1bc
begin
	ec50_lit = 2.2879 # [%-volume]
	ec50_vldl_voltg_mL = ec50_lit/100*vol_hepat_tot # [mL]
	ec50_vldl_mM = ec50_vldl_voltg_mL*rho_tg/tg_MW/vol_hepat_tot*1e6 # [mM], includes unit conversions for mol/mL --> mmol/L
end

# ╔═╡ 7077ec03-be71-4e99-a907-f13d8839ad4b
md"""
### TG Mass Balance
"""

# ╔═╡ 21e66a06-2102-48d5-ac21-cadd3423a435
begin
	emax_vldl_prod = 33.62629534 # [mmols/day], Adiels et al. 2005.
	vldl_flux_basal = emax_vldl_prod*tg_cyt_basal/(ec50_vldl_mM+tg_cyt_basal) # [mmols-TG/day], we assumed ER ~ Cytosolic concentration at basal
	tg_input_flux = vldl_flux_basal + chylo_flux_basal # [mmols-TG]/day
	tg_plasma_basal = 1.5385 # [mM], approximately a mean daily concentration
	vd_tg_plasma = 4.5 # [L], https://doi.org/10.1007/s11745-001-0696-6
end

# ╔═╡ ea93a794-3847-4828-8b47-bd72e9d3cb44
begin
	# Liver TG mass balances
	FA_TG = 3 # [mol FA/mol TG], conversion
	fa_liver_output_flux = vldl_flux_basal*FA_TG + liver_fat_ox_mmol # [mmol-FA/day]
	fa_liver_input_flux = fa_liver_output_flux # [mmol-FA/day], mass balance at basal
	frac_fa_dnl_basal = 0.05 # [fraction], Lambert et al. 2014.
	frac_fa_diet_basal = 0.05 # [fraction], Lambert et al. 2014.
	frac_fa_other = 1 - (frac_fa_dnl_basal + frac_fa_diet_basal) # [fraction]
	dnl_flux_basal = fa_liver_input_flux*frac_fa_dnl_basal # [mmol/day]
	# Dietary contribution is trickier because VLDL is originating from multiple sources. _Assume_ that chylos are diet and everything else is "fatty acid" derived (adipose storage)
	frac_tg_plasma_diet = chylo_flux_basal/tg_input_flux # [fraction], approximation
	kuptake_liver_tg = frac_fa_diet_basal*fa_liver_input_flux/(FA_TG*
		tg_plasma_basal*vd_tg_plasma*frac_tg_plasma_diet) # [1/day]
	
end

# ╔═╡ a7080cb5-6afc-4572-9c8f-adf0c70ef071
begin
	tot_tg_liver_uptake = kuptake_liver_tg*tg_plasma_basal*vd_tg_plasma # [mmols-TG/day]
	tg_lipase_clear = tg_input_flux - tot_tg_liver_uptake # [mmols/day]
	klipase_clear = tg_lipase_clear/(tg_plasma_basal*vd_tg_plasma) # [1/day], clearance by lipases in plasma
end

# ╔═╡ 42902b71-fd1b-424d-bfb3-dcf3071a8739
nefa_uptake_flux = fa_liver_input_flux - dnl_flux_basal - 3*tot_tg_liver_uptake # [mmols-FA/day]

# ╔═╡ e4fc0ae7-43f1-4cf8-917f-f3c2fc160076
md"""
### Remaining Mass Balances
We back-solved for the majority of the major fluxes at baseline. Now back-solve for the remaining rate constants needed to support the flux.
"""

# ╔═╡ 7c1c980e-b248-4932-86b8-f2fb878a7eb5
begin
	vd_er = vd_cyt*0.16/0.5 # [L], bionumbers.org
	fa_er_basal = fa_basal # [mM], assumption
	tg_synth_er_basal = vldl_flux_basal/vd_er # [mM-TG/day]
	ksynth_er_tg = vldl_flux_basal*FA_TG/(fa_er_basal^3*vd_er) # [1/(mM^2*days)]
	kuptake_er = vldl_flux_basal*FA_TG/(fa_basal*vd_cyt) # [1/days]

	er_uptake_flux = kuptake_er*fa_basal*vd_cyt # [mmols-FA-cyt/day]
	vd_er
end

# ╔═╡ 4c648d45-18a8-4225-a055-15817039c5ce
md"""
### Lipid Droplet Dynamics
In the basal or baseline state, the lipid droplet should be approximately in equilibrium, neither growing nor shrinking dramatically. Thus, estimating one arm (i.e., lipolysis or esterification) will determine the second arm. We base our estimate on estimates of liver TG being the primary contributor to VLDL, which indicates a need for a high cycling rate.
"""

# ╔═╡ bb78de7c-a91f-42f5-8658-6e0669b108ec
begin
	vldl_frac_ltg_cyt = 0.35 # [fraction], based off Fabrinni 2008 to approximate a residence time of 3 - 4 days
	lipolysis_flux_mmols = vldl_flux_basal*vldl_frac_ltg_cyt # [mmols-TG/day], estimation of a residence time
	klipo_cy_tg = lipolysis_flux_mmols/(tg_cyt_basal*vd_cyt) # [1/day]
	ksynth_cy_tg = lipolysis_flux_mmols*FA_TG/(fa_basal^3*vd_cyt) # [1/(mM^2*day]
	ester_cyt_flux = ksynth_cy_tg*fa_basal^3*vd_cyt
end

# ╔═╡ 2c97be58-f99b-4ebc-acf1-e315a08f1be1
begin
	liver_tg_residence_time = (tg_cyt_basal*vd_cyt)/(vldl_flux_basal*vldl_frac_ltg_cyt) # [days]
end

# ╔═╡ e94e45d1-3ee0-4987-b948-05c3835bdbe2
nefa_uptake_flux + dnl_flux_basal - er_uptake_flux - liver_fat_ox_mmol + 3*tot_tg_liver_uptake + 3*lipolysis_flux_mmols - ester_cyt_flux

# ╔═╡ 255b720e-1892-4d1f-b998-530209bca608
nefa_uptake_flux

# ╔═╡ a24ed6be-793e-46ae-8d9e-99ee70c6457d
dnl_flux_basal

# ╔═╡ a6fbeffb-48b3-4012-a333-2f347efa9353
er_uptake_flux

# ╔═╡ 717e1304-44d8-48f0-aaf7-5ba30b4fc9dd
liver_fat_ox_mmol

# ╔═╡ be038ffd-7198-4dbf-8437-413763d657a7
3*tot_tg_liver_uptake

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
InteractiveUtils = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
Markdown = "d6f4376e-aef5-505a-96c1-9c027394607a"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
"""

# ╔═╡ Cell order:
# ╠═4a97d273-8195-4493-b13d-2e41de123d5d
# ╠═68c18624-20a2-4085-85c0-5efb600ee505
# ╠═fb85dfb0-8e7d-4a1f-89b8-1ffab7549e89
# ╠═0bcd76cf-e186-4c73-b701-801c933c14b9
# ╠═dcabdd7a-f3da-41db-82f7-b0f6aba9ca02
# ╠═0feb0564-0a63-48ef-a748-3d3fbfe76d7d
# ╠═0244ae8f-b600-4087-82a3-3c104c51342a
# ╠═61be4a19-5772-4e31-8113-3d9cce2824d9
# ╠═7c06328e-6039-4077-b25e-d2285ba3c4ff
# ╠═2c37a8b3-3c8d-43ff-bdd3-4b71e5e1f9e1
# ╠═7fa8770d-f9b0-4534-ba89-476987209ab4
# ╠═1eec7d02-2720-44b7-9f15-fde740c0443e
# ╠═0bba348a-a829-416f-b2bf-88c588270641
# ╟─f4fdacf8-3c98-45b7-be71-dc2b2fc5e9c5
# ╠═407c2009-d18e-48b1-89d6-04ff4cafb1bc
# ╠═7077ec03-be71-4e99-a907-f13d8839ad4b
# ╠═21e66a06-2102-48d5-ac21-cadd3423a435
# ╠═ea93a794-3847-4828-8b47-bd72e9d3cb44
# ╠═a7080cb5-6afc-4572-9c8f-adf0c70ef071
# ╠═42902b71-fd1b-424d-bfb3-dcf3071a8739
# ╠═e4fc0ae7-43f1-4cf8-917f-f3c2fc160076
# ╠═7c1c980e-b248-4932-86b8-f2fb878a7eb5
# ╠═4c648d45-18a8-4225-a055-15817039c5ce
# ╠═bb78de7c-a91f-42f5-8658-6e0669b108ec
# ╠═2c97be58-f99b-4ebc-acf1-e315a08f1be1
# ╠═e94e45d1-3ee0-4987-b948-05c3835bdbe2
# ╠═255b720e-1892-4d1f-b998-530209bca608
# ╠═a24ed6be-793e-46ae-8d9e-99ee70c6457d
# ╠═a6fbeffb-48b3-4012-a333-2f347efa9353
# ╠═717e1304-44d8-48f0-aaf7-5ba30b4fc9dd
# ╠═be038ffd-7198-4dbf-8437-413763d657a7
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
