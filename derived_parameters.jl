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
	using XLSX

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

	md"""
	#### Load parameters.xlsx into a DataFrame
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
	rhoe_fa # [kcal/mmol-fat]
	
	md"""
	### Useful constants for calculations
	"""
end

# ╔═╡ 0bcd76cf-e186-4c73-b701-801c933c14b9
begin
	vol_frac_cytosol = 0.7 # [fraction], 				https://www.proteinatlas.org/humanproteome/cell/cytosol
	vol_hepat = 3.4e-9 # [cm^3/cell], https://en.wikipedia.org/wiki/Hepatocyte
	num_hepat = 139*10^6 # [cells/g liver], https://www.sciencedirect.com/science/article/pii/S088723330600124X?via%3Dihub
	mass_liver = 1500 # [g], https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/liver-weight
	vol_hepat_tot = vol_hepat*num_hepat*mass_liver # [mL]
	vd_cyt = vol_frac_cytosol * vol_hepat * num_hepat * mass_liver * 1/mL_to_L # [L]

	# Write the values back to the DataFrame:
	df[df.Parameter .== "vd_cyt", :TV] .= vd_cyt
	df[df.Parameter .== "vd_cyt", :LV] .= vd_cyt .* df[df.Parameter .== "vd_cyt", :LV_M]
	df[df.Parameter .== "vd_cyt", :HV] .= vd_cyt .* df[df.Parameter .== "vd_cyt", :HV_M]
	
	md"""
	### Cytosolic Volume (L)
	Determine the volume of the cytosol of a hepatocyte. This value can be floated for plausible patients, but reasonably not over orders-of-magnitude.
	"""
end

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
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
InteractiveUtils = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
Markdown = "d6f4376e-aef5-505a-96c1-9c027394607a"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
XLSX = "fdbf4ff8-1666-58a4-91e7-1b58723a45e0"

[compat]
DataFrames = "~1.3.4"
XLSX = "~0.7.10"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "b153278a25dd42c65abbf4e62344f9d22e59191b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.43.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "daa21eb85147f72e41f6352a57fccea377e310a9"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.4"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "cc1a8e22627f33c789ab60b36a9132ac050bbf75"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.12"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EzXML]]
deps = ["Printf", "XML2_jll"]
git-tree-sha1 = "0fa3b52a04a4e210aeb1626def9c90df3ae65268"
uuid = "8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615"
version = "1.1.0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[XLSX]]
deps = ["Dates", "EzXML", "Printf", "Tables", "ZipFile"]
git-tree-sha1 = "7fa8618da5c27fdab2ceebdff1da8918c8cd8b5d"
uuid = "fdbf4ff8-1666-58a4-91e7-1b58723a45e0"
version = "0.7.10"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[ZipFile]]
deps = ["Libdl", "Printf", "Zlib_jll"]
git-tree-sha1 = "3593e69e469d2111389a9bd06bac1f3d730ac6de"
uuid = "a5390f91-8eb1-5f08-bee0-b1d1ffed6cea"
version = "0.9.4"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═2cad0dea-172e-4779-96b1-259d8393612c
# ╠═4a97d273-8195-4493-b13d-2e41de123d5d
# ╠═8f830c81-e3ce-445f-a338-154d4a3f9038
# ╠═68c18624-20a2-4085-85c0-5efb600ee505
# ╠═0bcd76cf-e186-4c73-b701-801c933c14b9
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
