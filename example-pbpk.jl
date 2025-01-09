### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 8b19cc51-8ddf-4f4d-8e85-06953641d9f7
# ╠═╡ show_logs = false
using Pkg; Pkg.activate(".");

# ╔═╡ d242c44a-ce2e-11ef-3c3a-b366eaada965
using PlutoUI

# ╔═╡ 25751995-8656-4346-8dca-2df4c342b095
using DifferentialEquations, ModelingToolkit, Unitful

# ╔═╡ 1f85fbdc-7a64-4def-8e8f-2aafec12bdbb
TableOfContents()

# ╔═╡ 98a5a8f5-79ef-4cff-ab2e-592b5358a90f
md"# PBPK example"

# ╔═╡ 3179d612-eb4e-4007-b3c5-18bbc2d8e4ed
md"
Let's build a PBPK model for voriconazole
"

# ╔═╡ 9d8005dc-89bd-42e4-a0b5-abcb75d3e80b
# ╠═╡ show_logs = false
LocalResource(
	"figs/pbpk.png",
	:width => 400
)

# ╔═╡ e308fc4a-157d-486b-8d06-317cbad4283d
md"## Libraries"

# ╔═╡ e8c04c18-9ed6-4e06-a7a9-362b6bb55de2
md"## Model development"

# ╔═╡ 57966d02-b00e-4e65-b46d-43d9d343123d
PBPK = function(; name)
	@independent_variables t, [description = "time"]
	Dt = Differential(t)

	pars = @parameters begin
		# volumes (L); source: https://www.ncbi.nlm.nih.gov/pubmed/14506981
		Vad = 18.2, [description = "adipose volume"]
  		Vbo = 10.5, [description = "bone volume"] 
  		Vbr = 1.45, [description = "brain volume"]
  		VguWall = 0.65, [description = "gut volume"]
  		VguLumen = 0.35, [description = "gut volume"]
  		Vhe = 0.33, [description = "heart volume"]
  		Vki = 0.31, [description = "kidney volume"]
  		Vli = 1.8, [description = "liver volume"]
  		Vlu = 0.5, [description = "lung volume"]
  		Vmu = 29, [description = "muscle volume"]
  		Vsp = 0.15, [description = "spleen volume"]
  		Vbl = 5.6, [description = "blood volume"]

		# blood flows (L/h); Cardiac output = 6.5 (L/min); source: https://www.ncbi.nlm.nih.gov/pubmed/14506981 
		Qad = 0.05*6.5*60
  		Qbo = 0.05*6.5*60
  		Qbr = 0.12*6.5*60
  		Qgu = 0.15*6.5*60 
  		Qhe = 0.04*6.5*60
  		Qki = 0.19*6.5*60
  		Qmu = 0.17*6.5*60
  		Qsp = 0.03*6.5*60
  		Qha = 0.065*6.5*60  
  		Qlu = 6.5*60 

		# partition coefficients estimated by Poulin and Theil method https://jpharmsci.org/article/S0022-3549(16)30889-9/fulltext
  		Kpad = 9.89   # adipose:plasma
  		Kpbo = 7.91   # bone:plasma
  		Kpbr = 7.35   # brain:plasma
  		Kpgu = 5.82   # gut:plasma
  		Kphe = 1.95   # heart:plasma
  		Kpki = 2.9    # kidney:plasma
  		Kpli = 4.66   # liver:plasma
  		Kplu = 0.83   # lungs:plasma
  		Kpmu = 2.94   # muscle:plasma; optimized
  		Kpsp = 2.96   # spleen:plasma
  		Kpre = 4      # calculated as average of non adipose Kps
  		BP = 1        # blood:plasma ratio

		# other parameters
  		WEIGHT = 73
  		ka = 0.849   # absorption rate constant (/hr) 
  		fup = 0.42   # fraction of unbound drug in plasma
  
  		# in vitro hepatic clearance parameters http://dmd.aspetjournals.org/content/38/1/25.long
  		fumic = 0.711 # fraction of unbound drug in microsomes
  		MPPGL = 30.3  # adult mg microsomal protein per g liver (mg/g)
  		VmaxH = 40    # adult hepatic Vmax (pmol/min/mg)
  		KmH = 9.3     # adult hepatic Km (uM)
  
  		# renal clearance  https://link.springer.com/article/10.1007%2Fs40262-014-0181-y
  		CL_Ki = 0.096 # (L/hr) renal clearance
	end

	vars = @variables begin
		GUTLUMEN(t) 
		GUT(t) 
		ADIPOSE(t) 
		BRAIN(t) 
		HEART(t) 
		BONE(t) 
  		KIDNEY(t) 
		LIVER(t) 
		LUNG(t) 
		MUSCLE(t) 
		SPLEEN(t) 
		REST(t) 
  		ART(t) 
		VEN(t)
		CP(t)
	end

	# additional volume derivations
  	Vve = 0.705*Vbl  # venous blood
    Var = 0.295*Vbl  # arterial blood
    Vre = WEIGHT - (Vli+Vki+Vsp+Vhe+Vlu+Vbo+Vbr+Vmu+Vad+VguWall+Vbl)  # volume of rest of the body compartment
  
    # additional blood flow derivation
    Qli = Qgu + Qsp + Qha
    Qre = Qlu - (Qli + Qki + Qbo + Qhe + Qmu + Qad + Qbr)
  
    # intrinsic hepatic clearance calculation
  	CL_Li = ((VmaxH/KmH)*MPPGL*Vli*1000*60*1e-6) / fumic  # (L/hr) hepatic clearance

	# Calculation of tissue drug concentrations (mg/L)
  	Cadipose = ADIPOSE/Vad
    Cbone = BONE/Vbo
    Cbrain = BRAIN/Vbr 
    Cheart = HEART/Vhe 
    Ckidney = KIDNEY/Vki
    Cliver = LIVER/Vli 
    Clung = LUNG/Vlu 
    Cmuscle = MUSCLE/Vmu
    Cspleen = SPLEEN/Vsp
    Crest = REST/Vre
    Carterial = ART/Var
    Cvenous = VEN/Vve
    CgutLumen = GUTLUMEN/VguLumen
    Cgut = GUT/VguWall

	observed = [
		CP ~ Cvenous/BP
	]

	eqs = [
  		Dt(GUTLUMEN) ~ -ka*GUTLUMEN,
  		Dt(GUT) ~ ka*GUTLUMEN + Qgu*(Carterial - Cgut/(Kpgu/BP)),
  		Dt(ADIPOSE) ~ Qad*(Carterial - Cadipose/(Kpad/BP)),
  		Dt(BRAIN) ~ Qbr*(Carterial - Cbrain/(Kpbr/BP)),
  		Dt(HEART) ~ Qhe*(Carterial - Cheart/(Kphe/BP)),
  		Dt(KIDNEY) ~ Qki*(Carterial - Ckidney/(Kpki/BP)) - CL_Ki*(fup*Ckidney/(Kpki/BP)),
  		Dt(LIVER) ~ Qgu*(Cgut/(Kpgu/BP)) + Qsp*(Cspleen/(Kpsp/BP)) + Qha*(Carterial) - Qli*(Cliver/(Kpli/BP)) - CL_Li*(fup*Cliver/(Kpli/BP)),
  		Dt(LUNG) ~ Qlu*(Cvenous - Clung/(Kplu/BP)),
  		Dt(MUSCLE) ~ Qmu*(Carterial - Cmuscle/(Kpmu/BP)),
  		Dt(SPLEEN) ~ Qsp*(Carterial - Cspleen/(Kpsp/BP)),
  		Dt(BONE) ~ Qbo*(Carterial - Cbone/(Kpbo/BP)),
  		Dt(REST) ~ Qre*(Carterial - Crest/(Kpre/BP)),
  		Dt(VEN) ~ Qad*(Cadipose/(Kpad/BP)) + Qbr*(Cbrain/(Kpbr/BP)) +
    Qhe*(Cheart/(Kphe/BP)) + Qki*(Ckidney/(Kpki/BP)) + Qli*(Cliver/(Kpli/BP)) + 
    Qmu*(Cmuscle/(Kpmu/BP)) + Qbo*(Cbone/(Kpbo/BP)) + Qre*(Crest/(Kpre/BP)) - Qlu*Cvenous,
  		Dt(ART) ~ Qlu*(Clung/(Kplu/BP) - Carterial)
	]

	ODESystem(eqs, t, vars, pars; name=name, observed=observed)
end

# ╔═╡ 106b1b40-4458-453d-b5bd-931969c9471b
@mtkbuild pbpk = PBPK()

# ╔═╡ d63b718e-1847-4c92-90c7-6b507293582f
md"## Simulation"

# ╔═╡ Cell order:
# ╟─d242c44a-ce2e-11ef-3c3a-b366eaada965
# ╟─1f85fbdc-7a64-4def-8e8f-2aafec12bdbb
# ╟─98a5a8f5-79ef-4cff-ab2e-592b5358a90f
# ╟─3179d612-eb4e-4007-b3c5-18bbc2d8e4ed
# ╟─9d8005dc-89bd-42e4-a0b5-abcb75d3e80b
# ╟─e308fc4a-157d-486b-8d06-317cbad4283d
# ╠═8b19cc51-8ddf-4f4d-8e85-06953641d9f7
# ╠═25751995-8656-4346-8dca-2df4c342b095
# ╟─e8c04c18-9ed6-4e06-a7a9-362b6bb55de2
# ╠═57966d02-b00e-4e65-b46d-43d9d343123d
# ╠═106b1b40-4458-453d-b5bd-931969c9471b
# ╟─d63b718e-1847-4c92-90c7-6b507293582f
