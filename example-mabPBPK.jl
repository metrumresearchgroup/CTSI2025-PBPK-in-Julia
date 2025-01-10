### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 6f73203a-88b7-43cd-97bf-845384540559
# ╠═╡ show_logs = false
using Pkg; Pkg.activate(".");

# ╔═╡ 0bf4c192-cee6-11ef-1f91-1fd13e6da49e
using PlutoUI

# ╔═╡ 5c029511-b19e-4542-821b-90e55c892d1e
using DifferentialEquations, ModelingToolkit

# ╔═╡ 85db720a-2909-4073-957a-bb118babed7d
TableOfContents()

# ╔═╡ 434bca6d-cd0f-4a07-91bf-668a20df2699
# Reproduction of mouse model from Cao et al., 2014 “Incorporating Target-Mediated Drug Disposition in a Minimal Physiologically-Based Pharmacokinetic Model for Monoclonal Antibodies.” Journal of Pharmacokinetics and Pharmacodynamics 41 (4): 375–87.

#The model implemented here represents Model B in the publication, with the modification that the steady-state assumption was not used, and everything is expressed in concentrations instead.

# This version includes infusion dosing into the plasma compartment

function PBPK(; name)
	@independent_variables t 
	Dt = Differential(t)

	pars = @parameters begin
		s_per_h = 3600.0
		s_per_d = 3600.0*24.0
		KD_FAPa = 3.73e-11*(1e9) # DB:FAPalpha dissociation constant in cyno (M to nM)
	    k_on = 7.7e5*(1e-9)      # DB:FAPalpha binding rate constant in cyno
	    k_int = log(2)/(24.0 * s_per_h)
	    Thalf_FAPa = 24.0
	    k_syn = 0.3/s_per_h
		ISF = 0.579        # cyno total interstitial volume (L)
	    V_plasma = 0.0966  # cyno plasma volume (L)
	    V_lymph = 0.193    # cyno lymph volume (L)
	    L = 0.275*(1/s_per_d) #0.2*(1/s_per_d),  # cyno total lymph flow (L/d to L/s)
	    sig_tight = 0.945    # Sigma_1, vascular reflection coefficient for tight tissue (Table 1) 
	    sig_leaky = 0.687   # Sigma_2, vascular reflection coefficient for leaky tissue (Table 1) 
	    sig_L = 0.2        # lymphatic capillary reflection coefficient (unitless), predefined as in 'several previous PBPK models'
	    Kp = 0.8           # available fraction of ISF for mAb binding (unitless)
	    CL_p = 0.35*2.6*(10^-3)*(1/s_per_h) # clearance from plasma (mL/h/kg to L/s)
	    dose_mgkg = 1.0   # dose in mg/kg
		cyno_WT = 2.6   # cyno weight = 3 kg
		nmol_per_mol = 1.0e9
		DB_MW = 150000.0    # molecular weight
		dose = dose_mgkg*(cyno_WT)*(1/DB_MW)*(1/1e3)*(nmol_per_mol)
		infusion_dur = 30*60 # duration of infusion (30 min to s)
	    infusion_rate = dose/infusion_dur   # rate of infusion (nmol/s)
	    k_off = KD_FAPa*k_on
	    k_deg = log(2)/(Thalf_FAPa * s_per_h)
	    V_tight = 0.30108 #0.65 * ISF * Kp,
	    V_leaky = 0.16212 #0.35 * ISF * Kp,
		L_tight = 0.33 * L
    	L_leaky = 0.67 * L
	end

	vars = @variables begin
		C_plasma(t) = 0.0
		R_tight(t) = k_syn/k_deg
		C_tight(t) = 0.0
		AR_tight(t) = 0.0
		R_leaky(t) = k_syn/k_deg
		C_leaky(t) = 0.0
		AR_leaky(t) = 0.0
		C_lymph(t) = 0.0
	end
	
	eqs = [
		# Antibody concentration in plasma - Cp in 2014 paper, Eq (1)
	    Dt(C_plasma) ~ 1/V_plasma * (-L_leaky * (1 - sig_leaky) * C_plasma - L_tight * (1 - sig_tight) * C_plasma + L * C_lymph - CL_p * C_plasma + infusion_rate)
	    
	    # Receptor concentration in tight tissue compartment 
	    Dt(R_tight) ~ k_syn - R_tight*k_deg - k_on*C_tight*R_tight + k_off*AR_tight
	  
	    # Antibody concentration in tight tissue compartment
	    Dt(C_tight) ~ 1/V_tight * (L_tight * (1 - sig_tight) * C_plasma - L_tight * (1 - sig_L) * C_tight) - k_on*C_tight*R_tight + k_off*AR_tight
	
	    # Antibody:receptor complex concentration in leaky tissue compartment
	    Dt(AR_tight) ~ k_on*C_tight*R_tight - k_off*AR_tight - AR_tight*k_int
	
	    # Receptor concentration in leaky tissue compartment
	    Dt(R_leaky) ~ k_syn - R_leaky*k_deg - k_on*C_leaky*R_leaky + k_off*AR_leaky
	
	    # Antibody concentration in leaky tissue compartment
	    Dt(C_leaky) ~ 1/V_leaky * (L_leaky * (1 - sig_leaky) * C_plasma - L_leaky * (1 - sig_L) * C_leaky) - k_on*C_leaky*R_leaky + k_off*AR_leaky
	    
	    # Antibody:receptor complex concentration in leaky tissue compartment
	    Dt(AR_leaky) ~ k_on*C_leaky*R_leaky - k_off*AR_leaky - AR_leaky*k_int
	
	    # Antibody concentration in lymph
	    Dt(C_lymph) ~ 1/V_lymph * (L_tight * (1 - sig_L) * C_tight + L_leaky * (1 - sig_L) * C_leaky - L * C_lymph)
	]

	ODESystem(eqs, t, vars, pars; name=name)
end

# ╔═╡ 6228f0fb-7385-4de4-b681-1b03cf621644
@mtkbuild pbpk = PBPK()

# ╔═╡ b8fac44c-b674-4866-8b0b-dd3b648c4edc
infusion_rate_index = ModelingToolkit.parameter_index(pbpk, :infusion_rate).idx 

# ╔═╡ 6fa69291-8913-4f96-aad1-3b515294dc37
begin
	# Stop the infusion after the infusion duration ; 30 min
     affect!(integrator) = integrator.p[infusion_rate_index] = 0;
     cb = PresetTimeCallback([30*60], affect!);
end

# ╔═╡ e539ec22-55a1-4bdf-a165-87ed8973085d
begin
	tspan = (0.0, 504.5*3600.0)
    prob = ODEProblem(pbpk, [], tspan, [])
    sol = solve(prob, callback = cb, saveat = 60.0)
end

# ╔═╡ Cell order:
# ╟─0bf4c192-cee6-11ef-1f91-1fd13e6da49e
# ╟─85db720a-2909-4073-957a-bb118babed7d
# ╠═6f73203a-88b7-43cd-97bf-845384540559
# ╠═5c029511-b19e-4542-821b-90e55c892d1e
# ╠═434bca6d-cd0f-4a07-91bf-668a20df2699
# ╠═6228f0fb-7385-4de4-b681-1b03cf621644
# ╠═b8fac44c-b674-4866-8b0b-dd3b648c4edc
# ╠═6fa69291-8913-4f96-aad1-3b515294dc37
# ╠═e539ec22-55a1-4bdf-a165-87ed8973085d
