### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ c7030848-59db-42d9-b827-733f2afccd0f
# ╠═╡ show_logs = false
using Pkg; Pkg.activate(".");

# ╔═╡ 76902963-43fc-405d-9132-2ed8409fa56d
using PlutoUI

# ╔═╡ 01d1f036-99f2-4ef6-8ab3-f1e5180fd5c5
using DifferentialEquations, ModelingToolkit, Unitful, Plots;

# ╔═╡ 179bbe7a-16f1-466f-baa0-a132a4f916ce
TableOfContents()

# ╔═╡ c3f2054a-8df4-41ea-9062-e08fb4b1792b
md"
# PK example
"

# ╔═╡ a0e14355-0fbd-4b98-8790-e7ec50910a61
md"
Let's build a one-compartment PK model with oral absorption
"

# ╔═╡ 83701d28-60a8-4c11-8923-ef7be27d384d
# ╠═╡ show_logs = false
LocalResource(
	"figs/pk.png",
	:width => 400
)

# ╔═╡ 3a1c3888-1671-4bcc-bd0f-884d6e0e9365
md"
## Load libraries
"

# ╔═╡ c4a4b27e-d053-4f70-aa2a-5ab7538780ee
md"
## Build model
"

# ╔═╡ ac7b62aa-29ec-4b61-a01e-8a582b10e928
function PK(; name)
    
    @independent_variables t [description = "time", unit = u"hr"]
    Dt = Differential(t)
    
    pars = @parameters begin
        CL = 4.0, [description = "clearance", unit = u"L*hr^-1"]
        V = 35.0, [description = "volume", unit = u"L"]
        ka = 2.0, [description = "absorption rate constant", unit = u"hr^-1"]
    end

    vars = @variables begin
		(depot(t) = 320.0), [description = "Theophylline amount in depot compartment", unit = u"mg"]
        (cent(t) = 0.0), [description = "Theophylline amount in central compartment", unit = u"mg"]
		c_cent(t), [description = "Theophylline concentration in central compartment", unit = u"mg*L^-1"]
	end
	
    observed = [
        c_cent ~ cent/V
    ]

	eqs = [
		Dt(depot) ~ -ka * depot,
		Dt(cent) ~ ka * depot - (CL / V) * cent
	]

	ODESystem(eqs, t, vars, pars; name=name, observed=observed, checks = ~ModelingToolkit.CheckUnits)  # 100124 unit checking is broken at this point so will suppress for now
end

# ╔═╡ 3604e6c4-a29d-45ee-af90-ae6d480a185b
@mtkbuild pk = PK()

# ╔═╡ dd493723-8a28-41ca-9693-7ad328a4a618
ModelingToolkit.UnitfulUnitCheck.get_unit(pk.CL)

# ╔═╡ c1505fa8-b7fa-4f24-9a37-575cb87dcd89
ModelingToolkit.getdescription(pk.cent) 

# ╔═╡ b4c202a6-5561-46ab-98a8-ed157c6d9742
ModelingToolkit.unknowns(pk)

# ╔═╡ 5cf34c06-cc4a-4baf-8c9b-1017d006c971
ModelingToolkit.equations(pk) 

# ╔═╡ bbbe3a0c-7553-46c4-9543-ce8afba41cb5
ModelingToolkit.observed(pk)

# ╔═╡ 7ffd9484-d960-4d7a-b9ae-ca09d0ba5e62
md"
## Run simulation
"

# ╔═╡ df3ceb01-8b50-4771-89cb-7c54b77f5a03
md"
### Create ODE problem
"

# ╔═╡ 1a062802-ee3b-493c-b260-9a8c784a1f21
prob = ODEProblem(pk, [], (0.0, 25.0), [])

# ╔═╡ c49719ee-7f27-4eaf-aeb4-ee551e91427c
md"
### Solve
"

# ╔═╡ 445b8741-6716-4684-8c07-2e171f070abb
sol = solve(prob, Tsit5())

# ╔═╡ afe9cc47-1122-4376-944d-377bbd5fa0a4
md"
## Plot
"

# ╔═╡ 8166670f-bd07-4c79-99ef-ebe899177365
plot(sol)

# ╔═╡ 50a44f05-10ae-4f89-b6cc-68358efa126f
plot(sol, idxs = :cent)

# ╔═╡ Cell order:
# ╟─76902963-43fc-405d-9132-2ed8409fa56d
# ╟─179bbe7a-16f1-466f-baa0-a132a4f916ce
# ╟─c3f2054a-8df4-41ea-9062-e08fb4b1792b
# ╟─a0e14355-0fbd-4b98-8790-e7ec50910a61
# ╟─83701d28-60a8-4c11-8923-ef7be27d384d
# ╟─3a1c3888-1671-4bcc-bd0f-884d6e0e9365
# ╠═c7030848-59db-42d9-b827-733f2afccd0f
# ╠═01d1f036-99f2-4ef6-8ab3-f1e5180fd5c5
# ╟─c4a4b27e-d053-4f70-aa2a-5ab7538780ee
# ╠═ac7b62aa-29ec-4b61-a01e-8a582b10e928
# ╠═3604e6c4-a29d-45ee-af90-ae6d480a185b
# ╠═dd493723-8a28-41ca-9693-7ad328a4a618
# ╠═c1505fa8-b7fa-4f24-9a37-575cb87dcd89
# ╠═b4c202a6-5561-46ab-98a8-ed157c6d9742
# ╠═5cf34c06-cc4a-4baf-8c9b-1017d006c971
# ╠═bbbe3a0c-7553-46c4-9543-ce8afba41cb5
# ╟─7ffd9484-d960-4d7a-b9ae-ca09d0ba5e62
# ╟─df3ceb01-8b50-4771-89cb-7c54b77f5a03
# ╠═1a062802-ee3b-493c-b260-9a8c784a1f21
# ╟─c49719ee-7f27-4eaf-aeb4-ee551e91427c
# ╠═445b8741-6716-4684-8c07-2e171f070abb
# ╟─afe9cc47-1122-4376-944d-377bbd5fa0a4
# ╠═8166670f-bd07-4c79-99ef-ebe899177365
# ╠═50a44f05-10ae-4f89-b6cc-68358efa126f
