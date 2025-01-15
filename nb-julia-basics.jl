### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ ca4d8177-4680-4b6c-bde7-c209c501c74b
using PlutoUI

# ╔═╡ 195b30b8-2242-47c8-8006-1cad7c460883
TableOfContents()

# ╔═╡ dd5df476-ce4e-11ef-2e20-a90209692b8c
md"# What is Julia?"

# ╔═╡ 4e56ef31-9cf5-457c-ae45-2e1559defcd9
md"# Why Julia?"

# ╔═╡ 7516bd87-2218-4bf3-af57-c76619e6ac0c
md"
- Solves the two language problem
- Composability of the ecosystem
- SciML ecosystem
"

# ╔═╡ df3687ec-edff-471c-b64e-fac3af177f1a
LocalResource(
	"figs/sciml.png",
	:width => 600
)

# ╔═╡ f8bb3bee-cb49-43ed-86eb-fd47563b774f
md"https://docs.sciml.ai/Overview/stable/overview/"

# ╔═╡ 3aae6539-f566-4b30-9807-16d0d6acbdd2
md"# Syntax"

# ╔═╡ 84245768-6b83-4a54-a6c3-d6e735d0e6fb
md"## Commenting"

# ╔═╡ 1d84bfa1-f9da-446e-8f11-1198725df298
# This is a comment

# ╔═╡ 92fa5f26-8e77-4a62-accb-074879e236c4
#=
this is a multi-line comment
=#

# ╔═╡ 738a6cea-f68b-4450-8238-4ab8571b82ec
md"## Arithmetic operations"

# ╔═╡ db217b80-f814-417d-8ab7-0f0749e94750
# summation
3 + 2

# ╔═╡ 1fa73171-a500-4f7c-85c2-b4de4fd14c88
# subtraction
3 - 2

# ╔═╡ 017ece0a-eea0-4e9f-8d7c-d1d8acf3e358
# multiplication
3 * 2

# ╔═╡ 1cc66228-2fec-4127-9df3-55c451f783b4
# Division
3 / 2

# ╔═╡ e842ca5b-fe04-4622-93e8-ef6074911c1c
# exponentiation
3 ^ 2

# ╔═╡ 5e3c45ba-e50b-4bf7-937c-fe86ed44fb98
# remainder
3 % 2

# ╔═╡ eb011598-a820-49ca-9184-1e3fb4508cfe
# square root
sqrt(4.0)

# ╔═╡ b11bb123-99f0-4873-ae04-1f7c42b2656a
# log
log(3.0)

# ╔═╡ 6ecc2fd9-2726-412d-9f53-b043960bbdf3
# exponential
exp(3.0)

# ╔═╡ 6398ec02-35d6-438c-b0e9-6a708ef59e2d
md"## Logical operators"

# ╔═╡ 32b1bdd1-0693-4f25-9a7e-0a6049617435
# equal
2 == 2

# ╔═╡ 267812b4-8506-47f3-aeae-d75d19d6aa50
# not equal
2 != 2

# ╔═╡ 974b7833-573e-4263-811a-44309552159e
# greater than
3 > 2

# ╔═╡ 804deb4c-f27b-4724-9e8d-7faed48e0a1b
# greater than or equal
3 >= 2

# ╔═╡ 661aaeb9-9c8b-4638-9907-a44061f16945
# less than
1 < 2

# ╔═╡ 9d79e609-5951-4a70-9618-6ee532ef57a7
# less than or equal
1 <= 2

# ╔═╡ 24f28a55-6726-4771-9c92-767dc36bbc9c
# AND
1 < 2 && 2 < 3

# ╔═╡ e659f5e3-b111-4dde-bf72-751dbff95f4b
# OR
1 < 2 || 2 > 3

# ╔═╡ 0c777c3f-5e67-4bcf-bf89-885ee4cc6c19
md"## Variables and types"

# ╔═╡ 6c376b74-274d-4096-a4a2-d778798df06a
x = 2

# ╔═╡ a76e57b3-74d1-42aa-863d-2c4ece6cd909
typeof(x)

# ╔═╡ 36bb5baa-4a99-41c2-a98e-17fbeeb2f351
xx = 2.0

# ╔═╡ 32ec9234-0ab0-49ee-bf5b-5c2a30505d37
typeof(xx)

# ╔═╡ 544c350d-f69d-407d-8343-d98941c27076
# convert type
xxx = convert(Int, xx)

# ╔═╡ 1085b786-e45d-4269-b2ce-fc6b8d35a23d
typeof(xxx)

# ╔═╡ bf837599-09c9-4ad0-9976-345df752d6cc
# use unicode
α = 3.0

# ╔═╡ 0ffe61eb-90c5-4b58-b436-5205a6473efe
🐰 = 1

# ╔═╡ 32c555bc-2961-4f99-8299-2818511938e5
md"### Strings"

# ╔═╡ aca313af-36e4-43cd-820a-26ea113e032b
name = "ahmed "

# ╔═╡ fd9fb264-b052-4acf-b2c8-ac48af829e73
name2 = "elmokadem"

# ╔═╡ 87abf7ac-f745-4a27-8d7b-b47051f7975f
## concatenate
string(name, name2)

# ╔═╡ aeaf419a-3cfc-4788-959b-e5b2981bc8c9
name * name2

# ╔═╡ a30ed6d2-d863-45fb-8131-adeac7d02853
## interpolate
println("my name is $name")

# ╔═╡ 6f70466a-e65c-4b9c-85f9-8c0d18d22b2c
begin
	age = 20
	println("I am $age years old")
end

# ╔═╡ bb3faa9e-23ce-4ac6-b32b-83a1054cbfff
println("I am $(age^2) years old")

# ╔═╡ 46da2ed2-61be-42e3-8893-d8e0d8f63ea3
println("I am ", age, " years old")

# ╔═╡ 098dfe77-902a-4dc9-b646-d60a275eefe9
md"## Data structures"

# ╔═╡ d33ce047-19e3-44c3-ae07-9e608ee8b1a2
md"### Tuples"

# ╔═╡ e4a9603a-85ed-48b8-961d-6493fb30a826
## ordered but immutable
t = (4,5,6)

# ╔═╡ e48afea8-2f61-4cd2-b057-bcb49be949e2
t[2]

# ╔═╡ b559174e-af8b-4ccb-9027-f37ef4d24361
md"### Disctionaries"

# ╔═╡ 6e2a213e-43aa-4e81-a576-ac42db527e05
## unordered but mutable
d = Dict("a" => 1, "b" => 2, "c" => 3)

# ╔═╡ 61d413d6-cc2b-41f0-986d-a1b6f03f2dea
d["a"]

# ╔═╡ 791a60ec-f264-42e1-9bbf-9d7b7193885e
d["a"] = 4

# ╔═╡ a7df3ef8-3f41-4f6b-b320-56346b4f4afd
d

# ╔═╡ fa9f01ec-85e7-4f29-b5ae-294f3901138a
md"### Arrays"

# ╔═╡ 86e666a3-6c82-4be5-8e86-459a2005b800
md"#### Vectors"

# ╔═╡ f88935d9-6c50-4b86-8f66-83f64f441ae0
## ordered and mutable
v1 = [1,2,3]

# ╔═╡ 425d32c8-7fb0-4357-b646-4f3cea2ef88c
v2 = zeros(3)

# ╔═╡ 34ec342b-d81f-4c01-b5c4-c7d5e871d9f4
v3 = ones(3)

# ╔═╡ 813cbee0-391c-443f-8d24-4e9aa74a59d6
v4 = rand(3)

# ╔═╡ ad188ef1-c43f-4c8d-a930-7e7453d17c1a
v5 = rand(1:5, 3)

# ╔═╡ 66630bf4-068e-44ad-9dc6-4dd605de782c
v6 = [1,2,"hi"]

# ╔═╡ ebaca8c0-cdb4-46b4-8d2a-b5d4d46cf43b
push!(v6, "ahmed")

# ╔═╡ 1b60d6f2-01eb-4aca-bd20-c44441b60d71
pop!(v6)

# ╔═╡ 6af0f38e-9e1e-469f-a419-c2c7ece356f1
v6

# ╔═╡ d5242d23-28fa-4c30-a3ab-e83138430646
## indexing
begin
	v6[1]
	v6[2:end]
end

# ╔═╡ 0ac408b4-921a-414f-858b-38cfbf683d57
# operations
begin
	v6 .* v6
	v6 .== v6
end

# ╔═╡ ae43b9ae-2544-496a-a770-207291330a7c
md"#### Matrices"

# ╔═╡ d7c37f0d-13c7-406c-a202-c15cd719571b
md"#### Other"

# ╔═╡ d23010ed-3ec9-4d95-b9f9-88ac0b62ce9e
md"## Loops"

# ╔═╡ 4ac9e2a2-3cd9-4535-8b6d-8cdf5cf582aa
md"### While"

# ╔═╡ d7282fac-6185-49ac-a265-f79a705e4269
md"### For"

# ╔═╡ 9b988296-794c-47a8-8f2d-6c6acf12fa9b
md"## Conditionals"

# ╔═╡ 6f0109c1-def8-4ad2-8590-5ead4f963296
md"## Functions"

# ╔═╡ 1b695a57-76ca-4071-9066-81fc93dc8802
md"# Package Management"

# ╔═╡ ac7ddc10-86b4-4860-a678-0fa687bb1364
md"# Metaprogramming"

# ╔═╡ 2c643d24-09f8-4023-8b87-bd148f71c865
md"# Plotting"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.60"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.2"
manifest_format = "2.0"
project_hash = "8aa109ae420d50afa1101b40d1430cf3ec96e03e"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─ca4d8177-4680-4b6c-bde7-c209c501c74b
# ╟─195b30b8-2242-47c8-8006-1cad7c460883
# ╟─dd5df476-ce4e-11ef-2e20-a90209692b8c
# ╟─4e56ef31-9cf5-457c-ae45-2e1559defcd9
# ╟─7516bd87-2218-4bf3-af57-c76619e6ac0c
# ╟─df3687ec-edff-471c-b64e-fac3af177f1a
# ╟─f8bb3bee-cb49-43ed-86eb-fd47563b774f
# ╟─3aae6539-f566-4b30-9807-16d0d6acbdd2
# ╟─84245768-6b83-4a54-a6c3-d6e735d0e6fb
# ╠═1d84bfa1-f9da-446e-8f11-1198725df298
# ╠═92fa5f26-8e77-4a62-accb-074879e236c4
# ╟─738a6cea-f68b-4450-8238-4ab8571b82ec
# ╠═db217b80-f814-417d-8ab7-0f0749e94750
# ╠═1fa73171-a500-4f7c-85c2-b4de4fd14c88
# ╠═017ece0a-eea0-4e9f-8d7c-d1d8acf3e358
# ╠═1cc66228-2fec-4127-9df3-55c451f783b4
# ╠═e842ca5b-fe04-4622-93e8-ef6074911c1c
# ╠═5e3c45ba-e50b-4bf7-937c-fe86ed44fb98
# ╠═eb011598-a820-49ca-9184-1e3fb4508cfe
# ╠═b11bb123-99f0-4873-ae04-1f7c42b2656a
# ╠═6ecc2fd9-2726-412d-9f53-b043960bbdf3
# ╟─6398ec02-35d6-438c-b0e9-6a708ef59e2d
# ╠═32b1bdd1-0693-4f25-9a7e-0a6049617435
# ╠═267812b4-8506-47f3-aeae-d75d19d6aa50
# ╠═974b7833-573e-4263-811a-44309552159e
# ╠═804deb4c-f27b-4724-9e8d-7faed48e0a1b
# ╠═661aaeb9-9c8b-4638-9907-a44061f16945
# ╠═9d79e609-5951-4a70-9618-6ee532ef57a7
# ╠═24f28a55-6726-4771-9c92-767dc36bbc9c
# ╠═e659f5e3-b111-4dde-bf72-751dbff95f4b
# ╟─0c777c3f-5e67-4bcf-bf89-885ee4cc6c19
# ╠═6c376b74-274d-4096-a4a2-d778798df06a
# ╠═a76e57b3-74d1-42aa-863d-2c4ece6cd909
# ╠═36bb5baa-4a99-41c2-a98e-17fbeeb2f351
# ╠═32ec9234-0ab0-49ee-bf5b-5c2a30505d37
# ╠═544c350d-f69d-407d-8343-d98941c27076
# ╠═1085b786-e45d-4269-b2ce-fc6b8d35a23d
# ╠═bf837599-09c9-4ad0-9976-345df752d6cc
# ╠═0ffe61eb-90c5-4b58-b436-5205a6473efe
# ╟─32c555bc-2961-4f99-8299-2818511938e5
# ╠═aca313af-36e4-43cd-820a-26ea113e032b
# ╠═fd9fb264-b052-4acf-b2c8-ac48af829e73
# ╠═87abf7ac-f745-4a27-8d7b-b47051f7975f
# ╠═aeaf419a-3cfc-4788-959b-e5b2981bc8c9
# ╠═a30ed6d2-d863-45fb-8131-adeac7d02853
# ╠═6f70466a-e65c-4b9c-85f9-8c0d18d22b2c
# ╠═bb3faa9e-23ce-4ac6-b32b-83a1054cbfff
# ╠═46da2ed2-61be-42e3-8893-d8e0d8f63ea3
# ╟─098dfe77-902a-4dc9-b646-d60a275eefe9
# ╟─d33ce047-19e3-44c3-ae07-9e608ee8b1a2
# ╠═e4a9603a-85ed-48b8-961d-6493fb30a826
# ╠═e48afea8-2f61-4cd2-b057-bcb49be949e2
# ╟─b559174e-af8b-4ccb-9027-f37ef4d24361
# ╠═6e2a213e-43aa-4e81-a576-ac42db527e05
# ╠═61d413d6-cc2b-41f0-986d-a1b6f03f2dea
# ╠═791a60ec-f264-42e1-9bbf-9d7b7193885e
# ╠═a7df3ef8-3f41-4f6b-b320-56346b4f4afd
# ╟─fa9f01ec-85e7-4f29-b5ae-294f3901138a
# ╟─86e666a3-6c82-4be5-8e86-459a2005b800
# ╠═f88935d9-6c50-4b86-8f66-83f64f441ae0
# ╠═425d32c8-7fb0-4357-b646-4f3cea2ef88c
# ╠═34ec342b-d81f-4c01-b5c4-c7d5e871d9f4
# ╠═813cbee0-391c-443f-8d24-4e9aa74a59d6
# ╠═ad188ef1-c43f-4c8d-a930-7e7453d17c1a
# ╠═66630bf4-068e-44ad-9dc6-4dd605de782c
# ╠═ebaca8c0-cdb4-46b4-8d2a-b5d4d46cf43b
# ╠═1b60d6f2-01eb-4aca-bd20-c44441b60d71
# ╠═6af0f38e-9e1e-469f-a419-c2c7ece356f1
# ╠═d5242d23-28fa-4c30-a3ab-e83138430646
# ╠═0ac408b4-921a-414f-858b-38cfbf683d57
# ╟─ae43b9ae-2544-496a-a770-207291330a7c
# ╟─d7c37f0d-13c7-406c-a202-c15cd719571b
# ╟─d23010ed-3ec9-4d95-b9f9-88ac0b62ce9e
# ╟─4ac9e2a2-3cd9-4535-8b6d-8cdf5cf582aa
# ╟─d7282fac-6185-49ac-a265-f79a705e4269
# ╟─9b988296-794c-47a8-8f2d-6c6acf12fa9b
# ╟─6f0109c1-def8-4ad2-8590-5ead4f963296
# ╟─1b695a57-76ca-4071-9066-81fc93dc8802
# ╟─ac7ddc10-86b4-4860-a678-0fa687bb1364
# ╟─2c643d24-09f8-4023-8b87-bd148f71c865
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
