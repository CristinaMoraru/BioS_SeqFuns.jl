using BioS_SeqFuns
using Documenter

DocMeta.setdocmeta!(BioS_SeqFuns, :DocTestSetup, :(using BioS_SeqFuns); recursive=true)

makedocs(;
    modules=[BioS_SeqFuns],
    authors="Cristina Moraru",
    sitename="BioS_SeqFuns.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
