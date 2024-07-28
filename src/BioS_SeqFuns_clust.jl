export get_contig_pairs

function get_contig_pairs(inref::FnaP)
    contigs = get_contig_names(inref)

    dft = DataFrame()
    for i in eachindex(contigs)
        query = fill(contigs[i], length(contigs))
        ref = contigs

        dftmp = DataFrame(query = query, ref = ref)
        dft = vcat(dft, dftmp)
    end

    return dft
end