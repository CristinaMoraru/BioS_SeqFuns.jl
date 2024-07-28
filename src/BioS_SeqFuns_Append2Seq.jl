export read_appendix, append2fasta

function read_appendix(read_appendix_p::FnaP)
    appendix_seq = Vector{String}(undef, 1)
    appendix_id = Vector{String}(undef, 1)
    i = 0

    FastaReader(read_appendix_p.p) do FASTA
        for fr in FASTA
            i += 1
            appendix_seq[i] = fr[2]
            appendix_id[i] = fr[1]
        end
    end

    if length(appendix_id) > 1
        println("The appendix fasta file contains more than one sequence. Only the first sequence is used as appendix.")
    end

    return (appendix_id[1], appendix_seq[1])
end

function append2fasta(read_appendix_p::FnaP, read_fasta_p::FnaP, write_p::FnaP, position::Int64)
    appendix_id, appendix_seq = read_appendix(read_appendix_p)
    println("Appendix ID: $appendix_id")
    println("Appendix sequence: $appendix_seq")

    FastaWriter(write_p.p) do fw
        FastaReader(read_fasta_p.p) do FASTA
            for fr in FASTA
                if position == 5
                    writeentry(fw, "$(appendix_id)_$(fr[1])", "$(appendix_seq)$(fr[2])")
                elseif position == 3
                    writeentry(fw, "$(fr[1])_$(appendix_id)", "$(fr[2])$(appendix_seq)")
                end
            end
        end
    end

    return nothing
end