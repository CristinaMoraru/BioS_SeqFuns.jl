export fasta2mfasta, mfasta2fasta, mfasta2mfasta

function fasta2mfasta(fstlst::Vector{FnaP}, write_p::FnaP; addprefix::Bool=false)
    FastaWriter(write_p.p) do fw
        for fstp in fstlst
            FastaReader(fstp.p) do FASTA
                for fr in FASTA
                    if addprefix
                        id = getFileName(fstp.p)   
                        @suppress writeentry(fw, "$(id)_$(fr[1])", "$(fr[2])")
                    else
                        @suppress writeentry(fw, "$(fr[1])", "$(fr[2])")
                    end
                end
            end
        end
    end

    return nothing
end
    

function mfasta2fasta(mfasta_p::FnaP, write_p::String; prefix::String="", extension::String="fasta",
                        removeprev::Bool=false)
    if removeprev
        rm_mkpaths([write_p])
    else
        mkpath(write_p)
    end

    if prefix != ""
        prefix *= "_"
    end

    FastaReader(mfasta_p.p) do FASTA
        for fr in FASTA
            outp = "$(write_p)/$(prefix)$(fr[1]).$(extension)"

            FastaWriter(outp) do fw
                @suppress writeentry(fw, "$(prefix)$(fr[1])", "$(fr[2])")
            end

        end
    end 
    
    return nothing
end

function mfasta2mfasta(inp::FnaP, outp::FnaP; prefix::String="")
    
    namesdict = Dict{String, String}()
    i=0
    FastaWriter(outp.p) do fw
        FastaReader(inp.p) do FASTA
            for fr in FASTA
                i += 1
                namesdict["$(prefix)_$(i)"] = fr[1]
                @suppress writeentry(fw, "$(prefix)_$(i)", "$(fr[2])")
            end
        end
    end

    return namesdict
end