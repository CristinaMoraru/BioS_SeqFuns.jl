using BioS_Gen
using BioS_SeqFuns
using BioSequences
using FastaIO
using CSV
using DataFrames


in_p = "/home/cmoraru/MY_SPACE/BioinfSuite_4_GitHUB/BioS_SeqFuns.jl/src/tests/seq/s1.fna" |> FnaP
in_p = "/home/cmoraru/MY_SPACE/BioinfSuite_4_GitHUB/BioS_SeqFuns.jl/src/tests/seq/long.fna" |> FnaP


function detectContigShapeEnds(in_p::FnaP)
    df = DataFrame(contig_name = Vector{String}(), contig_shape = Vector{String}(), contig_end_repeat_type = Vector{String}(), contig_end_repeat_size = Vector{Integer}(), contig_trimmed_end = Vector{Int64}(), contig_full_end = Vector{Int64}())
    
    contig_length = get_contig_length(in_p)

        FastaReader(in_p.p) do FASTA
            iterator = 0

            for fr in FASTA
                iterator += 1
                
                max_repeat_len = Int(floor(contig_length[iterator]/2))
                
                cd = 0
                ci = 0
                dtr = 0
                itr = 0
                l = length(fr[2])
                
                for i in 20:max_repeat_len
                    try
                        if fr[2][1:i] == fr[2][(l-i+1):l]
                            cd = i
                            dtr += 1
                        elseif fr[2][1:i] == rc_string(fr[2][(l-i+1):l])
                            ci = i
                            itr += 1 
                        end
                    catch
                        break  # it will break when it finds the first N base
                    end
                end

                if cd == 0 && ci == 0  # if no repeats where found, cd and cd are zero
                    shape = "linear"
                    end_type = ""
                    cend = l
                    repeat_size = 0
                else                 
                    if cd > ci
                        shape = "circular"
                        end_type = "DTR"
                        
                        cend = (l-cd)
                        repeat_size = cd
                    elseif ci > cd
                        shape = "linear"
                        end_type = "ITR"

                        cend = l
                        repeat_size = ci
                    end
                end
                
                push!(df, Dict(:contig_name => fr[1], :contig_shape => shape, :contig_end_repeat_type => end_type, :contig_end_repeat_size => repeat_size,  :contig_trimmed_end => cend, :contig_full_end => l))
            end
        end

    return df
end

dr = detectContigShapeEnds(in_p)

