using BioS_Gen
using BioS_SeqFuns
using BioSequences
using FastaIO
using CSV
using DataFrames


#in_p = "/home/cmoraru/MY_SPACE/BioinfSuite_4_GitHUB/BioS_SeqFuns.jl/src/tests/seq/s1.fna" |> FnaP
#in_p = "/home/cmoraru/MY_SPACE/BioinfSuite_4_GitHUB/BioS_SeqFuns.jl/src/tests/seq/long.fna" |> FnaP
in_p = "/mnt/ASBRU/Projects/RESIST_phaseI_viruses/TEMP-TRANSFER/debug_dovipv09/ES22_IMP_S_C13_hybrid_min1000/N-02_checkV_Nonintegrated/N-02_00_aggregated_nonintegrated_virus_contigs.fna" |> FnaP

function find_dtrs(seq::String, min_repeat_len::Int, middle::Int, len::Int)

    tr = missing
    matches = findall("$(seq[1:min_repeat_len])", seq[middle:len])
    for match in matches
        pos2 = middle + match[1] - 1 
        pos1 = len - pos2 + 1

        if seq[pos2:len] == seq[1:pos1]
            tr = terminal_repeat("DTR", pos1, 1, pos1, pos2, len)
        end
    end

    if ismissing(tr)
        return missing
    else
        return tr
    end

end

function find_itrs(seq::String, min_repeat_len::Int, middle::Int, len::Int)
    # Initialize pos1 with a value smaller with one than the min_repeat_len
    pos1 = (min_repeat_len-1)
    for i in min_repeat_len:(middle-1)
        if seq[1:i] == rc_string(seq[(len-i+1):len])
            continue
        else
            #if the loop breaks at the first iteration, pos1 would have the same value as at initialization, meaning that there is no ITR
            # if the loop breaks after several iterations, pos1 would have the value of the last position where ITR was found
            pos1 = i-1
            break
        end
    end


    if pos1 > (min_repeat_len-1)
        tr = terminal_repeat("ITR", pos1, 1, pos1, (len-pos1+1), len)
    else
        tr = missing
    end

    return tr
end

"""
    detectContigShapeEnds(in_p::FnaP)
    Detects the presence of direct terminal repeats and inverted terminal repeats at the ends of the contigs. 
    If DTRs are present, it declares the contig of a circular shape. If no terminal repeats or ITRs are found, it declares the contig as linear.

"""

#=struct Fr_str
    name::String
    seq::String
end =#

#fr = Fr_str("ES22_IMP_S_C13_hybrid_87579_length_3858_cov_1", "TGAAGTTGCGCAGCTACTTCTTTACTATCCTTATCCCTGGTGGTGGTCACCTAACCCTCATGAGCTCGAGCTCCTTAAGGAGCTCGAGCTCATGAGGGTTAGGTGACCACCACCAGGGATAAGGATAGTAAAGAAGTAGCTGCGCAACTTCACCATAAAACAACATATTTTTGTAGAAGCTAAATGTATTAACTTCTACGAGTATACGTAGTAACACGAACACTCTTAGGTTGTCCGTCTCCTAAGCGATGGCGATGTGACGTGTCTCCTACTTTGTACCTCTTTCGCCAAGTCCTTCCCATTCCTCGCCTGCTGTTCCTCATCCTCGATTGTGCCTGTTGATGCCATTCATCCAGTTGTCTCTGATACTGACCTCTCGTTACCGCCATTATGACGTCACTGAACTGTTAATCGTAACACTCAGAGGTTGATTGGACAGAGGCCACTTCAACGTATCCTTCAGCTCCTTGAACTCCACAAGGAAGTTTAACTTGATCCGGATACGATACTTTAGGTCGTCACCGTATACCAAAGCGTCGTTTAGCCCAATCGCCTCGTTATAGGGCTTGATACGAATGATCAGACGCTCTGTCAACGATGGATTTGTTGTCACCAACGTCCAGTTTTCCACCTCGGAATCCAGTCGAACTTCTCGTTTGTAATCTCCAGGTCTGTACTCCCCGGAGAAAAAACAAGAATGACGTCCAGAAGCACCAACGTGATTCGTCGATTCGTAACTGTCCAGCGCCCCAGCTCCTGCTGTTCCCAGTTCTCGGTTTTCCGCACCCGTAACTTCACCAAGCTCCAGGTTACCACCTGCTTGGACCGCCAGGTAACGTCTATCCAAGTAGTGATATTTCACACCGTTCCACAATTGGTGGTCAGCGTTGGAACCTTCACCGGGAGGATGGGTGTCGTTACTGAACATTTCATATGCCCAAAACGGTTCTCCACTCATGTTCTCGATGTACACATTGTACTGGCATGAAACCACATGGTAGTAGTTGTAGATGCCAGCGTAGTAGTCAAACCATCGGGCACTACGAATGGATGATTCGTTCAGCTGGGACACGAATGTTGTTCCGGCTCCAATGTTCAAATCCACTTGCGCCATTTCCGTTCGATAATCATACGGACTCGTCATGCGATACACATGATCGCGGGAGTAATCGTTGCTTGAAGTAACGTCATACTGTTCAATCACAGGCAACGTTACAAACGTATAGGGTTCAGGGCCCCTCGTAACGTTGTACACATTGTCGACAGGGGTTTCCCGTAACCCAGCTGCGTTTCCGGAACCGTTTCCGTCTGCCATGTTTGCTCCTGTCTGCAAATTCGTGACTGACTCCATTTCCACATCGTCGGATTCATCAACTGGAATGCCTTTTGTGCGCATCTTTTTTTCCGCGGGAAGTTCTGGCTCTATACCACGTAGCCTTTTGTTGCCCGCCTTCTTGATGGGTGTAGTCATGGCGTCATCGTCCACGGATGAGATGAGTCCAGCGTCGTAGGCTGCCTTTTTTAACCCAAAAAAAGTCTTTGCAAGTGCCCCGCCATAATCTCCAGTAGTAAAGCTGCGATAGGCTGCAGCGTCGGCTTCGCTCCAATGTAGATACGGATTACCACCTCGATGGATAATCTTACCATAACCAATGTCATGGATCTTTGCAACAAAGTCATTGAAGTTCGATGGTGTCCCCTTATCGATCTTGTTACCAGGACCGAGATATTTGTATCCCGGAAGTGATAGGTCCCACATGTCTAATTCAAGACTGACTACACACGGGCACGTTTCTCACCAACTGGCTTACCGTAGTACGCCTCCTTGATCATCTTCTGGTATTCTTCTGGCGACTTGAGCTCTTCTCCCGTCGGGATGTCACGTAGATAGTCCCACGCCATGAAGATGCGGGTAGGGTCCGTGTTGCCTGCGGCCTTGAGTTGATTGATCTCGTACTTGGACAGCAGCTTGAGTCCGATCAGACGGAGTAGTTCACGAATGTTCACCTCCCAGAAACGGCGCTCCAACGCACGCAGCGTATCGTTGCGACCCTTGGTATCCTCGGGGAGAACGTCACTGATGGTGAAGTTGCTGGAGACAAGAACCGGACGATGGATGGGCATCGGAGACTTGTACTTCTGATCGATTGGGAACCCCGCCTCATCACAGATGGTCTTGAAGAACTGGGTACCGAGCGTTTCAAACGCGTTGTGGTCAACGTCACTGAGGAGGATCTGGTGTGGAACTTGGCGTCGAACAGGTCGAAGAAGCGGTTCATCAACGTCTTGTTGTAGTAGTTCGGGTAGATGACTTGGAGGATGGCCGATTTGCCGCTGCCAGGTGATCCTGTGATCCAGATGTTTGGGTTGCCCTCGGTCTTGTAGAAGTCGCGCTTTTGGCCAATCAGACTCTTGATTTTTTCACCATACTGCACAGCGGCGCGTGGATAGGTCTTGAAGGCCTCTTCGTCCTTTCCCTCTTCATAGAGACGACGCAGATCGATGATGACCTCGTTCATGGTGGCTTTCTTTTCTTGCGTGGACCTCTTGACGATGTCGTTGACAATCGTGCGGTCAGGAGGCAGCGTTCAATACTCGAAGAGCACGTTGTTCTCCGGGATCATCTTGGTAGCCGCCTTGATGTTGTGGGACTTCCATCCAGAGTCCGGGAGGTTCTTGTTCCGAGGGACCAGATAGTATCCGTAGCCGCTCATGATTTCCAGGTTCTTGATGATCGCAGACTTCGTGGTGCGGTTGTTGTACACCAGGCAGACGTGAACATGTTGGATCTTGTAATCGTCTTGGTGGGGGCTCTCACCAATCTCGACACCACCGACGAGCACGTAGCGGAACTTGTTGTTCTTGACCTCGTCTTTGACAGCGTTGACGAACCGGTCAATGTCAGTGTCCGTCGGCAGATTGAAGCGGGCGTCCCATTCGCGGTCCATCGTGTGCATTGGCACAAAACGAGAAGAAGACATGGTAGCTTGTGCACTTGGTAGTATTTTGCCGAATGACTTACTCAACCTGTTCTTCACTCTCGTCCAGAATGGTTTCTACCTCTGTGTCTGAGTCACTACCGAACAGCTCAGCTGCTTCACGCTCATCTTGCTCTTGCAGGTCCACAAACGCAACAGTGGAAGCAAAACGGTATCGCAAGCCCGGCTGCTCAACGAAGATCGTATGGAGCACGTTGGAGACCGCCACGTGGTCTAACTCACGCACGCGTGCCTCGACTGTCGCTTCTGCTAGCTGCTCGTTGGCTCGTTGAAGGTTTGTGTTCAGGGCGCGGATGGAGAAGTCGCGGACGGTGACGCTACGACGGAGACGGCGCACTTGGTCACGCAGAAGCTCGTTCTCACGACGCCAGTGGCTGTCGATGAGTTCCGAGATGATGTGCTTGGCCTCGGCGGTTGGGACGGGTTTCATGCTCTCCTTGCGATGGATCCATTTGCCCTTGAGCTGTGCTTTGGTCTTGGCCATGGCAGATTGAGAAGTGAACTCGGTTGGAATTGGCCGGGTGATTCTAGGCGCCGTCGGATAAATAAATCCAGTAATGAATCCAGTCGCCGTCGGATAAATGAATCTAGTTGTGAATCTGGCCGCCGTCGGACAAATAAATCCAGTAGTGAATCTAGATATATGGCCAGTTAGCGGGGTATAAAAAGGAGAATGAAAATTTTGAATCCTAACAGATCCTTGAAGATGAAGTTGCGCAGCTACTTCTTTACTATCCTTATCCCTGGTGGTGGTCACCTAACCCTCATGAGCTCGAGCTCCTTAAGGAGCTCGAGCTCATGAGGGTTAGGTGACCACCACCAGGGATAAGGATAGTAAAGAAGTAGCTGCGCAACTTCA")
#fr = Fr_str("end", "TGAAGTTGCGCAGCTACTTCTTTACTATCCTTATCCCTGGTGGTGGTCACCTAACCCTCATGAGCTCGAGCTCCTTAAGGAGCTCGAGCTCATGAGGGTTAGGTGACCACCACCAGGGATAAGGATAGTAAAGAAGTAGCTGCGCAACTTCA")
#function detectContigShapeEnds(in_p::FnaP)
    df = DataFrame(contig_name = Vector{String}(), contig_shape = Vector{String}(), contig_end_repeat_type = Vector{String}(), contig_shape_remarks = Vector{String}(),
                    contig_end_repeat_size = Vector{Integer}(), contig_trimmed_end = Vector{Int64}(), contig_full_end = Vector{Int64}())
    
    contig_length = get_contig_length(in_p)

        FastaReader(in_p.p) do FASTA
            iterator = 0

            for fr in FASTA
                iterator += 1     
  
                l = length(fr[2]) 
                 #=
                 if(iterator == 4755)
                    println(fr[1])
                    println(fr[2])
                    println(length(fr[2]))
                end =#

                if isodd(l)
                    middle = (l - Int(floor(l/2))) 
                else
                    middle = (l - Int(floor(l/2)) + 1) 
                end 
                
                # the Ns in sequences can lead to errors, therefore I'm wrapping this code in a try catch
                dtr = nothing
                try
                    dtr = find_dtrs(fr[2], 20, middle, l)
                catch
                    dtr = missing
                end

                itr = nothing
                try
                    itr = find_itrs(fr[2], 20, middle, l)
                catch
                    itr = missing
                end


                # choose the type of terminal repeat. If both dtr and itr are not missing, chose the one with the longest length
                if ismissing(dtr) && ismissing(itr)
                    tr = missing
                elseif !ismissing(dtr) && ismissing(itr)
                    tr = dtr
                elseif ismissing(dtr) && !ismissing(itr)
                    tr =itr
                elseif !ismissing(dtr) && !ismissing(itr)
                    if dtr.length > itr.length
                        tr = dtr
                    elseif dtr.length < itr.length
                        tr = itr
                    elseif dtr.length == itr.length    #changed here!!!!!!!!
                        tr = terminal_repeat(type = "DTRandITR",
                                              length = dtr.length,
                                              left_start = dtr.left_start,
                                              left_end = dtr.left_end,
                                              right_start = dtr.right_start,
                                              right_end = dtr.right_end)
                    end
                end


                if ismissing(tr) 
                    shape = "linear"
                    end_type = ""
                    cend = l
                    repeat_size = 0
                    shape_remarks = ""
                else                 
                    if tr.type == "DTR"
                        shape = "circular"
                        end_type = "DTR"
                        
                        cend = (tr.right_start  - 1)
                        repeat_size = tr.length

                        if repeat_size == l/2
                            shape_remarks = "duplicated sequence, DTR"
                        else
                            shape_remarks = ""
                        end

                    elseif tr.type == "ITR"
                        shape = "linear"
                        end_type = "ITR"

                        cend = l
                        repeat_size = tr.length

                        if repeat_size == l/2
                            shape_remarks = "duplicated sequence, ITR"
                        else
                            shape_remarks = ""
                        end
                    elseif tr.type == "DTRandITR"  #changed here
                        shape = "linear"
                        end_type = "DTRandITR"
                        cend = l
                        repeat_size = tr.length
                        shape_remarks = ""

                        if repeat_size == l/2
                            shape_remarks = "duplicated sequence, ITR"
                        else
                            shape_remarks = ""
                        end
                    end
                end
                
                push!(df, Dict(:contig_name => fr[1], :contig_shape => shape, :contig_end_repeat_type => end_type, :contig_shape_remarks => shape_remarks,
                                :contig_end_repeat_size => repeat_size, :contig_trimmed_end => cend, :contig_full_end => l))
            end
        end

    return df
#end