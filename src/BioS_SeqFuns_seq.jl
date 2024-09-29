export rc_string, gc_string, find_num_contigs, get_contig_names, get_contig_length
export contigSplit, contigLengthSel, contigNameSel, contigExtractRegions, renameContigs, detectContigShapeEnds, get_contig_names_4_contigs_with_seq
export contigFindIdent, contigDerepClusters

"""
   rc_string(String)
   It takes a DNA sequence of type String, converts it to type LongDNA, does reverse complement,
    and then converts it back to type String.
"""
rc_string = String ∘ reverse_complement! ∘ LongDNA{2}

"""
    gc_string(String)
    It takes a DNA sequence of type String, converts it to type LongDNA, and then it calculates and returns the GC content.
"""
function gc_string(seq::String)
    seq = LongDNA{2}(seq)
    GC = gc_content(seq)
    GC = round(GC, digits=2)
    GC = convert(Float16, GC)

    return GC
end


"""
    find_num_contigs()

    It returns the number of contigs in a fasta record.
"""
function find_num_contigs(fasta_p::FnaP)
    a= 0

    FastaReader(fasta_p.p) do FASTA
        for fr in FASTA
            a += 1
        end
    end

    return a
end

"""
    get_contig_names()
    It returns a vector where each element is the comment of a fasta record.
"""
function get_contig_names(fasta_p::FnaP)
    a = Vector{Union{Missing, String}}(missing, find_num_contigs(fasta_p))
    i=0
    FastaReader(fasta_p.p) do FASTA
        for fr in FASTA
            i += 1
            a[i] = fr[1]
        end
    end

    return a
end


function get_contig_names_4_contigs_with_seq(fasta_p::FnaP, seq::Vector{String})
    a = Vector{String}()
    #i=0
    FastaReader(fasta_p.p) do FASTA
        for fr in FASTA
            for j in eachindex(seq)
                if occursin(seq[j], fr[2])
                    #i += 1
                    push!(a, fr[1])
                end
            end
        end
    end

    return a
end

"""
    get_contig_length()
    It returns a vector where each element is the length of a fasta record.
"""
function get_contig_length(fasta_p::FnaP)
    a = Vector{Union{Missing, Int64}}(missing, find_num_contigs(fasta_p))
    i=0
    FastaReader(fasta_p.p) do FASTA
        for fr in FASTA
            i += 1
            a[i] = length(fr[2])
        end
    end

    return a
end



#region modify fasta files
"""
    contigSplit
    Split contigs longer than lengthTh into subcontigs of length lengthTh, 
        and writes all shorter contigs and the split ones to the output file.    
"""

function contigSplit(in_p::FnaP, out_p::FnaP, lengthTh::Int64)
    split_contigs = Dict{String, Vector{String}}()

    FastaWriter(out_p.p) do fw
        FastaReader(in_p.p) do FASTA
            for fr in FASTA
                if length(fr[2]) < lengthTh
                    writeentry(fw, "$(fr[1])", "$(fr[2])")
                else
                    no_subseqs = ceil(length(fr[2]) / lengthTh)
                    no_subseqs = convert(Int64, no_subseqs)

                    subcontings = Vector{String}(undef, no_subseqs)

                    for i in 1:(no_subseqs - 1)
                        subcontings[i] = "$(fr[1])__Subcontig$(i)"

                        start = (i-1)*lengthTh+1 # convert(Int64, ((i-1)*lengthTh+1))
                        stop = i*lengthTh #convert(Int64, (i*lengthTh))                        
                        subseq = fr[2][start:stop]

                        @suppress writeentry(fw, subcontings[i], "$(subseq)")
                    end

                    subcontings[no_subseqs] = "$(fr[1])__Subcontig$(no_subseqs)"

                    start = (no_subseqs-1)*lengthTh+1 #convert(Int64, ((no_subseqs-1)*lengthTh+1))
                    lastsubseq = fr[2][start:end]
                    @suppress writeentry(fw, "$(fr[1])__Subcontig$(no_subseqs)", "$(lastsubseq)")


                    split_contigs[fr[1]] = subcontings
                end
            end
        end
    end

    return split_contigs
end

"""
    contigLengthSel
    It selects contigs with length smaller or equal than maxLen, and writes them to the output file.
"""
function contigLengthSel(in_p::FnaP, out_p::FnaP, maxLen::Int64)
    FastaWriter(out_p.p) do fw
        FastaReader(in_p.p) do FASTA
            for fr in FASTA
                if length(fr[2]) <= maxLen
                    @suppress writeentry(fw, "$(fr[1])", "$(fr[2])")
                end
            end
        end
    end

    return nothing
end

"""
    contigLengthSel
    It selects contigs with length greater or equal than minLen, and writes them to the output file.
"""
function contigLengthSel(minLen::Int64, in_p::FnaP, out_p::FnaP)
    FastaWriter(out_p.p) do fw
        FastaReader(in_p.p) do FASTA
            for fr in FASTA
                if length(fr[2]) >= minLen
                    @suppress writeentry(fw, "$(fr[1])", "$(fr[2])")
                end
            end
        end
    end

    return nothing
end

"""
    contigLengthSel
    It selects contigs with length greater or equal than minLen and smaller or equal than maxLen, 
        and writes them to the output file.
"""
function contigLengthSel(minLen::Int64, in_p::FnaP, out_p::FnaP, maxLen::Int64)
    FastaWriter(out_p.p) do fw
        FastaReader(in_p.p) do FASTA
            for fr in FASTA
                if length(fr[2]) >= minLen && length(fr[2]) <= maxLen
                    @suppress writeentry(fw, "$(fr[1])", "$(fr[2])")
                end
            end
        end
    end

    return nothing
end

function contigNameSel(in_p::FnaP, out_p::FnaP, contnames) #I removed the type of the contig name, because dataframe columns can have weird types sometimes, the problem was PooledVector{String31, UInt32, Vector{UInt32}} ::Vector{T}) where T <: Union{String, InlineString}
    FastaWriter(out_p.p) do fw
        FastaReader(in_p.p) do FASTA
            for fr in FASTA
                if fr[1] in contnames
                    @suppress writeentry(fw, "$(fr[1])", "$(fr[2])")
                end
            end
        end
    end

    return nothing
end

function contigExtractRegions(in_p::FnaP, out_p::FnaP, contig_name, region_name, region_start, region_stop)
    FastaWriter(out_p.p) do fw
        FastaReader(in_p.p) do FASTA
            for c in eachindex(contig_name)
                for fr in FASTA
                    if fr[1] == contig_name[c]
                        region_seq = fr[2][region_start[c]:region_stop[c]]
                        @suppress writeentry(fw, region_name[c], region_seq)
                        break
                    end
                end
            end
        end
    end

    return nothing
end


function contigExtractRegions(in_p::FnaP, out_p::FnaP, contig_name, region_name, region_Lstart, region_Lstop, region_Rstart, region_Rstop)
    FastaWriter(out_p.p) do fw
        FastaReader(in_p.p) do FASTA
            for c in eachindex(contig_name)
                for fr in FASTA
                    if fr[1] == contig_name[c]
                        if ismissing(region_Rstart[c]) == false && ismissing(region_Rstop[c]) == false
                            #println()
                            region_seq = (fr[2][region_Lstart[c]:region_Lstop[c]]) * (fr[2][region_Rstart[c]:region_Rstop[c]])
                        else
                            region_seq = fr[2][region_Lstart[c]:region_Lstop[c]]
                        end
                        @suppress writeentry(fw, region_name[c], region_seq)
                        break
                    end
                end
            end
        end
    end

    return nothing
end

function renameContigs(in_p::FnaP, out_p::FnaP, namesdict::Dict{String, String})
    FastaWriter(out_p.p) do fw
        FastaReader(in_p.p) do FASTA
            for fr in FASTA
                for key in keys(namesdict)
                    if fr[1] == key
                        @suppress writeentry(fw, namesdict[key], "$(fr[2])")
                        break
                    end
                end
            end
        end
    end
    return nothing
end


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

function detectContigShapeEnds(in_p::FnaP)
    df = DataFrame(contig_name = Vector{String}(), contig_shape = Vector{String}(), contig_end_repeat_type = Vector{String}(), contig_shape_remarks = Vector{String}(),
                    contig_end_repeat_size = Vector{Integer}(), contig_trimmed_end = Vector{Int64}(), contig_full_end = Vector{Int64}())
    
    contig_length = get_contig_length(in_p)

        FastaReader(in_p.p) do FASTA
            iterator = 0

            for fr in FASTA
                iterator += 1           
                l = length(fr[2])

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
                    end
                end
                
                push!(df, Dict(:contig_name => fr[1], :contig_shape => shape, :contig_end_repeat_type => end_type, :contig_shape_remarks => shape_remarks,
                                :contig_end_repeat_size => repeat_size, :contig_trimmed_end => cend, :contig_full_end => l))
            end
        end

    return df
end


#endregion modify fasta files

#region cluster contigs
"""
    contigFindIdent(inref::FnaP)
    It detects which contigs are identical, based on their sequence hash, and it assigns them to the same cluster. It returns a dataframe.
"""

function contigFindIdent(inref::FnaP)
    l = find_num_contigs(inref)
    contig_names = fill("missing", l)
    contig_hash = Vector{Union{Missing, UInt64}}(missing, l)
    ident_contigs_dict = Dict{UInt64, Vector{String}}()

    # get hash
    i = 0
    FastaReader(inref.p) do FASTA
        for fr in FASTA
            i += 1
            contig_names[i] = fr[1]
            contig_hash[i] = hash(fr[2])
        end
    end

    # assign cluster
    for a in 1:l
        if !haskey(ident_contigs_dict, contig_hash[a])
            ident_contigs_dict[contig_hash[a]] = Vector{String}()
        end
        push!(ident_contigs_dict[contig_hash[a]], contig_names[a])
    end

    ident_contigs_dict2 =  Dict{String, Vector{String}}()

    for (key, value) in ident_contigs_dict
        ident_contigs_dict2[value[1]] = value[2:end]
    end
  
    contig_hash = nothing
    contig_names = nothing
    ident_contigs_dict = nothing

    return ident_contigs_dict2
end

"""
    contigDerepClusters(df::DataFrame, outfp::FnaP, col_cluster::Symbol, col_names::Symbol)
    It keeps only cluster representatives, and it writes them to the output file. It returns the corresponding dataframe.
"""
function contigDerepClusters(df::DataFrame, inref::FnaP, outfp::FnaP, col_cluster::Symbol, col_names::Symbol)
    dfderep = derepDf(df, col_cluster)
    contigNameSel(inref, outfp, dfderep[!, col_names])

    return dfderep
end
#endregion cluster contigs

