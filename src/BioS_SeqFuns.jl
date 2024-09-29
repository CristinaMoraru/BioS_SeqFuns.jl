module BioS_SeqFuns

using DataFrames
using BioSequences
using FASTX #.FASTQ: Reader
using FastaIO: FastaWriter, writeentry, FastaReader
using CSV
using Suppressor

#my packages
using BioS_Gen

include("BioS_SeqFuns_DataTypes.jl")
include("BioS_SeqFuns_seq.jl")
include("BioS_SeqFuns_Fastq2Fasta.jl")
include("BioS_SeqFuns_Append2Seq.jl")
include("BioS_SeqFuns_Fasta2MFasta.jl")
include("BioS_SeqFuns_clust.jl")


end # module BioS_SeqFuns
