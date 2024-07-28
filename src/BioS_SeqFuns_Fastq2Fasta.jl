#using FASTX
struct Read_orig_label
    read_orig_id::String
    read_orig_desc::String
end

function fastq2fasta(inP::String, outP::FnaP)
    num_reads = 0
    reads_Dict = Dict{String, Read_orig_label}()

    FastaWriter(outP.p) do fw
        reader = open(FASTQ.Reader, inP)
        for record in reader
            num_reads += 1
            read_id = string("read_", num_reads)
            reads_Dict[read_id] = Read_orig_label(FASTQ.identifier(record), FASTQ.description(record))
            writeentry(fw, read_id, FASTQ.sequence(record))
        end
        close(reader)
    end
    return reads_Dict
end