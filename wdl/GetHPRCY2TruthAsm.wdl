version 1.0

workflow GetHPRCY2TruthHaplotypes {
    input {
        File truth_asm_1
        File truth_asm_1_chromAlias
        File truth_asm_2
        File truth_asm_2_chromAlias
        String prefix
    }

    call get_truth_haplotypes_from_annotation {
        input:
            truth_asm_1 = truth_asm_1,
            truth_asm_chromoAlias_1 = truth_asm_1_chromAlias,
            truth_asm_2 = truth_asm_2,
            truth_asm_chromoAlias_2 = truth_asm_2_chromAlias,
            prefix = prefix
    }

    output {
        File fasta_file = get_truth_haplotypes_from_annotation.fasta_file
    }
}


task get_truth_haplotypes_from_annotation {
    input {
        File truth_asm_1
        File truth_asm_chromoAlias_1
        File truth_asm_2
        File truth_asm_chromoAlias_2
        String prefix
    }

    command <<<
        set -euxo pipefail

        python - --truth_fasta_1 ~{truth_asm_1} \
                 --truth_annotation_1 ~{truth_asm_chromoAlias_1} \
                 --truth_fasta_2 ~{truth_asm_2} \
                 --truth_annotation_2 ~{truth_asm_chromoAlias_2} \
                 --output_filename ~{prefix}.truth.fasta \
                 --output_header ~{prefix} \
                 <<-'EOF'
        import pysam
        import gzip
        import argparse

        def get_contig_info(filename):
            with open(filename, "r") as fp:
                data = fp.readlines()
                for line in data:
                    if line.startswith("#"):
                        continue
                    contig, chromo, genebank = line[:-1].split("\t")
                    if "chrM" in chromo:
                        return contig
            return

        def get_truth_seq(fastafile, chromoAliasfile):
            contig = get_contig_info(chromoAliasfile)
            fasta_records = pysam.Fastafile(fastafile)
            chromosomes = fasta_records.references
            truth_seq = []
            for chromosome in chromosomes:
                if contig == chromosome:
                    truth_seq.append(fasta_records.fetch(chromosome))
            return truth_seq

        def write_fasta(output_fasta_file, header, truth_seq):
            with open(output_fasta_file, "w") as f:
                for i, t_seq in enumerate(truth_seq):
                    f.write(">%s_Hap_%d\n" % (header, i+1))
                    char_length = 60
                    line_number = len(t_seq) // char_length
                    for i in range(line_number):
                        f.write(t_seq[i*char_length:(i+1)*char_length] + "\n")
                    f.write(t_seq[line_number*char_length:] + "\n")



        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--truth_fasta_1',
                                type=str)

            parser.add_argument('--truth_annotation_1',
                                type=str)

            parser.add_argument('--truth_fasta_2',
                                type=str)

            parser.add_argument('--truth_annotation_2',
                                type=str)

            parser.add_argument('--output_filename',
                                type=str)

            parser.add_argument('--output_header',
                            type=str)

            args = parser.parse_args()

            truth_seq_1 = get_truth_seq(args.truth_fasta_1, args.truth_annotation_1)
            truth_seq_2 = get_truth_seq(args.truth_fasta_2, args.truth_annotation_2)
            truth_seq = truth_seq_1 + truth_seq_2
            write_fasta(args.output_filename, args.output_header, truth_seq)

        if __name__ == "__main__":
            main()
        EOF
        
    >>>

    output {
        File fasta_file = "~{prefix}.truth.fasta"
        
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/slee/kage-lite:pr_29"
        memory: "8 GB"
        cpu: 2
        disks: "local-disk 100 SSD"
    }
}
