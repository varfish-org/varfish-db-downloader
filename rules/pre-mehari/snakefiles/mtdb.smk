def input_mitomap_tsv(wildcards):
    return {
        "ref": f"work/reference/{wildcards.genomebuild.lower()}/reference.fa",
        "txt": "rules/pre-mehari/tools/data/mtdb.tsv",
        "header": "rules/pre-mehari/header/mtdb.txt",
    }


rule result_GRChXX_mtdb_tsv:
    input:
        unpack(input_mitomap_tsv),
    output:
        tsv="output/pre-mehari/{genomebuild}/mtDB/{download_date}/MtDb.tsv",
        release_info="output/pre-mehari/{genomebuild}/mtDB/{download_date}/MtDb.release_info",
    run:
        import os
        import binning
        import csv

        from Bio import SeqIO

        if os.environ.get("CI", "false").lower() == "true":
            print("CI is true, skipping rule execution.")
            open(output.tsv, "w").close()
            open(output.release_info, "w").close()
            return

        if wildcards.genomebuild == "GRCh37":
            chrom = "MT"
        else:
            chrom = "chrM"


        with open(input.header, "r") as fh_header, open(input.txt, "r") as fh_input, open(
            input.ref, "r"
        ) as fasta:
            fasta_records = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
            mitochondrium = fasta_records[chrom].seq
            header = [fh_header.read().strip().split("\n")]
            result = []
            line_header = next(fh_input).strip("\n").split("\t")
            for row in fh_input:
                line = row.strip("\n").split("\t")
                record = dict(zip(line_header, line))
                alts = dict(A=0, C=0, G=0, T=0)
                start = int(record.get("Posn."))
                reference = record.get("Base")
                an = int(record[reference] or 0) + int(record["Gap"] or 0)
                # Sanity check ref base against mitochondrial genome and correct
                if not mitochondrium[start - 1] == reference:
                    if mitochondrium[start - 1] == "N":
                        print("Skipping position {} because ref is N".format(start))
                        continue
                    else:
                        print(
                            "Bases at position {} do not match, correcting: {} -> {}".format(
                                start, reference, mitochondrium[start - 1]
                            )
                        )
                        reference = mitochondrium[start - 1]
                        # Collect total number of alleles and all alternatice allele counts
                for alt in set(alts).difference(reference):
                    alts[alt] = int(record.get(alt) or 0)
                    an += alts[alt]
                    # Is an alternative registered, add it to the results list (normalization step)
                for alt, ac in alts.items():
                    if ac > 0:
                        _af = ac / an
                        # Warn if frequency is at 100%.
                        if _af == 1.0:
                            print("Frequency at position {} is 100%!".format(start))
                        result.append(
                            [
                                wildcards.genomebuild,  # release
                                chrom,  # chromosome
                                start,  # start
                                start,  # end
                                binning.assign_bin(start - 1, start),  # bin
                                reference,  # reference
                                alt,  # alternative
                                record.get("Location"),  # location
                                record.get("Codon"),  # codon
                                record.get("Position"),  # position
                                record.get("Amino Change"),  # aa_change
                                None
                                if not record.get("Syn?")
                                else record.get("Syn?") == "Yes",  # synonymous
                                ac,  # ac
                                an,  # an
                                _af,  # af
                            ]
                        )
            result.sort()

        with open(output.tsv, "w") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerows(header)
            writer.writerows(result)

        with open(output.release_info, "w") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerows(
                [
                    ["table", "version", "genomebuild", "null_value"],
                    ["MtDb", wildcards.download_date, wildcards.genomebuild, ""],
                ]
            )
