rule md5sum:
    input:
        "{file}"
    output:
        "{file}.md5"
    log:
        "logs/md5sum/{file}.log"
    shell:
        r"""
        md5sum {input} > {output} 2> {log}
        """