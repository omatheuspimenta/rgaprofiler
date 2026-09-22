// Minimal probe pipeline for the per-chunk resource-sizing and GPU-lock logic in
// conf/base.config. It has no containers and runs nothing heavy: each process just
// reports the cpus/memory that base.config's lazy directive closures resolved for it,
// so tests/bin/test_resource_sizing.sh can assert on them under mocked hosts.
// (`bin` next to this file is a symlink to the repo's bin/, so ${projectDir}/bin/...
// in base.config resolves the same way it does for the real pipeline.)

process CHUNK_MEDIUM {
    label 'process_medium_chunk'
    label 'process_gpu'
    input:
    path fasta
    output:
    stdout
    script:
    """
    echo "CHUNK_MEDIUM cpus=${task.cpus} mem_gb=${task.memory.toGiga()} use_gpu=${task.ext.use_gpu}"
    """
}

process CHUNK_HIGH {
    label 'process_high_chunk'
    input:
    path fasta
    output:
    stdout
    script:
    """
    echo "CHUNK_HIGH cpus=${task.cpus} mem_gb=${task.memory.toGiga()}"
    """
}

workflow {
    ch = Channel.fromPath(params.chunk)
    CHUNK_MEDIUM(ch).view()
    CHUNK_HIGH(ch).view()
}
