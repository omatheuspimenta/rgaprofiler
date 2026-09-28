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


// Named like the real module so conf/base.config's per-tool withName block applies too.
process SIGNALP6 {
    label 'process_medium_chunk'
    label 'process_gpu'
    input:
    path fasta
    output:
    stdout
    script:
    """
    echo "TOOL SIGNALP6 cpus=${task.cpus} mem_gb=${task.memory.toGiga()} time_h=${task.time.toHours()} vram_mb=${task.ext.gpu_vram_mb} longest=${task.ext.chunk.longest}"
    """
}

// Named like the real module so conf/base.config's per-tool withName block applies too.
process DEEPLOC2 {
    label 'process_medium_chunk'
    label 'process_gpu'
    input:
    path fasta
    output:
    stdout
    script:
    """
    echo "TOOL DEEPLOC2 cpus=${task.cpus} mem_gb=${task.memory.toGiga()} time_h=${task.time.toHours()} vram_mb=${task.ext.gpu_vram_mb} longest=${task.ext.chunk.longest}"
    """
}

// Named like the real module so conf/base.config's per-tool withName block applies too.
process DEEPTMHMM {
    label 'process_medium_chunk'
    label 'process_gpu'
    input:
    path fasta
    output:
    stdout
    script:
    """
    echo "TOOL DEEPTMHMM cpus=${task.cpus} mem_gb=${task.memory.toGiga()} time_h=${task.time.toHours()} vram_mb=${task.ext.gpu_vram_mb} longest=${task.ext.chunk.longest}"
    """
}

// Named like the real module so conf/base.config's per-tool withName block applies too.
process DEEPCOIL2 {
    label 'process_medium_chunk'
    label 'process_gpu'
    input:
    path fasta
    output:
    stdout
    script:
    """
    echo "TOOL DEEPCOIL2 cpus=${task.cpus} mem_gb=${task.memory.toGiga()} time_h=${task.time.toHours()} vram_mb=${task.ext.gpu_vram_mb} longest=${task.ext.chunk.longest}"
    """
}

workflow {
    ch = Channel.fromPath(params.chunk)
    CHUNK_MEDIUM(ch).view()
    CHUNK_HIGH(ch).view()
    SIGNALP6(ch).view()
    DEEPLOC2(ch).view()
    DEEPTMHMM(ch).view()
    DEEPCOIL2(ch).view()
}
